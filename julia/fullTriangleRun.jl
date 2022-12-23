using LinearAlgebra, Arpack, Random, SparseArrays
using BenchmarkTools, Printf, TimerOutputs

using LinearOperators, Krylov

include("HOLaGraF_lsmr.jl")
using .HOLaGraF_lsmr
include("generateDelauney.jl")

using ArnoldiMethod
using DataFrames, CSV

using Preconditioners, LinearOperators, LinearMaps, Krylov
using SuiteSparse
import LinearAlgebra.ldiv!
using IncompleteLU, LimitedLDLFactorizations, ILUZero
using ArnoldiMethod

function wrapper(G::NiceGraph, h0, α_st, α_fin, thrs::Thresh, h_ε, Es, Fs, inFun )
        
        while true

                (G.eps0 > 5.0) && break;

                if G.eps0 == h_ε
                        @printf "             alpha  stage started: \n";
                        G, thrs, track, p = alphaLevel(G, h0, α_st, α_fin, thrs, inFun; initial=true);
                        Fs = [Fs; track[end]];
                        Es = [Es G.e];
                        @printf " eps:  %f       ||||    functional :  %f  \n" G.eps0 track[end];
                        @printf "              alpha stage ended \n \n";
                else
                        G , thrs, track, p= alphaLevel(G, h0, α_st, α_fin, thrs, inFun);
                        Fs = [Fs; track[end]];
                        Es = [Es G.e];
                        @printf " <contrained>    eps = %f    ||||  functional :  %f  || %d steps \n" G.eps0 track[end] size(track, 1);
                end

               (Fs[end] < 1e-5) && break; 
                
                h0 = 0.1;
                ε_target = G.eps0 + h_ε; 
                ans, track, h_log, G, L1up, L1, L0, p = freeGradientTransition(ε_target, G, h0, thrs, p, inFun);
                G.e=ans;
                @printf  "  <free_grad>  eps_fin = %f  / %f        |||| functional:   %f      || steps: %d  \n"  G.eps0 ε_target track[end] size(track, 1);
     
        end
        return G, thrs, Es, Fs
end

function triangSampleRun(ν::Float64, N::Int)
        h0=1.0; n=N+4;
    
        edgs=zeros(Int32, 3*n-3-4, 2); trigs=zeros(Int32, 2*n-2-4, 3);
        points=(zeros(1, n), zeros(1, n));
        points, edgs, trigs = generateDelauney(N);

        indx=getIndx2Kill(edgs);
        edges2, trigs2=killEdge(indx, n, edgs, trigs);
        indx=getIndx2Kill(edges2);
        edges2, trigs2=killEdge(indx, n, edges2, trigs2);
        addNum=Int(3*round((ν*n*(n-1)/2)/3))-size(edges2, 1);
        if addNum<0
                @printf "DUDE, NOTHINS'S LEFT, CALM DOWN \n"
                return 0, 0, 0, 0, 0, 0
        end
        for rep2 in 1:addNum
                new_edge = getNewEdge(n, edges2);
                edges2, trigs2=addEdge(new_edge, n, edges2, trigs2);
        end

        w=rand(size(edges2, 1))*0.75 .+ 0.25;
        eps0=0; e=zeros(size(w, 1));
        G=NiceGraph(n, edges2, trigs2, w, eps0, e, points);

        L0, L1, L1up = getL0(G), getL1(G), getL1up(G);
        k_old, p, k = getK(L1), getP(L1), getK(L1up)

        m = size(L1up, 1);
        α = 0.001;
        C = cholesky( L1up + α * Diagonal(diag(L1up)) + 1e-8 * I(m) );
        op = opCholesky( L1up + α * Diagonal(diag(L1up)) + 1e-8 * I(m) );
        inFun = placeL1up(C);        # use preconditioner
        #inFun = placeL1up(I(m));    # do not use the preconditioner

        μ = eigs(L0, nev = 2, which = :SR)[1][2];
        thrs = Thresh( 0.75*μ, 100.0 );

        h0 = 1.0;
        h_ε=0.05;
        G.eps0 = h_ε;

        Es = Array{Float64}(undef, size(G.w, 1), 0);
        Es = [ Es G.e ];
        Fs = Vector{Float64}();

        @printf("\n sampled! \n");

        α_st, α_fin = 1.0, 100.0;
        
        times=@elapsed begin
                G, thrs, Es, Fs = wrapper(G, h0, α_st, α_fin, thrs, h_ε, Es, Fs, inFun);
                eps_MP, e_MP=doMaxPool(G);
        end

        return times, eps_MP, sum(e_MP .< -0.2), k, size(L1, 1), size(G.B2, 2)
end
    
Ns=[22, 28, 34, 40 ];
νs=[0.35; 0.5; ];
rep = 10;

mat = Array{Float64}(undef, 0, 8);
filename = "julia_lsmr_chol.csv";

for i in axes(νs, 1)
        for j in axes(Ns, 1)
                global Ns, νs, mat, rep

                ν=νs[i]; N=Ns[j];
                for ii in 1:rep
                        times, eps_ans, numE, k, numEdg, numTrig=triangSampleRun(ν, N);
                        mat = [mat; ν N times eps_ans numE k numEdg numTrig ]
                        df = DataFrame([ν N times eps_ans numE k numEdg numTrig]);
                        df = rename(df, ["nu", "N", "t", "eps", "num", "k", "numEdg", "numTrig"])

                        if isfile(filename)
                              old = CSV.read(filename, DataFrame);
                              new = [old; df];
                              CSV.write(filename, new);
                        else
                              CSV.write(filename, df);
                        end

                        @printf "number of vertices: %d  ||  sparsity:  %0.2f \n" N ν;
                        @printf "time: %0.2f     %d elim.edges   : total eps=%0.2f \n" times numE eps_ans;
                end
                @printf "\n";
        end
end









    