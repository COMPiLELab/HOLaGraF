using BenchmarkTools, Printf
using LinearAlgebra, Arpack

include("HOLaGraF_lsmr.jl")
using .HOLaGraF_lsmr

n = 8;
edges = readMatrix("julia/example/8.edges")
trigs = readMatrix("julia/example/8.trigs")
w = reshape([0.5; 0.95; 1.0; 1.0; 1.0; 1.1; 1.1; 1.1; 1.05; 0.31; 1.05; 0.9].^2, :, 1);
ε0 = 1e-8; e = -ones(size(w, 1)); e = e/norm(e, 2);

G=NiceGraph(n, edges, trigs, w, ε0, e, nothing);


L0 = getL0(G);
μ = eigs(L0, nev = 2, which = :SR)[1][2];
thrs = Thresh( 0.75*μ, 1.0 );

inFun = placeL1up(I(size(G.edges, 1)));

include("wrapper.jl")

h0 = 0.1;
h_ε=0.025;
G.eps0 = h_ε;

logSizes = Vector{Float64}(); logSteps = Vector{Float64}(); logTrack = Vector{Float64}();
logE = Array{Float64}(undef, size(G.w, 1), 0); logΛ = Array{Float64}(undef, size(G.w, 1), 0);

α_st, α_fin = 1.0, 100.0;

@time G, thrs, logSizes, logSteps, logTrack, logE, logΛ = wrapper(G, h0, α_st, α_fin, thrs, h_ε, logSizes, logSteps, logTrack, logE, logΛ, inFun );












