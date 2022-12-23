module HOLaGraF_lsmr

      export myEig,  readMatrix, NiceGraph, getW1, getW0, getW2, getL0, getL1, getL1up, getP, getK, getJ10, getJ12, F, getGrad, EulerStep, checkStep, innerLevel, alphaLevel, freeGradientTransition, doMaxPool, Thresh, placer, placeL1up

      using Random
      using DelimitedFiles
      using LinearAlgebra, Arpack
      using SparseArrays
      using LinearMaps
      using LinearOperators
      using Krylov
      using ArnoldiMethod

      using Printf

      include("NiceGraph.jl")
      include("thrs_struct.jl")

      function placer(Mat, C)
            m = size(Mat, 1);
            return Afun5(x) = Krylov.lsmr(Mat, x, M = C, N = I(m), ldiv = true)[1]
      end
      
      function placeL1up(C)
            return Afun6(Mat) = placer(Mat, C)
      end
      
      using SuiteSparse
      import LinearAlgebra.ldiv!
      using IncompleteLU, LimitedLDLFactorizations, ILUZero
      
      ldiv!(y::Vector{T}, F::SuiteSparse.CHOLMOD.Factor{T}, x::Vector{T}) where T = (y .= F \ x)
      ldiv!(F::SuiteSparse.CHOLMOD.Factor{Float64}, x::Vector{Float64}) = ( ldiv!(x, F, x) )
      ldiv!(y::Vector{T}, F::LinearOperator, x::Vector{T}) where T = (y .= F * x)
      ldiv!(F::LinearOperator, x::Vector{Float64}) = ( ldiv!(x, F, x) )
      ldiv!(F::Matrix{Float64}, x::Vector{Float64}) = ( ldiv!(x, F, x) )
      

      function myEig(L1up, inFun)
            Afun = inFun(L1up);
            m = size(L1up, 1);
            D = LinearMap(
                  Afun, Afun, m, m; ismutating = false, issymmetric = true
            );
            decomp,  = partialschur(D, nev=3, tol=1e-6, which=LM());
            λs_inv, X = partialeigen(decomp);
            λs = 1 ./ λs_inv;
            return real(λs[3]), real.(X[:, 3])
      end
      
      
      myinv(A) = sparse(pinv(Matrix(A)));
      Sym(A) = 0.5 * (A' + A);
 
      function readMatrix(filename::String)
            pre=readdlm(filename, Int);
            pre=sort(pre, dims = 2); 
            simplices=sortslices(pre, dims = 1);
            return simplices
      end

      getW1(G::NiceGraph) = Diagonal( vec(sqrt.(G.w) + G.eps0*G.e) );

      getW0(W1::Diagonal, G::NiceGraph; ρ=1.0) = Diagonal( vec( ρ*ones(size(G.B1, 1), 1) + abs.(G.B1)*diag(W1) ) );

      function getW2(W1::Diagonal, G::NiceGraph)
            (size(G.B2, 2) == 0) && return 0;
            weightedB2 = abs.(G.B2)'*W1;
            weightedB2[ weightedB2 .== 0 ] .= sum(W1);
            return Diagonal( vec( minimum(weightedB2; dims=2) ) )
      end

      function getL0(G::NiceGraph; ρ=1.0)
            W1 = getW1(G);
            W0 = getW0(W1, G; ρ=ρ);
            barB1 = myinv(W0) * G.B1 * W1;
            return barB1 * barB1'
      end
      
      function getL1(G::NiceGraph; ρ=1.0)
            W1 = getW1(G);
            W0 = getW0(W1, G; ρ=ρ);
            W2 = getW2(W1, G);

            barB1 = myinv(W0) * G.B1 * W1;
            barB2 = myinv(W1) * G.B2 * W2;

            return barB1' * barB1 + barB2 * barB2'
      end

      function getL1up(G::NiceGraph; ρ=1.0)
            W1 = getW1(G);
            W2 = getW2(W1, G);

            barB2 = myinv(W1) * G.B2 * W2;
            return barB2 * barB2'
      end

      getP(Mat; thr=1e-3) = sum( ( abs.(Mat) * ones(size(Mat, 1)) ) .< thr );

      getK(L1; thr=1e-6)= sum( (eigs(L1, nev=size(L1, 1)-1, which=:SR)[1]).< thr);


      getJ10(G::NiceGraph) = abs.(G.B1);

      function getJ12(W1::Diagonal, W2::Diagonal, G::NiceGraph)
            J12=spzeros(size(G.B2', 1), size(G.B2', 2));
            if length(G.trigs)==1
                    return J12
            end
            for i in axes(J12, 1)
                    ix=findfirst( abs.(diag(W1).-W2[i,i]).<1e-08);
                    J12[i, ix]=1;
            end
            return J12
      end

      function F(L1up, L0, thrs, inFun )
            vals, ~ =  myEig(L1up, inFun); λ = vals[1];
            μ = eigen( Symmetric(Matrix(L0)), 2:2).values[1];            
            return 0.5 * λ^2 + 0.5 * thrs.alph * (max(0, 1 - μ/thrs.mu )) ^2
      end

      function getGrad(L1up, L0, G::NiceGraph, thrs, inFun; normcor=true, ρ=1.0, thr0=1e-5)
            vals, vecs = myEig(L1up, inFun);
            λ = real(vals[1]); v = real.(vecs[:, 1]);
            vvT = v*v';

            #get all the weights
            W1 = getW1(G);
            W0 = getW0(W1, G; ρ=ρ);
            W2 = getW2(W1, G);

            #get jacobians
            J10 = getJ10(G);
            J12 = getJ12(W1, W2, G);

            # grad of L1
            term1 = 2 * Sym(- myinv(W1) * vvT * myinv(W1)*G.B2*W2^2*G.B2' * myinv(W1) );
            term1 = term1 + 2 * Diagonal( vec( J12' * diag( G.B2' * myinv(W1) * vvT * myinv(W1) * G.B2 * W2  ) ) ); 
            term1 = λ * term1;

            # get μ_2 and v_2 for L0
            res = eigen( Symmetric(Matrix(L0)), 2:2);
            μ = real(res.values[1]); v = real.( res.vectors[:, 1]);  vvT = v*v';

            # get penalization's gradient
            if thrs.mu>μ
                    term2 = - 2 * thrs.alph/thrs.mu*max(0, 1-μ/thrs.mu) * ( G.B1' * myinv(W0) * vvT * myinv(W0) * G.B1 * W1 - Diagonal( vec( J10' * diag( Sym(L0*vvT*myinv(W0))) ) ) );
                    term1=term1+term2;
            end

            # sum 2 terms
            grad = diag(term1);

            # norm correction
            matmask = ( .!(abs.(sqrt.(G.w)+G.e*G.eps0) .< thr0) );
            PE = G.e.*matmask;
            κ = (normcor ?  dot(grad.* matmask, PE)/dot(PE, PE) : 0);
            
            grad = grad .* matmask-κ*PE;

            return grad 
      end

      function EulerStep(e, ∇F, h::Float64, w, eps0::Float64; correction = true )
            e1 = e - h * ∇F; # step
            e1[ sqrt.(w) + eps0*e1 .< 0] = -1.0/eps0 * sqrt.(w[sqrt.(w) + eps0*e1 .< 0]); # non-negativity correction
            return e1/(  correction ?  norm(e1, 2) : 1.0 )# norm correction
      end

      function checkStep( G::NiceGraph, L0, L1, L1up, Fk::Float64, e1, p, thrs::Thresh, h, prevAccepted, inFun; ρ = 1.0,  β = 1.2 )
            e0 = G.e; # backup the inital perturbation
            G.e = e1;  #switch to new perturbation

            L0_1, L1_1, L1up_1 = getL0(G; ρ=ρ), getL1(G; ρ=ρ), getL1up(G; ρ=ρ);  # get updated laplacians 
            p1 = getP(L1_1); # get updated p
            (p1 != p) ? (  G.e = e0; return true, p1, h, L0_1, L1_1, L1up_1 ) : 0;  # accept if p changed
            Fk1=F(L1up_1, L0_1, thrs, inFun); # get new value for monotonicity
            ( Fk > Fk1 ) ?  ( G.e = e0; return true, p, h * ( ( prevAccepted ) ? β : 1), L0_1, L1_1,  L1up_1 )  :  (G.e=e0;  return false, p, h/β, L0, L1, L1up  );         
      end      

      function innerLevel( G::NiceGraph, h0::Float64,  thrs::Thresh, p, inFun; β = 1.2, thr_∇ = 1e-4 , thr_F = 1e-5, ρ = 1.0, max_iter = 50 )
            h = h0; prevAccepted = false; # intiialisation
            h_log = Vector{Float64}(); #remeber log of h's 
            track = Vector{Float64}(); # remember log of functionals
            e_log = Array{Float64}(undef, size(G.e, 1), 0); # log of the perturbations
            λ_log = Array{Float64}(undef, size(G.e, 1), 0); # log of the spectrum
            L0, L1, L1up = getL0(G; ρ=ρ), getL1(G; ρ=ρ), getL1up(G; ρ=ρ); track = [track; F(L1up, L0, thrs, inFun) ];  # first functional
            e_log = [e_log G.e]; #λ_log = [ λ_log eigen(Symmetric(Matrix(L1))).values ]; 
            h_log = [h_log; h0];
            while (track[end]>thr_F) && (size(track, 1) < max_iter )
                    ∇F = getGrad(L1up, L0, G, thrs, inFun);
                    ( maximum(abs.(∇F)) < thr_∇) ? break : 0;
                    e1 = zeros( size(G.e ) );
                    while  true
                            e1 = EulerStep(G.e, ∇F, h, G.w, G.eps0);
                            flag, p1, h, L0, L1, L1up = checkStep(G, L0, L1, L1up, track[end], e1, p, thrs, h, prevAccepted, inFun; β=β);
                            flag ?  (  (p1 == p) ?  prevAccepted=true : (prevAccepted=false; p1=p;)  ; G.e = e1;  break) : 0;
                    end
                    h_log = [h_log; h]; e_log = [e_log G.e]; 
                    track = [track; F(L1up, L0, thrs, inFun)];
            end

            return G.e, track, h_log, G, L1up, L1, L0, p, e_log, λ_log
      end

      function freeGradientTransition(ε_target, G::NiceGraph, h0::Float64,  thrs::Thresh, p,  inFun ; β = 1.2, thr_∇ = 1e-4 , thr_F = 1e-5, ρ = 1.0, max_iter = 100 )
            h = h0; prevAccepted = false; # intiialisation
            h_log = Vector{Float64}(); #remeber log of h's 
            track = Vector{Float64}(); # remember log of functionals
            e_log = Array{Float64}(undef, size(G.e, 1), 0); # log of the perturbations
            λ_log = Array{Float64}(undef, size(G.e, 1), 0); # log of the spectrum
            L0, L1, L1up = getL0(G; ρ=ρ), getL1(G; ρ=ρ), getL1up(G; ρ=ρ); track = [track; F(L1up, L0, thrs, inFun) ];  # first functional
            e_log = [e_log G.e]; #λ_log = [ λ_log eigen(Symmetric(Matrix(L1))).values ]; 
            h_log = [h_log; h0];
            while (track[end]>thr_F) && (size(track, 1) < max_iter ) && (G.eps0*norm(G.e, 2) < ε_target)
                    ∇F = getGrad(L1up, L0, G, thrs, inFun; normcor = false);
                    ( maximum(abs.(∇F)) < thr_∇) ? break : 0;
                    while  true
                            e1 = EulerStep(G.e, ∇F, h, G.w, G.eps0; correction = false);
                            flag, p1, h, L0, L1, L1up = checkStep(G, L0, L1, L1up, track[end], e1, p, thrs, h, prevAccepted, inFun; β=β);
                            flag ?  (  (p1 == p) ?  prevAccepted=true : (prevAccepted=false; p1=p;)  ; G.e = e1; break) : 0;
                    end
                    h_log = [h_log; h]; e_log = [e_log G.e]; 
                    track = [track; F(L1up, L0, thrs, inFun)];
            end
            G.eps0 = G.eps0 * norm(G.e, 2);
            G.e=G.e / norm(G.e, 2);
            return G.e, track, h_log, G, L1up, L1, L0, p, e_log, λ_log 
      end

      function alphaLevel(  G::NiceGraph,  h0::Float64, α_st::Float64, α_fin::Float64, thrs::Thresh, inFun  ; α_num=15, initial=false, ρ=1.0,  β = 1.2, thr_∇ = 1e-4 , thr_F = 1e-5, max_iter = 2000)
            if initial
                    thrs.alph = α_st;
                    ε0 = G.eps0;
                    G.e = ones(size(G.e)); G.e = G.e / norm(G.e, 2);
                    G.eps0 = 1e-6;
                    L0_1, L1_1, L1_1up = getL0(G; ρ=ρ), getL1(G; ρ=ρ), getL1up(G; ρ=ρ); p=getP(L1_1);
                    e0=getGrad(L1_1up, L0_1, G, thrs, inFun);
                    e0=-e0/norm(e0, 2); G.e=e0; G.eps0=ε0;
            end

            L0, L1, L1up = getL0(G; ρ=ρ), getL1(G; ρ=ρ), getL1up(G; ρ=ρ); 
            p = getP(L1);

            αs=range(sqrt(α_st), sqrt(α_fin), length=α_num).^2;
            ans, track, h_log, G,  L1up, L1, L0, p, e_log, λ_log = innerLevel(G, h0, thrs, p, inFun);
            
            Hs = Vector{Float64}(); Es = Array{Float64}(undef, size(G.e, 1), 0); Λs = Array{Float64}(undef, size(G.e, 1), 0);
            SZs = Vector{Float64}(); TrackS = Vector{Float64}();
            Hs = [Hs; h_log]; Es = [Es e_log]; Λs=[Λs λ_log]; SZs = [SZs; size(h_log, 1)]; TrackS = [TrackS; track ];
            i_st = (initial ? 2 : α_num+1 );
            for i in i_st : α_num
                    G.e = ans; thrs.alph = αs[i];
                    ans, track, h_log, G, L1up, L1, L0, p, e_log, λ_log = innerLevel(G, h0, thrs, p, inFun );
                    Hs = [Hs; h_log]; 
                    Es = [Es e_log]; 
                    Λs = [ Λs λ_log]; 
                    SZs = [ SZs; size(h_log, 1) ]; 
                    TrackS = [TrackS; track]
            end
            
            G.e=ans;
            return G, thrs, track, p, Hs, Es, Λs, SZs, TrackS
      end

      function doMaxPool(G::NiceGraph)
            e_new=zeros(size(G.e));
            e_new[G.e.<-0.2].=-1;
            e_new=e_new/norm(e_new, 2);
            eps_new=-sum(sqrt.(G.w).*e_new);
        
            return eps_new, e_new
      end      

end