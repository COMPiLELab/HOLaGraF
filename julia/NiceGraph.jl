#using LightGraphs
#using SimpleWeightedGraphs

mutable struct NiceGraph
    n::Integer
    edges
    trigs
    w
    eps0
    e
    B1
    B2
    points
    
   NiceGraph(n, edges, trigs, w, eps0, e, points=nothing)=new(
       n, edges, trigs, w, eps0, e, 
       B1fromEdges(n, edges),
       B2fromTrig(edges, trigs),
       getPositions(points, w, eps0, e, B1fromEdges(n, edges))
   )
end

function getPositions(points, w, eps0, e, B1)
    if isnothing(points)
        W=Diagonal( vec(sqrt.(w)+eps0*e) );
        Aw = Diagonal(diag(B1 * W * W * B1'))-B1 * W * W * B1';
        points=springLayout(Aw);
    else
        points=points;
    end
    return points
end

function B1fromEdges(n, edges)
    m = size(edges, 1);
    B1 = spzeros(n, m);
    
    for i in 1:m
        B1[edges[i, 1], i] = -1;
        B1[edges[i, 2], i] = 1;
    end
    return B1
end

function B2fromTrig(edges, trigs)
    
    m = size(edges, 1);
    if length(trigs)==1
        return spzeros(m, 1)
    end
    del = size(trigs, 1);
    B2 = spzeros(m, del);
    
    for i in 1:del
        B2[findfirst( all( [trigs[i, 1], trigs[i, 2]]' .== edges, dims=2 )[:, 1] ), i]=1;
        B2[findfirst( all( [trigs[i, 1], trigs[i, 3]]' .== edges, dims=2 )[:, 1] ), i]=-1;
        B2[findfirst( all( [trigs[i, 2], trigs[i, 3]]' .== edges, dims=2 )[:, 1] ), i]=1;
    end 
    return B2
end

function getAdjB1W(G::NiceGraph)
    W = getW(G);
    Aw = Diagonal(diag(G.B1 * W * W * G.B1'))-G.B1 * W * W * G.B1';
    return Aw
end

function getAdjB1(G::NiceGraph)
    A=Diagonal(diag(G.B1 * G.B1'))-G.B1 * G.B1';
    return A
end

function getW(G::NiceGraph)
    W=Diagonal(sqrt.(G.w)+G.eps0*G.e);
    return W
end


function springLayout(adj_matrix)

    nvg = size(adj_matrix, 1)
    locs_x=2*rand(nvg).-1.0;
    locs_y=2*rand(nvg).-1.0;
    C=2.0;
    MAXITER=100;
    INITTEMP=2.0;

    k = C * sqrt(4.0 / nvg)
    k² = k * k

    force_x = zeros(nvg)
    force_y = zeros(nvg)
    @inbounds for iter = 1:MAXITER
    for i = 1:nvg
        force_vec_x = 0.0
        force_vec_y = 0.0
        for j = 1:nvg
            i == j && continue
            d_x = locs_x[j] - locs_x[i]
            d_y = locs_y[j] - locs_y[i]
            dist²  = (d_x * d_x) + (d_y * d_y)
            dist = sqrt(dist²)

            if !( iszero(adj_matrix[i,j]) && iszero(adj_matrix[j,i]) )
    
                F_d = dist / k - k² / dist²
            else
                F_d = -k² / dist²
            end
            force_vec_x += F_d*d_x
            force_vec_y += F_d*d_y
        end
        force_x[i] = force_vec_x
        force_y[i] = force_vec_y
    end
    # Cool down
    temp = INITTEMP / iter
    # Now apply them, but limit to temperature
    for i = 1:nvg
        fx = force_x[i]
        fy = force_y[i]
        force_mag  = sqrt((fx * fx) + (fy * fy))
        scale      = min(force_mag, temp) / force_mag
        locs_x[i] += force_x[i] * scale
        locs_y[i] += force_y[i] * scale
        end
    end

    min_x, max_x = minimum(locs_x), maximum(locs_x)
    min_y, max_y = minimum(locs_y), maximum(locs_y)
    function scaler(z, a, b)
         2.0*((z - a)/(b - a)) - 1.0
    end
    map!(z -> scaler(z, min_x, max_x), locs_x, locs_x)
    map!(z -> scaler(z, min_y, max_y), locs_y, locs_y)

    return locs_x, locs_y
end