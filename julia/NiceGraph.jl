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

