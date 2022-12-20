
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