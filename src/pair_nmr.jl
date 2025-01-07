# using LinearAlgebra 

function G(nframes::Int64,maxlag::Int64,
           x1::Array{Float64,2},x2::Array{Float64,2},
           boxsize::Array{Float64,2},samplingvector::Array{Float64,1})

    (nrows,ncols) = size(x1)
    (nrows_box,ncols_box) = size(boxsize)

    if nframes != nrows_box
        return error("Inconsistency between number of frames and available box sizes")
    end 

    dx = x2 - x1
    phi_ij = zeros(Float64,nframes)

    for frame in 1:nframes

        currentbox = boxsize[frame,:]

        # println(size(dx[frame,:]), size(currentbox))

        dxcurrent = dx[frame,:] ./ currentbox

        wrap = round.(dxcurrent,RoundNearestTiesUp)
        dxcurrent = dxcurrent - wrap 

        dxcurrent = dxcurrent .* currentbox 
        rij = norm(dxcurrent)

        # if div(frame,100) == 0
           # println(rij)
        # end 

        uij = dxcurrent / rij

        # println(norm(uij))

        costheta = dot(uij,samplingvector)

        # if div(frame,100) == 0
           # println(costheta, uij)
        # end 

        phi_ij[frame] = (3 * costheta^2 - 1) / rij^3

    end 

    return acfNN(phi_ij,maxlag)

end
