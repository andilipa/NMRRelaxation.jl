# using MDTools
# using StatsBase
# using Combinatorics
# using Statistics
# include("/Users/dan/Julia/pair_NMR.jl")

function Inter(nsample::Int64, maxlag::Int64, sample_residues::Vector{Int64}, deltat::Float64,
               boxsize::Array{Float64,2}, nframes::Int64, sampling_vector::Vector{Vector{Float64}},
               MaxResidue::Int64, atoms_per_residue::Int64, base_h_indices::Vector{Int64}, ta)

    y = zeros(Float64,maxlag,nsample)
    nSourceHydrogens = length(base_h_indices)

    for i in 1:nsample

        sourceHydrogens = (sample_residues[i] - 1) * atoms_per_residue .+ base_h_indices;

        # For residue j, the target is 1:(j-1), and (j+1):MaxResidue
        targetResidues = [collect(1:(sample_residues[i]-1)); 
                        collect((sample_residues[i]+1):MaxResidue)]

        targetHydrogens = Vector{Int64}[];

        for j in targetResidues
            push!(targetHydrogens,(j-1) * atoms_per_residue .+ base_h_indices);
        end 
        targetHydrogens = reduce(vcat,targetHydrogens)

        tmp = zeros(maxlag)

        for k in sourceHydrogens

            x1 = ta[:,k].xyz;

            for l in targetHydrogens

                x2 = ta[:,l].xyz

                acf = [zeros(maxlag) for _ in 1:4]

                Threads.@threads for m in 1:4
                    sampleV = sampling_vector[m]
                    acf[m] = G(nframes,maxlag,x1,x2,boxsize,sampleV)
                end     

                for j in 1:4
                    tmp += acf[j]
                end     

            end 
        
        end 

        y[:,i] = tmp / ( 4 * nSourceHydrogens )
        # fout = "Inter$i.dat"
        # writedlm(fout,y)

    end 

    if nsample >= 5
        
        acf_out = zeros(Float64,maxlag,3)

        acf_out[:,1] = collect(0:(maxlag-1)) * deltat 

        acf_out[:,2] = mean(eachcol(y))

        acf_out[:,3] = std(eachcol(y)) / sqrt(nsample-1)

    else    

        acf_out = zeros(Float64,maxlag,2)

        acf_out[:,1] = collect(0:(maxlag-1)) * deltat

        acf_out[:,2] = mean(eachcol(y))

    end 

    return acf_out

end 
