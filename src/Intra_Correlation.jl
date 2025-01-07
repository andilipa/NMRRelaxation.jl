# using MDTools
# using StatsBase
# using Combinatorics
# using Statistics
# include("/Users/dan/Julia/pair_NMR.jl")

function Intra(nsample::Int64, maxlag::Int64, sample_residues::Vector{Int64}, deltat::Float64, 
    boxsize::Array{Float64,2}, nframes::Int64, sampling_vector::Vector{Vector{Float64}},  
    atoms_per_residue::Int64, base_h_indices::Vector{Int64}, ta)
    
    y = zeros(Float64,maxlag,nsample)

    for i in 1:nsample

        hydrogens = (sample_residues[i] - 1) * atoms_per_residue .+ base_h_indices;

        pairs = collect(combinations(hydrogens,2))

        number_of_hh_pairs = length(pairs)

        tmp = zeros(maxlag)

        for k in 1:number_of_hh_pairs

            h1 = pairs[k][1]   # Remember that pairs is a vector of vectors
            h2 = pairs[k][2] 

            x1 = ta[:,h1].xyz;
            x2 = ta[:,h2].xyz;

            acf = [zeros(maxlag) for _ in 1:4]

            Threads.@threads for j in 1:4
        
                sampleV = sampling_vector[j]

                acf[j] = G(nframes,maxlag,x1,x2,boxsize,sampleV)

            end     

            for j in 1:4
                tmp += acf[j]
            end     

        end 
    
        # Number of pairs is N * (N-1)/2. 
        # We need to normalize by number of hydrogens in the molecule, not number of pairs
        # Since we compute the autocorrelation over pairs, we need to 
        # correct by (N-1) =  (length(base_h_indices)-1)
        y[:,i] = tmp / ( 4 * number_of_hh_pairs) * (length(base_h_indices)-1)
    
        # fout = "Intra$i.dat"
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
