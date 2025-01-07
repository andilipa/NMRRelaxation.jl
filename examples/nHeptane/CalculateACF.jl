using NMRRelaxation
using StatsBase

# Load the molecule 
ta, boxsize = trajinfo("c7box.psf","c7box_pr.dcd");

# 
# System specific stuff 
#

# time in picoseconds; frames were saved every 100 fs = 0.1 ps
deltat = 0.1 
nframes = ta.nframe
maxlag  = div(nframes,2) 
#
# Need to create the base h-indices; so pick resid 1 and work with that. 
#
xx = ta["resid 1"];
base_h_indices = []
push!(base_h_indices, xx["name H11 or name H12 or name H13"].atomid)
push!(base_h_indices, xx["name H21 or name H22"].atomid)
push!(base_h_indices, xx["name H31 or name H32"].atomid)
push!(base_h_indices, xx["name H41 or name H42"].atomid)
push!(base_h_indices, xx["name H51 or name H52"].atomid)
push!(base_h_indices, xx["name H61 or name H62"].atomid)
push!(base_h_indices, xx["name H71 or name H72 or name H73"].atomid)
base_h_indices = reduce(vcat,base_h_indices);
nhydrogens = length(base_h_indices)
println("Number of protons is $nhydrogens\n")
atoms_per_residue = length(xx.atomid); 
#
# We will sample the four directions defined by a tetrahedron 
#
sampling_vector = [zeros(3) for _ in 1:4]
sampling_vector[1] = [0.0, 0.0, 1.0];
sampling_vector[2] = [√(8/9), 0, -1/3];
sampling_vector[3] = [-√(2/9), √(2/3), -1/3];
sampling_vector[4] = [-√(2/9), -√(2/3), -1/3];
#
# Let us sample 5 heptane molecules from the 261 residues 
#
MaxResidue=261
residues = 1:MaxResidue
nsample = 5
sample_residues = sample(residues, nsample; replace=false)
#
# Intramolecular correlation
#
acf = Intra(nsample, maxlag, sample_residues, deltat, boxsize, nframes, sampling_vector, 
            atoms_per_residue, base_h_indices, ta)
#
# Output the data
#
fout = "Intra.dat"
printfile(fout,acf) 

#
# Intermolecular correlation
#
acf = Inter(nsample, maxlag, sample_residues, deltat, boxsize, nframes, 
            sampling_vector, MaxResidue, atoms_per_residue, base_h_indices,ta)
#
# Output the data 
# 
fout = "Inter.dat"
printfile(fout,acf) 
