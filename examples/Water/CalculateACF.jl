using NMRRelaxation
using StatsBase

# Load the molecule 
ta, boxsize = trajinfo("w256.psf","w256_pr.dcd");
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
base_h_indices = xx["name H1 or name H2"].atomid;
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
# Let us sample 5 water molecules from the 256 residues 
#
MaxResidue=256
residues = 1:MaxResidue
nsample = 1
sample_residues = sample(residues, nsample; replace=false)
#
# End of system specific stuff 
#

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
