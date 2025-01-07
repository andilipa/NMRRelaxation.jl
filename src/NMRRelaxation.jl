module NMRRelaxation

using FFTW, LinearAlgebra, StatsBase, Combinatorics, Statistics, MDTools 

# greet() = print("Hello World!")

include("./TrajInfo.jl")
export trajinfo

include("./acfNN.jl")
export acfNN

include("./pair_NMR.jl") 
export G 

include("./printtofile.jl")
export printfile 

include("./Intra_Correlation.jl")
export Intra 

include("./Inter_Correlation.jl")
export Inter

end # module NMRRelaxation
