# using FFTW

"""
    acfNN(x,N)
    
    Compute the autocorrelation of data in vector x up to length N.
    The non-normalized autocorrelation is returned. 
    
"""
function acfNN(inputvector,maxlag)
    
    nd = length(inputvector)
    c  = zeros(Float64,nd)
    a  = vcat(inputvector,c) 

    denominator = nd .- collect(0:(nd-1))

    fvi = fft(a)
    acf = real( ifft(fvi .* conj(fvi)) )
    acf = acf[1:nd]
    acf = acf ./ denominator 

    return acf[1:maxlag]

end 
