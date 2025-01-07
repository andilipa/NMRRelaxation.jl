using Printf

function printfile(fout,acf)
    io = open(fout,"w");
    (nrow,ncol) = size(acf)
    if ncol == 3
        # We have error information
        for i in 1:nrow
            # Break if the uncertainty exceeds the signal
            if acf[i,2] < acf[i,3]
                break
            else    
                @printf(io,"%.2f  %.7f  %.7f\n",acf[i,1],acf[i,2],acf[i,3])
            end 
        end
    else    
        for i in 1:nrow
            @printf(io,"%.2f  %.7f\n",acf[i,1],acf[i,2])
        end 
    end  
end 
