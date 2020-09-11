function [ZCinv] = compute_ZCinvariants(ZCmoments, order)
%COMPUTE_ZCINVARIANTS
% Computes the Zernike moment based invariants, i.e. the norms of vectors with
% components of Z_nl^m with m being the running index.

n_invariants = round((((5+2)^2)/4));

ZCinv = zeros(n_invariants,1);

inv_count = 0;
inv_count2 = 0;

for n = 0:order
    
    sum_tmp = 0;
    
    for l = mod(n,2):2:n  
        
        for m=-l:l
            
            inv_count2 = inv_count2 +1;
            
            absM = abs(m);
            
            % The ZC_nlm moment
            mom = ZCmoments(n+1, l+1, absM+1);
            
            
            %conjugate if m negative
            if m<0
                mom = conj(mom);
                % take care of sign for odd m
                if mod(absM,2)
                    mom = -1*mom;
                end
            end
            
            sum_tmp = sum_tmp + norm(mom)^2; 
            % the C++ std:norm function gives the square of the L2 
            %(euclidian norm), which is the so called field norm
                                
        end
        
        inv_count = inv_count + 1;      
        ZCinv(inv_count,1) =  sqrt(sum_tmp);
        
    end
end

