function [zernikeMoments_list, zernikeMoments] = computeZCmom(order, chi_coeff_cell, chi_nlm_rst_cell, geo_moments)
% COMPUTE_ZCMOM
% Computes the Zernike moments of order n from geometric moments and the 
% object-independent chi-coefficients of the same order. This computation is data dependent
% and has to be performed for each new object and/or transformation.

% INPUT

% OUTPUT

zernikeMoments = zeros(order+1, order+1, order+1);

% geoMoms_test = zeros(order^5, 9);
zernikeMoments_list = zeros(order^3,5);
momcount = 0;

pi_factor = (3/(4*pi));

for n=0:order
    
    l0 = mod(n,2);
    
    for l = l0:2:n
        
        for m=0:l
            
            zm = complex(0);
            
            % get chi coeffs  
            chi_nlm_rst = chi_nlm_rst_cell{n+1,l+1,m+1};
            chi_values = chi_coeff_cell{n+1,l+1,m+1};
            nCoeffs = size(chi_nlm_rst,1);
            
            for i = 1:nCoeffs
                
                r = chi_nlm_rst(i,4)+1;
                s = chi_nlm_rst(i,5)+1;
                t = chi_nlm_rst(i,6)+1;
                
                zm = zm + conj(chi_values(i)) * geo_moments(r,s,t);
                
            end
            
            zm = zm * pi_factor;
            
            if n == 0 && l == 0 && m == 0
                nullMoment = real(zm);
            end
            
            zernikeMoments(n+1,l+1,m+1) = zm;
            
            momcount = momcount + 1;
            
            zernikeMoments_list(momcount,1) = n;
            zernikeMoments_list(momcount,2) = l;
            zernikeMoments_list(momcount,3) = m;
            
            zernikeMoments_list(momcount,4) = real(zm);
            zernikeMoments_list(momcount,5) = imag(zm);
            
            if m > 0
                momcount = momcount + 1;
                
                zernikeMoments_list(momcount,1) = n;
                zernikeMoments_list(momcount,2) = l;
                zernikeMoments_list(momcount,3) = -m;
                
                zm_minus = (-1)^m * conj(zm);
                
                zernikeMoments_list(momcount,4) = real(zm_minus);
                zernikeMoments_list(momcount,5) = imag(zm_minus);
                
            end
           
        end %m
    end % l
end % n

zernikeMoments_list(momcount+1:end,:) = [];

end

