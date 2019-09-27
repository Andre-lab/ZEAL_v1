function [chi_coeff, chi_coeff_cell, chi_nlm_rst_cell] = compute_chiCoeffs(order, UIhandle)

% Initialize waitbar in GUI
d = uiprogressdlg(UIhandle,'Title','Please Wait',...
        'Message','Computing ZC chi-coefficients ');
    
count_coeff = 0;
setStart = 1;
set_count = 0;

chi_coeff = zeros((order+1)^4,8);

chi_coeff_cell = cell(order+1, order+1, order+1);
chi_nlm_rst_cell = cell(order+1, order+1, order+1);

chi_ind_map = zeros((order+1)^3,6);

fprintf('\n Computing chi coefficients for order = %d', order);

n_count_vec = zeros(order+1,1);

for n = 0:order
    
    fprintf('\n Doing order %d', n);
    
    li=0;
    l0 = mod(n,2);
    
    n_count = 0;
       
    for l = l0:2:n % only even values of l since Zernike functions require that (n-l) is even
        li=li+1;
        
        for m = 0:l
            
            c_set_count=0;
            
            cs = c_lm(l,m) ;
            
            w = cs/ 2^(m);
            
            k=round( (n-l)/2 );
            
            for nu = 0:k
                
                qs = q_klnu(k,l,nu);
                w_Nu = w * qs;
                
                for alpha = 0:nu
                    
                    w_NuA = w_Nu * nchoosek(nu,alpha);
                    
                    for beta = 0:(nu-alpha)
                        
                        w_NuAB = w_NuA * nchoosek(nu-alpha, beta);
                        
                        for p = 0:m
                            
                            w_NuABP = w_NuAB * nchoosek(m,p);
                            
                            for mu = 0:floor((l-m)/2)
                                
                                w_NuABPMu = w_NuABP ...
                                    * nchoosek(l, mu) ...
                                    * nchoosek(l-mu, m+mu) ...
                                    / 2^(2 * mu);
                                
                                for q = 0:mu
                                    
                                    w_NuABPMuQ = w_NuABPMu * nchoosek(mu, q);
                                    
                                    % the sign
                                    if mod((m-p+mu),2)
                                        w_NuABPMuQ = -1 * w_NuABPMuQ;
                                    end
                                    
                                    rest = mod(p,4);
                                    
                                    switch rest
                                        case 0
                                            c = complex(w_NuABPMuQ, 0);
                                        case 1
                                            c = complex(0, w_NuABPMuQ);
                                        case 2
                                            c = complex(-1 * w_NuABPMuQ, 0);
                                        case 3
                                            c = complex(0, -1 * w_NuABPMuQ);
                                    end
                                    
                                    t_i = l - m + 2 * (nu - alpha - beta - mu);
                                    s_i = 2 * (mu - q + beta) + m - p;
                                    r_i = 2 * q + p + 2 * alpha;
                                    
                                    c_set_count = c_set_count +1;
                                    
                                    count_coeff = count_coeff + 1;
                                    n_count = n_count + 1;
                                    
                                    chi_coeff(count_coeff,1) = n;
                                    chi_coeff(count_coeff,2) = l;
                                    chi_coeff(count_coeff,3) = m;                                    
                                    
                                    chi_coeff(count_coeff,4) = r_i;
                                    chi_coeff(count_coeff,5) = s_i;
                                    chi_coeff(count_coeff,6) = t_i;
                                    
                                    chi_coeff(count_coeff,7) = real(c);
                                    chi_coeff(count_coeff,8) = imag(c);
                                                                        
                                end %q
                            end % mu
                        end % p
                    end % beta
                end % alpha
            end % nu
            
            if c_set_count > 0
                
                set_count = set_count + 1;
                
                chi_ind_map(set_count,1) = n;
                chi_ind_map(set_count,2) = l;
                chi_ind_map(set_count,3) = m;
                
                chi_ind_map(set_count, 4) = setStart;
                chi_ind_map(set_count, 5) = setStart + c_set_count - 1;
                
                setStart = chi_ind_map(set_count, 5) + 1;
            end
            
            sel_int = chi_ind_map(set_count, 4): chi_ind_map(set_count, 5);
            
            chi_coeff_cell{n+1,l+1,m+1} = complex( chi_coeff(sel_int, 7), chi_coeff(sel_int, 8) );
            chi_nlm_rst_cell{n+1,l+1,m+1} = chi_coeff(sel_int, 1:6);
            
        end % m
    end % l
    
    n_count_vec(n+1,1) = n_count;
    
    % update UI waitbar
    d.Value = (n+1)/(order+1);
    
end % n

chi_coeff(count_coeff+1:end,:) = [];
chi_ind_map(set_count+1:end,:) = [];

chi_ind_map(:, 6) = chi_ind_map(:, 5) - chi_ind_map(:, 4)+1;

close(d);

end

%%

function c_lm_val = c_lm(l,m)
% c_l^m = c_l^-m

c_lm_val = sqrt( (2 * l + 1) * factorial(l + m) * factorial(l - m) ) / factorial(l);


end


function q_klnu_val = q_klnu(k,l,nu)

% nominator of straight part
nom = nchoosek(2*k,k) * (nchoosek(k, nu)) * (nchoosek(2 * (k + l + nu) + 1, 2 * k));

if mod(k+nu,2)
    nom = -1*nom;
end

% denominator of straight part
den = 2^(2*k) * nchoosek(k+l+nu,k);

% nominator of sqrt part
n_sqrt = 2*l + 4*k + 3;

% denominator of sqrt part
d_sqrt = 3;

q_klnu_val =  nom / den * sqrt(n_sqrt / d_sqrt);

end
