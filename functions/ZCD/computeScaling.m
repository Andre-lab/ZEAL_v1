function [scaling_factor, max_extent] = computeScaling(scale_option, voxels, dim,  xCOG, yCOG, zCOG)
% compute_scaling
% To compute Zernike-Canterakis moments for a 3d object, the object has to 
% be scaled so that it can fit inside the unit-ball (where the ZC functions live). 

% This function gives that scaling factor by two different methods. Reproduction accuracy  
% falls off close to the ball surface (due to discretization effects?) and both methods make sure 
% the scaling avoids this. Novotni and Klein used method 2 as per their code (not described in paper)
% without any justification. Method 1 is used by Grandison, Roberts and Morris as per the paper,
% but no justification for why the factor of 0.6 was chosen. 


% Novotni, M., & Klein, R. (2003).
% 3D Zernike Descriptors for Content Based Shape Retrieval.
% Proceedings of the Eighth ACM Symposium on Solid Modeling and Applications, 
% 216–225. https://doi.org/10.1145/781636.781639


% Grandison, S., Roberts, C., & Morris, R. J. (2009). 
% The Application of 3D Zernike Moments for the Description of 
% “Model-Free” Molecular Structure, Functional Motion, and Structural Reliability. 
% Journal of Computational Biology, 16(3), 487–500. https://doi.org/10.1089/cmb.2008.0083


% INPUT
% <scale_option> : integer
    % 1 
    % scaling so that the maximum distance between a 
    % filled voxel and the grid centroid is 70 % of the unit ball radius. 

    % 2 (Used by Novotni & Klein)
     % scaling by the inverse of 2 times the average distance between filled
    % voxels and the grid centroid. 

    % 3 
    % scaling so that the maximum distance between a 
    % filled voxel and the grid centroid is 100 % of the unit ball radius. 

% <voxels> : 1D vector
% The flatted voxelized object, where x, y and z coordinates are mapped to
% the array elements as ((z-1) * dim + y-1) * dim + x-1+1

% <dim> : integer
% The side-length resolution of the cubic grid containing the voxelized
% object. 

% <xCOG>, <yCOG>, <zCOG> : float
% The center of gravity of the object in the x, y and z dimension
% respectively
                        
% OUTPUT
% <scaling_factor> : float
% The scaling factor for fitting object voxels into the unit ball. 
% <Rmax> : float
% Maximum distance from grid centroid
% <max_extent> : float
% For scaling option 1: The fraction of the unit ball radius where the largest distance from the
% centroid is 
% <voxel_resolution> : float
% The resolution (in Ångström) captured by a voxel

% Filip 2018-11-15

%% 

switch scale_option
    
    case 1
        Rmax = 0;
        for x = 1:dim
            for y = 1:dim
                for z = 1:dim
                    
                    
                    id = ((z-1) * dim + y-1) * dim + x-1+1;
                    
                    if voxels(id) > 0
                        
                        mx = x - xCOG;
                        my = y - yCOG;
                        mz = z - zCOG;
                        
                        temp = mx^2+my^2+mz^2;
                        
                        if temp > Rmax
                            Rmax = temp;
                        end
                        
                    end
                    
                end
            end
        end
        
        max_extent = 0.7;
        
        scaling_factor = max_extent*(1 / sqrt(Rmax));
%         Rmax = sqrt(Rmax);
        
        
    case 2
        temp=0;
        nvox = 0;
        for x = 1:dim
            for y = 1:dim
                for z = 1:dim
                    
                    id = ((z-1)* dim + y-1) * dim + x-1+1;
                    
                    if  voxels(id) > 0
                        
                        mx = x - xCOG;
                        my = y - yCOG;
                        mz = z - zCOG;
                        
                        temp = temp + mx^2+my^2+mz^2;
                        nvox=nvox+1;
                        
                    end
                end
            end
        end
                
        Rmean = sqrt(temp/nvox);       

        max_extent = 0.5;
        scaling_factor = 1/(2* Rmean);
        
    case 3
        
        Rmax = 0;
        for x = 1:dim
            for y = 1:dim
                for z = 1:dim
                    
                    
                    id = ((z-1) * dim + y-1) * dim + x-1+1;
                    
                    if voxels(id) > 0
                        
                        mx = x - xCOG;
                        my = y - yCOG;
                        mz = z - zCOG;
                        
                        temp = mx^2+my^2+mz^2;
                        
                        if temp > Rmax
                            Rmax = temp;
                        end
                        
                    end
                    
                end
            end
        end
        
        max_extent = 1;        
        scaling_factor = max_extent*(1 / sqrt(max(Rmax)));             

end

% voxel_resolution = Rmax / ( (dim/2));

end

