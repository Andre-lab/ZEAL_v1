function [elecDensityVolume, scaling_factor, voxel_res] = create_elecDensity(atom_list, grid_res, smear_factor, COM, scaling_factor, voxel_res)
%CREATE_ELECDENSITY Summary of this function goes here
%   Detailed explanation goes here


%smear_factor: fraction of the grid to splat atoms over (atoms are represented as Gaussians)
    
    
%% Step 1
padding = 0.6*smear_factor;

n_atoms = size(atom_list,1);

% center of mass
if isempty(COM)   
    COM = mean(atom_list(:,1:3));
end

% translate with COM at origo
XYZ_com = atom_list(:,1:3)-  (ones(n_atoms,1) * COM);


smear_range = round(smear_factor * grid_res/2);

if isempty(scaling_factor)
   % fprintf('\n Using fefault scaling');
    % find longest edge of box containing all points
    Rmax = max((max(XYZ_com)-min(XYZ_com))/2) + max(atom_list(:,6));
    
    % resolution of voxelgrid in Å
    voxel_res = Rmax / (grid_res/2) / (1-padding);
    scaling_factor = (1-padding) / (Rmax );
%else
    %fprintf('\n User-def scaling %2.2f', scaling_factor);
end

% scaled atoms and positional variance (B-factors)
sr_i = ceil(atom_list(:,6) / voxel_res); % scaled vdw
sr_B = ceil(atom_list(:,4) / voxel_res); % scaled positional var

% scale and translate atom coordinates with center of mass placed at origo
XYZ_comScaled =  XYZ_com * scaling_factor;

% pre-allocate memory for 3D grid
elecDensityVolume = zeros(grid_res, grid_res, grid_res);

%% Step 2
keep_on = true;
while keep_on
    breakflag=false;
    %     try
    for i = 1:n_atoms        
             
        sigma = sqrt(sr_B(i)) + sr_i(i); % scaled [position variance from B-factor] + scaled [vdw radius]
       
        xpos = round(0.5 * (XYZ_comScaled(i,1) + 1) * grid_res);
        ypos = round(0.5 * (XYZ_comScaled(i,2) + 1) * grid_res);
        zpos = round(0.5 * (XYZ_comScaled(i,3) + 1) * grid_res);
       
        for x_smear = -smear_range:smear_range
            for y_smear = -smear_range:smear_range
                for z_smear = -smear_range:smear_range
                    
                    d = (x_smear)^2 + (y_smear)^2 + (z_smear)^2;
                    
                    if d <= smear_range
                        if (xpos+x_smear)>0 && (xpos+x_smear)<grid_res && ...
                                (ypos+y_smear)>0 && (ypos+y_smear)<grid_res && ...
                                (zpos+z_smear)>0 && (zpos+z_smear)<grid_res
                            
                            val = atom_list(i,4) * exp(-d^2 / (2.0 * (sigma)^2));
                            
                            val = val / (sigma * 2 * pi)^0.2;
                            
                            elecDensityVolume(xpos+x_smear, ypos+y_smear, zpos+z_smear) = elecDensityVolume(xpos+x_smear, ypos+y_smear, zpos+z_smear) + val;
                            
                        else                            
                            breakflag=true;
                        end
                        
                    end
                    
                    if breakflag; break; end                                        
                end
                if breakflag; break; end
            end
            if breakflag; break; end
        end
        if breakflag; break; end
        
    end
    
    if breakflag
        
        keep_on = true;
        
%         fprintf('\n Could not fit SA solid on grid - rescaling again');
       
        scaling_factor = 0.95* scaling_factor;
        
        % rescale
        XYZ_comScaled =  XYZ_com * scaling_factor;
        
    else
        keep_on = false;
    end 
    
end



end

