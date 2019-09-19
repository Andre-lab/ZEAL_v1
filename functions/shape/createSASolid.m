function [SA_solid, voxel_res, scaling_factor ] = createSASolid(atom_list, grid_res, probe_radius, COM, scaling_factor, voxel_res)
% CREATE_SA_SOLID
% Create a binary volumetric representation of the solvent-accessible (SA)
% solid. Uses the coordinate extrema to find the minimal bounding cube
% around the object

% Step 1
% Scale and translate all atoms to fit inside a bounding box with side
% length L and with integer coordinates (voxels). The resolution of this
% grid is thus L^3 voxels. After this scaling, the vdW radii r_i and
% solvent probe radius r_p becomes sr_i and sr_p voxels.

% Step 2
% Create a solvent-accessible (SA) solid by assigning voxels a value of 1
% if they are within sr_i+sr_p for each atom, or 0 otherwise.

% INPUT
% atom_list     :   Nx6 matrix generated from parse_pdb_struct
% (see parse_pdb_struct.m for details)

% grid_res      :   Side length of the cubic grid (in grid intervals)

% probe_radius  :   Radius of the solvent-probe in Å

% padding       :   padding to add as the fraction of L (0-1)


% OUTPUT
% SA_solid      :   NxNxN binary matrix where N = grid_res
% Filled voxels = 1; Empty voxels = 0

%% Step 1
padding = 0.15;
n_atoms = size(atom_list,1);

% center of mass
% center of mass
if isempty(COM)
    COM = mean(atom_list(:,1:3));
    %fprintf('\n Using default COM');
    %else
    %   fprintf('\n User-def COM: %2.2f %2.2f %2.2f', COM(1), COM(2), COM(3));
end

% translate with COM at origo
XYZ_com = atom_list(:,1:3) - (ones(n_atoms,1) * COM);

if isempty(scaling_factor)
    % fprintf('\n Using default scaling');
    % find longest edge of box containing all points
    Rmax = max( (max(XYZ_com)-min(XYZ_com))/2 ) + probe_radius + max(atom_list(:,6));
    
    % resolution of voxelgrid in Å
    voxel_res = Rmax / (grid_res/2) / (1-padding);
    
    scaling_factor = (1-padding) / (Rmax );
    %else
    %  fprintf('\n User-def scaling: %3.3f', scaling_factor);
    
    
end

% scaled probe radius
sp_r = ceil(probe_radius / voxel_res);
sr_i = ceil(atom_list(:,6) / voxel_res);

% scale and translate atom coordinates with center of mass placed at origo
XYZ_comScaled =  XYZ_com * scaling_factor;

% pre-allocate memory for 3D grid
SA_solid = false(grid_res, grid_res, grid_res);

%% Step 2
keep_on = true;

while keep_on
    
    try
         atom_crd_grid = round((XYZ_comScaled+1)*(grid_res/2));
         
        for i = 1:n_atoms
            
            xpos = atom_crd_grid(i,1);
            ypos = atom_crd_grid(i,2);
            zpos = atom_crd_grid(i,3);
            
            Rip_vox = (sp_r+sr_i(i));
            Rip_vox2 = (Rip_vox^2);
            
            for x_fill = -Rip_vox:Rip_vox
                for y_fill = -Rip_vox:Rip_vox
                    for z_fill = -Rip_vox:Rip_vox
                        
                        d = (x_fill)^2 + (y_fill)^2 + (z_fill)^2;
                        
                        if d <= Rip_vox2
                            if (xpos+x_fill)>0 && (xpos+x_fill)<grid_res && ...
                                    (ypos+y_fill)>0 && (ypos+y_fill)<grid_res && ...
                                    (zpos+z_fill)>0 && (zpos+z_fill)<grid_res
                                
                                SA_solid(xpos+x_fill, ypos+y_fill, zpos+z_fill) = true;
                                
                            else
                                error('grid boundary crossing')
                            end
                            
                        end
                    end
                end
            end            
            
        end
        
        keep_on = false;
        
    catch
        
        keep_on = true;
        
        %         fprintf('\n Could not fit SA solid on grid - rescaling again');
        
        scaling_factor = 0.95 * scaling_factor;
        
        % rescale
        XYZ_comScaled =  XYZ_com * scaling_factor;
        
        
        
    end
    
    
end
