function [flat_voxed_mesh, voxed_mesh_list] = flattenGrid(meshgrid_3D, grid_res)
% GRID2BINVOX 
% Flatten 3D (cubic) voxel-grid to 1D vector.
% xyz coordinates are mapped to the vector as 
% element_id = dim*((z-1)*dim+(x-1)) + (y-1)+1;
% where dim is the side length in grid intervals of the cubic grid.

% INPUT
% meshgrid_3D : NxNxN matrix containing the voxelized object
                % N = grid_res
                
% grid_res    : Side length of the cubic grid (in grid intervals)

% OUTPUT
% flat_voxed_mesh       : Binary vector 

% voxed_mesh_list       : grid coordinates for filled voxels

%%

n_verts = grid_res^3;
voxed_mesh_list = zeros(n_verts,3);
v_count = 0;

flat_voxed_mesh = zeros(grid_res^3,1);

for x = 1:grid_res
    for y = 1:grid_res
        for z = 1:grid_res
            
            if meshgrid_3D(x,y,z)>0
                v_count = v_count + 1;
                
                voxed_mesh_list(v_count,1) = x;
                voxed_mesh_list(v_count,2) = y;
                voxed_mesh_list(v_count,3) = z;
                
%                 id = grid_res*((z-1)*grid_res+(y-1)) + (x-1)+1;
                id = grid_res*((z-1)*grid_res+(x-1)) + (y-1)+1;

                flat_voxed_mesh(id) = meshgrid_3D(x,y,z);
                
            end
            
        end
    end
end

voxed_mesh_list(v_count+1:end,:) = [];

end
