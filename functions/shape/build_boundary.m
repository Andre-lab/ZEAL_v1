function [SA_surf] = build_boundary(SA_solid)
%BUILD_BOUNDARY Find boundary (perimeter) of the SA soild. 
% This is the solvent-accessible surface

% Use a general technique that works for any dimensionality
% and any connectivity.

% INPUT
% SA_solid  :   NxNxN binary matrix containing the voxelized representation
                % of the solvent-accessible solid (see create_SA_solid.m)

% OUTPUT 
% SA_surf   :   NxNxN binary matrix containing the voxelized representation
                % of the solvent-accessible surface

%% 

% define connectivity matrix 
conn = conndef(3, 'minimal');
num_dims = max(3, ndims(conn));

% add padding to SA_solid 
b = padarray(SA_solid,ones(1,3),0,'both');

% erode the 3d image
b_eroded = imerode(b,conn);
p = b & ~b_eroded;
idx = cell(1,3);

for k = 1 : num_dims
    idx{k} = 2:(size(p,k) - 1);
end

SA_surf = p(idx{:});

end

