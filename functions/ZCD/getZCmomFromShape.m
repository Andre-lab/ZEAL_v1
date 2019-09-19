function [ZCmom_list, ZCmoms, ZCmom_cell, COM] = getZCmomFromShape(voxed_shape, order, grid_res, scale_option, chi_coeff_cell, chi_nlm_rst_cell)
%ZC_FUN Summary of this function goes here
%   Detailed explanation goes here


% flatten 3d grid to 1d grid
[flat_voxed_mesh, ~] = flattenGrid(voxed_shape, grid_res);

% Compute first order geometric moments
[~, geo_mom_COG] = computeGeoMoments(flat_voxed_mesh, grid_res , 0, 0, 0, 1, 1);

% 0'th order moments -> normalization
null_moment = geo_mom_COG(1,1,1);

% 1'st order moments -> center of gravity
xCOM = geo_mom_COG(2,1,1) / null_moment;
yCOM = geo_mom_COG(1,2,1) / null_moment;
zCOM = geo_mom_COG(1,1,2) / null_moment;

COM = [xCOM yCOM zCOM];

% Normalize to fit inside unit ball
[scale, ~] = computeScaling(scale_option, flat_voxed_mesh, grid_res, xCOM, yCOM, zCOM);

% cuts off object function for values outide unit ball
flat_voxed_mesh_norm = normalizeGrid(flat_voxed_mesh, grid_res , xCOM, yCOM, zCOM, scale);

% Compute geometric moments
% fprintf('\n Computing geometric moments');
[~, geo_mom] = computeGeoMoments(flat_voxed_mesh_norm, grid_res , xCOM, yCOM, zCOM, scale, order);

% Compute ZC moments
% fprintf('\n Computing ZC moments');
[ZCmom_list, ZCmom_cell] = computeZCmom(order, chi_coeff_cell, chi_nlm_rst_cell,  geo_mom);

ZCmoms = complex(ZCmom_list(:,4), ZCmom_list(:,5));

end

