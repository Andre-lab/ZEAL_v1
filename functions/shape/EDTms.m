function [mol_surf_thick, mol_surf_solid, mol_surf_mesh, scaling_factor] = EDTms(atom_list, grid_res, shape_settings, varargin)

if nargin > 3
    if nargin < 5
        smoothing = varargin{1}; % struct
    else
        smoothing = varargin{1}; % struct
        scaling = varargin{2}; % struct
    end
else
   
    smoothing.active = false;
    
    scaling.COM = [];
    scaling.factor = [];
    scaling.voxel_res = [];
end

%EDTMS Description

[sa_solid, vox_res, scaling_factor] = createSASolid(atom_list, grid_res, shape_settings.probe_radius, scaling.COM, scaling.factor, scaling.voxel_res);

% Find boundary -> SAS
sa_surf = build_boundary(sa_solid);

if strcmp(shape_settings.surf_type, 'MS')
    % Find interior
    M = logical(imfill(sa_surf,'holes'));
    
    % Do Euclidean distance transform (EDT) -> Euclidean distance map (EDM)
    EDM = round(bwdist(sa_surf,'euclidean'));
    
    % Change sign for interior voxels -> signed EDM (sEDM)
    sEDM=EDM;
    sEDM(M)=-sEDM(M);
    
    % scaled probe radius in voxel units
    sr_p = ceil(shape_settings.probe_radius/vox_res);
    
    % Create isosurface = molecular surface
    mol_surf = (sEDM == (-sr_p) );
    mol_surf_solid = (sEDM <= (-sr_p) );
    
    % Create thicker shells
    mol_surf_thick = mol_surf;
    
    if shape_settings.shell_thickness > 1
        for s = 1:shape_settings.shell_thickness-1
            mol_surf_thick = mol_surf_thick + (sEDM == (-sr_p - s) );
        end
    end
       
    
elseif strcmp(shape_settings.surf_type,'SA')
    
    % dilate surface (i.e. thicken the surface)    
    se = strel('cube', shape_settings.shell_thickness);
        
    mol_surf_thick = imdilate(sa_surf,se);
    mol_surf_solid = sa_solid;
else
    fprintf('\n Unrecognized surface type. \nChoose  between "SA" for solvent-accessible surface or \n "MS" for molecular (solvent-excluded surface)');
    mol_surf_thick = [];
    mol_surf_solid = [];
    return
end

if smoothing.active
    
    switch smoothing.filter_type
        case 1
            mol_surf_solid = smooth3(mol_surf_solid,'gaussian',smoothing.size,smoothing.var);
        case 2
            mol_surf_solid = smooth3(mol_surf_solid,'box',smoothing.size);
    end
end

if shape_settings.do_mesh
    %fprintf('\n Creating mesh representation of 3D image');
    gv = 1:grid_res;
    [xdata, ydata, zdata] = meshgrid(gv, gv, gv);
    mol_surf_mesh = isosurface(ydata,xdata,zdata,round(mol_surf_solid),0.01);
else
    mol_surf_mesh = [];
end

end

