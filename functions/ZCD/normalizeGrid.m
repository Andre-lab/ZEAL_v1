function [voxels] = normalizeGrid(voxels, dim, xCOG, yCOG, zCOG, scale)
%NORMALIZE_GRID Cuts of the object function for values outside the unit
%ball

radius = 1/scale;
sqrRadius = radius^2;

for x = 1:dim
    for y = 1:dim
        for z = 1:dim
            
            
            id = ((z-1) * dim + y-1) * dim + x-1+1;
            
            if voxels(id) > 0
                
                mx = x - xCOG;
                my = y - yCOG;
                mz = z - zCOG;
                
                sqrLen = mx^2+my^2+mz^2;
                
                if sqrLen > sqrRadius
                    voxels(id) = 0;
                end
                
            end
            
        end
    end
end




end

