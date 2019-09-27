function [geo_mom_list, geo_moments] = computeGeoMoments(voxels, dim , xCOM, yCOM, zCOM, scale, order)
%% compute_geoMoments
% Computes geometric moments M_rst up to order n for a voxelized object
% within a grid. The object, with center of mass at (x|y|z)COG, is scaled to fit within 
% the unit ball before M_rst is computed for each combination of indices,
% such that r,s,t>0 and r+s+t<n. The algorithm is based on the C++ code by 
% Novotni % Klein 
% 
% % Novotni, M., & Klein, R. (2003).
% 3D Zernike Descriptors for Content Based Shape Retrieval.
% Proceedings of the Eighth ACM Symposium on Solid Modeling and Applications, 
% 216–225. https://doi.org/10.1145/781636.781639

% INPUT 
% ---------------------
% <voxels> : 1D vector
% The flatted voxelized object, where x, y and z coordinates are mapped to
% the array elements as ((z-1) * dim + y-1) * dim + x-1+1

% <dim> : integer
% The side-length resolution of the cubic grid containing the voxelized
% object. 

% <xCOG>, <yCOG>, <zCOG> : float
% The center of gravity of the object in the x, y and z dimension
% respectively

% <scale> : float 
% The scaling factor to fit object into unit ball 

% <order> : integer
% The maximum order for geometric moments. 

% OUTPUT
% ---------------------
% <geo_mom_list> : [n x 4] matrix, float 
% The geometric moments for each combination of indices r, s, t.
% Column 1 = index r
% Column 2 = index s
% Column 3 = index t
% Column 4 = geo metric moment

% <geo_moments> : [(order+1) x (order+1) x (order +1)] matrix
% Each element contains the geometric moments for a given combinatino of
% indices r,s,t, i.g. M_{r=2,s=2,t=2} = geo_moments(r=2,s=2,t=2)

% Filip 2018-11-15

%% Scale object to fit into unit ball 
xDim = dim;
yDim = dim;
zDim = dim;

dimVec = [xDim yDim zDim];

min_v = zeros(3,1);

min_v(1) = (-xCOM+1)*scale;
min_v(2) = (-yCOM+1)*scale;
min_v(3) = (-zCOM+1)*scale;

samples = zeros(max(dimVec)+1,3);

for r=1:3
    for s=0:dimVec(r)
        samples(s+1,r)=min_v(r) + s*scale;
    end
end

%% Compute geometric moments

arrayDim = zDim;
layerDim = yDim*zDim;

diffArrayDim = zDim+1;
diffLayerDim = (yDim+1)*zDim;
diffGridDim = (xDim+1)*layerDim;

diffGrid = zeros(diffGridDim,1);
diffLayer = zeros(diffLayerDim,1);
diffArray = zeros(diffArrayDim,1);

layer = zeros(layerDim,1);
array = zeros(arrayDim,1);

geo_moments = zeros(order+1, order+1, order+1);
geo_mom_list = zeros((order+1)^3,4);

% generate the diff version of the voxel grid in x direction
iter = 1;
diffIter = 1;

for x = 0:layerDim-1
% for x = 0:0
    
    diffGrid = ComputeDiffFunction(voxels, diffGrid, iter, diffIter, xDim);
    
    iter = iter + xDim;
    diffIter = diffIter + xDim + 1;
    
end

count = 0;

for r=0:order
    
    % diffGrid
    diffIter = 1;
        
    for p=0:layerDim-1
           
        sampleIter = 1;
        [layer(p+1), diffGrid] = Multiply(diffGrid, samples(:,1), diffIter, sampleIter, xDim+1);
        
        diffIter = diffIter + xDim + 1;
        
    end
    
    % layer
    iter = 1;
    % diffLayer
    diffIter = 1;
    
    for y=0:arrayDim-1
        
        diffLayer = ComputeDiffFunction(layer, diffLayer, iter, diffIter, yDim);
        
        iter = iter + yDim;
        diffIter = diffIter + yDim + 1;
        
    end
    
    for s=0:order-r
        
        % diffLayer
        diffIter = 1;
        
        for p=0:arrayDim-1
            
            sampleIter = 1;
            [array(p+1), diffLayer] = Multiply(diffLayer, samples(:,2), diffIter, sampleIter, yDim+1);
            diffIter = diffIter + yDim + 1;
        end
        
        % array
        iter = 1;
        % diffarray
        diffIter = 1;
        
        diffArray = ComputeDiffFunction(array, diffArray, iter, diffIter, zDim);
        
        for t=0:order-r-s
            
            count = count+1;
            sampleIter = 1;
            [moment, diffArray] = Multiply(diffArray, samples(:,3), diffIter, sampleIter, zDim+1);
            geo_moments(r+1,s+1,t+1) = moment / ( (1+r)*(1+s)*(1+t) );
            geo_mom_list(count,:) = [r s t geo_moments(r+1,s+1,t+1)];
            
        end
        
    end % j
end %i

geo_mom_list(count+1:end,:) = [];

end


function [B] = ComputeDiffFunction(A, B, Aiter, Biter, dim)

B(Biter) = -1*A(Aiter);

for i = 1:dim-1
    
    B(Biter+i) = A(Aiter+i-1)-A(Aiter+i);
    
end

B(Biter+dim) = A(Aiter+(dim-1));

end



function [sum_val, diffGrid] = Multiply(diffGrid, sample, diffIter, sampleIter, dim)

sum_val = 0;

for i=0:dim-1
    
    diffGrid(diffIter+i) = diffGrid(diffIter+i) * sample(sampleIter+i);
    sum_val = sum_val + diffGrid(diffIter+i);
    
end

end


% vectorizing is slower than the loop
% function [sum_val, diffGrid] = Multiply_vec(diffGrid, sample, diffIter, sampleIter, dim)
% 
% i_vec = 0:dim-1;
% ii_vec = diffIter+i_vec;
% diffGrid(ii_vec) = diffGrid(ii_vec) .* sample(sampleIter+i_vec);
% 
% sum_val = sum(diffGrid(ii_vec));
% 
% 
% 
% end
% 







