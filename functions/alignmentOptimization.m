function [x,fval,exitflag,output,trials] = alignmentOptimization(ZCref, settings, atomlist_rot, options,lb,ub)
%OPTIM_FUN Summary of this function goes here
%   Detailed explanation goes here

[x,fval,exitflag,output,trials] = surrogateopt(@ZCcc_objfun,lb,ub,options);

    function [d] = ZCcc_objfun(x)
        %OBJ_FUN Summary of this function goes here       
        %R = euler2rotationMatrix(x(1), x(2), x(3), 'zyz');
        
        Ry = @(theta) [cos(theta) 0 -sin(theta); 0 1 0; sin(theta) 0 cos(theta)];
        Rz = @(theta) [cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1];
        
        R = Rz(x(3))*Ry(x(2))*Rz(x(1));
        
        atomlist_rot_new = atomlist_rot;
        atomlist_rot_new(:,1:3) = atomlist_rot_new(:,1:3) * R';
        
        rot_shape = createShape(atomlist_rot_new, settings);
                
        [~, ZC_mom_list, ~, ~] = getZCmomFromShape(rot_shape,  ZCref.order, ZCref.grid_res, ZCref.scale_option, ZCref.chi_coeff_cell, ZCref.chi_nlm_rst_cell);
        
        %d = -1*real(corr((ZC_mom_list), (ZCref.mom_list)));
        d = -1*real(dot(ZC_mom_list, ZCref.mom_list)) / (norm(ZC_mom_list)*norm(ZCref.mom_list));
        
    end

    function [shape] = createShape(atom_list, settings)
        %CREATESHAPE Summary of this function goes here
        %   Detailed explanation goes here
        switch settings.shape_op
            
            case 1 % solvent/molecular surface
                
                [shape, ~, ~, ~] = EDTms(atom_list, settings.grid_res, settings.shape);
                
            case 2 % density
                
                [shape, ~, ~] = create_elecDensity(atom_list, settings.grid_res, settings.shape.smear_factor, [], [], []);
                
                % normalize
                shape = shape/ max(shape(:));
        end
        
    end

end

