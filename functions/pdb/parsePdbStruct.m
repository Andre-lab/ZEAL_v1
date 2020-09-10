function [atom_list, varargout] = parsePdbStruct(pdb_data, varargin)
% PARSE_PDB_STRUCT
% Create list with atom coordinates together
% with their positional variance (from B-factors) and vdW radii

% INPUT
% pdb_data  :   structured array from pdbread() (Matlab bioinformatics toolbox )

% OUTPUT
% atom_list :   Nx6 matrix where N are the number of heavy atoms in the
% structrure
% column 1 - X coordinate
% column 2 - Y coordinate
% column 3 - Z coordinate
% column 4 - positional variance (Å)
% column 5 - element code: 1=H; 2=C; 3=N; 4=O; 5=S
% column 6 - vdW radii (Å)

% References
% Bondi, A. (1964). "Van der Waals Volumes and Radii". J. Phys. Chem. 68 (3): 441–451. doi:10.1021/j100785a001


%%

Xcrds = pdb_data.X';
Ycrds = pdb_data.Y';
Zcrds = pdb_data.Z';

atom_position_variance = (3* [pdb_data.betaFactor]) /  (8.0 * pi^2);
atom_sernumber = pdb_data.atomNum;

XYZ = [Xcrds Ycrds Zcrds];

if nargin>1
    centerCOMop=varargin{1};
else
    centerCOMop=true;
end

if centerCOMop
    COM = mean(XYZ);
else
    COM = [0 0 0];
end

n_atoms = size(Xcrds,1);


elements = {'H','C','N','O','S','X'};
element_codes = [1 2 3 4 5 6];
vdw_radii = [1.2 1.7 1.55 1.52 1.8];
avg_radius = mean(vdw_radii(2:end));

count = 0;
Hcount = 0;

atom_list = zeros(n_atoms,7);
hydrogen_list = zeros(n_atoms,1);

for i = 1:n_atoms
    
    element_i = pdb_data.element{i};
    if isempty(element_i) % no element annotation, check resname instead to get element        
        element_i = pdb_data.atomName{i}(1);
        % Take of cases such as 1H..; 2H..
        if ~isletter(element_i)
            element_i = pdb_data.atomName{i}(2);
        end
    end
    
    % ONLY HEAVY ATOMS
    try
        if (element_i=='C') || (element_i=='N')  || (element_i=='O') || (element_i=='S') || (element_i=='X')
            
            count = count + 1;
            
            atom_list(count,1) = Xcrds(i)-COM(1);
            atom_list(count,2) = Ycrds(i)-COM(2);
            atom_list(count,3) = Zcrds(i)-COM(3);
            atom_list(count,4) = atom_position_variance(i);
            atom_list(count,7) = atom_sernumber(i);
            
            element_code = element_codes(strcmp(element_i, elements));
            
            switch element_code
              
                case 2 % C
                    atom_list(count,6) = vdw_radii(2);
                    
                case 3 % N
                    atom_list(count,6) = vdw_radii(3);
                case 4 % O
                    atom_list(count,6) = vdw_radii(4);
                case 5 % S
                    atom_list(count,6) = vdw_radii(5);
                case 6
                    atom_list(count,6) = avg_radius;
            end
            
            atom_list(count,5) = element_code;
            
        elseif (element_i == 'H') || (element_i == 'D')
            
            Hcount = Hcount + 1;
            hydrogen_list(Hcount) = i;             
        end
        
    catch
        warning('Could not parse atom name.');
    end
end

atom_list(count+1:end,:) = [];
hydrogen_list(Hcount+1:end,:) = [];

varargout(1) = {hydrogen_list};

end

