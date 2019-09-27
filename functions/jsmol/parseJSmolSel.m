function [pdb_data_sel, atom_indices] = parseJSmolSel(JSmolSelFile, pdb_data)
%PARSEJSMOLSEL Summary of this function goes here
%   Detailed explanation goes here

JSmolStr = readJSmolSel(JSmolSelFile);
roiSel = JSmolStr{1};
roiSel(1:2)=[]; roiSel(end-1:end)=[];

atomIndIntCell = split(roiSel,' ');

nCells = length(atomIndIntCell);

atom_indices = [];

for n = 1:nCells
    
    n_range = split(atomIndIntCell{n},':');
   
    if length(n_range)>1
       
        n_ind1 = str2double( n_range{1});
        n_ind2 = str2double( n_range{2});
        
        n_atomInd_interval = n_ind1:n_ind2;        
        atom_indices = [atom_indices n_atomInd_interval];
        
    else
        
        atom_indices = [atom_indices str2double(n_range{1})];
        
    end
    
    
end

atom_indices = atom_indices + 1; % JSmol indices start from 0 (Matlab starts from 1)

% get pdb data from selection
pdb_data_sel = pdb_data;

pdb_data_sel.recordName = pdb_data.recordName(atom_indices);
pdb_data_sel.atomNum = pdb_data.atomNum(atom_indices);
pdb_data_sel.atomName = pdb_data.atomName(atom_indices);

pdb_data_sel.altLoc = pdb_data.altLoc(atom_indices);
pdb_data_sel.resName = pdb_data.resName(atom_indices);
pdb_data_sel.chainID = pdb_data.chainID(atom_indices);
pdb_data_sel.resNum = pdb_data.resNum(atom_indices);

pdb_data_sel.X = pdb_data.X(atom_indices);
pdb_data_sel.Y = pdb_data.Y(atom_indices);
pdb_data_sel.Z = pdb_data.Z(atom_indices);

pdb_data_sel.occupancy = pdb_data.occupancy(atom_indices);
pdb_data_sel.betaFactor = pdb_data.betaFactor(atom_indices);
pdb_data_sel.element = pdb_data.element(atom_indices);
pdb_data_sel.charge = pdb_data.charge(atom_indices);

end

function JSmolStr = readJSmolSel(filename)

dataLines = [1, Inf];

opts = delimitedTextImportOptions("NumVariables", 1);

opts.DataLines = dataLines;
opts.Delimiter = ",";

opts.PreserveVariableNames = true;
opts.VariableNames = "x__138147153154160_496499501503_528531536_564566_5715";
opts.VariableTypes = "string";

opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

opts = setvaropts(opts, "x__138147153154160_496499501503_528531536_564566_5715", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "x__138147153154160_496499501503_528531536_564566_5715", "EmptyFieldRule", "auto");

JSmolStr = readmatrix(filename, opts);

end

