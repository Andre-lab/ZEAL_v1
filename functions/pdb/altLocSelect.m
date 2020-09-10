function [pdbdata] = altLocSelect(pdbdata,altLoc)
%ALTLOCSELECT Summary of this function goes here
%   Detailed explanation goes here

% Find number of alternative-location identifiers in structure
% and convert the cell array to string array
altLocNames = string(unique(pdbdata.altLoc));

altLocNames(altLocNames==altLoc) = [];
altLocNames(altLocNames=='') = [];

removeList = [];
for n = 1:numel(altLocNames)
    
    removeList_n = find(strcmp(pdbdata.altLoc, altLocNames(n)));
    
    removeList = [removeList removeList_n];
    
end
    
nAtoms = numel(pdbdata.atomNum);
%fprintf('\n Removing %d alt loc entries in pdb', length(removeList));

fn = fieldnames(pdbdata);
for k=1:numel(fn)

    if numel(pdbdata.(fn{k})) == nAtoms
        pdbdata.(fn{k})(removeList) = [];
    end

end

end

