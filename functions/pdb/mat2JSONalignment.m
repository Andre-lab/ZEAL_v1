%% -- mat2PDB.m --
%
% this function creates a PDB from coordinate data. Represent all inputs as
% a structure field for it to be read. The output format is as given in
% online documentation (as of July 2012 when writing this program)
% http://www.wwpdb.org/documentation/format33/sect9.html#ATOM
%
% Make sure all data input is one-dimensional with the same length. If
% they are not the same length, the program ignores user input, states
% an error, and uses defaults. All string inputs should be in cell-format.
% Keep in mind the "element" and "charge" inputs are strings in
% cell-format, not numbers.
%
%
% -- required inputs (3) --
%
% input value        meaning
%
% input.X            orthagonal X coordinate data (angstroms)
% input.Y            orthagonal Y coordinate data (angstroms)
% input.Z            orthagonal Z coordinate data (angstroms)
%
% -- optional inputs (12): generates defaults when not user-specified --
%
% input value        meaning                           default value
%
% input.outfile      output file name                 "mat2PDB.pdb"
% input.recordName   output record name of atoms      "ATOM"
% input.atomNum      atom serial number                sequential number
% input.atomName     name of atoms                    "OW" (water oxygen)
% input.altLoc       alt. location indicator          " "
% input.resName      name of residue                  "SOL" (water)
%
% input.chainID      protein chain identifier         "A"
% input.resNum       residue sequence number           sequential number
% input.occupancy    occupancy factor                 "1.00"
% input.betaFactor   beta factor, temperature         "0.00"
% input.element      element symbol                   "O" (oxygen)
% input.charge       atomic charge                    " "
%
%
% -- example uses --
%
% % translates both X and Y coordinates of 3IJU.pdb by 5 angstroms
% PDBdata = pdb2mat('3IJU.pdb');
% PDBdata.X = PDBdata.X + 5;
% PDBdata.Y = PDBdata.Y + 5;
% mat2pdb(PDBdata)
%
% % make a PDB with 30 atoms in random places within a 10 angstrom box
% data.X = rand(1,20)*10;
% data.Y = rand(1,20)*10;
% data.Z = rand(1,20)*10;
% mat2pdb(data)
%
%

function JSONstr= mat2JSONalignment(fix, rot)

%% Fixed structure 

% review
% - XYZ coordinate data

input = fix;
% coordinate data is required! Checking XYZ input
if ~isfield(input, 'X') || ~isfield(input, 'Y') || ~isfield(input, 'Z')
    error('Field(s) for XYZ coordinates not found.');
end
X = input.X;
Y = input.Y;
Z = input.Z;
if length(X) ~= length(Y) || length(X) ~= length(Z)
    error('XYZ coordinates not of equal dimension.');
end

% - fields
fixReviewed = reviewInput(input);

input = rot;
% coordinate data is required! Checking XYZ input
if ~isfield(input, 'X') || ~isfield(input, 'Y') || ~isfield(input, 'Z')
    error('Field(s) for XYZ coordinates not found.');
end
X = input.X;
Y = input.Y;
Z = input.Z;
if length(X) ~= length(Y) || length(X) ~= length(Z)
    error('XYZ coordinates not of equal dimension.');
end

% - fields
rotReviewed = reviewInput(input);



%% create PDB

% fix
fixChainID = 'A';
mdlNum = 1;
JSONstr = sprintf('{"pdb":\n"');

mdlnumStr = sprintf('MODEL %8d\n', mdlNum);
JSONstr = [JSONstr mdlnumStr];

% output data
for n = 1:length(fixReviewed.atomNum)
    recName_n = cell2mat(fixReviewed.recordName(n));
    atomName_n = cell2mat(fixReviewed.atomName(n));
    altLoc_n = cell2mat(fixReviewed.altLoc(n));
    resName_n = cell2mat(fixReviewed.resName(n));
    %chainID_n = cell2mat(fixReviewed.chainID(n));
    chainID_n = fixChainID;
    element_n = cell2mat(fixReviewed.element(n));
    charge_n = cell2mat(fixReviewed.charge(n));
    
%     atomNum_n = cell2mat(fixReviewed.atomNum(n));    
    atomNum_n = (fixReviewed.atomNum(n));    

    resNum_n = fixReviewed.resNum(n);
    occupancy_n = fixReviewed.occupancy(n);
    betaFactor_n = fixReviewed.betaFactor(n);
    X_n = fixReviewed.X(n);
    Y_n = fixReviewed.Y(n);
    Z_n = fixReviewed.Z(n);
    
    % standard PDB output line
    pdbLineStr = sprintf('%-6s%5u%5s%1.1s%3s %1.1s%4u%12.3f%8.3f%8.3f%6.2f%6.2f%12s%2s\n', ...
        recName_n, atomNum_n, atomName_n, ...
        altLoc_n, resName_n, chainID_n, ...
        resNum_n, X_n, Y_n, Z_n, occupancy_n, betaFactor_n, ...
        element_n, charge_n);
    
    JSONstr = [JSONstr pdbLineStr];

end

JSONstr = [JSONstr sprintf('TER\n')];
JSONstr = [JSONstr sprintf('ENDMDL\n')];

% rot
rotChainID = 'B';
mdlNum = 2;

mdlnumStr = sprintf('MODEL %8d\n', mdlNum);
JSONstr = [JSONstr mdlnumStr];

% output data
for n = 1:length(rotReviewed.atomNum)
    recName_n = cell2mat(rotReviewed.recordName(n));
    atomName_n = cell2mat(rotReviewed.atomName(n));
    altLoc_n = cell2mat(rotReviewed.altLoc(n));
    resName_n = cell2mat(rotReviewed.resName(n));
    %chainID_n = cell2mat(rotReviewed.chainID(n));
    chainID_n = rotChainID;
    element_n = cell2mat(rotReviewed.element(n));
    charge_n = cell2mat(rotReviewed.charge(n));
   
    
    %atomNum_n = cell2mat(rotReviewed.atomNum(n));    
    atomNum_n = (rotReviewed.atomNum(n)); 
    resNum_n = rotReviewed.resNum(n);
    occupancy_n = rotReviewed.occupancy(n);
    betaFactor_n = rotReviewed.betaFactor(n);
    X_n = rotReviewed.X(n);
    Y_n = rotReviewed.Y(n);
    Z_n = rotReviewed.Z(n);
    
    % standard PDB output line
    pdbLineStr = sprintf('%-6s%5u%5s%1.1s%3s %1.1s%4u%12.3f%8.3f%8.3f%6.2f%6.2f%12s%2s\n', ...
        recName_n, atomNum_n, atomName_n, ...
        altLoc_n, resName_n, chainID_n, ...
        resNum_n, X_n, Y_n, Z_n, occupancy_n, betaFactor_n, ...
        element_n, charge_n);
    
    JSONstr = [JSONstr pdbLineStr];
    
end

JSONstr = [JSONstr sprintf('TER\n')];
JSONstr = [JSONstr sprintf('ENDMDL\n')];
JSONstr = [JSONstr sprintf('END"}')];

end


function inputReviewed = reviewInput(input)

%% review optional data inputs

% in case optional data  not given, fill in blanks
if ~isfield(input, 'outfile')
    input.outfile = 'mat2PDB.pdb';
end
if ~isfield(input, 'recordName')
    input.recordName = cell(1,length(X));
    input.recordName(1:end) = {'ATOM'};
end
if ~isfield(input, 'atomNum')
    input.atomNum = 1:length(X);
end
if ~isfield(input, 'atomName')
    input.atomName = cell(1,length(X));
    input.atomName(1:end) = {'OW'};
end
if ~isfield(input, 'altLoc')
    input.altLoc = cell(1,length(X));
    input.altLoc(1:end) = {' '};
end
if ~isfield(input, 'resName')
    input.resName = cell(1,length(X));
    input.resName(1:end) = {'SOL'};
end
if ~isfield(input, 'chainID')
    input.chainID = cell(1,length(X));
    input.chainID(1:end) = {'A'};
end
if ~isfield(input, 'resNum')
    input.resNum = 1:length(X);
end
if ~isfield(input, 'occupancy')
    input.occupancy = ones(1,length(X));
end
if ~isfield(input, 'betaFactor')
    input.betaFactor = zeros(1, length(X));
end
if ~isfield(input, 'element')
    input.element = cell(1,length(X));
    input.element(1:end) = {'O'};
end
if ~isfield(input, 'charge')
    input.charge = cell(1,length(X));
    input.charge(1:end) = {' '};
end


recordName = input.recordName;
atomNum    = input.atomNum;
atomName   = input.atomName;
altLoc     = input.altLoc;
resName    = input.resName;
chainID    = input.chainID;
resNum     = abs(input.resNum);
occupancy  = input.occupancy;
betaFactor = input.betaFactor;
element    = input.element;
charge     = input.charge;

X     = input.X;


%% remove faulty inputs

if length(recordName) ~= length(X)
    warning('recordName input is not the correct length!\n\tignoring user input\n');
    recordName = cell(1,length(X));
    recordName(1:end) = {'ATOM'};
end
if length(atomNum) ~= length(X)
    warning('atom serial number input is not the correct length!\n\tignoring user input\n');
    atomNum = 1:length(X);
end
if length(atomName) ~= length(X)
    warning('atom name input is not the correct length!\n\tignoring user input\n');
    atomName = cell(1,length(X));
    atomName(1:end) = {'OW'};
end
if length(altLoc) ~= length(X)
    warning('alternate location input is not the correct length!\n\tignoring user input\n');
    altLoc = cell(1,length(X));
    altLoc(1:end) = {' '};
end
if length(resName) ~= length(X)
    warning('residue name input is not the correct length!\n\tignoring user input\n');
    resName = cell(1,length(X));
    resName(1:end) = {'SOL'};
end
if length(chainID) ~= length(X)
    warning('chain ID input is not the correct length!\n\tignoring user input\n');
    chainID = cell(1,length(X));
    chainID(1:end) = {'A'};
end
if length(resNum) ~= length(X)
    warning('residue number input is not the correct length!\n\tignoring user input\n');
    resNum = 1:length(X);
end
if length(occupancy) ~= length(X)
    warning('occupancy input is not the correct length!\n\tignoring user input\n');
    occupancy = ones(1,length(X));
end
if length(betaFactor) ~= length(X)
    warning('beta factor input is not the correct length!\n\tignoring user input\n');
    betaFactor = zeros(1, length(X));
end
if length(element) ~= length(X)
    warning('element symbol input is not the correct length!\n\tignoring user input\n');
    element = cell(1,length(X));
    element(1:end) = {'O'};
end
if length(charge) ~= length(X)
    warning('charge input is not the correct length!\n\tignoring user input\n');
    charge = cell(1,length(X));
    charge(1:end) = {' '};
end

% fix atomName spacing
for n = 1:length(atomName)
    atomName(n) = {sprintf('%-3s',cell2mat(atomName(n)))};
end

inputReviewed.recordName = recordName;
inputReviewed.atomNum = atomNum;
inputReviewed.atomName = atomName;
inputReviewed.altLoc = altLoc;
inputReviewed.resName = resName;
inputReviewed.chainID = chainID;
inputReviewed.resNum = resNum;
inputReviewed.occupancy = occupancy;
inputReviewed.betaFactor = betaFactor;
inputReviewed.element = element;
inputReviewed.charge = charge;

inputReviewed.X = input.X;
inputReviewed.Y = input.Y;
inputReviewed.Z = input.Z;


end