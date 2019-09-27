%%  -- pdb2mat.m --
% This program is the most speedy way to read a PDB file that I could come
% up with. It's function is simple: give it a PDB file and out comes a
% matlab-friendly data structure. In cumbersomely large PDB's (such as those that
% include solvent), this can shave off a good amount of time relative to
% many programs. Unfortunately there is no easy way to hasten the slowest
% step, which is turning strings into doubles.

% The output format is as given in online documentation
% (as of July 2012 when writing this program)
% http://www.wwpdb.org/documentation/format33/sect9.html#ATOM
%
% It outputs 14 pieces total of information about the PDB.
%
% -- mandatory information (11) --
%
%     (the name of the PDB, this is the only input on the command line)
%
% recordName (the class or type of atom, such as ATOM, HETATM, SOL, etc)
% atomNum    (serial number of the atom)
% atomName   (elemental identification of the atom)
% altLoc     (alt. location indicator)
% resName    (name of the amino acid/residue)
%
% chainID    (protein chain identifier)
% resNum     (index number of the amino acid)
% X          (X position of atom)
% Y          (Y position of atom)
% Z          (Z position of atom)
%
% -- optional information (4) --
% These are extra data about the atoms. In PDBQT's they hold the partial
% charge, for CHARMM this is the chain name, and so on.
%
% occupancy
% betaFactor
% element
% charge
%

% Optional input
% pdb2mat(pdbStructure.pdb, chainID)
% where chainID is a stringspecifying which chain in the structure that
% should be eextracted by the reader.

% pdb2mat(pdbStructure.pdb, chainID, readHetatoms)
% where chainID is a stringspecifying which chain in the structure that
% should be eextracted by the reader.

% -- example usage: plot the atoms of 3IJU.pdb --
%
%
% PDBdata = pdb2mat('3IJU.pdb');               % read in data from PDB file
% plot3(PDBdata.X, PDBdata.Y, PDBdata.Z, '.'); % make a 3D plot of data
%
% -- example usage: translate the atoms of 3IJU.pdb by 10 angstroms in x direction --
%
% PDBdata = pdb2mat('3IJU.pdb');               % read in data from PDB file
% PDBdata.X = PDBdata.X + 10;                  % translate coordinates
% PDBdata.outfile = '3IJU_tran10angXdir.pdb';  % update file name
% mat2pdb(PDBdata);                            % output data in PDB format
%
%% --- HOW TO MAKE THIS CODE FASTER! --- >> COMMENT OUT WHAT YOU DON'T USE!!
%
% This program reads everything about the PDB by default. If you want a
% faster code for whatever reason, you can comment out the lines you don't
% need. Each numeric data removed (such at resNum, or betaFactor) speeds it
% up by 7-8%. Each string data removed (such as resName or atomName) speeds
% it up by 1-2%.

function [PDBdata, varargout] = pdb2mat(pdbfile, varargin)
%% -- OUTPUT --

if nargin>1
    chainOp = true;
    if isempty(varargin{1})
        chainOp = false;
    elseif ischar(varargin{1})
        chainIDsel = varargin{1};
        
    else
        error('EXPECTED STRING TYPE: The chainID argument has to be a string (for instance ''A'')\n');
    end
else
    chainOp = false;
end

if nargin>2
    
    if islogical(varargin{2})
        hetAtomOp = varargin{2};
    else
        error('EXPECTED LOGICAL TYPE: The readHetatoms-flag has to be logical. Use true/false to include/skip HETATOM records. If flag not set, then all HETATOM records will be inlcuded.');
    end
else
    hetAtomOp = true;
end



if exist(pdbfile,'file')>0
    PDBdata = readpdb(pdbfile);    
    varargout(1) = {true};
    varargout(2) = {false};
else
    try
        [filepath,name,ext] = fileparts(pdbfile);
        pdbid = name(1:4);
        %fprintf('\n Could not find file %s.%s in %s.\nTrying to download structure from RSCB PDB using pdb id = %s\n', name, ext, filepath, pdbid);
                
        matlabPDBdata = getpdb(name(1:4));
        
        PDBdata = parseMatlabPDBstruct(matlabPDBdata);
        varargout(1) = {false};
        varargout(2) = {true};
    catch
        varargout(1) = {false};
        varargout(2) = {false};
        PDBdata = [];
        fprintf('\n Could not download structure %s in the RCSB PDB.', pdbid);        
    end
end

    function PDBdata = parseMatlabPDBstruct(pdbStruct)
        
        nAtoms = size(pdbStruct.Model(1).Atom, 2);
        nHetAtoms = size(pdbStruct.Model(1).HeterogenAtom, 2);
        
        if ~hetAtomOp
            numLines = nAtoms;
        else
            numLines = nAtoms + nHetAtoms;        
        end
        
        recordName = cell(1,numLines);
        atomNum    = cell(1,numLines);
        atomName   = cell(1,numLines);
        altLoc     = cell(1,numLines);
        resName    = cell(1,numLines);
        
        chainID    = cell(1,numLines);
        resNum     = zeros(1,numLines);
        X          = zeros(1,numLines);
        Y          = zeros(1,numLines);
        Z          = zeros(1,numLines);
        
        occupancy  = zeros(1, length(recordName));
        betaFactor = zeros(1, length(recordName));
        element    = cell(1, length(recordName));
        charge     = cell(1, length(recordName));
        
        % get atom records
        m=1;
        for n = 1:nAtoms
            
            recordName(m) = {'ATOM  '};
            atomNum(m)    = {pdbStruct.Model(1).Atom(n).AtomSerNo};
            atomName(m)   = {pdbStruct.Model(1).Atom(n).AtomName};
            altLoc(m)     = {pdbStruct.Model(1).Atom(n).altLoc};
            resName(m)    = {pdbStruct.Model(1).Atom(n).resName};
            
            chainID(m)    = {pdbStruct.Model(1).Atom(n).chainID};
            resNum(m)     = pdbStruct.Model(1).Atom(n).resSeq;
            X(m)          = pdbStruct.Model(1).Atom(n).X;
            Y(m)          = pdbStruct.Model(1).Atom(n).Y;
            Z(m)          = pdbStruct.Model(1).Atom(n).Z;
            
            occupancy(m)  = pdbStruct.Model(1).Atom(n).occupancy;
            betaFactor(m) = pdbStruct.Model(1).Atom(n).tempFactor;
            element(m)    = {pdbStruct.Model(1).Atom(n).element};
            charge(m)    = {pdbStruct.Model(1).Atom(n).charge};
    
            m = m + 1;
        end        
        
        % get hetatom records
        if hetAtomOp
            for n = 1:nHetAtoms
                
                recordName(m) = {'HETATM'};
                atomNum(m)    = {pdbStruct.Model(1).HeterogenAtom(n).AtomSerNo};
                atomName(m)   = {pdbStruct.Model(1).HeterogenAtom(n).AtomName};
                altLoc(m)     = {pdbStruct.Model(1).HeterogenAtom(n).altLoc};
                resName(m)    = {pdbStruct.Model(1).HeterogenAtom(n).resName};
                
                chainID(m)    = {pdbStruct.Model(1).HeterogenAtom(n).chainID};
                resNum(m)     = pdbStruct.Model(1).HeterogenAtom(n).resSeq;
                X(m)          = pdbStruct.Model(1).HeterogenAtom(n).X;
                Y(m)          = pdbStruct.Model(1).HeterogenAtom(n).Y;
                Z(m)          = pdbStruct.Model(1).HeterogenAtom(n).Z;
                
                occupancy(m)  = pdbStruct.Model(1).Atom(n).occupancy;
                betaFactor(m) = pdbStruct.Model(1).Atom(n).tempFactor;
                element(m)    = {pdbStruct.Model(1).Atom(n).element};
                charge(m)    = {pdbStruct.Model(1).Atom(n).charge};
                
                m = m + 1;
                
            end
        end
        
        % trim exess
        if chainOp && ~hetAtomOp
            keepData = (strcmp(recordName,'ATOM  ') + strcmp(chainID,chainIDsel)) == 2;
%             fprintf('\n Keeping all ATOMS with chainID=%s, skipping all HETATOMS', chainIDsel);
        elseif hetAtomOp && ~chainOp
            keepData = logical(strcmp(recordName,'ATOM  ') + strcmp(recordName,'HETATM'));
%             fprintf('\n Keeping all ATOMS and HETATOMS');
        elseif chainOp && hetAtomOp
            keepData = logical( ((strcmp(recordName,'ATOM  ') + strcmp(chainID,chainIDsel))==2) + (( strcmp(chainID,chainIDsel) + strcmp(recordName,'HETATM'))==2) );
%             fprintf('\n Keeping all ATOMS and HETATOMS with chainID=%s', chainIDsel);
        else
            keepData = strcmp(recordName,'ATOM  ');
%             fprintf('\n Keeping all ATOMS, skipping all HETATOMS');
        end
        
        PDBdata.recordName = recordName(keepData);
        PDBdata.atomNum    = atomNum(keepData);
        PDBdata.atomName   = atomName(keepData);
        PDBdata.altLoc     = altLoc(keepData);
        PDBdata.resName    = resName(keepData);
        
        PDBdata.chainID    = chainID(keepData);
        PDBdata.resNum     = resNum(keepData);
        PDBdata.X          = X(keepData);
        PDBdata.Y          = Y(keepData);
        PDBdata.Z          = Z(keepData);
        
        PDBdata.occupancy  = occupancy(keepData);
        PDBdata.betaFactor = betaFactor(keepData);
        PDBdata.element    = element(keepData);
        PDBdata.charge     = charge(keepData);
        
    end


    function PDBdata = readpdb(pdbfile)
        % initialize file
        FileID = fopen(pdbfile);
        rawText = fread(FileID,inf,'*char');
        
        % parse lines by end-of-lines
        splitLines = strread(rawText, '%s', 'delimiter', '\n');
        
        % initialize variables
        numLines = length(splitLines);
        
        recordName = cell(1,numLines);
        atomNum    = cell(1,numLines);
        atomName   = cell(1,numLines);
        altLoc     = cell(1,numLines);
        resName    = cell(1,numLines);
        
        chainID    = cell(1,numLines);
        resNum     = cell(1,numLines);
        X          = cell(1,numLines);
        Y          = cell(1,numLines);
        Z          = cell(1,numLines);
        
        comment    = cell(1,numLines);
        
        % read each line
        m = 1;
        for n = 1:numLines
            
            thisLine = cell2mat(splitLines(n));
            
            if length(thisLine) > 53 && sum(isstrprop(thisLine(23:53), 'alpha')) == 0
                
                recordName(m) = {thisLine(1:6)};
                atomNum(m)    = {thisLine(7:11)};
                atomName(m)   = {thisLine(13:16)};
                altLoc(m)     = {thisLine(17)};
                resName(m)    = {thisLine(18:20)};
                
                chainID(m)    = {thisLine(22)};
                resNum(m)     = {thisLine(23:26)};
                X(m)          = {thisLine(31:38)};
                Y(m)          = {thisLine(39:46)};
                Z(m)          = {thisLine(47:54)};
                
                comment(m)            = {thisLine(55:end)};
                
                m = m + 1;
            end
            
        end
        
        % trim exess
        if chainOp && ~hetAtomOp
            keepData = (strcmp(recordName,'ATOM  ') + strcmp(chainID,chainIDsel)) == 2;
%             fprintf('\n Keeping all ATOMS with chainID=%s, skipping all HETATOMS', chainIDsel);
        elseif hetAtomOp && ~chainOp
            keepData = logical(strcmp(recordName,'ATOM  ') + strcmp(recordName,'HETATM'));
%             fprintf('\n Keeping all ATOMS and HETATOMS');
        elseif chainOp && hetAtomOp
            keepData = logical( ((strcmp(recordName,'ATOM  ') + strcmp(chainID,chainIDsel))==2) + (( strcmp(chainID,chainIDsel) + strcmp(recordName,'HETATM'))==2) );
%             fprintf('\n Keeping all ATOMS and HETATOMS with chainID=%s', chainIDsel);
        else
            keepData = strcmp(recordName,'ATOM  ');
%             fprintf('\n Keeping all ATOMS, skipping all HETATOMS');
        end
        
        recordName = recordName(keepData);
        atomNum    = atomNum(keepData);
        atomName   = atomName(keepData);
        altLoc     = altLoc(keepData);
        resName    = resName(keepData);
        
        chainID    = chainID(keepData);
        resNum     = resNum(keepData);
        X          = X(keepData);
        Y          = Y(keepData);
        Z          = Z(keepData);
        
        comment    = comment(keepData);
        
        % parse out "comment" section
        occupancy  = cell(1, length(recordName));
        betaFactor = cell(1, length(recordName));
        element    = cell(1, length(recordName));
        charge     = cell(1, length(recordName));
        
        % fix spacing
        for n = 1:length(recordName)
            thisLine = sprintf('%-26s',cell2mat(comment(n)));
            occupancy(n)  = {thisLine(1:6)};
            betaFactor(n) = {thisLine(7:12)};
            element(n)    = {thisLine(13:24)};
            charge(n)     = {thisLine(25:26)};
        end
        
        % reformat data for convenience
        PDBdata.recordName = strtrim(recordName);
        PDBdata.atomNum    = str2double(atomNum);
        PDBdata.atomName   = strtrim(atomName);
        PDBdata.altLoc     = altLoc;
        PDBdata.resName    = strtrim(resName);
        
        PDBdata.chainID    = chainID;
        PDBdata.resNum     = abs(str2double(resNum));
        PDBdata.X          = str2double(X);
        PDBdata.Y          = str2double(Y);
        PDBdata.Z          = str2double(Z);
        
        PDBdata.occupancy  = str2double(occupancy);
        PDBdata.betaFactor = str2double(betaFactor);
        PDBdata.element    = strtrim(element);
        PDBdata.charge     = strtrim(charge);
        
        % I commented these lines out, since they cause more problems than they
        % solve. They do clean up the output for certain situations.
        
        % if isnan(PDBdata.occupancy(1))
        %     PDBdata.occupancy = strtrim(PDBdata.occupancy);
        % end
        % if isnan(PDBdata.betaFactor(1))
        %     PDBdata.occupancy = strtrim(PDBdata.betaFactor);
        % end
        
        
        
        % close file
        fclose(FileID);
    end



end
