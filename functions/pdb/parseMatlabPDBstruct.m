
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