function ZEALexportPdb_all(filename, ZEALsettings, fixModel, rotModel)
%% ZEALexportPdb_all.m
% Save fixed and rotated model in one pdb file as separate models and different chains

[~, name, ~] = fileparts(filename);

fprintf('outputting PDB in file %s.pdb ...', name);
fid = fopen(filename, 'w');
if fid == -1
  error('Author:Function:OpenFile', 'Cannot open file: %s', filename);
end
% write settings used for the alignment to REMARK record
fprintf(fid, '%-12s ZEAL Protein shape alignment %s\n', 'REMARK' , datestr(now));
fprintf(fid, '%-12s \n', 'REMARK');

fixName = split(ZEALsettings.fix.name,'.');
fprintf(fid, '%-12s ZEAL Model 1 (fixed): %s CHAIN %s %s\n', 'REMARK' , fixName{1}, ZEALsettings.fix.chain);
rotName = split(ZEALsettings.rot.name,'.');
fprintf(fid, '%-12s ZEAL Model 2 (rotating): %s CHAIN %s %s\n', 'REMARK' , rotName{1}, ZEALsettings.rot.chain);

if ZEALsettings.globalSearchOp
    fprintf(fid, '%-12s ZEAL Correlation coefficient: %s %5.5f \n', 'REMARK' , ZEALsettings.zcCC);
    fprintf(fid, '%-12s \n', 'REMARK');
    fprintf(fid, '%-12s ZEAL SETTINGS\n', 'REMARK');
    fprintf(fid, '%-12s ZEAL Maximum order of Zernike-Canterakis moments: %d\n', 'REMARK', ZEALsettings.order);
    fprintf(fid, '%-12s ZEAL Grid resolution:\t%d\n', 'REMARK', ZEALsettings.gridRes);
    fprintf(fid, '%-12s ZEAL Shape function:\t%s\n', 'REMARK', ZEALsettings.shapeFcn);
    
    if strcmp(ZEALsettings.shapeFcn,'Molecular surface') || strcmp(ZEALsettings.shapeFcn,'Solvent-accessible surface')
        fprintf(fid, '%-12s ZEAL \t\t\t\t\t Probe radius: %2.2f Å\n', 'REMARK', ZEALsettings.probeRadius);
        fprintf(fid, '%-12s ZEAL \t\t\t\t\t Shell thickness: %d\n', 'REMARK', ZEALsettings.shellThickness);
    elseif strcmp(ZEALsettings.shapeFcn,'Density')
        fprintf(fid, '%-12s ZEAL \t\t\t\t Smear factor: %2.2f\n', 'REMARK', ZEALsettings.smearFactor);
    end
    
end
fprintf(fid, '%-12s\n', 'REMARK');

writeModel(fid, fixModel, 1,'A')

writeModel(fid, rotModel, 2,'B')

fprintf( fid, 'END\n');


% close file
fprintf('   %6.2f%%\n    done! closing file...\n', 100);
fclose(fid);

end

function writeModel(fid, model, modelNo, chainStr)

fprintf(fid, 'MODEL %8d\n', modelNo);

nAtoms = length(model.atomNum);

% output data
try
    for n = 1:nAtoms
        
        % fix atomName spacing
        model.atomName(n) = {sprintf('%-3s',cell2mat(model.atomName(n)))};
        
        % standard PDB output line
        fprintf( fid, '%-6s%5u%5s%1.1s%3s %1.1s%4i%12.3f%8.3f%8.3f%6.2f%6.2f%12s%2s\n', ...
            cell2mat(model.recordName(n)), model.atomNum(n), cell2mat(model.atomName(n)), ...
            cell2mat(model.altLoc(n)), cell2mat(model.resName(n)), chainStr, ...
            model.resNum(n), model.X(n), model.Y(n), model.Z(n), model.occupancy(n), model.betaFactor(n), ...
            cell2mat(model.element(n)), cell2mat(model.charge(n)));
    end
    
catch
    error('Failed to write PDB records.');
end

fprintf(fid, 'TER\n');
fprintf(fid, 'ENDMDL\n');

end

