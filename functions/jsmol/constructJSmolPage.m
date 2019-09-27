function [] = constructJSmolPage(JSmolPagePath, width, height, structureName, webGLOp)
%CONSTRUCTJSMOLPAGE Summary of this function goes here
%   Detailed explanation goes here

if webGLOp
    useStr = 'WEBGL HTML5';
else
    useStr = 'HTML5';
end

filename = fullfile(JSmolPagePath, 'JSmol_matlab.htm');
fid = fopen(filename, 'w');
if fid == -1
  error('Author:Function:OpenFile', 'Cannot open file: %s', filename);
end

scriptStr = sprintf('set antialiasDisplay; load ./tmp/%s; background {0.94, 0.94, 0.94};', structureName);

htmlCode = [...
        '<!DOCTYPE html>', ...
            '\n<html>', ...
            '\n<head>', ...
            '\n<meta charset="utf-8">', ...
            '\n<title>JSmol HTML5</title>', ...
            '\n<script type="text/javascript" src="JSmol.min.js"></script>', ...
            '\n<script type="text/javascript" src="JSmol.GLmol.min.js"></script>', ...
            '\n<script type="text/javascript">', ...
            '\nvar jmolApplet0; // set up in HTML table, below', ...
            '\nvar s = document.location.search;', ...
            '\nJmol._debugCode = (s.indexOf("debugcode") >= 0);', ...
            '\nvar Info = {', ...
            '\nwidth: %d,', ...
            '\nheight: %d,', ...
            '\ndebug: false,', ...
            '\ncolor: "0xFFFFFF",', ...
            '\naddSelectionOptions: false,', ...
            '\nuse: "%s",', ...
            '\nj2sPath: "j2s",', ...
            '\nscript: "%s",', ...
            '\ndisableJ2SLoadMonitor: false,', ...
            '\ndisableInitialConsole: false,', ...
            '\nallowJavaScript: false', ...
            '\n}', ...
            '\n\n', ...
            '\n$(document).ready(function() {', ...
            '\n$("#appdiv").html(Jmol.getAppletHtml("jmolApplet0", Info));', ...
            '\n})', ...
            '\n', ...
            '\nvar lastPrompt=0;', ...
            '\n</script>', ...
            '\n</head>', ...
            '\n<body>', ...
            '\n<div id="appdiv" style="width:%dpx"></div>', ...
            '\n</body>', ...
            '\n</html>', ...
            '\n' ...
            ];
        
fprintf(fid, htmlCode, width, height, useStr, scriptStr, width);

fclose(fid);

end

