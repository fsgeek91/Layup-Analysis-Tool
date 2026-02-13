function [isAvailable] = internal_checkToolbox(toolboxName)
%   Check for installed toolboxes.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 5.1.1 Copyright Louis Vallance 2026
%   Last modified 13-Feb-2026 11:01:37 UTC
%
    
    %%
    
% Check if the required toolbox is installed
versionData = ver;
isAvailable = any(strcmp(cellstr(char(versionData.Name)), toolboxName));

% Create appdata variable name
switch toolboxName
    case 'Image Processing Toolbox'
        tag = 'lat_noIPT';
    case 'Symbolic Math Toolbox'
        tag = 'lat_noSMT';
    case 'Statistics Toolbox'
        tag = 'lat_noSMLT';
    case 'Statistics and Machine Learning Toolbox'
        tag = 'lat_noSMLT';
    case 'Signal Processing Toolbox'
        tag = 'lat_noSPT';
    case 'Parallel Computing Toolbox'
        tag = 'lat_noPCT';
    otherwise
        isAvailable = -1.0;
        return
end

% Set the appdata flag accordingly
if isAvailable == 0.0
    setappdata(0, tag, 1.0)
else
    setappdata(0, tag, 0.0)
end
end