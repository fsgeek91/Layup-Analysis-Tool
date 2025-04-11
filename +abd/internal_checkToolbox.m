function [isAvailable] = internal_checkToolbox(toolboxName)
%   Check for installed toolboxes.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 3.0.5 Copyright Louis Vallance 2025
%   Last modified 11-Apr-2025 10:21:25 UTC
%
    
    %%
    
% Check if the required toolbox is installed
versionData = ver;
isAvailable = any(strcmp(cellstr(char(versionData.Name)), toolboxName));

% Create appdata variable name
switch toolboxName
    case 'Image Processing Toolbox'
        tag = 'qft_noIPT';
    case 'Symbolic Math Toolbox'
        tag = 'qft_noSMT';
    case 'Statistics Toolbox'
        tag = 'qft_noSMLT';
    case 'Statistics and Machine Learning Toolbox'
        tag = 'qft_noSMLT';
    case 'Signal Processing Toolbox'
        tag = 'qft_noSPT';
    case 'Parallel Computing Toolbox'
        tag = 'qft_noPCT';
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