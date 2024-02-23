function [versionData] = internal_getMATLABVersion()
%   Get the user's MATLAB version.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 3.0.1 Copyright Louis Vallance 2024
%   Last modified 23-Feb-2024 13:20:04 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%
% Initialize output
versionData = {[], []};

% Get the MATLAB release data
releaseData = char(matlabRelease.Release);

% Save the version year as string and double
versionData(1.0) = {releaseData(2.0:end)};
versionData(2.0) = {str2double(releaseData(2.0:5.0))};