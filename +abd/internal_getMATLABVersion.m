function [versionData] = internal_getMATLABVersion()
%   Get the user's MATLAB version.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 5.1.1 Copyright Louis Vallance 2026
%   Last modified 13-Feb-2026 11:01:37 UTC
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