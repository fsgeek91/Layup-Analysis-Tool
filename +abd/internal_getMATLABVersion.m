function [versionData] = internal_getMATLABVersion()
%   Get the user's MATLAB version.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 4.2.0 Copyright Louis Vallance 2025
%   Last modified 10-Jun-2025 08:28:19 UTC
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