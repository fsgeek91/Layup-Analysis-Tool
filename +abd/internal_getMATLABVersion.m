function [matlabVersion] = internal_getMATLABVersion()
%   Get the user's MATLAB version.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 2.6 Copyright Louis Vallance 2023
%   Last modified 17-May-2023 07:40:13 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%
% Initialize output
matlabVersion = {[], []};

% Get the MATLAB version year
versionNumber = ver('MATLAB');
versionYear = versionNumber.Release;

% Save the version year as string and double
matlabVersion(1.0) = {versionYear(3.0:7.0)};
matlabVersion(2.0) = {str2double(versionYear(3.0:6.0))};