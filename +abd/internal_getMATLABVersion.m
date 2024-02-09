function [matlabVersion] = internal_getMATLABVersion()
%   Get the user's MATLAB version.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 2.7 Copyright Louis Vallance 2024
%   Last modified 09-Feb-2024 09:10:19 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%
% Initialize output
matlabVersion = {[], []};

% Save the version year as string and double
matlabVersion(1.0) = {versionNumber.Release(3.0:7.0)};
matlabVersion(2.0) = {str2double(versionNumber.Release(3.0:6.0))};