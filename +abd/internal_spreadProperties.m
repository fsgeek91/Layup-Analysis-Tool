function [PROPERTY] = internal_spreadProperties(PROPERTY, nPlies, nSectionPoints)
%   Spread ply-wise properties over section points.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 2.4 Copyright Louis Vallance 2023
%   Last modified 10-May-2023 10:16:13 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%
PROPERTY_i = zeros(1.0, nSectionPoints);
index = 1.0;

for i = 1.0:nPlies
    % Update buffers to include section points
    PROPERTY_i(index:index + (nSectionPoints - 1.0)) = PROPERTY(i);

    index = index + nSectionPoints;
end

% Re-assign buffers to original variables
PROPERTY = PROPERTY_i;