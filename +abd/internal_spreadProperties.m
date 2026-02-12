function [PROPERTY] = internal_spreadProperties(PROPERTY, nPlies, nSectionPoints)
%   Spread ply-wise properties over section points.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 5.1.0 Copyright Louis Vallance 2026
%   Last modified 12-Feb-2026 12:33:07 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%
% Initialise output
PROPERTY_i = zeros(1.0, nSectionPoints);

% Initialise loop index
index = 1.0;

for i = 1.0:nPlies
    % Update buffers to include section points
    PROPERTY_i(index:index + (nSectionPoints - 1.0)) = PROPERTY(i);

    % Update the loop index
    index = index + nSectionPoints;
end

% Re-assign buffers to original variables
PROPERTY = PROPERTY_i;