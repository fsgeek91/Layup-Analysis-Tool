function [PROPERTY] = internal_spreadProperties(PROPERTY, nPlies, nSectionPoints)
%   Spread ply-wise properties over section points.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 3.0.4 Copyright Louis Vallance 2024
%   Last modified 23-Sep-2024 08:11:32 UTC
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