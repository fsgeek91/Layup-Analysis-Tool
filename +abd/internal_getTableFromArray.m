function [T] = internal_getTableFromArray(VAR, TYPE)
%   Convert numeric array to table for workspace output.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 4.1.0 Copyright Louis Vallance 2025
%   Last modified 06-Jun-2025 11:07:25 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%
% Set the variable names (section points numbers)
varNames = strcat("SP_", string(1.0:width(VAR)));

% Get the row names (components)
switch lower(TYPE)
    case 'stress_xy'
        rowNames = {'S_XX', 'S_YY', 'S_XY'};
    case 'stress_aligned'
        rowNames = {'S_11', 'S_22', 'S_12'};
    case 'strain_xy'
        rowNames = {'E_XX', 'E_YY', 'GAMMA_XY'};
    case 'strain_aligned'
        rowNames = {'E_11', 'E_22', 'GAMMA_12'};
    case 'therm_xy'
        rowNames = {'THERM_XX', 'THERM_YY', 'THERM_XY'};
    case 'therm_aligned'
        rowNames = {'THERM_11', 'THERM_22', 'THERM_12'};
    case 'moist_xy'
        rowNames = {'MOIST_XX', 'MOIST_YY', 'MOIST_XY'};
    case 'moist_aligned'
        rowNames = {'MOIST_11', 'MOIST_22', 'MOIST_12'};
    otherwise
end

% Generate the table
T = array2table(VAR, 'VariableNames', varNames, 'RowNames', rowNames);
end