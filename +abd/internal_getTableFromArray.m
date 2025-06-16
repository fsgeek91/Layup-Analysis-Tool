function [T] = internal_getTableFromArray(VAR, TYPE)
%   Convert numeric array to table for workspace output.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 4.2.0 Copyright Louis Vallance 2025
%   Last modified 10-Jun-2025 08:28:19 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%
% Initialise output
T = table();

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
        rowNames = {'E_THERM_XX', 'E_THERM_YY', 'E_THERM_XY'};
    case 'therm_aligned'
        rowNames = {'E_THERM_11', 'E_THERM_22', 'E_THERM_12'};
    case 'hydro_xy'
        rowNames = {'E_HYDRO_XX', 'E_HYDRO_YY', 'E_HYDRO_XY'};
    case 'hydro_aligned'
        rowNames = {'E_HYDRO_11', 'E_HYDRO_22', 'E_HYDRO_12'};
    case 'strain_midspan'
        rowNames = {'EXX_0', 'EYY_0', 'EXY_0', 'KAPPA_XX', 'KAPPA_YY', 'KAPPA_XY'};
        varNames = "MIDSPAN";
    case 'moduli_eq_tension'
        rowNames = {'MODULUS'};
        varNames = ["E_X_T", "E_Y_T", "G_XY_T", "NU_XY_T", "NU_YX_T"];
    case 'moduli_eq_bending'
        rowNames = {'MODULUS'};
        varNames = ["E_X_B", "E_Y_B", "G_XY_B", "NU_XY_B", "NU_YX_B"];
    case 'cfailure_stress'
        rowNames = {'MSTRS', 'TSAIH', 'TSAIW', 'AZZIT'};
        varNames = [strcat("SP_", string(1.0:width(VAR) - 1.0)), "SFAILRATIO"];
    case 'cfailure_strain'
        rowNames = {'MSTRN'};
        varNames = [strcat("SP_", string(1.0:width(VAR) - 1.0)), "SFAILRATIO"];
    case 'cfailure_hashin'
        rowNames = {'HSNFTCRT', 'HSNFCCRT', 'HSNMTCRT', 'HSNMCCRT'};
        varNames = [strcat("SP_", string(1.0:width(VAR) - 1.0)), "SFAILRATIO"];
    case 'cfailure_larc05'
        rowNames = {'LARPFCRT', 'LARMFCRT', 'LARKFCRT', 'LARSFCRT', 'LARTFCRT'};
        varNames = [strcat("SP_", string(1.0:width(VAR) - 1.0)), "SFAILRATIO"];
    otherwise
        % A valid variable type string was not provided, so RETURN
        return
end

% Generate the table
T = array2table(VAR, 'VariableNames', varNames, 'RowNames', rowNames);
end