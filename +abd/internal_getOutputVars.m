function [S] = internal_getOutputVars(ABD, Qij, Qt, E_midspan, E_ply_xy, E_ply_aligned, E_therm_xy, E_therm_aligned, E_hydro_xy, E_hydro_aligned, S_ply_xy, S_ply_aligned, EXT, EYT,...
    GXYT, NUXYT, NUYXT, EXB, EYB, GXYB, NUXYB, NUYXB, MSTRS, SFAILRATIO_STRESS, TSAIH, HOFFMAN, TSAIW, AZZIT, MSTRN, SFAILRATIO_STRAIN, HSNFTCRT, SFAILRATIO_HASHIN, HSNFCCRT,...
    HSNMTCRT, HSNMCCRT, LARPFCRT, SFAILRATIO_LARC05, LARMFCRT, LARKFCRT, LARSFCRT, LARTFCRT, UCRT, SFAILRATIO_UCRT, BEST_SEQUENCE, isStrengthOutput, outputLocation, settings,...
    noFailStress, noFailStrain, noHashin, noLaRC05, noUcrt)
%   Collect variables for output.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 5.1.1 Copyright Louis Vallance 2026
%   Last modified 13-Feb-2026 11:01:37 UTC
%
%#ok<*NASGU>

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%
%% Inverse ABD matrix
ABD_INV = inv(ABD);

%% Stiffness matrix
Q = struct('XY', Qij, 'PLY', Qt);

%% Midspan strains
if isempty(E_midspan) == false
    E_MIDSPAN = abd.internal_getTableFromArray(E_midspan, 'strain_midspan');
else
    E_MIDSPAN = [];
end

%% Strain tensors
if isempty(E_ply_xy) == false
    STRAIN = struct('TENSOR_XY', abd.internal_getTableFromArray(E_ply_xy, 'strain_xy'), 'TENSOR_PLY', abd.internal_getTableFromArray(E_ply_aligned, 'strain_aligned'),...
        'THERM_XY', abd.internal_getTableFromArray(E_therm_xy, 'therm_xy'), 'THERM_PLY',  abd.internal_getTableFromArray(E_therm_aligned, 'therm_aligned'),...
        'HYDRO_XY', abd.internal_getTableFromArray(E_hydro_xy, 'hydro_xy'), 'HYDRO_PLY', abd.internal_getTableFromArray(E_hydro_aligned, 'hydro_aligned'));
else
    STRAIN = [];
end

%% Stress tensors
if isempty(S_ply_xy) == false
    STRESS = struct('TENSOR_XY', abd.internal_getTableFromArray(S_ply_xy, 'stress_xy'), 'TENSOR_PLY', abd.internal_getTableFromArray(S_ply_aligned, 'stress_aligned'));
else
    STRESS = [];
end

%% Equivalent extension/bending moduli
if isempty(EXT) == false
    EQ_MODULI = {abd.internal_getTableFromArray([EXT, EYT, GXYT, NUXYT, NUYXT], 'moduli_eq_tension'), abd.internal_getTableFromArray([EXB, EYB, GXYB, NUXYB, NUYXB],...
        'moduli_eq_bending')};
else
    EQ_MODULI = [];
end

%% Failure/damage initiation
if isStrengthOutput == true
    % Add fail stress output to structure
    if noFailStress == false
        CFAILURE.STRESS = abd.internal_getTableFromArray([[MSTRS, SFAILRATIO_STRESS(1.0)]; [TSAIH, SFAILRATIO_STRESS(2.0)]; [HOFFMAN, SFAILRATIO_STRESS(3.0)];...
            [TSAIW, SFAILRATIO_STRESS(4.0)]; [AZZIT, SFAILRATIO_STRESS(5.0)]], 'cfailure_stress');
    end

    % Add fail strain output to structure
    if noFailStrain == false
        CFAILURE.STRAIN = abd.internal_getTableFromArray([MSTRN, SFAILRATIO_STRAIN], 'cfailure_strain');
    end

    % Add Hashin output to structure
    if noHashin == false
        CFAILURE.HASHIN = abd.internal_getTableFromArray([[HSNFTCRT, SFAILRATIO_HASHIN(1.0)]; [HSNFCCRT, SFAILRATIO_HASHIN(2.0)]; [HSNMTCRT, SFAILRATIO_HASHIN(3.0)];...
        [HSNMCCRT, SFAILRATIO_HASHIN(4.0)]], 'cfailure_hashin');
    end

    % Add LaRC05 output to structure
    if noLaRC05 == false
        CFAILURE.LARC05 = abd.internal_getTableFromArray([[LARPFCRT, SFAILRATIO_LARC05(1.0)]; [LARMFCRT, SFAILRATIO_LARC05(2.0)];...
        [LARKFCRT, SFAILRATIO_LARC05(3.0)]; [LARSFCRT, SFAILRATIO_LARC05(4.0)]; [LARTFCRT, SFAILRATIO_LARC05(5.0)]], 'cfailure_larc05');
    end

    % Add user-defined failure criterion output to structure
    if noUcrt == false
        CFAILURE.UCRT = abd.internal_getTableFromArray([UCRT, SFAILRATIO_UCRT], 'cfailure_ucrt'); %#ok<STRNU>
    end
else
    CFAILURE = [];
end

%% Stacking sequence optimiser
if isempty(BEST_SEQUENCE) == false
    % Add optimiser output to the structure
    OPT = struct('SEQUENCE', BEST_SEQUENCE{1.0}, 'CRITICAL_VALUE', BEST_SEQUENCE{2.0}, 'N_PERMUTATIONS', BEST_SEQUENCE{3.0}, 'WALL_CLOCK', BEST_SEQUENCE{4.0}, 'EXCEPTION',...
        BEST_SEQUENCE{5.0});

    % Get the optimised tensors
    TENSOR_OPTIMISED = BEST_SEQUENCE{6.0};

    % Add optimised tensors if they are available
    if isempty(TENSOR_OPTIMISED) == false
        % Convert arrays to tables
        TENSOR_TABLES = {abd.internal_getTableFromArray(TENSOR_OPTIMISED.STRESS_XY, 'stress_xy'), abd.internal_getTableFromArray(TENSOR_OPTIMISED.STRESS_PLY, 'stress_aligned'),...
            abd.internal_getTableFromArray(TENSOR_OPTIMISED.STRAIN_XY, 'strain_xy'), abd.internal_getTableFromArray(TENSOR_OPTIMISED.STRAIN_PLY, 'strain_aligned'),...
            TENSOR_OPTIMISED.SYMMETRIC_ABD};

        % Add cell contents to current structure
        OPT.TENSOR = struct('STRESS_XY', TENSOR_TABLES(1.0), 'STRESS_PLY', TENSOR_TABLES(2.0), 'STRAIN_XY', TENSOR_TABLES(3.0), 'STRAIN_PLY', TENSOR_TABLES(4.0), 'SYMMETRIC_ABD',...
            TENSOR_TABLES(5.0)); %#ok<STRNU>
    end
else
    OPT = [];
end

%% Save workspace variables to a MAT file in the output directory
% Get the non-empty variables
varNames = {'ABD', 'ABD_INV', 'Q', 'E_MIDSPAN', 'STRAIN', 'STRESS', 'EQ_MODULI', 'CFAILURE', 'OPT'};

% Initialize an empty cell array for non-empty variable names
nonEmptyVars = {};

% Loop over variable names to check if they are non-empty
for i = 1:length(varNames)
    if isempty(eval(varNames{i})) == false
        nonEmptyVars{end + 1.0} = varNames{i}; %#ok<AGROW>
    end
end

% Save the output
save([outputLocation, filesep, 'settings.mat'], 'settings');
save([outputLocation, filesep, 'output', '.mat'], nonEmptyVars{:});

%% Assign all non-empty variables to a structure (workspace output)
% Initialise the values list
values = cell(1.0, numel(nonEmptyVars));

for i = 1:numel(nonEmptyVars)
    % Assign the current variable into the CALLER workspace
    assignin('caller', nonEmptyVars{i}, eval(nonEmptyVars{i}));

    % Assign the variable's value to VALUES
    values{i} = evalin('caller', nonEmptyVars{i});
end

% Convert the CELL array to STRUCT
S = cell2struct(values, nonEmptyVars, 2.0);
end