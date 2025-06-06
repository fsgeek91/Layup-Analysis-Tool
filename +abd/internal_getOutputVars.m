function [ABD_INV, Q, E_MIDSPAN, E_PLY, S_PLY, EQ_MODULI, CFAILURE, BEST_SEQUENCE] = internal_getOutputVars(ABD, Qij, Qt, E_midspan, E_ply_xy, E_ply_aligned, E_therm_xy,...
    E_therm_aligned, E_moist_xy, E_moist_aligned, S_ply_xy, S_ply_aligned, EXT, EYT, GXYT, NUXYT, NUYXT, EXB, EYB, GXYB, NUXYB, NUYXB, MSTRS, SFAILRATIO_STRESS, TSAIH, TSAIW,...
    AZZIT, MSTRN, SFAILRATIO_STRAIN, HSNFTCRT, SFAILRATIO_HASHIN, HSNFCCRT, HSNMTCRT, HSNMCCRT, LARPFCRT, SFAILRATIO_LARC05, LARMFCRT, LARKFCRT, LARSFCRT, LARTFCRT, BEST_SEQUENCE)
%   Collect variables for output.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 4.1.0 Copyright Louis Vallance 2025
%   Last modified 06-Jun-2025 11:07:25 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%
% Inverse ABD matrix
ABD_INV = inv(ABD);

% Stiffness matrix
Q = struct('XY', Qij, 'PLY', Qt);

% Midspan 
E_MIDSPAN = abd.internal_getTableFromArray(E_midspan, 'strain_midplane');

% Strain tensors
E_PLY = {abd.internal_getTableFromArray(E_ply_xy, 'strain_xy'), abd.internal_getTableFromArray(E_ply_aligned, 'strain_aligned'),...
    abd.internal_getTableFromArray(E_therm_xy, 'therm_xy'), abd.internal_getTableFromArray(E_therm_aligned, 'therm_aligned'),...
    abd.internal_getTableFromArray(E_moist_xy, 'moist_xy'), abd.internal_getTableFromArray(E_moist_aligned, 'moist_aligned')};

% Stress tensors
S_PLY = {abd.internal_getTableFromArray(S_ply_xy, 'stress_xy'), abd.internal_getTableFromArray(S_ply_aligned, 'stress_aligned')};

% Equivalent extension/bending moduli
if isempty(EXT) == false
    EQ_MODULI = {abd.internal_getTableFromArray([EXT, EYT, GXYT, NUXYT, NUYXT], 'moduli_eq_tension'), abd.internal_getTableFromArray([EXB, EYB, GXYB, NUXYB, NUYXB],...
        'moduli_eq_bending')};
else
    EQ_MODULI = [];
end

% Failure/damage initiation
CFAILURE = struct('STRESS', abd.internal_getTableFromArray([[MSTRS, SFAILRATIO_STRESS(1.0)]; [TSAIH, SFAILRATIO_STRESS(2.0)]; [TSAIW, SFAILRATIO_STRESS(3.0)]; [AZZIT,...
    SFAILRATIO_STRESS(4.0)]], 'cfailure_stress'), 'STRAIN', abd.internal_getTableFromArray([MSTRN, SFAILRATIO_STRAIN], 'cfailure_strain'), 'HASHIN',...
    abd.internal_getTableFromArray([[HSNFTCRT, SFAILRATIO_HASHIN(1.0)]; [HSNFCCRT, SFAILRATIO_HASHIN(2.0)]; [HSNMTCRT, SFAILRATIO_HASHIN(3.0)]; [HSNMCCRT,...
    SFAILRATIO_HASHIN(4.0)]], 'cfailure_hashin'), 'LARC05', abd.internal_getTableFromArray([[LARPFCRT, SFAILRATIO_LARC05(1.0)]; [LARMFCRT, SFAILRATIO_LARC05(2.0)]; [LARKFCRT,...
    SFAILRATIO_LARC05(3.0)]; [LARSFCRT, SFAILRATIO_LARC05(4.0)]; [LARTFCRT, SFAILRATIO_LARC05(5.0)]], 'cfailure_larc05'));

% Stacking sequence optimiser
if isempty(BEST_SEQUENCE) == false
    TENSOR_OPTIMISED = BEST_SEQUENCE{6.0};
    BEST_SEQUENCE{6.0} = {abd.internal_getTableFromArray(TENSOR_OPTIMISED.STRESS_XY, 'stress_xy'), abd.internal_getTableFromArray(TENSOR_OPTIMISED.STRESS_PLY, 'stress_aligned'),...
        abd.internal_getTableFromArray(TENSOR_OPTIMISED.STRAIN_XY, 'strain_xy'), abd.internal_getTableFromArray(TENSOR_OPTIMISED.STRAIN_PLY, 'strain_aligned'),...
        TENSOR_OPTIMISED.SYMMETRIC_ABD};
    BEST_SEQUENCE = struct('SEQUENCE', BEST_SEQUENCE{1.0}, 'CRITICAL_VALUE', BEST_SEQUENCE{2.0}, 'N_PERMUTATIONS', BEST_SEQUENCE{3.0}, 'WALL_CLOCK', BEST_SEQUENCE{4.0}, 'EXCEPTION',...
        BEST_SEQUENCE{5.0}, 'TENSOR', struct('STRESS_XY', BEST_SEQUENCE{6.0}(1.0), 'STRESS_PLY', BEST_SEQUENCE{6.0}(2.0), 'STRAIN_XY', BEST_SEQUENCE{6.0}(3.0), 'STRAIN_PLY',...
        BEST_SEQUENCE{6.0}(4.0), 'SYMMETRIC_ABD', BEST_SEQUENCE{6.0}(5.0)));
end
end