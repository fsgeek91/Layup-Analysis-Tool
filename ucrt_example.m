function [UCRT] = ucrt_example(INFO, UCRT, MATERIAL_MECH, MATERIAL_FAIL, TENSORS)
%UCRT_EXAMPLE    Example routine for user-defined failure criterion.
%    To use this routine for a strength calculation with Layup Analysis
%    Tool, specify the following in the layup definition file:
%      OUTPUT_STRENGTH = {@ucrt_example, {'RESERVE' | 'VALUE'}};
%
%    Routine interface
%      [UCRT] = ucrt_example(INFO, UCRT, MATERIAL_MECH, MATERIAL_FAIL, TENSORS)
%
%    Variables to be provided to the routine
%      INFO [struct], Layup information
%        PLY_COUNT: Number of plies (layup)
%        SECTION_POINTS: Number of section points (per ply)
%        TOTAL_POINTS: Total number of sample points (layup)
%        FAILURE_PARAMETER: 1 (strength reserve factor); 2 (criterion value)
%      UCRT [1xTOTAL_POINTS double], Pre-initialised buffer for failure criterion values
%      MATERIAL_MECH [struct], Mechanical material properties
%      MATERIAL_FAIL [struct], Composite strength properties
%        STRESS: Stress-based failure criteria
%        STRAIN: Strain-based failure criteria
%        HASHIN: Hashin damage initiation criteria
%        LARC05: LaRC05 damage initiation criteria
%      TENSORS [struct], Stress/strain tensor data
%        E_MIDSPAN [PLY_COUNTx1 double], Midspan strains
%        E_PLY_XY [3xTOTAL_POINTS double], Ply strain in global (X-Y) coordinates
%        S_PLY_XY [3xTOTAL_POINTS double], Ply stress in global (X-Y) coordinates
%        E_PLY_ALIGNED [3xTOTAL_POINTS double], Ply strain in ply (1-2) coordinates
%        S_PLY_ALIGNED [3xTOTAL_POINTS double], Ply stress in ply (1-2) coordinates
%        E_THERM_XY [3xTOTAL_POINTS double], Stress-free ply strains due to thermal process in global (X-Y) coordinates
%        E_HYDRO_XY [3xTOTAL_POINTS double], Stress-free ply strains due to moisture process in global (X-Y) coordinates
%        E_THERM_ALIGNED [3xTOTAL_POINTS double], Stress-free ply strains due to thermal process in ply (1-2) coordinates
%        E_HYDRO_ALIGNED [3xTOTAL_POINTS double], Stress-free ply strains due to moisture process in ply (1-2) coordinates
%
%    NOTE: For a description of the mechanical material/composite strength
%    property identifieres, run the command >> help abd.main
%    and consult the sections "USE CASE I" and "USE CASE III",
%    respectively.
%
%    NOTE: The tensor variables are arranged as follows
%    (1, :) -> 11-component of stress/strain
%    (2, :) -> 22-component of stress/strain
%    (3, :) -> 12-component of stress/strain
%
%    NOTE: The value of FAILURE_PARAMETER is the requested failure
%    parameter. It is used for information purposes in case the author
%    wishes to provide strength reserve and criterion value output.
%
%    Variables returned from the routine
%      UCRT [1xTOTAL_POINTS double], Computed failure criterion values
%
%    FILE AUTO GNERATED BY LZ16-E5O-AUT ON 12-Feb-2026 09:37:52.
%
%    Insert code below:
%%
% Get the number of result points for UCRT output
N = INFO.TOTAL_POINTS;

% Get the stress tensor components
rows = num2cell(TENSORS.S_PLY_ALIGNED, 2.0);
[S11, S22, S12] = deal(rows{:});

% Get the stress-based failure parameters
STRESS = MATERIAL_FAIL.STRESS;
[XT, XC, YT, YC, S, C, B] = deal(STRESS.XT, STRESS.XC, STRESS.YT, STRESS.YC, STRESS.S, STRESS.C, STRESS.B);

% Tension-compression split (longitudinal)
X(S11 >= 0.0) = XT(S11 >= 0.0);
X(S11 < 0.0) = XC(S11 < 0.0);

% Tension-compression split (transverse)
Y(S22 >= 0.0) = YT(S22 >= 0.0);
Y(S22 < 0.0) = YC(S22 < 0.0);

% Compute individual safety factors in 11, 22 and 12 directions
[SF11, SF22, SF12] = deal(abs(S11./X), abs(S22./Y), abs(S12./S));

% Get criterion from overall maximum
UCRT = max([SF11', SF22', SF12'], [], 2.0)';
%__________________________________________________________________________
%%    END OF FILE
end
