function [S] = main(settings)
%ABD.MAIN    Analyse a user-defined composite layup.
%   [<output-structure>] = ABD.MAIN(<input-structure>) computes the ABD
%   matrix from an n-layer composite layup, performs a stress analysis
%   based on a user-defined load, and evaluates the layup strength from
%   stress and strain-based failure and damage initiation criteria.
%
%   A summary of the analysis is written to:
%     <output-folder>\<job-name>\summary.log
%
%   Detailed output is written to:
%     <output-folder>\<job-name>\output.mat
%
%   MATLAB figures are written to:
%     <output-folder>\<job-name>\figure\<figure-name>.fig
%
%   THE USER IS NOT REQUIRED TO RUN THIS FUNCTION. SEE USER_DEFINITIONS.M
%   FOR ANALYSIS TEMPLATE.
%
%   The following output is produced:
%   - A, B and D matrices (and their inverses)
%   - Induced midspan strains and curvatures
%   - X-Y stresses and strains based on specified forces and moments
%   - Ply stresses and strains based on specified forces and moments
%   - Equivalent extensional and bending moduli (symmetric layups only)
%   - Stress and strain-based failure criteria
%   - Stress-based damage initiation criteria
%   - Stacking sequence optimisation
%
%   Notes:
%   - Support for unidirectional fibre-reinforced composites only
%   - Strains are assumed to vary linearly through the layup thickness
%     (z-direction)
%   - Linear, orthotropic elasticity is assumed
%   - n: Number of plies defined in the layup
%
%   Units:
%   Force [N]
%   Length [mm]
%   Angle [degrees]
%   Moment [N.mm]
%   Stress [N/mm2]
%   Temperature [degC]
%   Thermal expansion [1/degC]
%   Moisture [%/100 moisture weight content change]
%   Moisture expansion [1/mm]
%
%   Theory References:
%   Python code (GitHub): https://tinyurl.com/48phdvv8
%   Coefficient of Hydroscopic Expansion: https://tinyurl.com/2p9b4mvp
%   Definition of Laminate Plane Element Behaviour:
%   https://tinyurl.com/suwrhfcs
%
%   Resources:
%   Interactive Composite Laminate Calculator: https://tinyurl.com/2nww8yft
%   Finding Stiffness Matrices A, B and D (eFunda):
%   https://tinyurl.com/4tvscjk7
%
%==========================================================================
%__________________________________________________________________________
%   USE CASE I - Stiffness calculation only (ABD matrix output):
%
%   [..] = ABD.MAIN(struct(...
%                          'material', {<1xn>},...
%                          'stacking_sequence', [<1xn>],...
%                          'ply_thickness', [<1xn>],...
%                          'symmetric_layup', [<1x1>],...
%                          'output_location', {<1x2>}),...
%                          ).
%
%   MATERIAL. A 1xn cell array specifying the mechanical material
%   properties:
%     {[E11, E22, G12, V12, A11, A22, B11, B22]}(1),
%     ...
%     {[E11, E22, G12, V12, A11, A22, B11, B22]}(n)
%
%     Note: MATERIAL(1) = Bottom; MATERIAL(n) = Top.
%
%     E11/22: Modulus of Elasticity (longitudinal/transverse)
%     G12: In-plane Shear Modulus
%     V12: In-plane Poisson's Ratio
%     A11/22: Coefficient of Thermal Expansion (longitudinal/transverse)
%     B11/22: Coefficient of Hydroscopic Expansion
%     (longitudinal/transverse)
%
%   STACKING_SEQUENCE. A 1xn array defining the layup stacking sequence
%   angles, THETA:
%     [THETA(1), ... , THETA(n)]
%
%     Note: THETA(1) = Bottom; THETA(n) = Top.
%
%   PLY_THICKNESS. A 1xn array containing the ply thickness values, T:
%     [T(1), ..., T(n)]
%
%     Note: T(1) = Bottom; T(n) = Top.
%
%   SYMMETRIC_LAYUP. A flag to make the calculated section symmetric:
%     true: Mirror the layup definition after the last ply in the stacking
%     sequence (positive-z side)
%     false: Do not mirror the layup
%
%   Note: When the symmetric layup definition is specified, material, layup
%   stacking and ply thickness definitions are automatically reflected on
%   the other side of the symmetry plane.
%
%   OUTPUT_LOCATION. A 1x2 cell array specifying the results location:
%     {'<parameter>', [true | false]}
%
%   OUTPUT_LOCATION(1) is a string specifying the results location:
%     DEFAULT: Save results under the output folder in the current working
%     directory
%     <path>: Specify the directory path directly
%
%   OUTPUT_LOCATION(2) is a flag to enable or disable opening of the
%   analysis summary file after the analysis has been completed:
%     true: Open analysis summary file after analysis
%     false: Do not open results file after analysis
%__________________________________________________________________________
%   USE CASE II - Stress analysis:
%
%   [..] = ABD.MAIN(struct(...
%                          'material', {<1xn>},...
%                          'stacking_sequence', [<1xn>],...
%                          'ply_thickness', [<1xn>],...
%                          'symmetric_layup', [<1x1>],...
%                          'section_points', ['DEFAULT' | sp],...
%                          'load_mech', [<2x3>],...
%                          'load_therm', [<1x1>],...
%                          'load_hydro', [<1x1>],...
%                          'output_ply', [<'param'> | [SP1,..., SPn]],...
%                          'output_figure', {<1x3>},...
%                          'output_location', {<1x2>}),...
%                          ).
%
%   SECTION_POINTS. The number of stress/strain section points per ply:
%     DEFAULT: Program controlled
%     sp: User-defined number of sections points
%
%   Note: Since the layup section is integrated once before the stress
%   analysis, section points are treated as sample points.
%
%   LOAD_MECH. A 2x3 array specifying the forces and moments:
%     [NXX, NYY, NXY; MXX, MYY, MXY]
%
%     NXX/YY: Direct forces along the x and y directions
%     NXY: Shear force in the x-y plane
%     MXX/YY: Moments about the x and y directions
%     MXY: Out-of-plane shear moment
%
%   LOAD_THERM. A 1x1 array specifying the thermal load, DELTA_T.
%
%   LOAD_HYDRO. A 1x1 array specifying the hydroscopic load, DELTA_M.
%
%   OUTPUT_PLY. The section output request:
%     DEFAULT: Program controlled
%       If:
%         SECTION_POINTS = 1: Midspans only
%       Else:
%         Top and bottom faces
%     TOP: Top faces only
%     MIDDLE: Midspans only
%     BOTTOM: Bottom faces only
%     ALL: All section points in all plies
%     ENVELOPEABSMAX: The largest (+ve/-ve) value for the layup
%     ENVELOPEMAX: The largest (+ve) value for the layup
%     ENVELOPEMIN: The largest (-ve) value for the layup
%     [SP1,..., SPn]: User-defined section point list
%
%   OUTPUT_FIGURE. A 1x3 cell array specifying settings for MATLAB figure
%   output:
%     {['DEFAULT' | 'SMOOTH'], ['POINTS'], ['SPLIT' | 'COMPACT']}
%
%   OUTPUT_FIGURE(1) is a parameter specifying the figure type:
%     []: Do not output MATLAB figures
%     DEFAULT: Create MATLAB figures from the raw stress/strain data
%     SMOOTH: Use the built-in MATLAB SMOOTHDATA function to smooth the
%     output at the ply boundaries
%
%   OUTPUT_FIGURE(2) is a parameter specifying the visualisation of the
%   section points:
%     []: Do not plot section points
%     POINTS: Plot the section points over the stress/strain MATLAB
%     figures. If a strength analysis is performed with OUTPUT_STRENGTH,
%     then the section points are coloured green or red for points which
%     are safe or unsafe, respectively
%
%     Note: A section point is coloured red if it failed according to at
%     least on failure/damage initiation criterion)
%
%   OUTPUT_FIGURE(3) is a parameter specifying the figure layout:
%     SPLIT: Create a single MATLAB figure (separate plot for each tensor
%     component)
%     COMPACT: Create a single MATLAB figure (one plot for all tensor
%     components)
%     DETACHED: Create a separate MATLAB figure for each tensor component
%     (one plot for each MATLAB figure)
%__________________________________________________________________________
%   USE CASE III - Strength analysis:
%
%   [..] = ABD.MAIN(struct(...
%                          'material', {<1xn>},...
%                          'fail_stress', {<1xn>},...
%                          'fail_strain', {<1xn>},...
%                          'hashin', {<1xn>},...
%                          'larc05', {<1xn>},...
%                          'stacking_sequence', [<1xn>],...
%                          'ply_thickness', [<1xn>],...
%                          'symmetric_layup', [<1x1>],...
%                          'section_points', ['DEFAULT' | sp],...
%                          'load_mech', [<2x3>],...
%                          'load_therm', [<1x1>],...
%                          'load_hydro', [<1x1>],...
%                          'output_ply', [<'param'> | [SP1,..., SPn]],...
%                          'output_figure', {<1x3>},...
%                          'output_strength', {<1x2>},...
%                          'output_location', {<1x2>}),...
%                          ).
%
%   FAIL_STRESS. A 1xn cell array specifying the strength properties for
%   stress-based failure criteria:
%     {[XT, XC, YT, YC, S, C, B]}(1),
%     ...
%     {[XT, XC, YT, YC, S, C, B]}(n)
%
%     Note: FAIL_STRESS(1) = Bottom; FAIL_STRESS(n) = Top.
%
%     XT/C: Tensile/compressive stress limit (longitudinal)
%     YT/C: Tensile/compressive stress limit (transverse)
%     S: Shear strength in the XY-plane
%     C: Cross-product coefficient
%     B: Biaxial stress limit
%
%     Note: If B = 0, the coupling term is computed from C.
%
%   FAIL_STRAIN. A 1xn cell array specifying the strength properties for
%   strain-based failure criteria:
%     {[XET, XEC, YET, YEC, SE]}(1),
%     ...
%     {[XET, XEC, YET, YEC, SE]}(n)
%
%     Note: FAIL_STRAIN(1) = Bottom; FAIL_STRAIN(n) = Top.
%
%     XET/C: Tensile/compressive strain limit (longitudinal)
%     YET/C: Tensile/compressive strain limit (transverse)
%     SE: Shear strain limit in the XY-plane
%
%   HASHIN. A 1xn cell array specifying the strength properties for the
%   Hashin damage initiation criteria:
%     {[ALPHA, XHT, XHC, YHT, YHC, SHX, SHY]}(1),
%     ...
%     {[ALPHA, XHT, XHC, YHT, YHC, SHX, SHY]}(n)
%
%     Note: HASHIN(1) = Bottom; HASHIN(n) = Top.
%
%     ALPHA: Shear influence parameter;
%     XHT/C: Lamina tensile/compressive strength (longitudinal)
%     YHT/C: Lamina tensile/compressive strength (transverse)
%     SHX/Y: Lamina in-plane/transverse shear strength
%
%   LARC05. A 1xn cell array specifying the strength properties for the
%   LaRC05 damage initiation criteria:
%     {[XLT, XLC, YLT, YLC, SLX, SLY, GL12, NL, NT, A0, PHI0]}(1),
%     ...
%     {[XLT, XLC, YLT, YLC, SLX, SLY, GL12, NL, NT, A0, PHI0]}(n)
%
%     Note: LARC05(1) = Bottom; LARC05(n) = Top.
%
%     XLT/C: Tensile/compressive strength (longitudinal)
%     YLT/C: Tensile/compressive strength (transverse)
%     SLX/Y: In-plane/transverse shear strength
%     GL12: Shear modulus in the 12-plane
%     NL/T: Longitudinal/transverse shear friction coefficient
%     A0: Fracture plane angle for pure compression
%     PHI0: Misalignment angle at failure for pure compression
%
%     Note: The parameters YLC, SLY, NL, NT, A0 and PHI0 are derived if a
%     value of -1 is specified. Unspecified criteria are left empty ( [] ).
%
%   OUTPUT_STRENGTH. A 1x2 cell array specifying settings for the strength
%   assessment:
%     {[true | false], '<parameter>'}
%
%   OUTPUT_STRENGTH(1) is a flag to enable or disable the strength
%   assessment:
%     true: Enable strength assessment
%     false: Disable strength assessment
%
%   OUTPUT_STRENGTH(2) is the failure assessment parameter:
%     RESERVE: Strength reserve factor, R. For stress-based and
%     strain-based failure criteria, the inverse of the strength reserve
%     factor [1/R] is the scaling factor by which the load matrix must be
%     multiplied to hit the failure surface.
%     VALUE: Criterion value, V. The value of the index obtained directly
%     from the failure criterion.
%   
%   Note: The strength assessment is performed on the basis of the
%   available material data given by FAIL_STRESS, FAIL_STRAIN, HASHIN and
%   LARC05.
%
%   Note: The Tsai-Hill, Tsai-Wu and Azzi-Tsai-Hill failure criteria can be
%   expressed in terms of the strength reserve factor or the criterion
%   value; for all other failure criteria, the criterion value is identical
%   to the strength reserve factor.
%
%   Note: For Hashin's theory and LaRC05, R is not evaluated; output for
%   these criteria is quoted as the damage initiation criterion index.
%__________________________________________________________________________
%   USE CASE IV - Stacking sequence optimisation:
%
%   [..] = ABD.MAIN(struct(...
%                          'material', {<1xn>},...
%                          'fail_stress', {<1xn>},...
%                          'fail_strain', {<1xn>},...
%                          'hashin', {<1xn>},...
%                          'larc05', {<1xn>},...
%                          'stacking_sequence', [<1xn>],...
%                          'ply_thickness', [<1xn>],...
%                          'symmetric_layup', [<1x1>],...
%                          'section_points', ['DEFAULT' | sp],...
%                          'load_mech', [<2x3>],...
%                          'load_therm', [<1x1>],...
%                          'load_hydro', [<1x1>],...
%                          'output_ply', [<'param'> | [SP1,..., SPn]],...
%                          'output_figure', {<1x3>},...
%                          'output_strength', {<1x2>},...
%                          'output_optimised', {<1x4>},...
%                          'optimiser_settings', {<1x3>},...
%                          'output_location', {<1x2>}),...
%                          ).
%
%   OUTPUT_OPTIMISED. A 1x4 cell array specifying settings for the stacking
%   sequence optimiser:
%     {'<criterion>', '<parameter>', '<fun>', STEP}
%
%   OUTPUT_OPTIMISED(1) is the failure criterion for the optimisation:
%     <criterion>: Mstrs (Maximum stress); Tsaih (Tsai-Hill);
%     Tsaiw (Tsai-Wu); Azzit (Azzi-Tsai-Hill); Mstrn (Maximum strain);
%     Hashin; LaRC05
%
%   OUTPUT_OPTIMISED(2) is the failure assessment parameter:
%     <parameter>: Reserve (strength reserve factor); Value (criterion
%     value)
%
%   OUTPUT_OPTIMISED(3) is the objective function:
%     <fun>: MinMax (minimise the maximum criterion value); MinMean
%     (minimise the average criterion value)distrubte 
%
%   OUTPUT_OPTIMISED(4) is the angular step size for the stacking sequence
%   permutations, STEP.
%
%   OPTIMISER_SETTINGS. A 1x3 cell array specifying the solver settings for
%   the stacking sequence optimiser:
%     {'<solver>', '<chunk-size>', '<constant>'}
%
%   OPTIMISER_SETTINGS(1) is the solver type:
%     FULL MATRIX: This method computes all stacking sequence combinations
%     before the start of the optimisation. This method is not recommended.
%     For large stacking sequences, the optimisation may exit with an
%     error.
%     MIXED-RADIX: This method uses index-based generation to compute the
%     stacking sequence combinations dynamically. This method has a much
%     lower memory cost than the FULL MATRIX method and can handle very
%     large stacking sequences. It is the recommended method in most cases.
%     CHUNKS: This method uses chunking and worker looping to distribute
%     the stacking sequence combinations between workers. It scales better
%     than the MIXED-RADIX method for large problems, but requires a chunk
%     size to be specified. This method provides no benefit over the
%     MIXED-RADIX method unless a parallel pool is connected.
%
%   OPTIMISER_SETTINGS(2) is the chunk size when the CHUNKS method is
%   specified:
%     DEFAULT: Compute the chunk size automatically based on a heuristic
%     algorithm
%     s: User-specified
%
%   OPTIMISER_SETTINGS(3) is the tuning constant for the chunk size when
%   the CHUNKS method is specified:
%     DEFAULT: Program controlled
%     k: User-defined (typically in the range 2-10)
%
%   Note: The stacking sequence optimiser solver settings are intended for
%   advanced users only. The use of non-standard settings may result in an
%   increase in the analysis time, or the analysis may crash.
%__________________________________________________________________________
%   Optional output arguments:
%
%   [S] = ABD.MAIN(..).
%
%   S. A structure containing the same output as 'output_variables.mat'.
%
%   ABD. A 6x6 matrix of the computed ABD matrix.
%
%   ABD_INV. A 6x6 matrix of the inverse ABD matrix.
%
%   Q. A structure containing the stiffness matrices in X-Y coordinates and
%   ply coordinates.
%
%   Q.XY is a 3x3xn array of the stiffness matrices in X-Y coordinates,
%   where n is the number of plies in the layup.
%
%   Q.PLY is a 3x3xn array of the stiffness matrices in ply coordinates,
%   where n is the number of plies in the layup.
%
%   E_MIDSPAN. A 6x1 table of the midspan strain and curvature components
%   (EXX_0, EYY_0, EXY_0, KAPPA_XX, KAPPA_YY, KAPPA_XY) induced in the
%   laminate. These strains represent the deflections of the laminate about
%   the neutral axis
%
%   STRAIN. A 1x6 structure containing the ply strain components (EXX/E11,
%   EYY/E22, GAMMA_XY/GAMMA_12) for all section points, where n (below) is
%   the total number of section points in the layup.
%
%   STRAIN.TENSOR_XY is a 3xn table of the ply strains in global (X-Y)
%   coordinates.
%
%   STRAIN.TENSOR_PLY is a 3xn table of the ply strains in ply (1-2)
%   coordinates.
%
%   STRAIN.THERM_XY is a 3xn table of the stress-free ply strains due to
%   thermal process in global (X-Y) coordinates.
%
%   STRAIN.THERM_PLY is a 3xn table of the stress-free ply strains due to
%   thermal process in ply (1-2) coordinates.
%
%   STRAIN.HYDRO_XY is a 3xn table of the stress-free ply strains due to
%   moisture process in global (X-Y) coordinates.
%
%   STRAIN.HYDRO_PLY is a 3xn table of the stress-free ply strains due to
%   moisture process in ply (1-2) coordinates.
%
%   Note: For stress-free thermal/hydroscopic strains, contractions have
%   positive values.
%
%   S_PLY. A structure containing the ply stress components (SXX/S11,
%   SYY/S22, SXY/S12) for all section points, where n (below) is the total
%   number of section points in the layup.
%
%   STRESS.TENSOR_XY is a 3xn table of the ply stresses in global (X-Y)
%   coordinates.
%
%   STRESS.TENSOR_PLY is a 3xn table of the ply stresses in ply (1-2)
%   coordinates.
%
%   EQ_MODULI. A 1x2 cell array of the equivalent moduli in tension and
%   bending.
%
%   EQ_MODULI(1) = EXT, EYT, GXYT, NUXYT, NUYXT.
%
%   EQ_MODULI(2) = EXB, EYB, GXYB, NUXYB, NUYXB.
%
%   Note: The equivalent moduli are only calculated for symmetric
%   laminate stacking sequences.
%
%   CFAILURE. A structure containing the failure/damage initiation measure
%   components for all section points.
%
%   CFAILURE.STRESS is a 4xN table of the stress-based failure measure
%   components for all (N-1) section points. The Nth column contains the
%   value of SFAILRATIO for each failure measure component:
%       MSTRS: Maximum stress theory failure measure
%       TSAIH: Tsai-Hill theory failure measure (reserve/value)
%       TSAIW: Tsai-Wu theory failure measure (reserve/value)
%       AZZIT: Azzi-Tsai-Hill theory failure measure (reserve/value)
%
%   CFAILURE.STRAIN is a 1xN table of the strain-based failure measure
%   component for all (N-1) section points. The Nth column contains the
%   value of SFAILRATIO for the failure measure component:
%       MSTRN: Maximum strain theory failure measure
%
%   CFAILURE.HASHIN is a 4xN table of the Hashin damage initiation criteria
%   for all (N-1) section points. The Nth column contains the value of
%   SFAILRATIO for each damage initiation criterion:
%       HSNFTCRT: Hashin’s fibre tensile damage initiation criterion
%       HSNFCCRT: Hashin’s fibre compression damage initiation criterion
%       HSNMTCRT: Hashin’s matrix tensile damage initiation criterion
%       HSNMCCRT: Hashin’s matrix compression damage initiation criterion
%
%   CFAILURE.LARC05 is a 4xN table of the LaCC05 damage initiation criteria
%   for all (N-1) section points. The Nth column contains the value of
%   SFAILRATIO for each damage initiation criterion:
%       LARPFCRT: LaRC05 polymer failure measure
%       LARMFCRT: LaRC05 matrix failure measure
%       LARKFCRT: LaRC05 fibre kinking failure measure
%       LARSFCRT: LaRC05 fibre splitting failure measure
%       LARTFCRT: LaRC05 fibre tensile failure measure
%
%   OPT. A structure containing the results of the stacking sequence
%   optimisation.
%
%   OPT.SEQUENCE is a 1xn array of the optimum stacking sequence,
%   where n is the number of plies in the layup.
%
%   OPT.CRITICAL_VALUE is the critical failure criterion value.
%
%   OPT.N_PERMUTATIONS is the total number of stacking permutations
%   considered by the optimiser.
%
%   OPT.WALL_CLOCK is the analysis time (seconds).
%
%   OPT.EXCEPTION is an exception object returned by the optimiser.
%
%   OPT.TENSOR is a structure containing the stress and strain tensors
%   corresponding to the optimum stacking sequence.
%
%   OPT.TENSOR.STRESS_XY is a 3xn table of the post-optimisation ply
%   stresses in global (X-Y) coordinates.
%
%   OPT.TENSOR.STRESS_PLY is a 3xn table of the post-optimisation ply
%   stresses in ply (1-2) coordinates.
%
%   OPT.TENSOR.STRAIN_XY is a 3xn table of the post-optimisation ply
%   strains in global (X-Y) coordinates.
%
%   OPT.TENSOR.STRAIN_PLY is a 3xn table of the post-optimisation ply
%   strains in ply (1-2) coordinates.
%
%   OPT.TENSOR.SYMMETRIC_ABD is a flag indicating if the optimised
%   stacking sequence is symmetric.
%
%   MATLAB figures:
%
%   EP. MATLAB figure of E_PLY in X-Y coordinates and ply coordinates for
%   all section points.
%
%   SP. MATLAB figure of S_PLY in X-Y coordinates and ply coordinates for
%   all section points.
%
%   CB, Cumulative Density Function (CFD) plot of optimiser criterion for
%   all stacking permutations.
%
%   See also examples, user_definitions.
%
%==========================================================================
%
%   With contributions from: Matt Fig and Jan Simon.
%
%   Permission to user the above author's work is granted under the BSD and
%   CC by-nc-sa 4.0 licenses, where applicable. Third-party source code is
%   clearly indicated in its own subfolder.
%
%   Layup Analysis Tool 4.2.1 Copyright Louis Vallance 2025
%   Last modified 17-Jun-2025 14:50:26 UTC

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%
clc

%% INITIALISE OUTPUT
S = [];

%% GET USER INPUTS FROM VARARGIN
[enableTensor, printTensor, materialDataMechanical, materialDataFailStress, materialDataFailStrain, materialDataHashin, materialDataLaRC05, theta, t_ply, symmetricPly,...
    SECTION_POINTS, OUTPUT_PLY, OUTPUT_FIGURE, OUTPUT_STRENGTH, OUTPUT_OPTIMISED, OPTIMISER_SETTINGS, OUTPUT_LOCATION, Nxx, Nyy, Nxy, Mxx, Myy, Mxy, deltaT, deltaM, jobName,...
    jobDescription, settings, error] =...
    ...
    abd.internal_initialise(settings);

% An error occurred, so RETURN
if error == true
    return
end

%% MIRROR THE LAYUP DEFINITION (IF APPLICABLE)
[t_ply, theta, nPlies, error] =...
    ...
    abd.internal_mirror(symmetricPly, t_ply, theta);

% An error occurred, so RETURN
if error == true
    return
end

%% PROCESS OUTPUT_FIGURE
[error, OUTPUT_FIGURE] =...
    ...
    abd.internal_plot.getSettings(OUTPUT_FIGURE);

% An error occurred, so RETURN
if error == true
    return
end

%% PROCESS OUTPUT_STRENGTH
[error, OUTPUT_STRENGTH] =...
    ...
    abd.internal_strength.getSettings(OUTPUT_STRENGTH);

% An error occurred, so RETURN
if error == true
    return
end

%% GET MATERIAL DATA (MECHANICAL)
[error, ~, E11, E22, G12, V12, A11, A22, B11, B22] =...
    ...
    abd.internal_getMaterial(materialDataMechanical, nPlies, symmetricPly, 1.0, 'MATERIAL');

% An error occurred, so RETURN
if error == true
    return
end

%% GET MATERIAL DATA (STRENGTH)
if OUTPUT_STRENGTH{1.0} == true
    % Get fail stress properties
    [error, noFailStress, XT, XC, YT, YC, S, C, B] =...
        ...
        abd.internal_getMaterial(materialDataFailStress, nPlies, symmetricPly, 2.0, 'FAIL_STRESS');

    % An error occurred, so RETURN
    if error == true
        return
    end

    % Correct the sign (if applicable)
    XT = abd.internal_correctSign(XT, 1.0);
    XC = abd.internal_correctSign(XC, -1.0);
    YT = abd.internal_correctSign(YT, 1.0);
    YC = abd.internal_correctSign(YC, -1.0);
    S = abd.internal_correctSign(S, 1.0);
    
    % Get fail strain properties
    [error, noFailStrain, XET, XEC, YET, YEC, SE] =...
        ...
        abd.internal_getMaterial(materialDataFailStrain, nPlies, symmetricPly, 3.0, 'FAIL_STRAIN');

    % An error occurred, so RETURN
    if error == true
        return
    end

    % Correct the sign (if applicable)
    XET = abd.internal_correctSign(XET, 1.0);
    XEC = abd.internal_correctSign(XEC, -1.0);
    YET = abd.internal_correctSign(YET, 1.0);
    YEC = abd.internal_correctSign(YEC, -1.0);
    SE = abd.internal_correctSign(SE, 1.0);

    % Get Hashin properties
    [error, noHashin, ALPHA, XHT, XHC, YHT, YHC, SHX, SHY] =...
        ...
        abd.internal_getMaterial(materialDataHashin, nPlies, symmetricPly, 2.0, 'HASHIN');

    % An error occurred, so RETURN
    if error == true
        return
    end

    % Correct the sign (if applicable)
    XHT = abd.internal_correctSign(XHT, 1.0);
    XHC = abd.internal_correctSign(XHC, 1.0);
    YHT = abd.internal_correctSign(YHT, 1.0);
    YHC = abd.internal_correctSign(YHC, 1.0);
    SHX = abd.internal_correctSign(SHX, 1.0);
    SHY = abd.internal_correctSign(SHY, 1.0);

    % Get LaRC05 properties
    [error, noLaRC05, XLT, XLC, YLT, YLC, SLX, SLY, GL12, NL, NT, A0, PHI0] =...
        ...
        abd.internal_getMaterial(materialDataLaRC05, nPlies, symmetricPly, 4.0, 'LARC05');

    % An error occurred, so RETURN
    if error == true
        return
    end

    % Correct the sign (if applicable)
    XLT = abd.internal_correctSign(XLT, 1.0);
    XLC = abd.internal_correctSign(XLC, 1.0);
    YLT = abd.internal_correctSign(YLT, 1.0);
    YLC = abd.internal_correctSign(YLC, 2.0);
    SLX = abd.internal_correctSign(SLX, 1.0);
    SLY = abd.internal_correctSign(SLY, 2.0);
    GL12 = abd.internal_correctSign(GL12, 1.0);

    if (noFailStress == true) && (noFailStrain == true) && (noHashin == true) && (noLaRC05 == true)
        %{
            Since the strength calculation has been requested, at least
            one of FAIL_STRESS, FAIL_STRAIN, HASHIN or LARC05 must be
            defined for the layup!
        %}
        fprintf('[ERROR] The strength calculation requires at least\nFAIL_STRESS, FAIL_STRAIN, HASHIN or LARC05 material properties\n');
        return
    end

    if noLaRC05 == false
        % Check if the symbolic math toolbox is available for LaRC05
        symsAvailable = abd.internal_checkToolbox('Symbolic Math Toolbox');
    else
        symsAvailable = [];
    end
else
    noFailStress = true;
    noFailStrain = true;
    noHashin = true;
    noLaRC05 = true;
end

%% SET THE NEAR-ZERO TOLERANCE VALUE
tolerance = 1e-6;

%% GET THICKNESS FRACTIONS
[z, t] = abd.internal_getThickness(nPlies, t_ply, tolerance);

%% PROCESS SECTION_POINTS
[error, z_points, theta_points, nPlies_points, A11_points, A22_points, B11_points, B22_points, plyBuffer, thickness, SECTION_POINTS, OUTPUT_PLY] =...
    ...
    abd.internal_getSectionPoints(SECTION_POINTS, 'SECTION_POINTS', nPlies, theta, z, A11, A22, B11, B22, tolerance, OUTPUT_PLY);

% An error occurred, so RETURN
if error == true
    return
end

%% INITIALISE VARIABLES
BEST_SEQUENCE = [];

%% PROCESS OUTPUT_PLY
[error, OUTPUT_PLY_POINTS, plyBuffer, OUTPUT_ENVELOPE, ENVELOPE_MODE, outputApproximate, plyBuffer_sfailratio] =...
    ...
    abd.internal_getOutputPoints(OUTPUT_PLY, z, z_points, nPlies, nPlies_points, plyBuffer, SECTION_POINTS, tolerance, enableTensor);

% An error occurred, so RETURN
if error == true
    return
end

%% GET OPTIMISER SETTINGS
if isempty(OUTPUT_OPTIMISED{1.0}) == false
    [error, OUTPUT_OPTIMISED] =...
        ...
        abd.internal_optimise.getSettings(OUTPUT_OPTIMISED, noFailStress, noFailStrain, noHashin, noLaRC05, OUTPUT_STRENGTH{1.0});

    % An error occurred, so RETURN
    if error == true
        return
    end
end

%% COMPUTE REDUCED STIFFNESS TERMS
[Q11, Q22, Q66, Q12, Qij] =...
    ...
    abd.internal_getReducedQ(nPlies, E11, E22, V12, G12);

%% COMPUTE TRANSFORMED REDUCED STIFFNESS MATRIX COMPONENTS
[Q11t, Q12t, Q16t, Q22t, Q26t, Q66t, Qt] =...
    ...
    abd.internal_getTransformedQ(nPlies, theta, Q11, Q12, Q66, Q22);

%% GET EFFECTIVE THEMAL AND MOISTURE EXPANSION COEFFICIENTS FOR EACH PLY
[axx, ayy, axy, bxx, byy, bxy] =...
    ...
    abd.internal_getThermoHydro(theta_points, A11_points, A22_points, B11_points, B22_points);

%% COMPUTE A, B and D MATRICES
[ABD, ABD_INV, Qijt, NxxT, NyyT, NxyT, MxxT, MyyT, MxyT, NxxM, NyyM, NxyM, MxxM, MyyM, MxyM] =...
    ...
    abd.internal_getABD(nPlies, Q11t, Q12t, Q16t, Q22t, Q26t, Q66t, z, nargin, deltaT, deltaM, axx, ayy, axy, bxx, byy, bxy, SECTION_POINTS);

%% COMPUTE TENSOR QUANTITIES
if enableTensor == true
    [E_midspan, E_ply_xy, S_ply_xy, E_ply_aligned, S_ply_aligned, E_therm_xy, E_hydro_xy, E_therm_aligned, E_hydro_aligned] =...
        ...
        abd.internal_getTensor(ABD, Nxx, NxxT, NxxM, Nyy, NyyT, NyyM, Nxy, NxyT, NxyM, Mxx, MxxT, MxxM, Myy, MyyT, MyyM, Mxy, MxyT, MxyM, nPlies_points, z_points, theta_points,...
        Qijt, deltaT, deltaM, axx, ayy, axy, bxx, byy, bxy, tolerance);

    if any(E_midspan) == false
        printTensor = -1.0;
    end
else
    % Initialise values to default
    printTensor = false;
    E_midspan = [];
    E_ply_xy = [];  E_ply_aligned = [];
    S_ply_xy = [];  S_ply_aligned = [];
    E_therm_xy = [];    E_therm_aligned = [];
    E_hydro_xy = [];    E_hydro_aligned = [];
end

%% DETERMINE IF ABD MATRIX IS SYMMETRIC
symmetricAbd = abd.internal_getSymmetry(ABD, tolerance);

%% GET THE EQUIVALENT MODULI
if symmetricAbd == true
    [EXT, EYT, GXYT, NUXYT, NUYXT, EXB, EYB, GXYB, NUXYB, NUYXB] =...
        ...
    abd.internal_getModuli(t, ABD_INV);
else
    % Equivalent moduli are not calculated for unsymmetric layups
    EXT = []; EYT = []; GXYT = []; NUXYT = []; NUYXT = []; EXB = [];
    EYB = []; GXYB = []; NUXYB = []; NUYXB = [];
end

%% ROUND SMALL ABD VALUES TO ZERO
ABD(abs(ABD) < tolerance) = 0.0;

%% PERFORM STRENGTH CALCULATION ON PLY STRESSES
% Preallocate failed section points colour buffer (for plotting)
SP_COLOUR_BUFFER = repmat([0.0, 0.0, 0.0], [nPlies_points, 1.0]);

% Initialise chunk size and number of chunks
CHUNK_SIZE = [];
N_CHUNKS = [];
EXECUTION_MODE = [];

% Initialise failure criteria component buffers
[MSTRS, TSAIH, TSAIW, AZZIT, MSTRN, HSNFTCRT, HSNFCCRT, HSNMTCRT, HSNMCCRT, LARPFCRT, LARMFCRT, LARKFCRT, LARSFCRT, LARTFCRT] = abd.internal_strength.init(nPlies_points);

if (OUTPUT_STRENGTH{1.0} == true) && (printTensor == 1.0)
    [MSTRS, TSAIH, TSAIW, AZZIT, MSTRN, HSNFTCRT, HSNFCCRT, HSNMTCRT, HSNMCCRT, LARPFCRT, LARMFCRT, LARKFCRT, LARSFCRT, LARTFCRT, XT, XC, YT, YC, S, C, B, E11, E22, G12, V12, XET,...
        XEC, YET, YEC, SE, ALPHA, XHT, XHC, YHT, YHC, SHX, SHY, XLT, XLC, YLT, YLC, SLX, SLY, GL12, NL, NT, A0, PHI0, S1, S2, S3] =...
        ...
        abd.internal_strength.main(noFailStress, noFailStrain, noHashin, noLaRC05, symsAvailable, XT, XC, YT, YC, S, C, B, E11, E22, G12, V12, XET, XEC, YET, YEC, SE, ALPHA, XHT,...
        XHC, YHT, YHC, SHX, SHY, XLT, XLC, YLT, YLC, SLX, SLY, GL12, NL, NT, A0, PHI0, S_ply_aligned, nPlies, nPlies_points, SECTION_POINTS, OUTPUT_STRENGTH{2.0}, MSTRS, TSAIH,...
        TSAIW, AZZIT, MSTRN, HSNFTCRT, HSNFCCRT, HSNMTCRT, HSNMCCRT, LARPFCRT, LARMFCRT, LARKFCRT, LARSFCRT, LARTFCRT);

    if OUTPUT_OPTIMISED{1.0} == true
        %% FIND THE OPTIMUM STACKING SEQUENCE
        [BEST_SEQUENCE, CRITERION_BUFFER, ~, CHUNK_SIZE, N_CHUNKS, EXECUTION_MODE] =...
            ...
            abd.internal_optimise.main(OUTPUT_OPTIMISED, nargin, nPlies, nPlies_points, SECTION_POINTS, z, z_points, Q11, Q22, Q66, Q12, A11_points, A22_points, B11_points,...
            B22_points, tolerance, XT, XC, YT, YC, S, C, B, XET, XEC, YET, YEC, SE, ALPHA, XHT, XHC, YHT, YHC, SHX, SHY, XLT, XLC, YLT, YLC, SLX, SLY, GL12, NL, NT, A0, PHI0,...
            deltaT, deltaM, Nxx, Nyy, Nxy, Mxx, Myy, Mxy, E11, E22, V12, G12, symsAvailable, S1, S2, S3, SECTION_POINTS, OPTIMISER_SETTINGS);
    else
        CRITERION_BUFFER = [];
    end

    % Set failed section points colour buffer (for plotting)
    SP_COLOUR_BUFFER = abd.internal_strength.getFailedSpBuffer(MSTRS, TSAIH, TSAIW, AZZIT, MSTRN, HSNFTCRT, HSNFCCRT, HSNMTCRT, HSNMCCRT, LARPFCRT, LARMFCRT, LARKFCRT, LARSFCRT,...
        LARTFCRT, SP_COLOUR_BUFFER, nPlies_points);
else
    % Initialise values to default
    CRITERION_BUFFER = [];

    % Suppress strength output
    OUTPUT_STRENGTH{1.0} = false;
end

%% CREATE OUTPUT DIRECTORY
% Create the root output folder if it does not already exist
if exist(OUTPUT_LOCATION{1.0}, 'dir') ~= 7.0
    mkdir(OUTPUT_LOCATION{1.0})
end

% Add the root output folder to the MATLAB path
addpath(OUTPUT_LOCATION{1.0})

% Get the date string for the output folder
dateString = char(datetime('now'));
dateStringFile = dateString;

% Replace unsupported characters for file name string
for i = 1:length(dateStringFile)
    if (strcmpi(dateStringFile(i), ':') == true) || (strcmpi(dateStringFile(i), ' ') == true)
        dateStringFile(i) = '_';
    end
end

% Construct the output location path
outputLocation = [OUTPUT_LOCATION{1.0}, filesep, jobName];

% Create the folder if it does not already exist
if exist(outputLocation, 'dir') ~= 7.0
    mkdir(outputLocation)
end

%% PLOT STRAINS AND STRESSES IN A MATLAB FIGURE
if (isempty(OUTPUT_FIGURE{1.0}) == false) && (printTensor == 1.0) && (nPlies_points > 1.0)
    % Construct the output location path
    outputLocationFigure = [OUTPUT_LOCATION{1.0}, filesep, jobName, filesep, 'figure'];

    % Create the MATLAB figure folder if it does not already exist
    if exist(outputLocationFigure, 'dir') ~= 7.0
        mkdir(outputLocationFigure)
    end

    % Create the MATLAB figures
    abd.internal_plot.main(OUTPUT_FIGURE{1.0}, OUTPUT_FIGURE{2.0}, OUTPUT_FIGURE{3.0}, outputLocationFigure, nPlies, E_ply_xy, S_ply_xy, E_ply_aligned, S_ply_aligned, z, z_points,...
        CRITERION_BUFFER, OUTPUT_OPTIMISED, SP_COLOUR_BUFFER)
end

%% WRITE RESULTS TO A TEXT FILE
[SFAILRATIO_STRESS, SFAILRATIO_STRAIN, SFAILRATIO_HASHIN, SFAILRATIO_LARC05] =...
    abd.internal_outputToFile(dateString, outputLocation, OUTPUT_STRENGTH, nPlies, t_ply, theta, enableTensor, printTensor, S_ply_aligned, S_ply_xy, E_ply_aligned, E_ply_xy,...
    E_therm_xy, E_hydro_xy, E_therm_aligned, E_hydro_aligned, ABD, symmetricAbd, EXT, EYT, GXYT, NUXYT, NUYXT, EXB, EYB, GXYB, NUXYB, NUYXB, MSTRS, TSAIH, TSAIW, AZZIT, MSTRN,...
    HSNFTCRT, HSNFCCRT, HSNMTCRT, HSNMCCRT, LARPFCRT, LARMFCRT, LARKFCRT, LARSFCRT, LARTFCRT, noFailStress, noFailStrain, noHashin, noLaRC05, SECTION_POINTS, OUTPUT_PLY_POINTS,...
    plyBuffer, thickness, OUTPUT_ENVELOPE, ENVELOPE_MODE, outputApproximate, BEST_SEQUENCE, OUTPUT_OPTIMISED, OUTPUT_FIGURE{1.0}, plyBuffer_sfailratio, axx, ayy, axy, bxx, byy,...
    bxy, E_midspan, OUTPUT_PLY, z_points, OPTIMISER_SETTINGS, CHUNK_SIZE, N_CHUNKS, EXECUTION_MODE, jobName, jobDescription);

%% COLLECT OUTPUT
[S] = abd.internal_getOutputVars(ABD, Qij, Qt, E_midspan, E_ply_xy, E_ply_aligned, E_therm_xy, E_therm_aligned, E_hydro_xy, E_hydro_aligned, S_ply_xy, S_ply_aligned, EXT, EYT, GXYT,...
    NUXYT, NUYXT, EXB, EYB, GXYB, NUXYB, NUYXB, MSTRS, SFAILRATIO_STRESS, TSAIH, TSAIW, AZZIT, MSTRN, SFAILRATIO_STRAIN, HSNFTCRT, SFAILRATIO_HASHIN, HSNFCCRT, HSNMTCRT, HSNMCCRT,...
    LARPFCRT, SFAILRATIO_LARC05, LARMFCRT, LARKFCRT, LARSFCRT, LARTFCRT, BEST_SEQUENCE, OUTPUT_STRENGTH{1.0}, outputLocation, settings, noFailStress, noFailStrain, noHashin, noLaRC05);

%% Add the output location to the MATLAB path
addpath(genpath(outputLocation));

%% Open the results file now (if applicable)
if OUTPUT_LOCATION{2.0} == true
    try
        open([outputLocation, filesep, 'summary.log'])
    catch
        % Do nothing
    end
end
end