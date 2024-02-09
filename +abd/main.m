function [varargout] = main(varargin)
%ABD.MAIN    Analyse a user-defined composite layup.
%   [VARARGOUT] = ABD.MAIN(VARARGIN) computes the ABD matrix from an
%   n-layer composite layup definition, and evaluates the layup strength
%   based on a range of failure and damage initiation criteria.
%
%   THE USER IS NOT REQUIRED TO RUN THIS FUNCTION.
%
%   The following output is produced:
%
%   - A, B and D matrices (and their inverses)
%   - Induced midplane strains and curvatures
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
%
%   ABD matrix calculation only:
%
%   [..] = ABD.MAIN(MATERIAL, LAYUP, OUTPUT_DEF).
%
%   MATERIAL. A 1x4 cell array specifying mechanical and strength material
%   properties.
%
%   MATERIAL(1) is a 1xn cell array specifying the mechanical material
%   properties, E11, E22, G12, V12, A11, A22, B11 and B22 for each ply.
%
%   MATERIAL(2) is a 1xn cell array specifying the strength properties for
%   stress-based failure criteria, XT, XC, YT, YC, S, C and B for each ply.
%
%   MATERIAL(3) is a 1xn cell array specifying the strength properties for
%   strain-based failure criteria, XET, XEC, YET, YEC and SE for each ply.
%
%   MATERIAL(4) is a 1xn cell array specifying the strength properties for
%   the Hashin damage initiation criteria ALPHA, XHT, XHC, YHT, YHC, SHX
%   and SHY for each ply.
%
%   Note: n = 1 for constant material properties; for ply-wise material
%   properties, specify n sets of properties corresponding to the n ply
%   definitions given by LAYUP.
%
%   LAYUP. A 1x4 cell specifying the layup orientations, laminate thickness
%   values, stacking symmetry and the number of section points.
%
%   LAYUP(1) is a 1xn array defining the layup stacking sequence,
%   STACKING_SEQUENCE.
%
%   LAYUP(2) is a 1xn array containing the ply thickness values,
%   PLY_THICKNESS: n = 1 for a constant thickness laminate; for variable
%   thickness laminates, specify the thickness of each ply, where 
%   n = length(STACKING_SEQUENCE).
%
%   LAYUP(3) is a flag to make the calculated section symmetric,
%   SYMMETRIC_LAYUP. When SYMMETRIC_LAYUP = true, the layup definition is
%   mirrored after the last ply in the stacking sequence (the positive-z
%   side).
%
%   Note: When the symmetric layup definition is specified, material, layup
%   stacking and ply thickness definitions are automatically reflected on
%   the other side of the symmetry plane.
%
%   LAYUP(4) is an integer specifying the number of stress/strain section
%   points per ply, SECTION_POINTS. Since the layup section is integrated
%   once before the stress analysis, section points are treated as sample
%   points.
%
%   OUTPUT_DEF. A 1x5 cell array specifying the ply output location,
%   MATLAB figures, strength calculation, stacking sequence optimisation
%   and the results output location.
%
%   OUTPUT_DEF(1) is the section output request, OUTPUT_PLY. When
%   OUTPUT_PLY is a string, it specifies the output location of each ply:
%
%     DEFAULT: Top and bottom faces
%     TOP: Top faces only
%     MIDDLE: Midspans only
%     BOTTOM: Bottom faces only
%     ALL: All section points in all plies
%     ENVELOPEABSMAX: The largest (+ve/-ve) value for the layup
%     ENVELOPEMAX: The largest (+ve) value for the layup
%     ENVELOPEMIN: The largest (-ve) value for the layup
%
%   When OUTPUT_PLY is a 1xn array, it specifies a user-defined section
%   point list.
%
%   OUTPUT_DEF(2) is a 1x2 cell array specifying settings for MATLAB figure
%   output, OUTPUT_FIGURE. OUTPUT_FIGURE(1) is a parameter specifying the
%   figure type ([], 'DEFAULT' or 'SMOOTH'). The parameter 'DEFAULT'
%   creates figures from the raw stress/strain data; the parameter 'SMOOTH'
%   uses the built-in SMOOTHDATA function to smooth the output at the ply
%   boundaries; an empty assignment disables MATLAB figure output;
%   OUTPUT_FIGURE(2) is a parameter specifying the figure layout ('SPLIT'
%   or 'COMPACT'). The parameter 'SPLIT' creates a separate plot for each
%   tensor component; the parameter 'COMPACT' overlays each tensor
%   component in a single plot.
%
%   OUTPUT_DEF(3) is a 1x2 cell array specifying settings for the strength
%   assessment, OUTPUT_STRENGTH. OUTPUT_STRENGTH(1) is a flag to enable or
%   disable the strength assessment; OUTPUT_STRENGTH(2) is the failure
%   assessment parameter:
%   
%     RESERVE: Strength reserve factor, R. For stress-based and
%     strain-based failure criteria, the inverse of the strength reserve
%     factor [1/R] is the scaling factor by which the load matrix must be
%     multiplied to hit the failure surface.
%     VALUE: Criterion value, V. The value of the index obtained directly
%     from the failure criterion.
%   
%   The strength calculation requires strength properties defined by
%   FAIL_STRESS, FAIL_STRAIN and HASHIN.
%
%   Note: The Tsai-Hill, Tsai-Wu and Azzi-Tsai-Hill failure criteria can be
%   expressed in terms of the strength reserve factor or the criterion
%   value; for all other failure criteria, the criterion value is identical
%   to the strength reserve factor.
%
%   Note: For Hashin's theory, R is not evaluated; output for these
%   criteria is quoted as the damage initiation criterion index.
%
%   OUTPUT_DEF(4) is a 1x4 cell array specifying settings for the stacking
%   sequence optimiser, OUTPUT_OPTIMISED. OUTPUT_OPTIMISED(1) is the
%   failure criterion for the optimisation ('MSTRS', 'TSAIH', 'TSAIW',
%   'AZZIT', 'MSTRN' or 'HASHIN'); OUTPUT_OPTIMISED(2) is the failure
%   assessment parameter ('RESERVE' or 'VALUE'); OUTPUT_OPTIMISED(3) is the
%   objective function ('MINMAX' or 'MINMEAN'); OUTPUT_OPTIMISED(4) is the
%   angular step size for the stacking sequence permutations.
%
%   OUTPUT_DEF(5) is a string specifying the results location,
%   OUTPUT_LOCATION. Use 'DEFAULT' to save results under a new folder in
%   the current working directory, or specify the directory directly.
%
%   Specify the load matrix:
%
%   [..] = ABD.MAIN(.., LOAD).
%
%   LOAD. A 1x6 array specifying the applied load, N11, N22, N12, M11, M22
%   and M12.
%
%   Include thermal/hydroscopic loads:
%
%   [..] = ABD.MAIN(.., THERM_HYDRO).
%
%   THERM_HYDRO. A 1x2 array specifying the thermal and hydroscopic load,
%   DELTA_T and DELTA_M, respectively.
%
%   Optional output arguments:
%
%   [ABD, ABD_INV, EI, EP, SP, EMTB] = ABD.MAIN(..).
%
%   ABD. A 6x6 matrix of the computed ABD matrix.
%
%   ABD_INV. A 6x6 matrix of the inverse ABD matrix.
%
%   EI. A 6x1 array of the midplane strains and curvatures induced in the
%   laminate. These strains represent the deflections of the laminate
%   about the neutral axis
%
%   E_PLY. A 1x6 cell array of the ply strains for all section points,
%   where n (below) is the total number of section points in the layup.
%
%   E_PLY(1) is a 3xn array of the ply strains in X-Y coordinates.
%
%   E_PLY(2) is a 3xn array of the ply strains in ply coordinates.
%
%   E_PLY(3) is a 3xn array of the stress-free ply strains due to thermal
%   process in X-Y coordinates.
%
%   E_PLY(4) is a 3xn array of the stress-free ply strains due to thermal
%   process in ply coordinates.
%
%   E_PLY(5) is a 3xn array of the stress-free ply strains due to moisture
%   process in X-Y coordinates.
%
%   E_PLY(6) is a 3xn array of the stress-free ply strains due to moisture
%   process in ply coordinates.
%
%   Note: For stress-free thermal/moisture strains, contractions have
%   positive values.
%
%   S_PLY. A 1x2 cell array of the ply stresses for all section points,
%   where n (below) is the total number of section points in the layup.
%
%   S_PLY(1) is a 3xn array of the ply stresses in X-Y coordinates.
%
%   S_PLY(2) is a 3xn array of the ply stresses in ply coordinates.
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
%   CFAILURE. A structure of the failure/damage initiation measure
%   components for all section points.
%
%   Failure/damage initiation analysis output variable identifiers:
%       - MSTRS, Maximum stress theory failure measure
%       - TSAIH, Tsai-Hill theory failure measure (reserve/value)
%       - TSAIW, Tsai-Wu theory failure measure (reserve/value)
%       - AZZIT, Azzi-Tsai-Hill theory failure measure (reserve/value)
%       - MSTRN, Maximum strain theory failure measure
%       - HSNFTCRT, Hashin’s fibre tensile damage initiation criterion
%       - HSNFCCRT, Hashin’s fibre compression damage initiation criterion
%       - HSNMTCRT, Hashin’s matrix tensile damage initiation criterion
%       - HSNMCCRT, Hashin’s matrix compression damage initiation criterion
%       - SFAILRATIO, The section failure ratio across all plies [%/100]
%
%   OPT_SEQ. A 1x6 cell of the results of the stacking sequence
%   optimisation.
%
%   OPT_SEQ(1) is a 1xn array of the optimum stacking sequence, where n is
%   the number of plies in the layup.
%
%   OPT_SEQ(2) is the critical failure criterion value.
%
%   OPT_SEQ(3) is the total number of stacking permutations considered by
%   the optimiser.
%
%   OPT_SEQ(4) is the analysis time (in seconds).
%
%   OPT_SEQ(5) is an exception returned in case of optimisation failure.
%
%   OPT_SEQ(6) is a structure of the stress and strain tensors
%   corresponding to the optimum stacking sequence.
%
%   Note: Replace unrequested outputs with an empty ( [] ) assignment.
%
%   MATLAB figures:
%
%   EP. MATLAB figure of E_PLY in X-Y coordinates and ply coordinates
%   for all section points.
%
%   SP. MATLAB figure of S_PLY in X-Y coordinates and ply coordinates
%   for all section points.
%
%   CB, Optimiser criterion for all stacking permutations.
%
%   See also abd.user_definitions.
%
%==========================================================================
%
%   With contributions from: Matt Fig.
%
%   Permission to user the above author's work is granted under the BSD and
%   CC by-nc-sa 4.0 licenses, where applicable. Third-party source code is
%   clearly indicated in its own subfolder.
%
%   Layup Analysis Tool 2.7.1 Copyright Louis Vallance 2024
%   Last modified 09-Feb-2024 09:10:19 UTC

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%
clc

%% INITIALISE VARARGOUT
varargout{1.0} = [];
varargout{2.0} = [];
varargout{3.0} = [];
varargout{4.0} = [];
varargout{5.0} = [];
varargout{6.0} = [];
varargout{7.0} = [];
varargout{8.0} = [];

%% GET USER INPUTS FROM VARARGIN
[enableTensor, printTensor, materialDataMechanical, materialDataFailStress, materialDataFailStrain, materialDataHashin, theta, t_ply, symmetricPly, SECTION_POINTS, OUTPUT_PLY,...
    OUTPUT_FIGURE, OUTPUT_STRENGTH, OUTPUT_OPTIMISED, OUTPUT_LOCATION, Nxx, Nyy, Nxy, Mxx, Myy, Mxy, deltaT, deltaM, error] =...
    ...
    abd.internal_initialise(nargin, varargin);

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
    abd.internal_getMaterial(materialDataMechanical, nPlies, symmetricPly, 1.0, 'MECHANICAL');

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

    if (noFailStress == true) && (noFailStrain == true) && (noHashin == true)
        %{
            Since the strength calculation has been requested, at least one
            of FAIL_STRESS, FAIL_STRAIN or HASHIN must be defined for the
            layup!
        %}
        fprintf('[ERROR] The strength calculation requires at least\nFAIL_STRESS, FAIL_STRAIN or HASHIN material properties\n');
        return
    end
else
    noFailStress = true;
    noFailStrain = true;
    noHashin = true;
end

%% SET THE NEAR-ZERO TOLERANCE VALUE
tolerance = 1e-6;

%% GET THICKNESS FRACTIONS
[z, t] = abd.internal_getThickness(nPlies, t_ply, tolerance);

%% PROCESS SECTION_POINTS
[error, z_points, theta_points, nPlies_points, A11_points, A22_points,...
    B11_points, B22_points, plyBuffer, thickness] =...
    ...
    abd.internal_getSectionPoints(SECTION_POINTS, 'SECTION_POINTS', nPlies, theta, z, A11, A22, B11, B22, tolerance);

% An error occurred, so RETURN
if error == true
    return
end

%% INITIALISE VARIABLES
BEST_SEQUENCE = [];
axx = zeros(1.0, nPlies_points);
ayy = axx;
axy = axx;
bxx = axx;
byy = axx;
bxy = axx;

%% PROCESS OUTPUT_PLY
[error, OUTPUT_PLY_POINTS, plyBuffer, OUTPUT_ENVELOPE, ENVELOPE_MODE, outputApproximate, plyBuffer_sfailratio] =...
    ...
    abd.internal_getOutputPoints(OUTPUT_PLY, z, z_points, nPlies,nPlies_points, plyBuffer, SECTION_POINTS, tolerance);

% An error occurred, so RETURN
if error == true
    return
end

%% GET OPTIMISER SETTINGS
if isempty(OUTPUT_OPTIMISED{1.0}) == false
    [error, OUTPUT_OPTIMISED] =...
        ...
        abd.internal_optimise.getSettings(OUTPUT_OPTIMISED, noFailStress, noFailStrain, noHashin, OUTPUT_STRENGTH{1.0});

    % An error occurred, so RETURN
    if error == true
        return
    end
end

%% COMPUTE REDUCED STIFFNESS TERMS
[Q11, Q22, Q66, Q12] =...
    ...
    abd.internal_getReducedQ(E11, E22, V12, G12);

%% COMPUTE TRANSFORMED REDUCED STIFFNESS MATRIX COMPONENTS
[Q11t, Q12t, Q16t, Q22t, Q26t, Q66t] =...
    ...
    abd.internal_getTransformedQ(theta, Q11, Q12, Q66, Q22);

%% GET EFFECTIVE THEMAL AND MOISTURE EXPANSION COEFFICIENTS FOR EACH PLY
if nargin == 5.0
    [axx, ayy, axy, bxx, byy, bxy] =...
        ...
        abd.internal_getThermoHydro(theta_points, A11_points, A22_points, B11_points, B22_points);
end

%% COMPUTE A, B and D MATRICES
[ABD, ABD_INV, Qijt, NxxT, NyyT, NxyT, MxxT, MyyT, MxyT, NxxM, NyyM, NxyM, MxxM, MyyM, MxyM] =...
    ...
    abd.internal_getABD(nPlies, Q11t, Q12t, Q16t, Q22t, Q26t, Q66t, z, nargin, deltaT, deltaM, axx, ayy, axy, bxx, byy, bxy, SECTION_POINTS);

%% COMPUTE TENSOR QUANTITIES
if enableTensor == true
    [E_midplane, E_ply_xy, S_ply_xy, E_ply_aligned, S_ply_aligned, E_therm_xy, E_moist_xy, E_therm_aligned, E_moist_aligned] =...
        ...
        abd.internal_getTensor(ABD, Nxx, NxxT, NxxM, Nyy, NyyT, NyyM, Nxy, NxyT, NxyM, Mxx, MxxT, MxxM, Myy, MyyT, MyyM, Mxy, MxyT, MxyM, nPlies_points, z_points, theta_points,...
        Qijt, deltaT, deltaM, axx, ayy, axy, bxx, byy, bxy, tolerance);

    if any(E_midplane) == false
        printTensor = -1.0;
    end
else
    % Initialise values to default
    printTensor = 0.0;
    E_midplane = [];
    E_ply_xy = [];
    S_ply_xy = [];
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
if (OUTPUT_STRENGTH{1.0} == true) && (printTensor == 1.0)
    [MSTRS, TSAIH, TSAIW, AZZIT, MSTRN, HSNFTCRT, HSNFCCRT, HSNMTCRT, HSNMCCRT, XT, XC, YT, YC, S, C, B, E11, E22, G12, V12, XET, XEC, YET, YEC, SE, ALPHA, XHT, XHC, YHT, YHC,...
        SHX, SHY] =...
        ...
        abd.internal_strength.main(noFailStress, noFailStrain, noHashin, XT, XC, YT, YC, S, C, B, E11, E22, G12, V12, XET, XEC, YET, YEC, SE, ALPHA, XHT, XHC, YHT, YHC, SHX, SHY,...
        S_ply_aligned, nPlies, nPlies_points, SECTION_POINTS, OUTPUT_STRENGTH{2.0});

    if OUTPUT_OPTIMISED{1.0} == true
        %% FIND THE OPTIMUM STACKING SEQUENCE
        [BEST_SEQUENCE, CRITERION_BUFFER, MIN_CRITERION] =...
            ...
            abd.internal_optimise.main(OUTPUT_OPTIMISED, nargin, nPlies, nPlies_points, SECTION_POINTS, z, z_points, Q11, Q22, Q66, Q12, A11_points, A22_points, B11_points,...
            B22_points, tolerance, XT, XC, YT, YC, S, C, B, XET, XEC, YET, YEC, SE, ALPHA, XHT, XHC, YHT, YHC, SHX, SHY, deltaT, deltaM, Nxx, Nyy, Nxy, Mxx, Myy, Mxy, E11, E22,...
            V12, G12);
    else
        CRITERION_BUFFER = [];
    end
else
    % Initialise values to default
    MSTRS = [];
    TSAIH = [];
    TSAIW = [];
    AZZIT = [];
    MSTRN = [];
    HSNFTCRT = [];
    HSNFCCRT = [];
    HSNMTCRT = [];
    HSNMCCRT = [];
    CRITERION_BUFFER = [];

    % Suppress strength output
    OUTPUT_STRENGTH{1.0} = false;
end

%% OUTPUT TO VARARGOUT
varargout{1.0} = ABD;
varargout{2.0} = inv(ABD);
varargout{3.0} = E_midplane;
varargout{4.0} = {E_ply_xy, E_ply_aligned, E_therm_xy, E_therm_aligned, E_moist_xy, E_moist_aligned};
varargout{5.0} = {S_ply_xy, S_ply_aligned};
varargout{6.0} = {[EXT, EYT, GXYT, NUXYT, NUYXT],...
                  [EXB, EYB, GXYB, NUXYB, NUYXB]};
varargout{7.0} = struct('MSTRS', MSTRS, 'TSAIH', TSAIH, 'TSAIW', TSAIW, 'AZZIT', AZZIT, 'MSTRN', MSTRN, 'HSNFTCRT', HSNFTCRT, 'HSNFCCRT', HSNFCCRT, 'HSNMTCRT', HSNMTCRT,...
    'HSNMCCRT', HSNMCCRT);
varargout{8.0} = BEST_SEQUENCE;

%% CREATE OUTPUT DIRECTORY
% Create the root output folder if it does not already exist
if exist(OUTPUT_LOCATION, 'dir') ~= 7.0
    mkdir(OUTPUT_LOCATION)
end

% Add the root output folder to the MATLAB path
addpath(OUTPUT_LOCATION)

% Get the date string for the output folder
dateString = char(datetime('now'));
for i = 1:length(dateString)
    if (strcmpi(dateString(i), ':') == 1.0) || (strcmpi(dateString(i), ' ') == 1.0)
        dateString(i) = '_';
    end
end

% Construct the output location path
outputLocation = [OUTPUT_LOCATION, [filesep, 'abd_results_', dateString]];

% Create the folder if it does not already exist
if exist(outputLocation, 'dir') ~= 7.0
    mkdir(outputLocation)
end

%% PLOT STRAINS AND STRESSES IN A MATLAB FIGURE
if (isempty(OUTPUT_FIGURE{1.0}) == false) && (printTensor == 1.0) && (nPlies_points > 1.0)
    abd.internal_plot.main(OUTPUT_FIGURE{1.0}, OUTPUT_FIGURE{2.0}, outputLocation, nPlies, E_ply_xy, S_ply_xy, E_ply_aligned, S_ply_aligned, z, z_points, CRITERION_BUFFER,...
        OUTPUT_OPTIMISED)
end

%% WRITE RESULTS TO A TEXT FILE
abd.internal_outputToFile(dateString, outputLocation, OUTPUT_STRENGTH, nPlies, t_ply, theta, enableTensor, printTensor, S_ply_aligned, S_ply_xy, E_ply_aligned, E_ply_xy,...
    E_therm_xy, E_moist_xy, E_therm_aligned, E_moist_aligned, ABD, symmetricAbd, EXT, EYT, GXYT, NUXYT, NUYXT, EXB, EYB, GXYB, NUXYB, NUYXB, MSTRS, TSAIH, TSAIW, AZZIT, MSTRN,...
    HSNFTCRT, HSNFCCRT, HSNMTCRT, HSNMCCRT, noFailStress, noFailStrain, noHashin, SECTION_POINTS, OUTPUT_PLY_POINTS, plyBuffer, thickness, OUTPUT_ENVELOPE, ENVELOPE_MODE,...
    outputApproximate, BEST_SEQUENCE, OUTPUT_OPTIMISED, OUTPUT_FIGURE{1.0}, plyBuffer_sfailratio, axx, ayy, axy, bxx, byy, bxy, E_midplane)

%% Add the output location to the MATLAB path
addpath(genpath(outputLocation));
end