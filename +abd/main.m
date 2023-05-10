function [varargout] = main(varargin)
%ABD.MAIN    Analyse a user-defined composite layup.
%   [VARARGOUT] = ABD.MAIN(VARARGIN) computes the ABD matrix based on an
%   N-layer composite layup and plane stress material properties. The
%   following output is produced:
%
%   - A, B and D matrices (and their inverses)
%   - Induced midplane strains and curvatures
%   - X-Y stresses and strains based on specified forces and moments
%   - Ply stresses and strains based on specified forces and moments
%   - Equivalent extensional and bending moduli (symmetric layups only)
%   - Stress and strain-based failure criteria
%
%   Notes:
%   - Support for unidirectional fibre-reinforced composites only
%   - Strains are assumed to vary linearly through the layup thickness
%   (z-direction)
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
%   properties E11, E22, G12, V12, A11, A22, B11 and B22 for each ply.
%
%   MATERIAL(2) is a 1xn cell array specifying the strength properties for
%   stress-based failure criteria XT, XC, YT, YC, S, C and B for each ply.
%
%   MATERIAL(3) is a 1xn cell array specifying the strength properties for
%   strain-based failure criteria XET, XEC, YET, YEC and SE for each ply.
%
%   MATERIAL(4) is a 1xn cell array specifying the strength properties for
%   the Hashin damage initiation criteria ALPHA, XHT, XHC, YHT, YHC, SHX
%   and SHY for each ply.
%
%   Note: n = 1 for constant material properties; for ply-wise material
%   properties, specify n sets of properties corresponding to the n ply
%   definitions given by LAYUP.
%
%   LAYUP. A 1x4 cell specifying the layup orientation, the laminate
%   thickness, stacking symmetry and the number of section points.
%
%   LAYUP(1) is a 1xn array defining the layup stacking sequence,
%   STACKING_SEQUENCE.
%
%   LAYUP(2) is a 1xn array containing the ply thickness values,
%   PLY_THICKNESS: n = 1 for a constant thickness laminate; for variables
%   thickness laminates, specify the thickness of each ply such that
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
%   LAYUP(4) is an integer specifying the number of section points (per
%   ply) for stress/strain computation, SECTION_POINTS.
%
%   OUTPUT_DEF. A 1x5 cell array specifying the ply output location, MATLAB
%   figures, strength calculation, stacking sequence optimisation and the
%   results output location.
%
%   OUTPUT_DEF(1) is the section output request, OUTPUT_PLY. When
%   OUTPUT_PLY is a string, it specifies the output location of each ply:
%
%     DEFAULT: Top and bottom faces
%     TOP: Top faces only
%     MIDDLE: Midspan only
%     BOTTOM: Bottom faces only
%     ALL: All section points in all plies
%     ENVELOPEABSMAX: The largest (+ve/-ve) value for the layup
%     ENVELOPEMAX: The largest (+ve) value for the layup
%     ENVELOPEMIN: The largest (-ve) value for the layup
%
%   When OUTPUT_PLY is a 1xn array, it specifies a user-defined section
%   point list.
%
%   OUTPUT_DEF(2) is a flag to request MATLAB figures of stress/strain
%   output for the layup, OUTPUT_FIGURE. The parameter 'DEFAULT' creates
%   figures from the raw stress/strain data; the parameter 'SMOOTH' uses
%   the built-in SMOOTHDATA function to smooth the output at the ply
%   boundaries.
%
%   OUTPUT_DEF(3) is a flag to request the strength calculation for each
%   ply, OUTPUT_STRENGTH. The strength calculation requires strength
%   properties defined by FAIL_STRESS and FAIL_STRAIN.
%
%   OUTPUT_DEF(4) is a 1x4 cell array specifying settings for the stacking
%   sequence optimiser. OUTPUT_OPTIMISED(1) is the failure criterion for
%   the optimisation ('MSTRS', 'TSAIH', 'TSAIW', 'AZZIT', 'MSTRN' or
%   'HASHIN'); OUTPUT_OPTIMISED(2) is the failure assessment parameter
%   ('RESERVE' or 'VALUE'); OUTPUT_OPTIMISED(3) is the objective function
%   ('MINMAX' or 'MINMEAN'); OUTPUT_OPTIMISED(4) is the angular step size
%   for the stacking sequence permutations.
%
%   OUTPUT_DEF(5) is a string specifying the results location,
%   OUTPUT_LOCATION. Use 'DEFAULT' to save results under a new folder in
%   the current working directory, or specify the directory directly.
%
%   Specify the load matrix:
%
%   [..] = ABD.MAIN(.., LOAD).
%
%   LOAD. A 1x6 array specifying the applied load N11, N22, N12, M11, M22
%   and M12.
%
%   Include thermal/hydroscopic loads:
%
%   [..] = ABD.MAIN(.., THERM_HYDRO).
%
%   THERM_HYDRO. A 1x2 array specifying the thermal and hydroscopic load
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
%   E_PLY(2) is a 3xn array of the ply strains in the ply directions.
%
%   E_PLY(3) is a 3xn array of the stress-free ply strains due to thermal
%   process in X-Y coordinates.
%
%   E_PLY(4) is a 3xn array of the stress-free ply strains due to thermal
%   process in the ply directions.
%
%   E_PLY(5) is a 3xn array of the stress-free ply strains due to moisture
%   process in X-Y coordinates.
%
%   E_PLY(6) is a 3xn array of the stress-free ply strains due to moisture
%   process in the ply directions.
%
%   Note: For stress-free thermal/moisture strains, contractions have
%   positive values.
%
%   S_PLY. A 1x2 cell array of the ply stresses for all section points,
%   where n (below) is the total number of section points in the layup.
%
%   S_PLY(1) is a 3xn array of the ply stresses in X-Y coordinates.
%
%   S_PLY(2) is a 3xn array of the ply stresses in the ply directions.
%
%   EQ_MODULI. A 1x2 cell array of the equivalent moduli in tension and
%   bending.
%
%   EQ_MODULI(1) = EXT, EYT, GXYT, NUXYT, NUYXT.
%
%   EQ_MODULI(2) = EXB, EYB, GXYB, NUXYB, NUYXB.
%
%   Note: The equivalent moduli are only calculated for symmetric laminate
%   stacking sequences.
%
%   CFAILURE. A structure of the failure/damage initiation measure
%   components.
%
%   Failure/damage initiation analysis output variable identifiers:
%       - MSTRS, Maximum stress theory failure measure
%       - TSAIH, Tsai-Hill theory failure measure
%       - TSAIW, Tsai-Wu theory failure measure
%       - AZZIT, Azzi-Tsai-Hill theory failure measure
%       - MSTRN, Maximum strain theory failure measure
%       - HSNFTCRT, Hashin’s fibre tensile damage initiation criterion
%       - HSNFCCRT, Hashin’s fibre compression damage initiation criterion
%       - HSNMTCRT, Hashin’s matrix tensile damage initiation criterion
%       - HSNMCCRT, Hashin’s matrix compression damage initiation criterion
%       - SFAILRATIO, The section failure ratio across all plies [%/100]
%
%   Note: Failure criteria Tsai-Hill, Tsai-Wu and Azzi-Tsai-Hill are
%   expressed as the strength reserve factor.
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
%   EP. MATLAB figure of E_PLY in X-Y coordinates and the ply directions
%   for all section points.
%
%   SP. MATLAB figure of S_PLY in X-Y coordinates and the ply directions
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
%   Layup Analysis Tool 2.4 Copyright Louis Vallance 2023
%   Last modified 10-May-2023 10:16:13 UTC

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
[enableTensor, printTensor, materialDataMechanical,...
    materialDataFailStress, materialDataFailStrain, materialDataHashin,...
    theta, t_ply, symmetricPly, SECTION_POINTS, OUTPUT_PLY,...
    OUTPUT_FIGURE, OUTPUT_STRENGTH, OUTPUT_OPTIMISED, OUTPUT_LOCATION,...
    Nxx, Nyy, Nxy, Mxx, Myy, Mxy, deltaT, deltaM] =...
    ...
    abd.internal_initialise(nargin, varargin);

%% MIRROR THE LAYUP DEFINITION (IF APPLICABLE)
[t_ply, theta, nPlies, error] =...
    ...
    abd.internal_mirror(symmetricPly, t_ply, theta);

% An error occurred, so RETURN
if error == true
    return
end

%% GET MATERIAL DATA (MECHANICAL)
[error, ~, E11, E22, G12, V12, A11, A22, B11, B22] =...
    ...
    abd.internal_getMaterial(materialDataMechanical, nPlies,...
    symmetricPly, 1.0, 'MECHANICAL');

% An error occurred, so RETURN
if error == true
    return
end

%% GET MATERIAL DATA (STRENGTH)
if OUTPUT_STRENGTH == true
    % Get fail stress properties
    [error, noFailStress, XT, XC, YT, YC, S, C, B] =...
        ...
        abd.internal_getMaterial(materialDataFailStress, nPlies,...
        symmetricPly, 2.0, 'FAIL_STRESS');

    % An error occurred, so RETURN
    if error == true
        return
    end

    % Get fail strain properties
    [error, noFailStrain, XET, XEC, YET, YEC, SE] =...
        ...
        abd.internal_getMaterial(materialDataFailStrain, nPlies,...
        symmetricPly, 3.0, 'FAIL_STRAIN');

    % An error occurred, so RETURN
    if error == true
        return
    end

    % Get Hashin properties
    [error, noHashin, ALPHA, XHT, XHC, YHT, YHC, SHX, SHY] =...
        ...
        abd.internal_getMaterial(materialDataHashin, nPlies,...
        symmetricPly, 2.0, 'HASHIN');

    % An error occurred, so RETURN
    if error == true
        return
    end

    if (noFailStress == true) && (noFailStrain == true) && (noHashin == true)
        %{
            Since the strength calculation has been requested, at least one
            of FAIL_STRESS, FAIL_STRAIN or HASHIN must be defined for the
            layup!
        %}
        fprintf(['[ABD ERROR] The strength calculation requies at leas',...
            't FAIL_STRESS, FAIL_STRAIN or HASHIN material properties\n']);
        return
    end
else
    noFailStress = true;
    noFailStrain = true;
    noHashin = true;
end

%% GET THICKNESS FRACTIONS
% Buffer for thickness fractions
z = zeros(1.0, nPlies + 1.0);

% Check for invalid thickness values
if any(t_ply <= 0.0) == true
    fprintf(['[ABD ERROR] Zero or negative ply thickness values are no',...
        't allowed\n']);
    return
end

% Get the total ply thickness
t = sum(t_ply);

% Compute thickness fraction values
z(1.0) = -t/2.0;
for i = 2.0:nPlies + 1.0
    z(i) = z(i - 1.0) + t_ply(i-1);
end

%% SET THE NEAR-ZERO TOLERANCE VALUE
tolerance = 1e-6;

%% PROCESS SECTION_POINTS
[error, z_points, theta_points, nPlies_points, A11_points, A22_points,...
    B11_points, B22_points, plyBuffer, thickness] =...
    ...
    abd.internal_getSectionPoints(SECTION_POINTS, 'SECTION_POINTS',...
    nPlies, theta, z, A11, A22, B11, B22, tolerance);

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
[error, OUTPUT_PLY_POINTS, plyBuffer, OUTPUT_ENVELOPE, ENVELOPE_MODE,...
    outputApproximate, plyBuffer_sfailratio] =...
    ...
    abd.internal_getOutputPoints(OUTPUT_PLY, z, z_points, nPlies,...
    nPlies_points, plyBuffer, SECTION_POINTS, tolerance);

% An error occurred, so RETURN
if error == true
    return
end

%% GET OPTIMISER SETTINGS
if isempty(OUTPUT_OPTIMISED{1.0}) == false
    [error, OUTPUT_OPTIMISED] =...
        ...
        abd.internal_optimise.getSettings(OUTPUT_OPTIMISED, noFailStress,...
        noFailStrain, noHashin, OUTPUT_STRENGTH);

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
        abd.internal_getThermoHydro(theta_points, A11_points,...
        A22_points, B11_points, B22_points);
end

%% COMPUTE A, B and D MATRICES
[ABD, ABD_INV, Qijt, NxxT, NyyT, NxyT, MxxT, MyyT, MxyT, NxxM, NyyM,...
    NxyM, MxxM, MyyM, MxyM] =...
    ...
    abd.internal_getABD(nPlies, Q11t, Q12t, Q16t, Q22t, Q26t, Q66t, z,...
    nargin, deltaT, deltaM, axx, ayy, axy, bxx, byy, bxy, SECTION_POINTS);

%% COMPUTE TENSOR QUANTITIES
if enableTensor == true
    [E_midplane, E_ply_xy, S_ply_xy, E_ply_aligned, S_ply_aligned,...
        E_therm_xy, E_moist_xy, E_therm_aligned, E_moist_aligned] =...
        ...
        abd.internal_getTensor(ABD, Nxx, NxxT, NxxM, Nyy, NyyT, NyyM,...
        Nxy, NxyT, NxyM, Mxx, MxxT, MxxM, Myy, MyyT, MyyM, Mxy, MxyT,...
        MxyM, nPlies_points, z_points, theta_points, Qijt, deltaT,...
        deltaM, axx, ayy, axy, bxx, byy, bxy, tolerance);

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
if (OUTPUT_STRENGTH == true) && (printTensor == 1.0)
    [MSTRS, TSAIH, TSAIW, AZZIT, MSTRN, HSNFTCRT, HSNFCCRT, HSNMTCRT,...
        HSNMCCRT, XT, XC, YT, YC, S, C, B, E11, E22, G12, V12, XET, XEC,...
        YET, YEC, SE, ALPHA, XHT, XHC, YHT, YHC, SHX, SHY] =...
        ...
        abd.internal_strength.main(noFailStress, noFailStrain, noHashin,...
        XT, XC, YT, YC, S, C, B, E11, E22, G12, V12, XET, XEC, YET, YEC,...
        SE, ALPHA, XHT, XHC, YHT, YHC, SHX, SHY, S_ply_aligned, nPlies,...
        nPlies_points, SECTION_POINTS);

    if OUTPUT_OPTIMISED{1.0} == true
        %% FIND THE OPTIMUM STACKING SEQUENCE
        [BEST_SEQUENCE, CRITERION_BUFFER, MIN_CRITERION] =...
            ...
            abd.internal_optimise.main(OUTPUT_OPTIMISED, nargin, nPlies,...
            nPlies_points, SECTION_POINTS, z, z_points, Q11, Q22, Q66,...
            Q12, A11_points, A22_points, B11_points, B22_points,...
            tolerance, XT, XC, YT, YC, S, C, B, XET, XEC, YET, YEC, SE,...
            ALPHA, XHT, XHC, YHT, YHC, SHX, SHY, deltaT, deltaM, Nxx,...
            Nyy, Nxy, Mxx, Myy, Mxy, E11, E22, V12,...
            G12);
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
    OUTPUT_STRENGTH = false;
end

%% OUTPUT TO VARARGOUT
varargout{1.0} = ABD;
varargout{2.0} = inv(ABD);
varargout{3.0} = E_midplane;
varargout{4.0} = {E_ply_xy, E_ply_aligned, E_therm_xy, E_therm_aligned,...
    E_moist_xy, E_moist_aligned};
varargout{5.0} = {S_ply_xy, S_ply_aligned};
varargout{6.0} = {[EXT, EYT, GXYT, NUXYT, NUYXT],...
                  [EXB, EYB, GXYB, NUXYB, NUYXB]};
varargout{7.0} = struct('MSTRS', MSTRS, 'TSAIH', TSAIH, 'TSAIW', TSAIW,...
    'AZZIT', AZZIT, 'MSTRN', MSTRN, 'HSNFTCRT', HSNFTCRT, 'HSNFCCRT',...
    HSNFCCRT, 'HSNMTCRT', HSNMTCRT, 'HSNMCCRT', HSNMCCRT);
varargout{8.0} = BEST_SEQUENCE;

%% CREATE OUTPUT DIRECTORY
% Get the date string for the output folder
dateString = char(datetime('now'));
for i = 1:length(dateString)
    if (strcmpi(dateString(i), ':') == 1.0) ||...
            (strcmpi(dateString(i), ' ') == 1.0)
        dateString(i) = '_';
    end
end

outputLocation = [OUTPUT_LOCATION, ['\abd_results_', dateString]];
if exist(outputLocation, 'dir') == false
    mkdir(outputLocation)
end

%% PLOT STRAINS AND STRESSES IN A MATLAB FIGURE
if (isempty(OUTPUT_FIGURE) == false) && (printTensor == 1.0) && (nPlies_points > 1.0)
    abd.internal_plot(OUTPUT_FIGURE, outputLocation, nPlies, E_ply_xy,...
        S_ply_xy, E_ply_aligned, S_ply_aligned, z, z_points,...
        CRITERION_BUFFER, OUTPUT_OPTIMISED)
end

%% WRITE RESULTS TO A TEXT FILE
abd.internal_outputToFile(dateString, outputLocation, OUTPUT_STRENGTH,...
    nPlies, t_ply, theta, enableTensor, printTensor, S_ply_aligned,...
    S_ply_xy, E_ply_aligned, E_ply_xy, E_therm_xy, E_moist_xy,...
    E_therm_aligned, E_moist_aligned, ABD, symmetricAbd, EXT, EYT,...
    GXYT, NUXYT, NUYXT, EXB, EYB, GXYB, NUXYB, NUYXB, MSTRS, TSAIH,...
    TSAIW, AZZIT, MSTRN, HSNFTCRT, HSNFCCRT, HSNMTCRT, HSNMCCRT,...
    noFailStress, noFailStrain, noHashin, SECTION_POINTS,...
    OUTPUT_PLY_POINTS, plyBuffer, thickness, OUTPUT_ENVELOPE,...
    ENVELOPE_MODE, outputApproximate, BEST_SEQUENCE, OUTPUT_OPTIMISED,...
    OUTPUT_FIGURE, plyBuffer_sfailratio)

%% Add the output location to the MATLAB path
addpath(genpath(outputLocation));
end