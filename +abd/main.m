function [S] = main(settings)
%ABD.MAIN    Analyse a user-defined composite layup.
%   [<output-structure>] = ABD.MAIN(<input-structure>) computes the ABD
%   matrix from an n-layer composite layup, performs a stress analysis
%   based on a user-defined load, and evaluates the layup strength from
%   stress and strain-based failure and damage initiation criteria.
%
%   For the full documentation, type the following command:
%   >> help lat
%
%==========================================================================
%
%   With contributions from: Matt Fig and Jan Simon.
%
%   Permission to user the above author's work is granted under the BSD and
%   CC by-nc-sa 4.0 licenses, where applicable. Third-party source code is
%   clearly indicated in its own subfolder.
%
%   Layup Analysis Tool 5.1.5 Copyright Louis Vallance 2026
%   Last modified 04-Sep-2026 12:21:27 UTC

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
% Flag indicating if strength output is requested
isStrengthOutput = (isa(OUTPUT_STRENGTH{1.0}, 'function_handle') == true) || ((islogical(OUTPUT_STRENGTH{1.0}) == true) && (OUTPUT_STRENGTH{1.0} == true));

if isStrengthOutput == true
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
        fprintf('[ERROR] The strength calculation requires material properties for at least\none failure/damage initiation criterion: FAIL_STRESS, FAIL_STRAIN, HASHIN\nor LARC05\n');
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
[z, t, error] = abd.internal_getThickness(nPlies, t_ply, tolerance);

% An error occurred, so RETURN
if error == true
    return
end

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
        abd.internal_optimise.getSettings(OUTPUT_OPTIMISED, noFailStress, noFailStrain, noHashin, noLaRC05, isStrengthOutput, isa(OUTPUT_STRENGTH{1.0}, 'function_handle'));

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
[ABD, ABD_INV, Qijt, NxxT, NyyT, NxyT, MxxT, MyyT, MxyT, NxxM, NyyM, NxyM, MxxM, MyyM, MxyM, error] =...
    ...
    abd.internal_getABD(nPlies, Q11t, Q12t, Q16t, Q22t, Q26t, Q66t, z, nargin, deltaT, deltaM, axx, ayy, axy, bxx, byy, bxy, SECTION_POINTS);

% An error occurred, so RETURN
if error == true
    return
end

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
UCRT_MException = [];
noUcrt = true;

% Initialise failure criteria component buffers
[MSTRS, TSAIH, HOFFMAN, TSAIW, AZZIT, MSTRN, HSNFTCRT, HSNFCCRT, HSNMTCRT, HSNMCCRT, LARPFCRT, LARMFCRT, LARKFCRT, LARSFCRT, LARTFCRT, UCRT] = abd.internal_strength.init(nPlies_points);

if (isStrengthOutput == true) && (printTensor == 1.0)
    % Get the tensor data into a structure
    TENSORS = struct('E_MIDSPAN', E_midspan, 'E_PLY_XY', E_ply_xy, 'S_PLY_XY', S_ply_xy, 'E_PLY_ALIGNED', E_ply_aligned, 'S_PLY_ALIGNED', S_ply_aligned, 'E_THERM_XY', E_therm_xy,...
        'E_HYDRO_XY', E_hydro_xy, 'E_THERM_ALIGNED', E_therm_aligned, 'E_HYDRO_ALIGNED', E_hydro_aligned);

    [MSTRS, TSAIH, HOFFMAN, TSAIW, AZZIT, MSTRN, HSNFTCRT, HSNFCCRT, HSNMTCRT, HSNMCCRT, LARPFCRT, LARMFCRT, LARKFCRT, LARSFCRT, LARTFCRT, UCRT, XT, XC, YT, YC, S, C, B, E11, E22,...
        G12, V12, XET, XEC, YET, YEC, SE, ALPHA, XHT, XHC, YHT, YHC, SHX, SHY, XLT, XLC, YLT, YLC, SLX, SLY, GL12, NL, NT, A0, PHI0, S1, S2, S3, UCRT_MException] =...
        ...
        abd.internal_strength.main(noFailStress, noFailStrain, noHashin, noLaRC05, symsAvailable, XT, XC, YT, YC, S, C, B, E11, E22, G12, V12, axx, ayy, axy, bxx, byy, bxy, XET,...
        XEC, YET, YEC, SE, ALPHA, XHT, XHC, YHT, YHC, SHX, SHY, XLT, XLC, YLT, YLC, SLX, SLY, GL12, NL, NT, A0, PHI0, TENSORS, nPlies, nPlies_points, SECTION_POINTS,...
        OUTPUT_STRENGTH{1.0}, OUTPUT_STRENGTH{2.0}, OUTPUT_STRENGTH{3.0}, MSTRS, TSAIH, HOFFMAN, TSAIW, AZZIT, MSTRN, HSNFTCRT, HSNFCCRT, HSNMTCRT, HSNMCCRT, LARPFCRT, LARMFCRT,...
        LARKFCRT, LARSFCRT, LARTFCRT, UCRT, 1.0);

    % Update UCRT flag based on outcome of routine
    noUcrt = ((isa(OUTPUT_STRENGTH{1.0}, 'function_handle') == true) & (all(UCRT == -1.0) == false)) == false;

    if (isempty(OUTPUT_OPTIMISED{1.0}) == false) && (OUTPUT_OPTIMISED{1.0} == true) && (isempty(UCRT_MException) == true)
        % Flag indicating if Parallel Computing Toolbox is available
        pctAvail = abd.internal_checkToolbox('Parallel Computing Toolbox');

        %% FIND THE OPTIMUM STACKING SEQUENCE
        [BEST_SEQUENCE, CRITERION_BUFFER, ~, CHUNK_SIZE, N_CHUNKS, EXECUTION_MODE] =...
            ...
            abd.internal_optimise.main(OUTPUT_OPTIMISED, nargin, nPlies, nPlies_points, SECTION_POINTS, z, z_points, Q11, Q22, Q66, Q12, A11_points, A22_points, B11_points,...
            B22_points, tolerance, XT, XC, YT, YC, S, C, B, XET, XEC, YET, YEC, SE, ALPHA, XHT, XHC, YHT, YHC, SHX, SHY, XLT, XLC, YLT, YLC, SLX, SLY, GL12, NL, NT, A0, PHI0,...
            deltaT, deltaM, Nxx, Nyy, Nxy, Mxx, Myy, Mxy, E11, E22, V12, G12, symsAvailable, S1, S2, S3, SECTION_POINTS, OUTPUT_STRENGTH{1.0}, OPTIMISER_SETTINGS, pctAvail,...
            OUTPUT_STRENGTH{3.0});
    else
        CRITERION_BUFFER = [];
    end

    % Set failed section points colour buffer (for plotting)
    SP_COLOUR_BUFFER = abd.internal_strength.getFailedSpBuffer(MSTRS, TSAIH, HOFFMAN, TSAIW, AZZIT, MSTRN, HSNFTCRT, HSNFCCRT, HSNMTCRT, HSNMCCRT, LARPFCRT, LARMFCRT, LARKFCRT,...
        LARSFCRT, LARTFCRT, UCRT, SP_COLOUR_BUFFER, nPlies_points);
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
[SFAILRATIO_STRESS, SFAILRATIO_STRAIN, SFAILRATIO_HASHIN, SFAILRATIO_LARC05, SFAILRATIO_UCRT] =...
    abd.internal_outputToFile(dateString, outputLocation, OUTPUT_STRENGTH, nPlies, t_ply, theta, enableTensor, printTensor, S_ply_aligned, S_ply_xy, E_ply_aligned, E_ply_xy,...
    E_therm_xy, E_hydro_xy, E_therm_aligned, E_hydro_aligned, ABD, symmetricAbd, EXT, EYT, GXYT, NUXYT, NUYXT, EXB, EYB, GXYB, NUXYB, NUYXB, MSTRS, TSAIH, HOFFMAN, TSAIW, AZZIT,...
    MSTRN, HSNFTCRT, HSNFCCRT, HSNMTCRT, HSNMCCRT, LARPFCRT, LARMFCRT, LARKFCRT, LARSFCRT, LARTFCRT, UCRT, noFailStress, noFailStrain, noHashin, noLaRC05, noUcrt, SECTION_POINTS,...
    OUTPUT_PLY_POINTS, plyBuffer, thickness, OUTPUT_ENVELOPE, ENVELOPE_MODE, outputApproximate, BEST_SEQUENCE, OUTPUT_OPTIMISED, OUTPUT_FIGURE{1.0}, plyBuffer_sfailratio, axx,...
    ayy, axy, bxx, byy, bxy, E_midspan, OUTPUT_PLY, z_points, OPTIMISER_SETTINGS, CHUNK_SIZE, N_CHUNKS, EXECUTION_MODE, jobName, jobDescription, isStrengthOutput, UCRT_MException);

%% COLLECT OUTPUT
[S] = abd.internal_getOutputVars(ABD, Qij, Qt, E_midspan, E_ply_xy, E_ply_aligned, E_therm_xy, E_therm_aligned, E_hydro_xy, E_hydro_aligned, S_ply_xy, S_ply_aligned, EXT, EYT, GXYT,...
    NUXYT, NUYXT, EXB, EYB, GXYB, NUXYB, NUYXB, MSTRS, SFAILRATIO_STRESS, TSAIH, HOFFMAN, TSAIW, AZZIT, MSTRN, SFAILRATIO_STRAIN, HSNFTCRT, SFAILRATIO_HASHIN, HSNFCCRT, HSNMTCRT,...
    HSNMCCRT, LARPFCRT, SFAILRATIO_LARC05, LARMFCRT, LARKFCRT, LARSFCRT, LARTFCRT, UCRT, SFAILRATIO_UCRT, BEST_SEQUENCE, isStrengthOutput, outputLocation, settings, noFailStress,...
    noFailStrain, noHashin, noLaRC05, noUcrt);

%% Add the output location to the MATLAB path
addpath(genpath(outputLocation));

%% Open the results file now (if applicable)
if OUTPUT_LOCATION{2.0} == true
    try
        open([outputLocation, filesep, 'summary.log'])
    catch MException
        % Save the MATLAB exception object
        abd.internal_saveMExceptionObj(MException)
    end
end
end