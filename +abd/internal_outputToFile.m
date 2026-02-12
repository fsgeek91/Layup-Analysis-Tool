function [SFAILRATIO_STRESS, SFAILRATIO_STRAIN, SFAILRATIO_HASHIN, SFAILRATIO_LARC05, SFAILRATIO_UCRT] =...
    internal_outputToFile(dateString, outputLocation, outputStrength, nPlies, t_ply, theta, enableTensor, printTensor, S_ply_aligned, S_ply_xy, E_ply_aligned, E_ply_xy,...
    E_therm_xy, E_hydro_xy, E_therm_aligned, E_hydro_aligned, ABD, symmetricAbd, EXT, EYT, GXYT, NUXYT, NUYXT, EXB, EYB, GXYB, NUXYB, NUYXB, MSTRS, TSAIH, HOFFMAN, TSAIW, AZZIT,...
    MSTRN, HSNFTCRT, HSNFCCRT, HSNMTCRT, HSNMCCRT, LARPFCRT, LARMFCRT, LARKFCRT, LARSFCRT, LARTFCRT, UCRT, noFailStress, noFailStrain, noHashin, noLaRC05, noUcrt, SECTION_POINTS,...
    outputPoints, plyBuffer, thickness, OUTPUT_ENVELOPE, ENVELOPE_MODE, outputApproximate, BEST_SEQUENCE, OUTPUT_OPTIMISED, OUTPUT_FIGURE, plyBuffer_sfailratio, axx, ayy, axy,...
    bxx, byy, bxy, E_midspan, OUTPUT_PLY, z_points, OPTIMISER_SETTINGS, CHUNK_SIZE, N_CHUNKS, EXECUTION_MODE, JOB_NAME, JOB_DESCRIPTION, isStrengthOutput, UCRT_MException)
%   Write results output to a text file.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 5.1.0 Copyright Louis Vallance 2026
%   Last modified 12-Feb-2026 12:33:07 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%

%% Save the MATLAB Exception object (if applicable)
% Flag indicating that UCRT was run and failed
ucrtFail = (isempty(UCRT_MException) == false) && (exist('UCRT_MException', 'var') == 1.0);

if ucrtFail == true
    save([outputLocation, filesep, UCRT_MException.stack(1.0).name, '_exception.mat'], 'UCRT_MException');
end

%% Open the results file and print the header
fid = fopen([outputLocation, filesep, 'summary', '.log'], 'w+');

% Get the user's machine name
[~, hostname] = system('hostname');

% Print header
fprintf(fid, '***************************************************************************\n');
fprintf(fid, '*   For questions, comments or suggestions, please contact the author:    *\n*   help.qft@gmail.com                                                    *\n');
fprintf(fid, '*                                                                         *\n');
fprintf(fid, '*   File Exchange: 128914-layup-analysis-tool                             *\n');
fprintf(fid, '*   GitHub: https://github.com/fsgeek91/Layup-Analysis-Tool/releases      *\n');
fprintf(fid, '***************************************************************************\n\n');
fprintf(fid, 'Layup Analysis Tool 5.1.0 on machine %s\nMATLAB version %s on %s\n\n', hostname(1.0:end - 1.0), version, computer);
fprintf(fid, 'Copyright Louis Vallance 2026\nLast modified 12-Feb-2026 12:33:07 UTC\n\n');
fprintf(fid, 'ANALYSIS RESULTS GENERATED ON %s\n\n', upper(dateString));
fprintf(fid, 'Job name:  %s\n', JOB_NAME);
if isempty(JOB_DESCRIPTION) == false
    fprintf(fid, 'Job description: %s\n\n', JOB_DESCRIPTION);
else
    fprintf(fid, 'Job description: NONE\n\n');
end

% Print the units and CSYS conventions
if printTensor == true
    fprintf(fid, 'Length units: [mm]; Stress units: [N/mm2]; Strain units: [mm/mm]\n[xx, yy, xy] -> Global (x-y) CSYS\n[11, 22, 12] -> Layup (longitudinal-transverse) CSYS\n\n');
end

%% Initialise variables
SFAILRATIO_STRESS = -1.0*ones(1.0, 4.0);
SFAILRATIO_STRAIN = -1.0;
SFAILRATIO_HASHIN = -1.0*ones(1.0, 4.0);
SFAILRATIO_LARC05 = -1.0*ones(1.0, 5.0);
SFAILRATIO_UCRT = -1.0;

%% Print layup summary
fprintf(fid, 'Composite layup summary (all section points):\n');

if (enableTensor == true) && (printTensor == true)
    % Print layup summary header
    fprintf(fid, 'PLY    THICKNESS    ORIENTATION    MAX. FIBRE    MAX. TRANSVERSE    MAX. SHEAR    \n');
    fprintf(fid, '                                   STRESS        STRESS             STRESS        \n');

    % Initialise the section point index
    spIndex = 1.0;

    for i = 1.0:nPlies
        % Get the stresses for all section points of the current ply
        Si = S_ply_aligned(:, spIndex:spIndex + (SECTION_POINTS - 1.0));

        % Update the section point index
        spIndex = spIndex + SECTION_POINTS;

        % Get the numerically largest stresses over the ply
        S1iMax = abd.internal_getAbsMax(Si(1.0, :), 1.0);
        S2iMax = abd.internal_getAbsMax(Si(2.0, :), 1.0);
        S3iMax = abd.internal_getAbsMax(Si(3.0, :), 1.0);

        % Print information for the current ply
        fprintf(fid, '%-7.0f%-13g%-15g%-14g%-19g%-14g\n', i, t_ply(i), theta(i), S1iMax, S2iMax, S3iMax);

        % Draw symmetry plane (if applicable)
        if (symmetricAbd == true) && (i == nPlies/2.0)
            fprintf(fid, '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - SYM\n');
        end
    end
else
    % Print layup summary header
    fprintf(fid, 'PLY   THICKNESS    ORIENTATION\n');

    for i = 1.0:nPlies
        % Print information for the current ply
        fprintf(fid, '%-6.0f%-13g%-11g\n', i, t_ply(i), theta(i));

        % Draw symmetry plane (if applicable)
        if (symmetricAbd == true) && (i == nPlies/2.0)
            fprintf(fid, '- - - - - - - - - - - - - - - - SYM\n');
        end
    end
end

fprintf(fid, '\n===========================================================================\n');

%% Print ABD matrix
fprintf(fid, '\nA, B and D matrices:\n');

for i = 1.0:6.0
    % Print matrix components
    fprintf(fid, '%-14.5g%-14.5g%-14.5g | %-14.5g%-14.5g%-14.5g\n', ABD(i, 1.0:6.0));

    % Print separator
    if i == 3.0
        fprintf(fid, '---------------------------------------------------------------------------\n');
    end
end

if printTensor == true
    fprintf(fid, '\nNote: The layup section is integrated once before the stress analysis.\n');
end

fprintf(fid, '\n===========================================================================\n');

%% Print equivalent moduli (if applicable)
if symmetricAbd == true
    % Equivalent extensional moduli
    fprintf(fid, '\nEquivalent extensional moduli:\n');
    fprintf(fid, 'EXT = %g\nEYT = %g\nGXYT = %g\nNUXYT = %g\nNUYXT = %g\n\n', EXT, EYT, GXYT, NUXYT, NUYXT);

    % Equivalent bending moduli
    fprintf(fid, 'Equivalent bending moduli:\n');
    fprintf(fid, 'EXB = %g\nEYB = %g\nGXYB = %g\nNUXYB = %g\nNUYXB = %g\n', EXB, EYB, GXYB, NUXYB, NUYXB);
else
    % Print message about symmetric stacking sequences
    fprintf(fid, '\nNote: The equivalent extensional and bending moduli are only printed for\nsymmetric laminate stacking sequences.\n');
end

fprintf(fid, '\n===========================================================================\n');

%% Print effective thermal and hydroscopic properties
if any(any([axx; ayy; axy; bxx; byy; bxy])) == true
    %{
        Print the effective thermal and hydroscopic properties only if they
        are not all ZERO
    %}
    fprintf(fid, '\nEffective thermal/hydroscopic constants:\n');
    fprintf(fid, 'PLY   axx           ayy           axy           bxx           byy           bxy           \n');
    for i = 1.0:nPlies
        %{
            Get the section point for the current ply. Since the thermal
            and hydroscopic properties are constant over each ply, it is
            safe to accept the first value in each ply
        %}
        s = outputPoints(find(plyBuffer == i, 1.0));

        if isempty(s) == true
            % There is no data for the current ply
            fprintf(fid, '%-6.0fNO RESULTS\n', i);
        else
            % Print the data
            fprintf(fid, '%-6.0f%-14g%-14g%-14g%-14g%-14g%-14g\n', i, axx(s), ayy(s), axy(s), bxx(s), byy(s), bxy(s));
        end
    end

    fprintf(fid, '\n===========================================================================\n');
end

%% Print stress/strain tensors
if (isempty(OUTPUT_FIGURE) == false) && (printTensor == true) && (isscalar(z_points) == true)
    %{
        Inform the user if there is only one total section point for output
        and MATLAB figures were requested
    %}
    fprintf(fid, '\nNote: Insufficient section points for MATLAB figure output\n');
end

if printTensor == true
    % Set the ply location string
    if OUTPUT_ENVELOPE == 1.0
        % Envelope plot
        plyOutputString = '1 (ENVELOPE)';
    else
        if ischar(OUTPUT_PLY) == true
            % Location
            if (mod(SECTION_POINTS, 2.0) == 0.0) && (strcmpi(OUTPUT_PLY, 'MIDDLE') == true)
                plyOutputString = [sprintf('%.0f', length(outputPoints)/nPlies), ' (MIDDLE - APPROXIMATE)'];
            else
                if strcmpi(OUTPUT_PLY, 'DEFAULT') == true
                    plyOutputString = [sprintf('%.0f', length(outputPoints)/nPlies), ' (TOP AND BOTTOM)'];
                else
                    plyOutputString = [sprintf('%.0f', length(outputPoints)/nPlies), ' (', upper(OUTPUT_PLY), ')'];
                end
            end
        else
            % User-sppecified
            plyOutputString = [sprintf('%.0f', length(OUTPUT_PLY)), ' (USER-SPECIFIED)'];
        end
    end

    % Set the stress/strain tensor header
    header = sprintf('Stress/strain tensor calculation summary for user-defined stacking sequence\nSection points per ply: %.0f\nOutput locations: %s',...
        SECTION_POINTS, plyOutputString);

    % Print the stress/strain tensor
    abd.internal_printTensor(fid, OUTPUT_ENVELOPE, ENVELOPE_MODE, S_ply_xy, S_ply_aligned, E_ply_xy, E_ply_aligned, E_therm_xy, E_hydro_xy, E_therm_aligned, E_hydro_aligned,...
        nPlies, outputPoints, plyBuffer, symmetricAbd, outputApproximate, thickness, header, OUTPUT_PLY, z_points, SECTION_POINTS)

    %% Print the midspan strain
    fprintf(fid, '\nMidspan strains and curvatures:\nExx_0         Eyy_0         Exy_0         Kxx           Kyy           Kxy\n%-14g%-14g%-14g%-14g%-14g%-14g\n',...
        E_midspan(1.0), E_midspan(2.0), E_midspan(3.0), E_midspan(4.0), E_midspan(5.0), E_midspan(6.0));
elseif printTensor == -1.0
    % Print message about zero load
    fprintf(fid, '\nNote: There is zero load in the layup. Stress/strain tensor information has\nnot been printed.\n');
    fprintf(fid, '\n===========================================================================\n');
end

%% Print critical ply summary
if isStrengthOutput == true
    fprintf(fid, '\nFAILURE CRITERIA ASSESSMENT RESULTS\n');
    fprintf(fid, '\nCritical ply summary (all criteria):\n');
    fprintf(fid, 'CRITERION     PLY           SYMMETRIC?\n');

    if noFailStress == false
        % Get the critical ply and the symmetry condition
        [MAX_MSTRS, MAX_MSTRS_DATA, MSTRS_SYM, MAX_TSAIH, MAX_TSAIH_DATA, TSAIH_SYM, MAX_HOFFMAN, MAX_HOFFMAN_DATA, HOFFMAN_SYM, MAX_TSAIW, MAX_TSAIW_DATA, TSAIW_SYM, MAX_AZZIT,...
            MAX_AZZIT_DATA, AZZIT_SYM, SFAILRATIO_STRESS, FAILED_PLY_MSTRS] =...
            ...
            abd.internal_getCriticalPly([MSTRS', TSAIH', HOFFMAN', TSAIW', AZZIT'], symmetricAbd, plyBuffer_sfailratio, nPlies);

        % Print the result
        fprintf(fid, 'MSTRS         %-14.0f%-9s\nTSAIH         %-14.0f%-9s\nHOFFMAN       %-14.0f%-9s\nTSAIW         %-14.0f%-9s\nAZZIT         %-14.0f%-9s\n', MAX_MSTRS,...
            MSTRS_SYM, MAX_TSAIH, TSAIH_SYM, MAX_HOFFMAN, HOFFMAN_SYM, MAX_TSAIW, TSAIW_SYM, MAX_AZZIT, AZZIT_SYM);

        % Get worst criterion values (per ply)
        MAX_MSTRS_VAL = MAX_MSTRS_DATA(:, 1.0);
        MAX_TSAIH_VAL = MAX_TSAIH_DATA(:, 1.0);
        MAX_HOFFMAN_VAL = MAX_HOFFMAN_DATA(:, 1.0);
        MAX_TSAIW_VAL = MAX_TSAIW_DATA(:, 1.0);
        MAX_AZZIT_VAL = MAX_AZZIT_DATA(:, 1.0);

        % Get worst section point per criterion value (per ply)
        MAX_MSTRS_SP = MAX_MSTRS_DATA(:, 2.0);
        MAX_TSAIH_SP = MAX_TSAIH_DATA(:, 2.0);
        MAX_HOFFMAN_SP = MAX_HOFFMAN_DATA(:, 2.0);
        MAX_TSAIW_SP = MAX_TSAIW_DATA(:, 2.0);
        MAX_AZZIT_SP = MAX_AZZIT_DATA(:, 2.0);
    end

    if noFailStrain == false
        % Get the critical ply and the symmetry condition
        [MAX_MSTRN, MAX_MSTRN_DATA, MSTRN_SYM, SFAILRATIO_STRAIN, FAILED_PLY_MSTRN] =...
            ...
            abd.internal_getCriticalPly(MSTRN', symmetricAbd, plyBuffer_sfailratio, nPlies);

        % Print the result
        fprintf(fid, 'MSTRN         %-14.0f%-9s\n', MAX_MSTRN, MSTRN_SYM);

        % Get worst criterion values (per ply)
        MAX_MSTRN_VAL = MAX_MSTRN_DATA(:, 1.0);

        % Get worst section point per criterion value (per ply)
        MAX_MSTRN_SP = MAX_MSTRN_DATA(:, 2.0);
    end

    if noHashin == false
        % Get the critical ply and the symmetry condition
        [MAX_HSNFTCRT, MAX_HSNFTCRT_DATA, HSNFTCRT_SYM, MAX_HSNFCCRT, MAX_HSNFCCRT_DATA, HSNFCCRT_SYM, MAX_HSNMTCRT, MAX_HSNMTCRT_DATA, HSNMTCRT_SYM, MAX_HSNMCCRT,...
            MAX_HSNMCCRT_DATA, HSNMCCRT_SYM, SFAILRATIO_HASHIN, FAILED_PLY_HSN] =...
            ...
            abd.internal_getCriticalPly([HSNFTCRT', HSNFCCRT', HSNMTCRT', HSNMCCRT'], symmetricAbd, plyBuffer_sfailratio, nPlies);

        % Print the result
        fprintf(fid, 'HSNFTCRT      %-14.0f%-9s\nHSNFCCRT      %-14.0f%-9s\nHSNMTCRT      %-14.0f%-9s\nHSNMCCRT      %-14.0f%-9s\n', MAX_HSNFTCRT, HSNFTCRT_SYM, MAX_HSNFCCRT,...
            HSNFCCRT_SYM, MAX_HSNMTCRT, HSNMTCRT_SYM, MAX_HSNMCCRT, HSNMCCRT_SYM);

        % Get worst criterion values (per ply)
        MAX_HSNFTCRT_VAL = MAX_HSNFTCRT_DATA(:, 1.0);
        MAX_HSNFCCRT_VAL = MAX_HSNFCCRT_DATA(:, 1.0);
        MAX_HSNMTCRT_VAL = MAX_HSNMTCRT_DATA(:, 1.0);
        MAX_HSNMCCRT_VAL = MAX_HSNMCCRT_DATA(:, 1.0);

        % Get worst section point per criterion value (per ply)
        MAX_HSNFTCRT_SP = MAX_HSNFTCRT_DATA(:, 2.0);
        MAX_HSNFCCRT_SP = MAX_HSNFCCRT_DATA(:, 2.0);
        MAX_HSNMTCRT_SP = MAX_HSNMTCRT_DATA(:, 2.0);
        MAX_HSNMCCRT_SP = MAX_HSNMCCRT_DATA(:, 2.0);
    end

    if noLaRC05 == false
        % Get the critical ply and the symmetry condition
        [MAX_LARPFCRT, MAX_LARPFCRT_DATA, LARPFCRT_SYM, MAX_LARMFCRT, MAX_LARMFCRT_DATA, LARMFCRT_SYM, MAX_LARKFCRT, MAX_LARKFCRT_DATA, LARKFCRT_SYM, MAX_LARSFCRT,...
            MAX_LARSFCRT_DATA, LARSFCRT_SYM, MAX_LARTFCRT, MAX_LARTFCRT_DATA, LARTFCRT_SYM, SFAILRATIO_LARC05, FAILED_PLY_LARC05] =...
            ...
            abd.internal_getCriticalPly([LARPFCRT', LARMFCRT', LARKFCRT', LARSFCRT', LARTFCRT'], symmetricAbd, plyBuffer_sfailratio, nPlies);

        % Print the result
        fprintf(fid, 'LARPFCRT      %-14.0f%-9s\nLARMFCRT      %-14.0f%-9s\nLARKFCRT      %-14.0f%-9s\nLARSFCRT      %-14.0f%-9s\nLARTFCRT      %-14.0f%-9s\n', MAX_LARPFCRT,...
            LARPFCRT_SYM, MAX_LARMFCRT, LARMFCRT_SYM, MAX_LARKFCRT, LARKFCRT_SYM, MAX_LARSFCRT, LARSFCRT_SYM, MAX_LARTFCRT, LARTFCRT_SYM);

        % Get worst criterion values (per ply)
        MAX_LARPFCRT_VAL = MAX_LARPFCRT_DATA(:, 1.0);
        MAX_LARMFCRT_VAL = MAX_LARMFCRT_DATA(:, 1.0);
        MAX_LARKFCRT_VAL = MAX_LARKFCRT_DATA(:, 1.0);
        MAX_LARSFCRT_VAL = MAX_LARSFCRT_DATA(:, 1.0);
        MAX_LARTFCRT_VAL = MAX_LARTFCRT_DATA(:, 1.0);

        % Get worst section point per criterion value (per ply)
        MAX_LARPFCRT_SP = MAX_LARPFCRT_DATA(:, 2.0);
        MAX_LARMFCRT_SP = MAX_LARMFCRT_DATA(:, 2.0);
        MAX_LARKFCRT_SP = MAX_LARKFCRT_DATA(:, 2.0);
        MAX_LARSFCRT_SP = MAX_LARSFCRT_DATA(:, 2.0);
        MAX_LARTFCRT_SP = MAX_LARTFCRT_DATA(:, 2.0);
    end

    if noUcrt == false
        % Get the critical ply and the symmetry condition
        [MAX_UCRT, MAX_UCRT_DATA, UCRT_SYM, SFAILRATIO_UCRT, FAILED_PLY_UCRT] =...
            ...
            abd.internal_getCriticalPly(UCRT', symmetricAbd, plyBuffer_sfailratio, nPlies);

        % Print the result
        fprintf(fid, 'UCRT          %-14.0f%-9s\n', MAX_UCRT, UCRT_SYM);

        % Get worst criterion values (per ply)
        MAX_UCRT_VAL = MAX_UCRT_DATA(:, 1.0);

        % Get worst section point per criterion value (per ply)
        MAX_UCRT_SP = MAX_UCRT_DATA(:, 2.0);
    end

    fprintf(fid, '\n===========================================================================\n');
end

%% Print results of failure criteria analysis (stress-based)
if (isStrengthOutput == true) && (noFailStress == false)
    % Get the parameter string
    if outputStrength{2.0} == 1.0
        parameter = '(R)';
    else
        parameter = '(V)';
    end

    % Get maximum criterion values
    FAIL_STRESS_ALL = [MAX_MSTRS_VAL, MAX_TSAIH_VAL, MAX_HOFFMAN_VAL, MAX_TSAIW_VAL, MAX_AZZIT_VAL];
    FAIL_STRESS_ALL_MAX = max(FAIL_STRESS_ALL, [], 2.0);

    % Print table header
    fprintf(fid, '\nAssessment summary for stress-based failure criteria\nOutput location: <Worst section point>\n');
    fprintf(fid, 'PLY           MSTRS(V)      @SP   TSAIH%s      @SP   HOFFMAN%s    @SP   TSAIW%s      @SP   AZZIT%s      @SP   (WORST)       STATUS\n', parameter, parameter,...
        parameter, parameter);

    % Print ply-wise results
    for i = 1.0:nPlies
        if FAILED_PLY_MSTRS(i) == true
            % All section points in the current ply have failed
            fprintf(fid, '%-14.0f%-14g%-6.0f%-14g%-6.0f%-14g%-6.0f%-14g%-6.0f%-14g%-6.0f%-14g%-6s\n', i, MAX_MSTRS_VAL(i), MAX_MSTRS_SP(i), MAX_TSAIH_VAL(i), MAX_TSAIH_SP(i),...
                MAX_HOFFMAN_VAL(i), MAX_HOFFMAN_SP(i), MAX_TSAIW_VAL(i), MAX_TSAIW_SP(i), MAX_AZZIT_VAL(i), MAX_AZZIT_SP(i), FAIL_STRESS_ALL_MAX(i), 'FAILED');
        elseif FAIL_STRESS_ALL_MAX(i) >= 1.0
            % At least one section point in the current ply have failed
            fprintf(fid, '%-14.0f%-14g%-6.0f%-14g%-6.0f%-14g%-6.0f%-14g%-6.0f%-14g%-6.0f%-14g%-6s\n', i, MAX_MSTRS_VAL(i), MAX_MSTRS_SP(i), MAX_TSAIH_VAL(i), MAX_TSAIH_SP(i),...
                MAX_HOFFMAN_VAL(i), MAX_HOFFMAN_SP(i), MAX_TSAIW_VAL(i), MAX_TSAIW_SP(i), MAX_AZZIT_VAL(i), MAX_AZZIT_SP(i), FAIL_STRESS_ALL_MAX(i), 'UNSAFE');
        elseif MAX_MSTRS_VAL(i) == -1.0
            % There is no data for the current ply
            fprintf(fid, '%-14.0fNO RESULTS\n', i);
        else
            % No section points in the current ply have failed
            fprintf(fid, '%-14.0f%-14g%-6.0f%-14g%-6.0f%-14g%-6.0f%-14g%-6.0f%-14g%-6.0f%-14g%-6s\n', i, MAX_MSTRS_VAL(i), MAX_MSTRS_SP(i), MAX_TSAIH_VAL(i), MAX_TSAIH_SP(i),...
                MAX_HOFFMAN_VAL(i), MAX_HOFFMAN_SP(i), MAX_TSAIW_VAL(i), MAX_TSAIW_SP(i), MAX_AZZIT_VAL(i), MAX_AZZIT_SP(i), FAIL_STRESS_ALL_MAX(i), 'SAFE');
        end
    end

    % Print SFAILRATIO
    %{
        Note: The value of SFAILRATIO considers results over ALL section
        points
    %}
    fprintf(fid, 'SFAILRATIO    %-20g%-20g%-20g%-20g%-20g\n', SFAILRATIO_STRESS);
end

%% Print results of failure criteria analysis (strain-based)
if (isStrengthOutput == true) && (noFailStrain == false)
    % Print table header
    fprintf(fid, '\nAssessment summary for strain-based failure criteria\nOutput location: <Worst section point>\n');
    fprintf(fid, 'PLY           MSTRN(V)      @SP   STATUS\n');

    % Print ply-wise results
    for i = 1.0:nPlies
        if FAILED_PLY_MSTRN(i) == true
            % All section points in the current ply have failed
            fprintf(fid, '%-14.0f%-14g%-6.0f%-6s\n', i, MAX_MSTRN_VAL(i), MAX_MSTRN_SP(i), 'FAILED');
        elseif MSTRN(i) >= 1.0
            % At least one section point in the current ply have failed
            fprintf(fid, '%-14.0f%-14g%-6.0f%-6s\n', i, MAX_MSTRN_VAL(i), MAX_MSTRN_SP(i), 'UNSAFE');
        elseif MAX_MSTRN_VAL(i) == -1.0
            % There is no data for the current ply
            fprintf(fid, '%-14.0fNO RESULTS\n', i);
        else
            % No section points in the current ply have failed
            fprintf(fid, '%-14.0f%-14g%-6.0f%-6s\n', i, MAX_MSTRN_VAL(i), MAX_MSTRN_SP(i), 'SAFE');
        end
    end

    % Print SFAILRATIO
    %{
        Note: The value of SFAILRATIO considers results over ALL section
        points
    %}
    fprintf(fid, 'SFAILRATIO    %-14g\n', SFAILRATIO_STRAIN);
end

%% Print parameter descriptor
if ((isStrengthOutput == true) && (noFailStress == false)) || ((isStrengthOutput == true) && (noFailStrain == false))
    fprintf(fid, '\n(R): Strength reserve factor; (V): Criterion value\n');
end

%% Print summary of failure criteria assessments
if (SFAILRATIO_STRESS(1.0) ~= -1.0) || (SFAILRATIO_STRAIN(1.0) ~= -1.0)
    if any([SFAILRATIO_STRESS, SFAILRATIO_STRAIN] == 1.0) == true
        fprintf(fid, '\nEVERY PLY IN THE LAYUP WILL FAIL BASED ON THE EVALUATED CRITERIA\n');
    elseif any([SFAILRATIO_STRESS, SFAILRATIO_STRAIN] > 0.0) == true
        fprintf(fid, '\nAT LEAST ONE PLY IN THE LAYUP WILL FAIL BASED ON THE EVALUATED CRITERIA\n');
    else
        fprintf(fid, '\nLAYUP WILL NOT FAIL BASED ON THE EVALUATED CRITERIA\n');
    end
end

%% Print results of damage initiation criteria analysis (HASHIN)
if (isStrengthOutput == true) && (noHashin == false)
    % Get maximum criterion values
    HASHIN_ALL = [MAX_HSNFTCRT_VAL, MAX_HSNFCCRT_VAL, MAX_HSNMTCRT_VAL, MAX_HSNMCCRT_VAL];
    HASHIN_ALL_MAX = max(HASHIN_ALL, [], 2.0);

    % Print table header
    fprintf(fid, '\nAssessment summary for Hashin''s damage initiation criteria\nOutput location: <Worst section point>\n');
    fprintf(fid, 'PLY           HSNFTCRT      @SP   HSNFCCRT      @SP   HSNMTCRT      @SP   HSNMCCRT      @SP   (WORST)       STATUS\n');

    % Print ply-wise results
    for i = 1.0:nPlies
        if FAILED_PLY_HSN(i) == true
            % All section points in the current ply have failed
            fprintf(fid, '%-14.0f%-14g%-6.0f%-14g%-6.0f%-14g%-6.0f%-14g%-6.0f%-14g%-6s\n', i, MAX_HSNFTCRT_VAL(i), MAX_HSNFTCRT_SP(i), MAX_HSNFCCRT_VAL(i), MAX_HSNFCCRT_SP(i),...
                MAX_HSNMTCRT_VAL(i), MAX_HSNMTCRT_SP(i), MAX_HSNMCCRT_VAL(i), MAX_HSNMCCRT_SP(i), HASHIN_ALL_MAX(i), 'FAILED');
        elseif HASHIN_ALL_MAX(i) >= 1.0
            % At least one section point in the current ply have failed
            fprintf(fid, '%-14.0f%-14g%-6.0f%-14g%-6.0f%-14g%-6.0f%-14g%-6.0f%-14g%-6s\n', i, MAX_HSNFTCRT_VAL(i), MAX_HSNFTCRT_SP(i), MAX_HSNFCCRT_VAL(i), MAX_HSNFCCRT_SP(i),...
                MAX_HSNMTCRT_VAL(i), MAX_HSNMTCRT_SP(i), MAX_HSNMCCRT_VAL(i), MAX_HSNMCCRT_SP(i), HASHIN_ALL_MAX(i), 'UNSAFE');
        elseif MAX_HSNFTCRT_VAL(i) == -1.0
            % There is no data for the current ply
            fprintf(fid, '%-14.0fNO RESULTS\n', i);
        else
            % No section points in the current ply have failed
            fprintf(fid, '%-14.0f%-14g%-6.0f%-14g%-6.0f%-14g%-6.0f%-14g%-6.0f%-14g%-6s\n', i, MAX_HSNFTCRT_VAL(i), MAX_HSNFTCRT_SP(i), MAX_HSNFCCRT_VAL(i), MAX_HSNFCCRT_SP(i),...
                MAX_HSNMTCRT_VAL(i), MAX_HSNMTCRT_SP(i), MAX_HSNMCCRT_VAL(i), MAX_HSNMCCRT_SP(i), HASHIN_ALL_MAX(i), 'SAFE');
        end
    end

    % Print SFAILRATIO
    %{
        Note: The value of SFAILRATIO considers results over ALL section
        points
    %}
    fprintf(fid, 'SFAILRATIO    %-20g%-20g%-20g%-20g\n', SFAILRATIO_HASHIN);
end

% Print summary of assessment
if SFAILRATIO_HASHIN(1.0) ~= -1.0
    if any(SFAILRATIO_HASHIN == 1.0) == true
        fprintf(fid, '\nEVERY PLY IN THE LAYUP WILL BE DAMAGED BASED ON THE EVALUATED CRITERIA\n');
    elseif any(SFAILRATIO_HASHIN > 0.0) == true
        fprintf(fid, '\nAT LEAST ONE PLY IN THE LAYUP WILL BE DAMAGED BASED ON THE EVALUATED CRITERIA\n');
    else
        fprintf(fid, '\nLAYUP WILL NOT BE DAMAGED BASED ON THE EVALUATED CRITERIA\n');
    end
end

%% Print results of damage initiation criteria analysis (LARC05)
if (isStrengthOutput == true) && (noLaRC05 == false)
    % Get maximum criterion values
    LARC05_ALL = [MAX_LARPFCRT_VAL, MAX_LARMFCRT_VAL, MAX_LARKFCRT_VAL, MAX_LARSFCRT_VAL, MAX_LARTFCRT_VAL];
    LARC05_ALL_MAX = max(LARC05_ALL, [], 2.0);

    % Print table header
    fprintf(fid, '\nAssessment summary for LaRC05 damage initiation criteria\nOutput location: <Worst section point>\n');
    fprintf(fid, 'PLY           LARPFCRT      @SP   LARMFCRT      @SP   LARKFCRT      @SP   LARSFCRT      @SP   LARTFCRT      @SP   (WORST)       STATUS\n');

    % Print ply-wise results
    for i = 1.0:nPlies
        if FAILED_PLY_LARC05(i) == true
            % All section points in the current ply have failed
            fprintf(fid, '%-14.0f%-14g%-6.0f%-14g%-6.0f%-14g%-6.0f%-14g%-6.0f%-14g%-6.0f%-14g%-6s\n', i, MAX_LARPFCRT_VAL(i), MAX_LARPFCRT_SP(i), MAX_LARMFCRT_VAL(i),...
            MAX_LARMFCRT_SP(i), MAX_LARKFCRT_VAL(i), MAX_LARKFCRT_SP(i), MAX_LARSFCRT_VAL(i), MAX_LARSFCRT_SP(i), MAX_LARTFCRT_VAL(i), MAX_LARTFCRT_SP(i), LARC05_ALL_MAX(i),...
            'FAILED');
        elseif LARC05_ALL_MAX(i) >= 1.0
            % At least one section point in the current ply have failed
            fprintf(fid, '%-14.0f%-14g%-6.0f%-14g%-6.0f%-14g%-6.0f%-14g%-6.0f%-14g%-6.0f%-14g%-6s\n', i, MAX_LARPFCRT_VAL(i), MAX_LARPFCRT_SP(i), MAX_LARMFCRT_VAL(i),...
                MAX_LARMFCRT_SP(i), MAX_LARKFCRT_VAL(i), MAX_LARKFCRT_SP(i), MAX_LARSFCRT_VAL(i), MAX_LARSFCRT_SP(i), MAX_LARTFCRT_VAL(i), MAX_LARTFCRT_SP(i), LARC05_ALL_MAX(i),...
                'UNSAFE');
        elseif MAX_LARPFCRT_VAL(i) == -1.0
            % There is no data for the current ply
            fprintf(fid, '%-14.0fNO RESULTS\n', i);
        else
            % No section points in the current ply have failed
            fprintf(fid, '%-14.0f%-14g%-6.0f%-14g%-6.0f%-14g%-6.0f%-14g%-6.0f%-14g%-6.0f%-14g%-6s\n', i, MAX_LARPFCRT_VAL(i), MAX_LARPFCRT_SP(i), MAX_LARMFCRT_VAL(i),...
                MAX_LARMFCRT_SP(i), MAX_LARKFCRT_VAL(i), MAX_LARKFCRT_SP(i), MAX_LARSFCRT_VAL(i), MAX_LARSFCRT_SP(i), MAX_LARTFCRT_VAL(i), MAX_LARTFCRT_SP(i), LARC05_ALL_MAX(i), 'SAFE');
        end
    end

    % Print SFAILRATIO
    %{
        Note: The value of SFAILRATIO considers results over ALL section
        points
    %}
    fprintf(fid, 'SFAILRATIO    %-20g%-20g%-20g%-20g%-20g\n', SFAILRATIO_LARC05);
end

% Print summary of assessment
if SFAILRATIO_LARC05(1.0) ~= -1.0
    if any(SFAILRATIO_LARC05 == 1.0) == true
        fprintf(fid, '\nEVERY PLY IN THE LAYUP WILL BE DAMAGED BASED ON THE EVALUATED CRITERIA\n');
    elseif any(SFAILRATIO_LARC05 > 0.0) == true
        fprintf(fid, '\nAT LEAST ONE PLY IN THE LAYUP WILL BE DAMAGED BASED ON THE EVALUATED CRITERIA\n');
    else
        fprintf(fid, '\nLAYUP WILL NOT BE DAMAGED BASED ON THE EVALUATED CRITERIA\n');
    end
end

%% Print results of damage initiation criterion analysis (UCRT)
if (isStrengthOutput == true) && ((isa(outputStrength{1.0}, 'function_handle') == true) && (all(UCRT == -1.0) == false))
    % Print table header
    fprintf(fid, '\nAssessment summary for user-defined damage initiation criterion\nUser routine: @%s\nOutput location: <Worst section point>\n', char(outputStrength{1.0}));
    fprintf(fid, 'PLY           UCRT          @SP   STATUS\n');

    % Print ply-wise results
    for i = 1.0:nPlies
        if FAILED_PLY_UCRT(i) == true
            % All section points in the current ply have failed
            fprintf(fid, '%-14.0f%-14g%-6.0f%-6s\n', i, MAX_UCRT_VAL(i), MAX_UCRT_SP(i), 'FAILED');
        elseif MAX_UCRT_VAL(i) >= 1.0
            % At least one section point in the current ply have failed
            fprintf(fid, '%-14.0f%-14g%-6.0f%-6s\n', i, MAX_UCRT_VAL(i), MAX_UCRT_SP(i), 'UNSAFE');
        elseif MAX_UCRT_VAL(i) == -1.0
            % There is no data for the current ply
            fprintf(fid, '%-14.0fNO RESULTS\n', i);
        else
            % No section points in the current ply have failed
            fprintf(fid, '%-14.0f%-14g%-6.0f%-6s\n', i, MAX_UCRT_VAL(i), MAX_UCRT_SP(i), 'SAFE');
        end
    end

    % Print SFAILRATIO
    %{
        Note: The value of SFAILRATIO considers results over ALL section
        points
    %}
    fprintf(fid, 'SFAILRATIO    %-20g\n', SFAILRATIO_UCRT);
end

% Print summary of assessment
if SFAILRATIO_UCRT(1.0) ~= -1.0
    if any(SFAILRATIO_UCRT == 1.0) == true
        fprintf(fid, '\nEVERY PLY IN THE LAYUP WILL BE DAMAGED BASED ON THE EVALUATED CRITERION\n');
    elseif any(SFAILRATIO_UCRT > 0.0) == true
        fprintf(fid, '\nAT LEAST ONE PLY IN THE LAYUP WILL BE DAMAGED BASED ON THE EVALUATED CRITERION\n');
    else
        fprintf(fid, '\nLAYUP WILL NOT BE DAMAGED BASED ON THE EVALUATED CRITERION\n');
    end
end

%% Print failure assessment summary
if (isStrengthOutput == true) && (any(~[noFailStress, noFailStrain, noHashin]) == true)
    fprintf(fid, '\nNotes about failure/damage initiation assessment output:\n\t');
    fprintf(fid, '- The assessment is performed at every section point in the layup,\n\t  regardless of the setting of OUTPUT_PLY\n\t');
    fprintf(fid, '- The assessment criteria report the worst section point for each ply\n\t');
    fprintf(fid, ['- A section point is considered to have failed when failure or damage\n\t  is reported according to at least one of the evaluated failure\n\t  indexes or damage',...
        ' initiation criteria for the selected strength\n\t  assessment\n\t']);
    fprintf(fid, '- A ply is considered to have failed when all of the section points in\n\t  the ply have failed\n\t');
    fprintf(fid, '- SFAILRATIO is the section failure ratio across all the plies\n\t  (NUMBER_OF_FAILED_PLIES/TOTAL_NUMBER_OF_PLIES)\n\t');
    fprintf(fid, ['- The plies are marked with the STATUS flags as follows:\n\t\tSAFE: No section points in the ply have failed\n\t\tFAILED: All of the section points in the ply h',...
        'ave failed\n\t\tUNSAFE: The ply has not failed, but at least one section point in\n\t\tthe ply has failed\n\t']);
    fprintf(fid, '- The worst section point value for the ply may be greater than 1\n\t  without the ply''s status being marked as FAILED\n');

    fprintf(fid, '\n===========================================================================\n');
end

%% Print optimisation results summary
if (isempty(BEST_SEQUENCE) == false) && (isempty(BEST_SEQUENCE{5.0}) == false)
    % Get the exception object
    exception = BEST_SEQUENCE{5.0};

    % Print message about no optimisation output
    fprintf(fid, '\nException: Stacking optimisation was not performed.\n\tidentifier: %s\n\tmessage: %s\n\n', exception.identifier, exception.message);

    % Advise the user to select a different optimization method
    if OPTIMISER_SETTINGS{1.0} == 1.0
        fprintf(fid, '\nWarning: The FULL MATRIX method is not recommended. Use MIXED-RADIX or\nCHUNKS instead.\n\n');
    end

    % Get the error stack object
    stack = exception.stack;

    % Print contents of the error stack object
    fprintf(fid, 'Error stack info:');
    for i = 1:length(stack)
        fprintf(fid, '\n\t<LEVEL %.0f>\n\tfile: %s\n\tname: %s\n\tline: %.0f', i, stack(i).file, stack(i).name, stack(i).line);
    end
    fprintf(fid, '\n\nRequired number of iterations for optimisation: %g', BEST_SEQUENCE{3.0});

    fprintf(fid, '\n\n===========================================================================\n');
elseif isempty(BEST_SEQUENCE) == false
    % Print the settings header
    fprintf(fid, '\nStacking sequence optimisation settings:\n');

    % Get the criterion selected for optimisation
    optiCriterion = lower(OUTPUT_OPTIMISED{2.0});

    % Print the criterion
    switch optiCriterion
        case 'mstrs'
            criterionString = 'Maximum stress';
        case 'tsaih'
            criterionString = 'Tsai-Hill';
        case 'hoffman'
            criterionString = 'Hoffman';
        case 'tsaiw'
            criterionString = 'Tsai-Wu';
        case 'azzit'
            criterionString = 'Azzi-Tsai-Hill';
        case 'mstrn'
            criterionString = 'Mean strain';
        case 'hashin'
            criterionString = 'Hashin (worst criterion)';
        case 'larc05'
            criterionString = 'LaRC05 (worst criterion)';
        case 'ucrt'
            criterionString = 'User-defined';
        otherwise
            % This condition should never be reached!
    end

    % Modify the criterion string (if applicable)
    if (strcmpi(optiCriterion, 'tsaih') == true) || (strcmpi(optiCriterion, 'hoffman') == true) || (strcmpi(optiCriterion, 'tsaiw') == true) ||...
            (strcmpi(optiCriterion, 'azzit') == true) || (strcmpi(optiCriterion, 'ucrt') == true)
        if OUTPUT_OPTIMISED{3.0} == 1.0
            % Criterion uses strength reserve factor
            criterionString = [criterionString, ' (strength reserve factor'];
        else
            % Criterion uses computed value
            criterionString = [criterionString, ' (criterion value'];
        end

        % Append additional info in case criterion is UCRT
        if strcmpi(OUTPUT_OPTIMISED{2.0}, 'ucrt') == true
            criterionString = [criterionString, ' - requested)'];
        else
            criterionString = [criterionString, ')'];
        end
    end
    fprintf(fid, 'Criterion: %s', criterionString);

    % Print the objective function
    if OUTPUT_OPTIMISED{4.0} == 1.0
        fprintf(fid, '\nObjective function: MinMax (minimise the maximum criterion value)');
    else
        fprintf(fid, '\nObjective function: MinMean (minimise the average criterion value)');
    end

    % Print the precision
    fprintf(fid, '\nPrecision: %g degrees', OUTPUT_OPTIMISED{5.0}(2.0) - OUTPUT_OPTIMISED{5.0}(1.0));

    % Print the optimisation method
    switch OPTIMISER_SETTINGS{1.0}
        case 1.0
            fprintf(fid, '\nMethod: Full matrix');
            fprintf(fid, '\n\tWarning: This method is not recommended. For improved scalability, use\n\tMIXED-RADIX or CHUNKS instead');
        case 2.0
            fprintf(fid, '\nMethod: Index-based generation (mixed-radix)');
        case 3.0
            fprintf(fid, '\nMethod: Chunking + worker looping');

            % Print the chunk size and number of chunks
            if isempty(CHUNK_SIZE) == false
                fprintf(fid, '\n\tChunk size: %.0f', CHUNK_SIZE);
                fprintf(fid, '\n\tNumber of chunks: %.0f', N_CHUNKS);
            end
        otherwise
            fprintf(fid, '\nMethod: UNKNOWN');
    end

    % Print the execution mode
    fprintf(fid, '\nExecution mode: %s', EXECUTION_MODE);

    if (OPTIMISER_SETTINGS{1.0} == 3.0) && (strcmpi(EXECUTION_MODE, 'SERIAL') == true)
        fprintf(fid, ['\n\tNote: The chunking + worker looping method is not recommended for\n\tserial execution. Either connect a parallel pool or use index-based\n\tgeneration (',...
            'mixed-radix) instead']);
    end

    % Print the results header
    fprintf(fid, '\n\nStacking sequence optimisation results summary:');

    % Print the formatted stacking sequence string
    seqStr = BEST_SEQUENCE{1.0};
    seqStrFormatted = '[';
    for i = 1:length(seqStr)
        seqStrFormatted = [seqStrFormatted, sprintf('%g, ', seqStr(i))]; %#ok<AGROW>
    end
    seqStrFormatted(end - 1.0:end) = [];
    seqStrFormatted = [seqStrFormatted, ']'];

    % Print the optimised stacking sequence
    fprintf(fid, '\nOptimised stacking sequence: %s', sprintf('%s', seqStrFormatted));

    % Print the critical value
    fprintf(fid, '\nCritical value: %g', BEST_SEQUENCE{2.0});

    % Print other calculation data
    str = sprintf('%.0f', BEST_SEQUENCE{3.0});
    str_with_commas = regexprep(fliplr(str), '(\d{3})(?=\d)', '$1,');
    str_with_commas = fliplr(str_with_commas);

    fprintf(fid, '\n(Checked %s stacking sequence permutations in %g seconds)\n', str_with_commas, BEST_SEQUENCE{4.0});

    % Print the optimised stress/strain tensor
    abd.internal_printTensor(fid, OUTPUT_ENVELOPE, ENVELOPE_MODE, BEST_SEQUENCE{6.0}.STRESS_XY, BEST_SEQUENCE{6.0}.STRESS_PLY, BEST_SEQUENCE{6.0}.STRAIN_XY,...
        BEST_SEQUENCE{6.0}.STRAIN_PLY, [], [], [], [], nPlies, outputPoints, plyBuffer, BEST_SEQUENCE{6.0}.SYMMETRIC_ABD, outputApproximate, thickness,...
        sprintf('Stress/strain calculation summary for optimised stacking sequence\nSection points per ply: <From layup definition>\nOutput locations: <From layup definition>'),...
        OUTPUT_PLY, z_points, SECTION_POINTS)

    % Print post-optimisation layup summary
    fprintf(fid, '\nPost-optimisation composite layup summary (all section points):\n');

    % Print layup summary header
    fprintf(fid, 'PLY    THICKNESS    ORIENTATION    MAX. FIBRE    (%% Change)       MAX. TRANSVERSE    (%% Change)       MAX. SHEAR    (%% Change)       \n');
    fprintf(fid, '                                   STRESS                         STRESS                              STRESS        \n');

    % Initialise the section point index
    spIndex = 1.0;

    for i = 1.0:nPlies
        % Get the stresses for all section points of the current ply
        Si = S_ply_aligned(:, spIndex:spIndex + (SECTION_POINTS - 1.0));

        %{
            Get the post-optimisation stresses for all section points of
            the current ply
        %}
        Si_opt = BEST_SEQUENCE{6.0}.STRESS_PLY(:, spIndex:spIndex + (SECTION_POINTS - 1.0));

        % Update the section point index
        spIndex = spIndex + SECTION_POINTS;

        % Get the numerically largest stresses over the ply (pre-optimisation)
        S1iMax = abd.internal_getAbsMax(Si(1.0, :), 1.0);
        S2iMax = abd.internal_getAbsMax(Si(2.0, :), 1.0);
        S3iMax = abd.internal_getAbsMax(Si(3.0, :), 1.0);

        % Get the numerically largest stresses over the ply (post-optimisation)
        S1_opt_iMax = abd.internal_getAbsMax(Si_opt(1.0, :), 1.0);
        S2_opt_iMax = abd.internal_getAbsMax(Si_opt(2.0, :), 1.0);
        S3_opt_iMax = abd.internal_getAbsMax(Si_opt(3.0, :), 1.0);

        % Get the % stress change (symmetric formula)
        if (S1iMax == 0.0) && (S1_opt_iMax == 0.0)
            S1_opt_iMax_reduction = 0.0;
        else
            S1_opt_iMax_reduction = ((S1_opt_iMax - S1iMax) / max(abs(S1iMax), abs(S1_opt_iMax))) * 100.0;
        end
        if (S2iMax == 0.0) && (S2_opt_iMax == 0.0)
            S2_opt_iMax_reduction = 0.0;
        else
            S2_opt_iMax_reduction = ((S2_opt_iMax - S2iMax) / max(abs(S2iMax), abs(S2_opt_iMax))) * 100.0;
        end
        if (S3iMax == 0.0) && (S3_opt_iMax == 0.0)
            S3_opt_iMax_reduction = 0.0;
        else
            S3_opt_iMax_reduction = ((S3_opt_iMax - S3iMax) / max(abs(S3iMax), abs(S3_opt_iMax))) * 100.0;
        end

        % Format % stress change into string
        if S1_opt_iMax_reduction >= 0.0
            S1_opt_iMax_reduction = sprintf('+%g', S1_opt_iMax_reduction);
        else
            S1_opt_iMax_reduction = sprintf('%g', S1_opt_iMax_reduction);
        end
        if S2_opt_iMax_reduction >= 0.0
            S2_opt_iMax_reduction = sprintf('+%g', S2_opt_iMax_reduction);
        else
            S2_opt_iMax_reduction = sprintf('%g', S2_opt_iMax_reduction);
        end
        if S3_opt_iMax_reduction >= 0.0
            S3_opt_iMax_reduction = sprintf('+%g', S3_opt_iMax_reduction);
        else
            S3_opt_iMax_reduction = sprintf('%g', S3_opt_iMax_reduction);
        end

        % Print information for the current ply
        fprintf(fid, '%-7.0f%-13g%-15g%-14g%-17s%-19g%-17s%-14g%-17s\n', i, t_ply(i), BEST_SEQUENCE{1.0}(i), S1_opt_iMax, S1_opt_iMax_reduction, S2_opt_iMax,...
            S2_opt_iMax_reduction, S3_opt_iMax, S3_opt_iMax_reduction);

        % Draw symmetry plane (if applicable)
        if (BEST_SEQUENCE{6.0}.SYMMETRIC_ABD == true) && (i == nPlies/2.0)
            fprintf(fid, '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - SYM\n');
        end
    end

    fprintf(fid, '\n===========================================================================\n');
elseif ((isempty(OUTPUT_OPTIMISED{1.0}) == false) && (OUTPUT_OPTIMISED{1.0} == true))
    % Print message about no optimisation output
    if (strcmpi(OUTPUT_OPTIMISED{2.0}, 'ucrt') == true) && (ucrtFail == true)
        fprintf(fid, ['\nError: Stacking sequence optimisation results are unavailable because the\nuser-defined failure criterion failed upstream. Please select a different\nfail',...
            'ure/damage initiation criterion.\n']);
    elseif printTensor ~= -1.0
        fprintf(fid, '\nNote: Stacking sequence optimisation results are unavailable. The strength\ncalculation must first be enabled with OUTPUT_STRENGTH = {true, <param>}.\n');
    else
        fprintf(fid, '\nNote: Stacking sequence optimisation results are unavailable. A load matrix\ndefinition is required using NXX/NYY/NXY and MXX/MYY/MXY.\n');
    end
    fprintf(fid, '\n===========================================================================\n');
end

%% Footer
% Print the file footer
fprintf(fid, '\nEND OF FILE');

% Close the output file
fclose(fid);

% Print notice to MATLAB command window
fprintf('[NOTICE] Results have been written to:\n''%s''\n', [outputLocation,  filesep, 'analysis_results.txt']);
end