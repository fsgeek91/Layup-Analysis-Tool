function [] = internal_outputToFile(dateString, outputLocation,...
    outputStrength, nPlies, t_ply, theta, enableTensor, printTensor,...
    S_ply_aligned, S_ply_xy, E_ply_aligned, E_ply_xy, ABD, symmetricAbd,...
    EXT, EYT, GXYT, NUXYT, NUYXT, EXB, EYB, GXYB, NUXYB, NUYXB, MSTRS,...
    TSAIH, TSAIW, AZZIT, MSTRN, noFailStress, noFailStrain,...
    nSectionPoints, outputPoints, plyBuffer, thickness, OUTPUT_ENVELOPE,...
    ENVELOPE_MODE, outputApproximate, BEST_SEQUENCE, OUTPUT_OPTIMISED,...
    OUTPUT_FIGURE)
%   Write results output to a text file.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 2.2 Copyright Louis Vallance 2023
%   Last modified 09-May-2023 07:31:07 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%

%% Open the results file and print the header
fid = fopen([outputLocation, '\analysis_results.txt'], 'w+');

% Print header
fprintf(fid, 'Layup Analysis Tool\n');
fprintf(fid, 'ANALYSIS RESULTS GENERATED ON %s\n\n', upper(dateString));

% Print the units
if printTensor == 1.0
    fprintf(fid, 'Stress units - [N/mm2]; Strain units - [mm/mm]\n\n');
end

%% Print layup summary
fprintf(fid, 'Composite layup summary:\n');

if (enableTensor == 1.0) && (printTensor == 1.0)
    % Print layup summary header
    fprintf(fid, ['PLY   THICKNESS    ORIENTATION    MAX. FIBRE    MAX',...
        '. TRANSVERSE    MAX. SHEAR    \n']);
    fprintf(fid, ['                                  STRESS        STR',...
        'ESS             STRESS        \n']);

    % Initialise the section point index
    spIndex = 1.0;

    for i = 1.0:nPlies
        % Get the stresses for all section points of the current ply
        Si = S_ply_aligned(:, spIndex:spIndex + (nSectionPoints - 1.0));

        % Update the section point index
        spIndex = spIndex + nSectionPoints;

        % Get the numerically largest stresses over the ply
        S1iMax = abd.internal_getAbsMax(Si(1.0, :), 1.0);
        S2iMax = abd.internal_getAbsMax(Si(2.0, :), 1.0);
        S3iMax = abd.internal_getAbsMax(Si(3.0, :), 1.0);

        % Print information for the current ply
        fprintf(fid, '%-6.0f%-13g%-15g%-14g%-19g%-14g\n', i, t_ply(i),...
            theta(i), S1iMax, S2iMax, S3iMax);

        % Draw symmetry plane (if applicable)
        if (symmetricAbd == true) && (i == nPlies/2.0)
            fprintf(fid, ['- - - - - - - - - - - - - - - - - - - - - -',...
                ' - - - - - - - - - - - - - - - - - SYM\n']);
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

fprintf(fid, ['\n=====================================================',...
    '======================\n']);

%% Print ABD matrix
fprintf(fid, '\nA, B and D matrices:\n');

for i = 1.0:6.0
    % Print matrix components
    fprintf(fid, '%-14.5g%-14.5g%-14.5g | %-14.5g%-14.5g%-14.5g\n',...
        ABD(i, 1.0:6.0));

    % Print separator
    if i == 3.0
        fprintf(fid, ['-----------------------------------------------',...
            '----------------------------\n']);
    end
end

if printTensor == true
    fprintf(fid, ['\nNote: The layup section is integrated once before',...
        ' the stress analysis\n']);
end

fprintf(fid, ['\n=====================================================',...
    '======================\n']);

%% Print equivalent moduli (if applicable)
if symmetricAbd == true
    % Equivalent extensional moduli
    fprintf(fid, '\nEquivalent extensional moduli:\n');
    fprintf(fid, ['EXT = %g\nEYT = %g\nGXYT = %g\nNUXYT = %g\nNUYXT = ',...
        '%g\n\n'], EXT, EYT, GXYT, NUXYT, NUYXT);

    % Equivalent bending moduli
    fprintf(fid, 'Equivalent bending moduli:\n');
    fprintf(fid, ['EXB = %g\nEYB = %g\nGXYB = %g\nNUXYB = %g\nNUYXB = ',...
        '%g\n'], EXB, EYB, GXYB, NUXYB, NUYXB);
else
    % Print message about symmetric stacking sequences
    fprintf(fid, ['\nNote: The equivalent extensional and bending modu',...
        'li are only printed for\nsymmetric laminate stacking sequence',...
        's.\n']);
end

fprintf(fid, ['\n=====================================================',...
    '======================\n']);

%% Print stress/strain tensors
if (isempty(OUTPUT_FIGURE) == false) && (printTensor == 1.0) && (nSectionPoints == 1.0)
    %{
        Inform the user if there is only one total section point for output
        and MATLAB figures were requested
    %}
    fprintf(fid, ['\nNote: Insufficient section points for MATLAB figu',...
        're output\n']);
end

if printTensor == 1.0
    abd.internal_printTensor(fid, OUTPUT_ENVELOPE, ENVELOPE_MODE,...
        S_ply_xy, S_ply_aligned, E_ply_xy, E_ply_aligned, nPlies,...
        outputPoints, plyBuffer, symmetricAbd, outputApproximate,...
        thickness, sprintf(['Stress/strain calculation summary for use',...
        'r-defined stacking sequence\nSection points per ply: %.0f'],...
        nSectionPoints))
elseif printTensor == -1.0
    % Print message about zero load
    fprintf(fid, ['\nNote: There is zero load in the layup. Stress/str',...
        'ain tensor information has\nnot been printed.\n']);
end

fprintf(fid, ['\n=====================================================',...
    '======================\n']);

%% Print critical ply summary
if outputStrength == 1.0
    fprintf(fid, '\nFAILURE CRITERIA ASSESSMENT RESULTS\n');
    fprintf(fid, '\nCritical ply summary (all criteria):\n');
    fprintf(fid, 'CRITERION     PLY           SYMMETRIC?\n');

    if noFailStress == false
        % Get the critical ply and the symmetry condition
        [MAX_MSTRS, MAX_MSTRS_VAL, MSTRS_SYM, MAX_TSAIH, MAX_TSAIH_VAL,...
            TSAIH_SYM, MAX_TSAIW, MAX_TSAIW_VAL, TSAIW_SYM, MAX_AZZIT,...
            MAX_AZZIT_VAL, AZZIT_SYM] =...
            ...
            abd.internal_getCriticalPly([MSTRS', TSAIH', TSAIW',...
            AZZIT'], symmetricAbd, outputPoints, plyBuffer, nPlies);

        % Print the result
        fprintf(fid, ['MSTRS         %-14.0f%-9s\nTSAIH         %-14.0',...
            'f%-9s\nTSAIW         %-14.0f%-9s\nAZZIT         %-14.0f%-',...
            '9s\n'],...
            MAX_MSTRS, MSTRS_SYM, MAX_TSAIH, TSAIH_SYM, MAX_TSAIW,...
            TSAIW_SYM, MAX_AZZIT, AZZIT_SYM);
    end

    if noFailStrain == false
        % Get the critical ply and the symmetry condition
        [MAX_MSTRN, MAX_MSTRN_VAL, MSTRN_SYM] =...
            ...
            abd.internal_getCriticalPly(MSTRN', symmetricAbd,...
            outputPoints, plyBuffer, nPlies);

        % Print the result
        fprintf(fid, 'MSTRN         %-14.0f%-9s\n',...
            MAX_MSTRN, MSTRN_SYM);
    end

    fprintf(fid, ['\n=================================================',...
        '==========================\n']);
end

%% Print results of failure criteria analysis (MAXIMUM STRESS)
if (outputStrength == 1.0) && (noFailStress == false)
    FAIL_STRESS_ALL = [MAX_MSTRS_VAL, MAX_TSAIH_VAL, MAX_TSAIW_VAL,...
        MAX_AZZIT_VAL];
    FAIL_STRESS_ALL_MAX = max(FAIL_STRESS_ALL, [], 2.0);

    fprintf(fid, ['\nFailure analysis summary (stress-based criteria):',...
        '\n']);
    fprintf(fid, ['PLY           MSTRS         TSAIH         TSAIW    ',...
        '     AZZIT         (WORST)       STATUS\n']);
    for i = 1.0:nPlies
        if FAIL_STRESS_ALL_MAX(i) >= 1.0
            fprintf(fid, '%-14.0f%-14g%-14g%-14g%-14g%-14g%-6s\n',...
                i, MAX_MSTRS_VAL(i), MAX_TSAIH_VAL(i), MAX_TSAIW_VAL(i),...
                MAX_AZZIT_VAL(i), FAIL_STRESS_ALL_MAX(i), 'UNSAFE');
        elseif MAX_MSTRS_VAL(i) == -1.0
            % There is no data for the current ply
            fprintf(fid, '%-14.0fNO RESULTS\n', i);
        else
            fprintf(fid, '%-14.0f%-14g%-14g%-14g%-14g%-14g%-6s\n',...
                i, MAX_MSTRS_VAL(i), MAX_TSAIH_VAL(i), MAX_TSAIW_VAL(i),...
                MAX_AZZIT_VAL(i), FAIL_STRESS_ALL_MAX(i), 'SAFE');
        end
    end

    % Print SFAILRATIO
    %{
        Note: The value of SFAILRATIO considers results over ALL section
        points
    %}
    SFAILRATIO = zeros(1.0, 4.0);
    for i = 1.0:4.0
        FAIL_STRESS_CRITERION = FAIL_STRESS_ALL(:, i);
        nFailedPlies = length(FAIL_STRESS_CRITERION(FAIL_STRESS_CRITERION >= 1.0));
        SFAILRATIO(i) = nFailedPlies/nPlies;
    end
    fprintf(fid, 'SFAILRATIO    %-14g%-14g%-14g%-14g\n', SFAILRATIO);
end

%% Print results of failure criteria analysis (MAXIMUM STRAIN)
if (outputStrength == 1.0) && (noFailStrain == false)
    fprintf(fid, ['\nFailure analysis summary (strain-based criteria):',...
        '\n']);
    fprintf(fid, 'PLY           MSTRN         STATUS\n');
    for i = 1.0:nPlies
        if MSTRN(i) >= 1.0
            fprintf(fid, '%-14.0f%-14g%-6s\n', i, MAX_MSTRN_VAL(i),...
                'UNSAFE');
        elseif MAX_MSTRN_VAL(i) == -1.0
            % There is no data for the current ply
            fprintf(fid, '%-14.0fNO RESULTS\n', i);
        else
            fprintf(fid, '%-14.0f%-14g%-6s\n', i, MAX_MSTRN_VAL(i),...
                'SAFE');
        end
    end

    % Print SFAILRATIO
    nFailedPlies = length(MAX_MSTRN_VAL(MAX_MSTRN_VAL >= 1.0));
    SFAILRATIO = nFailedPlies/nPlies;
    fprintf(fid, 'SFAILRATIO    %-14g\n', SFAILRATIO);
end

%% Print failure assessment summary
if (outputStrength == 1.0) && (any(~[noFailStress, noFailStrain]) == true)
    if any([MSTRS, TSAIH, TSAIW, AZZIT, MSTRN] >= 1.0) == true
        fprintf(fid, '\nLAYUP IS UNSAFE BASED ON EVALUATED CRITERIA\n');
    else
        fprintf(fid, '\nLAYUP IS SAFE BASED ON EVALUATED CRITERIA\n');
    end

    fprintf(fid, ['\nNote: The failure assessment considers ALL sectio',...
        'n points in the layup\n']);

    fprintf(fid, ['\n=================================================',...
        '==========================\n']);
end

%% Print optimisation results summary
if (isempty(BEST_SEQUENCE) == false) && (isempty(BEST_SEQUENCE{5.0}) == false)
    % Get the exception object
    exception = BEST_SEQUENCE{5.0};

    % Print message about no optimisation output
    fprintf(fid, ['\nException: Stacking optimisation was not performe',...
        'd.\n\tidentifier: %s\n\tmessage: %s\n'], exception.identifier,...
        exception.message);

    % Get the error stack object
    stack = exception.stack;

    % Print contents of the error stack object
    fprintf(fid, 'Error stack info:');
    for i = 1:length(stack)
        fprintf(fid, ['\n\t<LEVEL %.0f>\n\tfile: %s\n\tname: %s\n\tlin',...
            'e: %.0f'], i, stack(i).file, stack(i).name, stack(i).line);
    end

    fprintf(fid, ['\n\n===============================================',...
        '============================\n']);
elseif isempty(BEST_SEQUENCE) == false
    % Print the header
    fprintf(fid, '\nStacking sequence optimisation summary:\n');

    % Print the criterion
    switch lower(OUTPUT_OPTIMISED{2.0})
        case 'mstrs'
            criterionString = 'Maximum stress';
        case 'tsaih'
            criterionString = 'Tsai-Hill';
        case 'tsaiw'
            criterionString = 'Tsai-Wu';
        case 'aazit'
            criterionString = 'Azzi-Tsai-Hill';
        case 'mstrn'
            criterionString = 'Mean strain';
    end
    if OUTPUT_OPTIMISED{3.0} == 1.0
        criterionString = [criterionString, ' (reserve)'];
    else
        criterionString = [criterionString, ' (value)'];
    end
    fprintf(fid, 'Criterion: %s', criterionString);

    % Print the objective function
    if OUTPUT_OPTIMISED{4.0} == 1.0
        fprintf(fid, '\nObjective function: MinMax');
    else
        fprintf(fid, '\nObjective function: MinMean');
    end

    % Print the optimised stacking sequence
    fprintf(fid, '\nOptimised stacking sequence: %s',...
        sprintf('%g ', BEST_SEQUENCE{1.0}));

    % Print the critical value
    fprintf(fid, '\nCritical value: %g', BEST_SEQUENCE{2.0});

    % Print other calculation data
    fprintf(fid, ['\n(Checked %.0f stacking permutations in %g seconds',...
        ')\n'], BEST_SEQUENCE{3.0}, BEST_SEQUENCE{4.0});

    % Print the optimised stress/strain tensor
    abd.internal_printTensor(fid, OUTPUT_ENVELOPE, ENVELOPE_MODE,...
        BEST_SEQUENCE{6.0}.STRESS_XY, BEST_SEQUENCE{6.0}.STRESS_PLY,...
        BEST_SEQUENCE{6.0}.STRAIN_XY, BEST_SEQUENCE{6.0}.STRAIN_PLY,...
        nPlies, outputPoints, plyBuffer,...
        BEST_SEQUENCE{6.0}.SYMMETRIC_ABD, outputApproximate, thickness,...
        sprintf(['Stress/strain calculation summary for optimised stac',...
        'king sequence\nSection points per ply: From layup definition']))

    fprintf(fid, ['\n=================================================',...
        '==========================\n']);
elseif OUTPUT_OPTIMISED{1.0} == true && printTensor ~= -1.0
    % Print message about no optimisation output
    fprintf(fid, ['\nNote: Stacking optimisation was not performed. En',...
        'able the strength\ncalculation with OUTPUT_STRENGTH = true\n']);

    fprintf(fid, ['\n=================================================',...
        '==========================\n']);
end

%% Footer
% Print the file footer
fprintf(fid, '\nEND OF FILE');

% Close the output file
fclose(fid);

% Print notice to MATLAB command window
fprintf('[ABD] Results have been written to:\n''%s''\n',...
    [outputLocation, '\analysis_results.txt']);
end