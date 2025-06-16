classdef internal_optimise < handle
%   Find the optimum stacking sequence based on the load and selected
%   failure measure.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 4.2.0 Copyright Louis Vallance 2025
%   Last modified 10-Jun-2025 08:28:19 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%

    methods(Static = true, Access = public)
        %% RUN THE OPTIMISER
        function [BEST_SEQUENCE, CRITERION_BUFFER, MIN_CRITERION, CHUNK_SIZE, N_CHUNKS, EXECUTION_MODE] =...
                main(OUTPUT_OPTIMISED, nargin, nPlies, nPlies_points, nSectionPoints, z, z_points, Q11, Q22, Q66, Q12, A11_points, A22_points, B11_points, B22_points, tolerance,...
                XT, XC, YT, YC, S, C12, B12, XET, XEC, YET, YEC, SE, ALPHA, XHT, XHC, YHT, YHC, SHX, SHY, XLT, XLC, YLT, YLC, SLX, SLY, GL12, NL, NT, A0, PHI0, deltaT, deltaM,...
                Nxx, Nyy, Nxy, Mxx, Myy, Mxy, E11, E22, V12, G12, symsAvailable, S1, S2, S3, SECTION_POINTS, OPTIMISER_SETTINGS)
            % Initialise output
            %{
                BEST_SEQUENCE(1) = Optimum stacking sequence
                BEST_SEQUENCE(2) = Critical value
                BEST_SEQUENCE(3) = Number of permutations
                BEST_SEQUENCE(4) = Analysis time (s)
                BEST_SEQUENCE(5) = Exception
                BEST_SEQUENCE(6) = Best tensor structure
            %}
            BEST_SEQUENCE = cell(1.0, 6.0);
            CRITERION_BUFFER = [];
            MIN_CRITERION = [];
            CHUNK_SIZE = [];
            N_CHUNKS = [];
            EXECUTION_MODE = 'N/A';

            % Get data from OUTPUT_OPTIMISED
            [enabled, failureCriterion, parameter, objective, thetaAll] = deal(OUTPUT_OPTIMISED{1.0}, OUTPUT_OPTIMISED{2.0}, OUTPUT_OPTIMISED{3.0}, OUTPUT_OPTIMISED{4.0},...
                OUTPUT_OPTIMISED{5.0});

            if enabled == false
                % Do not run optimisation if it was disabled by the user
                return
            else
                % Inform the user that optimisation has started
                %fprintf('[NOTICE] Stacking sequence optimisation started. Please be patient...\n');
            end

            % Get the number of ply orientations to be considered
            numAngles = numel(thetaAll);

            % Get the number of permutations
            nPermutations = length(thetaAll)^nPlies;

            % Warn the user if a large number of iteration will be required
            if nPermutations > 1e6
                fprintf('[WARNING] More than one million iterations are required by the stacking\nsequence optimiser\n');
                fprintf('-> The analysis may take a long time to complete\n');
                fprintf('-> The required number of iterations is given by ANGLE_PRECISION^NUMBER_OF_PLIES\n');
            end

            % Set dummy variable
            dummy = zeros(1.0, nPlies);
            
            % Get the MATLAB parallel pool object
            poolObj = gcp("nocreate");

            % Get the execution mode
            if isempty(poolObj) == true
                EXECUTION_MODE = 'Serial';
            elseif poolObj.Connected == true
                EXECUTION_MODE = 'Parallel';
            else
                EXECUTION_MODE = 'Unknown';
            end

            % Save number of analysed permutations
            BEST_SEQUENCE{3.0} = nPermutations;

            % Start the timer
            timer = tic;

            switch OPTIMISER_SETTINGS{1.0}
                case 1.0  % FULL MATRIX
                    % Get the stacking permutation matrix
                    try
                        indexPermutations =...
                            ...
                            fig.combinator(numAngles, nPlies, 'p', 'r');
                    catch seqException
                        % A problem occurred while getting the combinations
                        BEST_SEQUENCE{5.0} = seqException;
                        return
                    end

                    % Get the angles
                    anglePermutations = thetaAll(indexPermutations);

                    % Transpose ANGLEPERMUTATIONS in case there is only one ply
                    if (width(indexPermutations) == 1.0) && (nPlies == 1.0)
                        anglePermutations = anglePermutations';
                    end

                    % Run the stacking sequence optimiser
                    CRITERION_BUFFER = abd.internal_optimise.parfor1(nPermutations, anglePermutations, nSectionPoints, nPlies, z, dummy, tolerance, Q11, Q12, Q66, Q22, A11_points,...
                        A22_points, B11_points, B22_points, nargin, deltaT, deltaM, Nxx, Nyy, Nxy, Mxx, Myy, Mxy, nPlies_points, z_points, failureCriterion, XT, XC, YT, YC, S,...
                        parameter, C12, B12, E11, E22, V12, G12, XET, XEC, YET, YEC, SE, ALPHA, XHT, XHC, YHT, YHC, SHX, SHY, symsAvailable, S1, S2, S3, GL12, XLT, XLC, YLT, YLC,...
                        SLX, SLY, A0, PHI0, NL, NT, SECTION_POINTS, objective);
                case 2.0 % MIXED-RADIX REPRESENTATION
                    % Run the stacking sequence optimiser
                    try
                        CRITERION_BUFFER = abd.internal_optimise.parfor2(nPermutations, thetaAll, numAngles, nSectionPoints, nPlies, z, dummy, tolerance, Q11, Q12, Q66, Q22,...
                            A11_points, A22_points, B11_points, B22_points, nargin, deltaT, deltaM, Nxx, Nyy, Nxy, Mxx, Myy, Mxy, nPlies_points, z_points, failureCriterion, XT,...
                            XC, YT, YC, S, parameter, C12, B12, E11, E22, V12, G12, XET, XEC, YET, YEC, SE, ALPHA, XHT, XHC, YHT, YHC, SHX, SHY, symsAvailable, S1, S2, S3, GL12,...
                            XLT, XLC, YLT, YLC, SLX, SLY, A0, PHI0, NL, NT, SECTION_POINTS, objective);
                    catch seqException
                        % A problem occurred while running the optimisation
                        BEST_SEQUENCE{5.0} = seqException;
                        return
                    end
                    
                case 3.0 % CHUNKS
                    if strcmpi(OPTIMISER_SETTINGS{2.0}, 'DEFAULT') == true
                        if isempty(poolObj) == false
                            % Parallel mode: Optimise chunk size for worker load

                            % Get the tuning constant
                            if strcmpi(OPTIMISER_SETTINGS{3.0}, 'DEFAULT') == true
                                % Default value
                                k = 5.0;
                            else
                                % User-defined
                                k = OPTIMISER_SETTINGS{3.0};
                            end

                            % Set the chunk size
                            CHUNK_SIZE = ceil(nPermutations/(poolObj.NumWorkers*k));
                        else
                            %{
                                Serial mode: Use smaller chunks for
                                frequent saves. Aim for about 40 total
                                chunks
                            %}
                            CHUNK_SIZE = ceil(nPermutations/40.0);
                        end
                    else
                        % User-defined chunk size
                        CHUNK_SIZE = OPTIMISER_SETTINGS{2.0};
                    end

                    % Get the number of chunks
                    N_CHUNKS = ceil(nPermutations/CHUNK_SIZE);

                    % Run the stacking sequence optimiser
                    try
                        CRITERION_BUFFER = abd.internal_optimise.parfor3(N_CHUNKS, CHUNK_SIZE, nPermutations, thetaAll, numAngles, nSectionPoints, nPlies, z, dummy, tolerance,...
                            Q11, Q12, Q66, Q22, A11_points, A22_points, B11_points, B22_points, nargin, deltaT, deltaM, Nxx, Nyy, Nxy, Mxx, Myy, Mxy, nPlies_points, z_points,...
                            failureCriterion, XT, XC, YT, YC, S, parameter, C12, B12, E11, E22, V12, G12, XET, XEC, YET, YEC, SE, ALPHA, XHT, XHC, YHT, YHC, SHX, SHY,...
                            symsAvailable, S1, S2, S3, GL12, XLT, XLC, YLT, YLC, SLX, SLY, A0, PHI0, NL, NT, SECTION_POINTS, objective);
                    catch seqException
                        % A problem occurred while running the optimisation
                        BEST_SEQUENCE{5.0} = seqException;
                        return
                    end
                    
                otherwise
            end

            % Stop the timer
            BEST_SEQUENCE{4.0} = toc(timer);

            % Get the permutation with the minimum value of CRITERION
            MIN_CRITERION = find(CRITERION_BUFFER == min(CRITERION_BUFFER), 1.0);

            %{
                Get the stacking sequence corresponding to the optimal
                permutation
            %}
            switch OPTIMISER_SETTINGS{1.0}
                case 1.0
                    BEST_SEQUENCE{1.0} = anglePermutations(MIN_CRITERION, :);
                otherwise
                    BEST_SEQUENCE{1.0} = abd.internal_optimise.indexToStacking(MIN_CRITERION, thetaAll, numAngles, nPlies);
            end

            % Save the critical value
            BEST_SEQUENCE{2.0} = CRITERION_BUFFER(MIN_CRITERION);

            % Get the optimised stress/strain tensors
            [E_ply_xy, S_ply_xy, E_ply_aligned, S_ply_aligned, symmetricAbd] =...
                ...
                abd.internal_optimise.getOptiStressStrain(BEST_SEQUENCE{1.0}, nSectionPoints, nPlies, z, dummy, tolerance, Q11, Q12, Q66, Q22, A11_points, A22_points, B11_points,...
                B22_points, nargin, deltaT, deltaM, Nxx, Nyy, Nxy, Mxx, Myy, Mxy, nPlies_points, z_points);

            % Collect output from stress/strain analysis
            BEST_SEQUENCE{6.0} = struct('STRESS_XY', S_ply_xy, 'STRESS_PLY', S_ply_aligned, 'STRAIN_XY', E_ply_xy, 'STRAIN_PLY', E_ply_aligned, 'SYMMETRIC_ABD', symmetricAbd);
        end

        %% GET THE OPTIMISED STRESS/STRAIN TENSORS
        function [E_ply_xy, S_ply_xy, E_ply_aligned, S_ply_aligned, symmetricAbd] =...
                getOptiStressStrain(theta, nSectionPoints, nPlies, z, dummy, tolerance, Q11, Q12, Q66, Q22, A11_points, A22_points, B11_points, B22_points, nargin, deltaT, deltaM,...
                Nxx, Nyy, Nxy, Mxx, Myy, Mxy, nPlies_points, z_points)
            % Get the values of theta over the section points
            [~, ~, theta_points, ~, ~, ~, ~, ~, ~, ~] =...
                ...
                abd.internal_getSectionPoints(nSectionPoints, '', nPlies, theta, z, dummy, dummy, dummy, dummy, tolerance, []);

            % Compute transformed reduced stiffness matrix components
            [Q11t, Q12t, Q16t, Q22t, Q26t, Q66t] =...
                ...
                abd.internal_getTransformedQ(nPlies, theta, Q11, Q12, Q66, Q22);

            %{
                Get effective thermal and moisture expansion coefficients
                for each ply
            %}
            [axx, ayy, axy, bxx, byy, bxy] =...
                ...
                abd.internal_getThermoHydro(theta_points, A11_points, A22_points, B11_points, B22_points);

            % Compute A, B and D matrices
            [ABD, ~, Qijt, NxxT, NyyT, NxyT, MxxT, MyyT, MxyT, NxxM, NyyM, NxyM, MxxM, MyyM, MxyM] =...
                ...
                abd.internal_getABD(nPlies, Q11t, Q12t, Q16t, Q22t, Q26t, Q66t, z, nargin, deltaT, deltaM, axx, ayy, axy, bxx, byy, bxy, nSectionPoints);

            % Compute tensor quantities
            [~, E_ply_xy, S_ply_xy, E_ply_aligned, S_ply_aligned, ~, ~, ~, ~] =...
                ...
                abd.internal_getTensor(ABD, Nxx, NxxT, NxxM, Nyy, NyyT, NyyM, Nxy, NxyT, NxyM, Mxx, MxxT, MxxM, Myy, MyyT, MyyM, Mxy, MxyT, MxyM, nPlies_points, z_points,...
                theta_points, Qijt, deltaT, deltaM, axx, ayy, axy, bxx, byy, bxy, tolerance);

            % Determine if ABD matrix is symmetric
            symmetricAbd = abd.internal_getSymmetry(ABD, tolerance);
        end

        %% GET DATA FROM OUTPUT_OPTIMISED
        function [error, output] = getSettings(OUTPUT_OPTIMISED, noFailStress, noFailStrain, noHashin, noLaRC05, OUTPUT_STRENGTH)
            % Initialise output
            error = false;
            output = cell(1.0, 4.0);

            if iscell(OUTPUT_OPTIMISED) == false
                % Convert to cell if necessary
                OUTPUT_OPTIMISED = {OUTPUT_OPTIMISED};
            end

            if (all(cellfun(@isempty, OUTPUT_OPTIMISED)) == true) || (length(OUTPUT_OPTIMISED) ~= 4.0)
                % Incorrect number of arguments
                fprintf('[ERROR] The setting OUTPUT_OPTIMISED requires four\narguments: {''<criterion>'', ''<parameter>'', ''<fun>'', theta}\n');

                % Reset the error flag and RETURN
                error = true;
                return
            end

            % Process the first argument
            argument = OUTPUT_OPTIMISED{1.0};

            % Check validity of the argument
            if isempty(argument) == true
                output{1.0} = false;
            elseif ischar(argument) == false
                % Incorrect argument type
                fprintf('[ERROR] OUTPUT_OPTIMISED(1) must be a string\n');

                    % Reset the error flag and RETURN
                    error = true;
                    return
            elseif ischar(argument) == true
                if (strcmpi(argument, 'mstrs') == false) && (strcmpi(argument, 'tsaih') == false) && (strcmpi(argument, 'tsaiw') == false) &&...
                        (strcmpi(argument, 'azzit') == false) && (strcmpi(argument, 'mstrn') == false) && (strcmpi(argument, 'hashin') == false) &&...
                        (strcmpi(argument, 'larc05') == false)
                    % Unregognised parameter
                    if (strcmpi(argument, 'hsnftcrt') == true) || (strcmpi(argument, 'hsnfccrt') == true) || (strcmpi(argument, 'hsnmtcrt') == true) ||...
                            (strcmpi(argument, 'hsnmccrt') == true)
                        fprintf('[ERROR] Parameter ''HSNFTCRT'', ''HSNFCCRT'', ''HSNMTCRT'' or ''HSNMCCRT'' is used for OUTPUT_OPTIMISED(1). Specify ''HASHIN'' instead\n');
                    elseif (strcmpi(argument, 'larpfcrt') == true) || (strcmpi(argument, 'larmfcrt') == true) || (strcmpi(argument, 'larkfcrt') == true) ||...
                            (strcmpi(argument, 'larsfcrt') == true) || (strcmpi(argument, 'lartfcrt') == true)
                        fprintf('[ERROR] Parameter ''LARPFCRT'', ''LARMFCRT'', ''LARKFCRT''m ''LARSFCRT'' or ''LARTFCRT'' is used for OUTPUT_OPTIMISED(1). Specify ''LARC05'' instead\n');
                    else
                        fprintf('[ERROR] OUTPUT_OPTIMISED(1) must be one of the following parameters: MSTRS, TSAIH, TSAIW, AZZIT, MSTRN, HASHIN or LARC05\n');
                    end

                    % Reset the error flag and RETURN
                    error = true;
                    return
                elseif (noFailStress == true) && (OUTPUT_STRENGTH == true) && (strcmpi(argument, 'mstrs') == true || strcmpi(argument, 'tsaih') == true ||...
                        strcmpi(argument, 'tsaiw') == true || strcmpi(argument, 'azzit') == true)
                    % Insufficient material data
                    fprintf('[ERROR] Requested a stress-based criterion for optimisation, but FAIL_STRESS properties are not available\n');

                    % Reset the error flag and RETURN
                    error = true;
                    return
                elseif (noFailStrain == true) && (OUTPUT_STRENGTH == true) && (strcmpi(argument, 'mstrn') == true)
                    % Insufficient material data
                    fprintf('[ERROR] Requested a strain-based criterion for optimisation, but FAIL_STRAIN properties are not available\n');

                    % Reset the error flag and RETURN
                    error = true;
                    return
                elseif (noHashin == true) && (OUTPUT_STRENGTH == true) && (strcmpi(argument, 'hashin') == true)
                    % Insufficient material data
                    fprintf('[ERROR] Requested a Hashin criterion for optimisation, but HASHIN properties are not available\n');

                    % Reset the error flag and RETURN
                    error = true;
                    return
                elseif (noLaRC05 == true) && (OUTPUT_STRENGTH == true) && (strcmpi(argument, 'larc05') == true)
                    % Insufficient material data
                    fprintf('[ERROR] Requested a LaRC05 criterion for optimisation, but LARC05 properties are not available\n');

                    % Reset the error flag and RETURN
                    error = true;
                    return
                else
                    % Everything is OK
                    output{1.0} = true;
                    output{2.0} = argument;
                end
            else
                % Everything is OK
                output{1.0} = true;
                output{2.0} = argument;
            end

            % Process the second argument
            argument = OUTPUT_OPTIMISED{2.0};

            % Get the parameter type
            switch lower(argument)
                case 'reserve'
                    output{3.0} = 1.0;
                case 'value'
                    output{3.0} = 2.0;
                otherwise
                    output{3.0} = 1.0;
            end

            % Process the third argument
            argument = OUTPUT_OPTIMISED{3.0};

            switch lower(argument)
                case 'minmax'
                    output{4.0} = 1.0;
                case 'minmean'
                    output{4.0} = 2.0;
                otherwise
                    output{4.0} = 1.0;
            end

            % Process the fourth argument
            argument = OUTPUT_OPTIMISED{4.0};

            if (argument <= 0.0) || (argument > 90.0)
                fprintf('[ERROR] Invalid value of OUTPUT_OPTIMISED(4). The angular step size must be in the range {0 < theta <= 90}\n');

                % Reset the error flag and RETURN
                error = true;
                return
            else
                output{5.0} = linspace(0.0, 90.0, 1.0 + 90.0/argument);
            end
        end

        %% INDEX-TO-SEQUENCE CONVERTER
        function seq = indexToStacking(index, angles, numAngles, nPlies)
            seq = zeros(1.0, nPlies);
            idx = index - 1.0;

            for i = nPlies:-1.0:1.0
                angleIndex = mod(idx, numAngles) + 1.0;
                seq(i) = angles(angleIndex);
                idx = floor(idx / numAngles);
            end
        end

        %% PARFOR LOOP FOR METHOD 1: FULL MATRIX
        function [CRITERION_BUFFER] = parfor1(nPermutations, anglePermutations, nSectionPoints, nPlies, z, dummy, tolerance, Q11, Q12, Q66, Q22, A11_points, A22_points,...
                B11_points, B22_points, nargin, deltaT, deltaM, Nxx, Nyy, Nxy, Mxx, Myy, Mxy, nPlies_points, z_points, failureCriterion, XT, XC, YT, YC, S, parameter, C12, B12,...
                E11, E22, V12, G12, XET, XEC, YET, YEC, SE, ALPHA, XHT, XHC, YHT, YHC, SHX, SHY, symsAvailable, S1, S2, S3, GL12, XLT, XLC, YLT, YLC, SLX, SLY, A0, PHI0, NL, NT,...
                SECTION_POINTS, objective)
            parfor i = 1.0:nPermutations                
                % Get the current stacking order
                theta = anglePermutations(i, :);

                % Get the values of theta over the section points
                [~, ~, theta_points, ~, ~, ~, ~, ~, ~, ~] =...
                    ...
                    abd.internal_getSectionPoints(nSectionPoints, '', nPlies, theta, z, dummy, dummy, dummy, dummy, tolerance, []);

                % Compute transformed reduced stiffness matrix components
                [Q11t, Q12t, Q16t, Q22t, Q26t, Q66t] =...
                    ...
                    abd.internal_getTransformedQ(nPlies, theta, Q11, Q12, Q66, Q22);

                %{
                    Get effective thermal and moisture expansion
                    coefficients for each ply
                %}
                [axx, ayy, axy, bxx, byy, bxy] =...
                    ...
                    abd.internal_getThermoHydro(theta_points, A11_points, A22_points, B11_points, B22_points);

                % Compute A, B and D matrices
                [ABD, ~, Qijt, NxxT, NyyT, NxyT, MxxT, MyyT, MxyT, NxxM, NyyM, NxyM, MxxM, MyyM, MxyM] =...
                    ...
                    abd.internal_getABD(nPlies, Q11t, Q12t, Q16t, Q22t, Q26t, Q66t, z, nargin, deltaT, deltaM, axx, ayy, axy, bxx, byy, bxy, nSectionPoints);

                % Compute tensor quantities
                [~, ~, ~, ~, stress, ~, ~, ~, ~] =...
                    ...
                    abd.internal_getTensor(ABD, Nxx, NxxT, NxxM, Nyy, NyyT, NyyM, Nxy, NxyT, NxyM, Mxx, MxxT, MxxM, Myy, MyyT, MyyM, Mxy, MxyT, MxyM, nPlies_points, z_points,...
                    theta_points, Qijt, deltaT, deltaM, axx, ayy, axy, bxx, byy, bxy, tolerance);

                % Perform strength calculation on ply stresses
                switch lower(failureCriterion)
                    case 'mstrs' % Maximum stress
                        CRITERION =...
                            ...
                            abd.internal_strength.getMstrs(stress, XT, XC, YT, YC, S);
                    case 'tsaih' % Tsai-Hill
                        CRITERION =...
                            ...
                            abd.internal_strength.getTsaih(parameter, stress, XT, XC, YT, YC, S);
                    case 'tsaiw' % Tsai-Wu
                        CRITERION =...
                            ...
                            abd.internal_strength.getTsaiw(parameter, stress, XT, XC, YT, YC, S, C12, B12);
                    case 'azzit' % Azzi-Tsai-Hill
                        CRITERION =...
                            ...
                            abd.internal_strength.getAzzit(parameter, nPlies_points, stress, XT, XC, YT, YC, S);
                    case 'mstrn' % Maximum strain
                        CRITERION =...
                            ...
                            abd.internal_strength.getMstrn(nPlies_points, stress, E11, E22, V12, G12, XET, XEC, YET, YEC, SE);
                    case 'hashin' % Hashin
                        [HSNFTCRT, HSNFCCRT, HSNMTCRT, HSNMCCRT] =...
                            ...
                            abd.internal_strength.getHashin(nPlies_points, stress, ALPHA, XHT, XHC, YHT, YHC, SHX, SHY);

                        % Get the worst criterion of all four calculations
                        CRITERION = max([HSNFTCRT; HSNFCCRT; HSNMTCRT; HSNMCCRT], [], 1.0);
                    case 'larc05' % LaRC05
                        [LARPFCRT, LARMFCRT, LARKFCRT, LARSFCRT, LARTFCRT] = ...
                            ...
                            abd.internal_getLaRC05(nPlies_points, stress, symsAvailable, S1, S2, S3, GL12, XLT, XLC, YLT, YLC, SLX, SLY, A0, PHI0, NL, NT, SECTION_POINTS);

                        %{
                            Reset large LARMFCT values to one to prevent
                            very large values from breaking the CB MATLAB
                            figure
                        %}
                        LARMFCRT(LARMFCRT > 1.0) = 1.0;

                        % Get the worst criterion of all four calculations
                        CRITERION = max([LARPFCRT, LARMFCRT, LARKFCRT, LARSFCRT, LARTFCRT], [], 1.0);
                    otherwise
                        % Default to Tsai-Hill
                        CRITERION =...
                            ...
                            abd.internal_strength.getTsaih(parameter, stress, XT, XC, YT, YC, S);
                end

                % Get the worst value of CRITERION based on objective
                if objective == 1.0
                    % Objective function: MinMax
                    CRITERION_BUFFER(i) = max(CRITERION);
                else
                    % Objective function: MinMean
                    CRITERION_BUFFER(i) = mean(CRITERION);
                end
            end
        end

        %% PARFOR LOOP FOR METHOD 2: MIXED-RADIX REPRESENTATION
        function [CRITERION_BUFFER] = parfor2(nPermutations, thetaAll, numAngles, nSectionPoints, nPlies, z, dummy, tolerance, Q11, Q12, Q66, Q22, A11_points, A22_points,...
                B11_points, B22_points, nargin, deltaT, deltaM, Nxx, Nyy, Nxy, Mxx, Myy, Mxy, nPlies_points, z_points, failureCriterion, XT, XC, YT, YC, S, parameter, C12, B12,...
                E11, E22, V12, G12, XET, XEC, YET, YEC, SE, ALPHA, XHT, XHC, YHT, YHC, SHX, SHY, symsAvailable, S1, S2, S3, GL12, XLT, XLC, YLT, YLC, SLX, SLY, A0, PHI0, NL, NT,...
                SECTION_POINTS, objective)
            parfor i = 1.0:nPermutations
                % Get the current stacking order
                theta = abd.internal_optimise.indexToStacking(i, thetaAll, numAngles, nPlies);

                % Get the values of theta over the section points
                [~, ~, theta_points, ~, ~, ~, ~, ~, ~, ~] =...
                    ...
                    abd.internal_getSectionPoints(nSectionPoints, '', nPlies, theta, z, dummy, dummy, dummy, dummy, tolerance, []);

                % Compute transformed reduced stiffness matrix components
                [Q11t, Q12t, Q16t, Q22t, Q26t, Q66t] =...
                    ...
                    abd.internal_getTransformedQ(nPlies, theta, Q11, Q12, Q66, Q22);

                %{
                    Get effective thermal and moisture expansion
                    coefficients for each ply
                %}
                [axx, ayy, axy, bxx, byy, bxy] =...
                    ...
                    abd.internal_getThermoHydro(theta_points, A11_points, A22_points, B11_points, B22_points);

                % Compute A, B and D matrices
                [ABD, ~, Qijt, NxxT, NyyT, NxyT, MxxT, MyyT, MxyT, NxxM, NyyM, NxyM, MxxM, MyyM, MxyM] =...
                    ...
                    abd.internal_getABD(nPlies, Q11t, Q12t, Q16t, Q22t, Q26t, Q66t, z, nargin, deltaT, deltaM, axx, ayy, axy, bxx, byy, bxy, nSectionPoints);

                % Compute tensor quantities
                [~, ~, ~, ~, stress, ~, ~, ~, ~] =...
                    ...
                    abd.internal_getTensor(ABD, Nxx, NxxT, NxxM, Nyy, NyyT, NyyM, Nxy, NxyT, NxyM, Mxx, MxxT, MxxM, Myy, MyyT, MyyM, Mxy, MxyT, MxyM, nPlies_points, z_points,...
                    theta_points, Qijt, deltaT, deltaM, axx, ayy, axy, bxx, byy, bxy, tolerance);

                % Perform strength calculation on ply stresses
                switch lower(failureCriterion)
                    case 'mstrs' % Maximum stress
                        CRITERION =...
                            ...
                            abd.internal_strength.getMstrs(stress, XT, XC, YT, YC, S);
                    case 'tsaih' % Tsai-Hill
                        CRITERION =...
                            ...
                            abd.internal_strength.getTsaih(parameter, stress, XT, XC, YT, YC, S);
                    case 'tsaiw' % Tsai-Wu
                        CRITERION =...
                            ...
                            abd.internal_strength.getTsaiw(parameter, stress, XT, XC, YT, YC, S, C12, B12);
                    case 'azzit' % Azzi-Tsai-Hill
                        CRITERION =...
                            ...
                            abd.internal_strength.getAzzit(parameter, nPlies_points, stress, XT, XC, YT, YC, S);
                    case 'mstrn' % Maximum strain
                        CRITERION =...
                            ...
                            abd.internal_strength.getMstrn(nPlies_points, stress, E11, E22, V12, G12, XET, XEC, YET, YEC, SE);
                    case 'hashin' % Hashin
                        [HSNFTCRT, HSNFCCRT, HSNMTCRT, HSNMCCRT] =...
                            ...
                            abd.internal_strength.getHashin(nPlies_points, stress, ALPHA, XHT, XHC, YHT, YHC, SHX, SHY);

                        % Get the worst criterion of all four calculations
                        CRITERION = max([HSNFTCRT; HSNFCCRT; HSNMTCRT; HSNMCCRT], [], 1.0);
                    case 'larc05' % LaRC05
                        [LARPFCRT, LARMFCRT, LARKFCRT, LARSFCRT, LARTFCRT] = ...
                            ...
                            abd.internal_getLaRC05(nPlies_points, stress, symsAvailable, S1, S2, S3, GL12, XLT, XLC, YLT, YLC, SLX, SLY, A0, PHI0, NL, NT, SECTION_POINTS);

                        %{
                            Reset large LARMFCT values to one to prevent
                            very large values from breaking the CB MATLAB
                            figure
                        %}
                        LARMFCRT(LARMFCRT > 1.0) = 1.0;

                        % Get the worst criterion of all four calculations
                        CRITERION = max([LARPFCRT, LARMFCRT, LARKFCRT, LARSFCRT, LARTFCRT], [], 1.0);
                    otherwise
                        % Default to Tsai-Hill
                        CRITERION =...
                            ...
                            abd.internal_strength.getTsaih(parameter, stress, XT, XC, YT, YC, S);
                end

                % Get the worst value of CRITERION based on objective
                if objective == 1.0
                    % Objective function: MinMax
                    CRITERION_BUFFER(i) = max(CRITERION);
                else
                    % Objective function: MinMean
                    CRITERION_BUFFER(i) = mean(CRITERION);
                end
            end
        end

        %% PARFOR LOOP FOR METHOD 3: CHUNKS
        function [CRITERION_BUFFER] = parfor3(numChunks, CHUNK_SIZE, nPermutations, thetaAll, numAngles, nSectionPoints, nPlies, z, dummy, tolerance, Q11, Q12, Q66, Q22,...
                A11_points, A22_points, B11_points, B22_points, nargin, deltaT, deltaM, Nxx, Nyy, Nxy, Mxx, Myy, Mxy, nPlies_points, z_points, failureCriterion, XT, XC, YT, YC, S,...
                parameter, C12, B12, E11, E22, V12, G12, XET, XEC, YET, YEC, SE, ALPHA, XHT, XHC, YHT, YHC, SHX, SHY, symsAvailable, S1, S2, S3, GL12, XLT, XLC, YLT, YLC, SLX,...
                SLY, A0, PHI0, NL, NT, SECTION_POINTS, objective)
            % Initialise the outer buffer for the PARFOR loop
            OUTER_BUFFER = cell(numChunks, 1.0);

            parfor chunkId = 1.0:numChunks
                % Get iterator parameters
                startIdx = (chunkId - 1.0) * CHUNK_SIZE + 1.0;
                endIdx = min(chunkId * CHUNK_SIZE, nPermutations);
                numItems = endIdx - startIdx + 1.0;

                % Initialise the inner buffer for the FOR loop
                INNER_BUFFER = cell(numItems, 1.0);

                for i = startIdx:endIdx
                    % Get the current stacking order
                    theta = abd.internal_optimise.indexToStacking(i, thetaAll, numAngles, nPlies);

                    % Get the values of theta over the section points
                    [~, ~, theta_points, ~, ~, ~, ~, ~, ~, ~] =...
                        ...
                        abd.internal_getSectionPoints(nSectionPoints, '', nPlies, theta, z, dummy, dummy, dummy, dummy, tolerance, []);

                    % Compute transformed reduced stiffness matrix components
                    [Q11t, Q12t, Q16t, Q22t, Q26t, Q66t] =...
                        ...
                        abd.internal_getTransformedQ(nPlies, theta, Q11, Q12, Q66, Q22);

                    %{
                        Get effective thermal and moisture expansion
                        coefficients for each ply
                    %}
                    [axx, ayy, axy, bxx, byy, bxy] =...
                        ...
                        abd.internal_getThermoHydro(theta_points, A11_points, A22_points, B11_points, B22_points);

                    % Compute A, B and D matrices
                    [ABD, ~, Qijt, NxxT, NyyT, NxyT, MxxT, MyyT, MxyT, NxxM, NyyM, NxyM, MxxM, MyyM, MxyM] =...
                        ...
                        abd.internal_getABD(nPlies, Q11t, Q12t, Q16t, Q22t, Q26t, Q66t, z, nargin, deltaT, deltaM, axx, ayy, axy, bxx, byy, bxy, nSectionPoints);

                    % Compute tensor quantities
                    [~, ~, ~, ~, stress, ~, ~, ~, ~] =...
                        ...
                        abd.internal_getTensor(ABD, Nxx, NxxT, NxxM, Nyy, NyyT, NyyM, Nxy, NxyT, NxyM, Mxx, MxxT, MxxM, Myy, MyyT, MyyM, Mxy, MxyT, MxyM, nPlies_points, z_points,...
                        theta_points, Qijt, deltaT, deltaM, axx, ayy, axy, bxx, byy, bxy, tolerance);

                    % Perform strength calculation on ply stresses
                    switch lower(failureCriterion)
                        case 'mstrs' % Maximum stress
                            CRITERION =...
                                ...
                                abd.internal_strength.getMstrs(stress, XT, XC, YT, YC, S);
                        case 'tsaih' % Tsai-Hill
                            CRITERION =...
                                ...
                                abd.internal_strength.getTsaih(parameter, stress, XT, XC, YT, YC, S);
                        case 'tsaiw' % Tsai-Wu
                            CRITERION =...
                                ...
                                abd.internal_strength.getTsaiw(parameter, stress, XT, XC, YT, YC, S, C12, B12);
                        case 'azzit' % Azzi-Tsai-Hill
                            CRITERION =...
                                ...
                                abd.internal_strength.getAzzit(parameter, nPlies_points, stress, XT, XC, YT, YC, S);
                        case 'mstrn' % Maximum strain
                            CRITERION =...
                                ...
                                abd.internal_strength.getMstrn(nPlies_points, stress, E11, E22, V12, G12, XET, XEC, YET, YEC, SE);
                        case 'hashin' % Hashin
                            [HSNFTCRT, HSNFCCRT, HSNMTCRT, HSNMCCRT] =...
                                ...
                                abd.internal_strength.getHashin(nPlies_points, stress, ALPHA, XHT, XHC, YHT, YHC, SHX, SHY);

                            % Get the worst criterion of all four calculations
                            CRITERION = max([HSNFTCRT; HSNFCCRT; HSNMTCRT; HSNMCCRT], [], 1.0);
                        case 'larc05' % LaRC05
                            [LARPFCRT, LARMFCRT, LARKFCRT, LARSFCRT, LARTFCRT] = ...
                                ...
                                abd.internal_getLaRC05(nPlies_points, stress, symsAvailable, S1, S2, S3, GL12, XLT, XLC, YLT, YLC, SLX, SLY, A0, PHI0, NL, NT, SECTION_POINTS);

                            %{
                                Reset large LARMFCT values to one to
                                prevent very large values from breaking the
                                CB MATLAB figure
                            %}
                            LARMFCRT(LARMFCRT > 1.0) = 1.0;

                            % Get the worst criterion of all four calculations
                            CRITERION = max([LARPFCRT, LARMFCRT, LARKFCRT, LARSFCRT, LARTFCRT], [], 1.0);
                        otherwise
                            % Default to Tsai-Hill
                            CRITERION =...
                                ...
                                abd.internal_strength.getTsaih(parameter, stress, XT, XC, YT, YC, S);
                    end

                    % Get the worst value of CRITERION based on objective
                    if objective == 1.0
                        % Objective function: MinMax
                        INNER_BUFFER{i} = max(CRITERION);
                    else
                        % Objective function: MinMean
                        INNER_BUFFER{i} = mean(CRITERION);
                    end
                end

                % Update the outer loop buffer
                OUTER_BUFFER{chunkId} = INNER_BUFFER(startIdx:endIdx);
            end

            % Collect outer buffers into function output
            CRITERION_BUFFER = cell2mat(vertcat(OUTER_BUFFER{:})');
        end
    end
end
