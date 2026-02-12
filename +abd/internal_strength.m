classdef internal_strength < handle
%   Perform strength calculations based on the ply stresses.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 5.0.0 Copyright Louis Vallance 2026
%   Last modified 11-Feb-2026 08:06:52 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%

    methods(Static = true, Access = public)
        %% MAIN FUNCTION FOR STRENGTH CALCULATION
        function [MSTRS, TSAIH, HOFFMAN, TSAIW, AZZIT, MSTRN, HSNFTCRT, HSNFCCRT, HSNMTCRT, HSNMCCRT, LARPFCRT, LARMFCRT, LARKFCRT, LARSFCRT, LARTFCRT, UCRT, XT, XC, YT, YC, S,...
                C12, B12, E11, E22, G12, V12, XET, XEC, YET, YEC, SE, ALPHA, XHT, XHC, YHT, YHC, SHX, SHY, XLT, XLC, YLT, YLC, SLX, SLY, GL12, NL, NT, A0, PHI0, S1, S2, S3,...
                UCRT_MException]...
                =...
                main(noFailStress, noFailStrain, noHashin, noLaRC05, symsAvailable, XT, XC, YT, YC, S, C12, B12, E11, E22, G12, V12, AXX, AYY, AXY, BXX, BYY, BXY, XET, XEC, YET,...
                YEC, SE, ALPHA, XHT, XHC, YHT, YHC, SHX, SHY, XLT, XLC, YLT, YLC, SLX, SLY, GL12, NL, NT, A0, PHI0, TENSORS, nPlies, nPlies_points, SECTION_POINTS, fcnHandle,...
                parameter, MSTRS, TSAIH, HOFFMAN, TSAIW, AZZIT, MSTRN, HSNFTCRT, HSNFCCRT, HSNMTCRT, HSNMCCRT, LARPFCRT, LARMFCRT, LARKFCRT, LARSFCRT, LARTFCRT, UCRT)
            % Initialise the output
            S1 = [];    S2 = [];    S3 = [];    UCRT_MException = [];

            % Flag for user-defined failure criterion
            userCrit = true;

            % Get the stress components for criteria
            stress = TENSORS.S_PLY_ALIGNED;

            if userCrit == true
                E11 = abd.internal_spreadProperties(E11, nPlies, SECTION_POINTS);
                E22 = abd.internal_spreadProperties(E22, nPlies, SECTION_POINTS);
                G12 = abd.internal_spreadProperties(G12, nPlies, SECTION_POINTS);
                V12 = abd.internal_spreadProperties(V12, nPlies, SECTION_POINTS);
            end

            % Spread material data over section points
            if noFailStress == false
                XT = abd.internal_spreadProperties(XT, nPlies, SECTION_POINTS);
                XC = abd.internal_spreadProperties(XC, nPlies, SECTION_POINTS);
                YT = abd.internal_spreadProperties(YT, nPlies, SECTION_POINTS);
                YC = abd.internal_spreadProperties(YC, nPlies, SECTION_POINTS);
                S = abd.internal_spreadProperties(S, nPlies, SECTION_POINTS);
                C12 = abd.internal_spreadProperties(C12, nPlies, SECTION_POINTS);
                B12 = abd.internal_spreadProperties(B12, nPlies, SECTION_POINTS);
            end

            % Spread material data over section points
            if noFailStrain == false
                E11 = abd.internal_spreadProperties(E11, nPlies, SECTION_POINTS);
                E22 = abd.internal_spreadProperties(E22, nPlies, SECTION_POINTS);
                G12 = abd.internal_spreadProperties(G12, nPlies, SECTION_POINTS);
                V12 = abd.internal_spreadProperties(V12, nPlies, SECTION_POINTS);
                XET = abd.internal_spreadProperties(XET, nPlies, SECTION_POINTS);
                XEC = abd.internal_spreadProperties(XEC, nPlies, SECTION_POINTS);
                YET = abd.internal_spreadProperties(YET, nPlies, SECTION_POINTS);
                YEC = abd.internal_spreadProperties(YEC, nPlies, SECTION_POINTS);
                SE = abd.internal_spreadProperties(SE, nPlies, SECTION_POINTS);
            end

            % Spread material data over section points
            if noHashin == false
                ALPHA = abd.internal_spreadProperties(ALPHA, nPlies, SECTION_POINTS);
                XHT = abd.internal_spreadProperties(XHT, nPlies, SECTION_POINTS);
                XHC = abd.internal_spreadProperties(XHC, nPlies, SECTION_POINTS);
                YHT = abd.internal_spreadProperties(YHT, nPlies, SECTION_POINTS);
                YHC = abd.internal_spreadProperties(YHC, nPlies, SECTION_POINTS);
                SHX = abd.internal_spreadProperties(SHX, nPlies, SECTION_POINTS);
                SHY = abd.internal_spreadProperties(SHY, nPlies, SECTION_POINTS);
            end

            % Spread material data over section points
            if noLaRC05 == false
                XLT = abd.internal_spreadProperties(XLT, nPlies, SECTION_POINTS);
                XLC = abd.internal_spreadProperties(XLC, nPlies, SECTION_POINTS);
                YLT = abd.internal_spreadProperties(YLT, nPlies, SECTION_POINTS);
                YLC = abd.internal_spreadProperties(YLC, nPlies, SECTION_POINTS);
                SLX = abd.internal_spreadProperties(SLX, nPlies, SECTION_POINTS);
                SLY = abd.internal_spreadProperties(SLY, nPlies, SECTION_POINTS);
                GL12 = abd.internal_spreadProperties(GL12, nPlies, SECTION_POINTS);
                NL = abd.internal_spreadProperties(NL, nPlies, SECTION_POINTS);
                NT = abd.internal_spreadProperties(NT, nPlies, SECTION_POINTS);
                A0 = abd.internal_spreadProperties(A0, nPlies, SECTION_POINTS);
                PHI0 = abd.internal_spreadProperties(PHI0, nPlies, SECTION_POINTS);

                % GET PRINCIPAL STRESS TERMS FOR LARC05
                % Get the individual tensor components from S_PLY_ALIGNED
                S11 = stress(1.0, :);
                S22 = stress(2.0, :);
                S12 = stress(3.0, :);

                % % Use Eigenvalues
                % S = [S11(1), S12(1), 0; 0, S22(1), 0; 0, 0, 0];
                % I = eig(S);
                % S1 = max(I);
                % S2 = median(I);
                % S3 = min(I);

                % Get the two in-plane principal stress components
                S1 = 0.5.*(S11 + S22) + sqrt((0.5.*(S11 - S22)).^2.0 + S12.^2.0);
                S2 = 0.5.*(S11 + S22) - sqrt((0.5.*(S11 - S22)).^2.0 + S12.^2.0);
                S3 = zeros(1.0, nPlies_points);
            end

            if noFailStress == false
                % Failure calculation: MSTRS
                MSTRS =...
                    ...
                    abd.internal_strength.getMstrs(stress, XT, XC, YT, YC, S);

                % Failure calculation: TSAIH
                TSAIH =...
                    ...
                    abd.internal_strength.getTsaih(parameter, stress, XT, XC, YT, YC, S);

                % Failure calculation: HOFFMAN
                HOFFMAN =...
                    ...
                    abd.internal_strength.getHoffman(parameter, stress, XT, XC, YT, YC, S, zeros(1.0, nPlies_points), linspace(-1.0, -1.0, nPlies_points));

                % Failure calculation: TSAIW
                TSAIW =...
                    ...
                    abd.internal_strength.getTsaiw(parameter, stress, XT, XC, YT, YC, S, C12, B12);

                % Failure calculation: AZZIT
                AZZIT =...
                    ...
                    abd.internal_strength.getAzzit(parameter, nPlies_points, stress, XT, XC, YT, YC, S);
            end

            % Failure calculation: MSTRN
            if noFailStrain == false
                MSTRN = ...
                    ...
                    abd.internal_strength.getMstrn(nPlies_points, stress, E11, E22, V12, G12, XET, XEC, YET, YEC, SE);
            end

            % Failure calculation: HASHIN
            if noHashin == false
                [HSNFTCRT, HSNFCCRT, HSNMTCRT, HSNMCCRT] = ...
                    ...
                    abd.internal_strength.getHashin(nPlies_points, stress, ALPHA, XHT, XHC, YHT, YHC, SHX, SHY);
            end

            % Failure calculation: LARC05
            if noLaRC05 == false
                [LARPFCRT, LARMFCRT, LARKFCRT, LARSFCRT, LARTFCRT] = ...
                    ...
                    abd.internal_getLaRC05(nPlies_points, stress, symsAvailable, S1, S2, S3, GL12, XLT, XLC, YLT, YLC, SLX, SLY, A0, PHI0, NL, NT, SECTION_POINTS);
            end

            % Failure calculation: UCRT
            if isa(fcnHandle, 'function_handle') == true
                INFO = struct('PLY_COUNT', nPlies, 'SECTION_POINTS', SECTION_POINTS, 'TOTAL_POINTS', nPlies_points, 'FAILURE_PARAMETER', parameter);
                MATERIAL_MECH = struct('E11', E11, 'E22', E22, 'G12', G12, 'V12', V12, 'AXX', AXX, 'AYY', AYY, 'AXY', AXY, 'BXX', BXX, 'BYY', BYY, 'BXY', BXY);
                MATERIAL_FAIL = struct('STRESS', struct('XT', XT, 'XC', XC, 'YT', YT, 'YC', YC, 'S', S, 'C12', C12, 'B12', B12),...
                    'STRAIN', struct('XET', XET, 'XEC', XEC, 'YET', YET, 'YEC', YEC, 'SE', SE),...
                    'HASHIN', struct('ALPHA', ALPHA, 'XHT', XHT, 'XHC', XHC, 'YHT', YHT, 'YHC', YHC, 'SHX', SHX, 'SHY', SHY),...
                    'LARC05', struct('XLT', XLT, 'XLC', XLC, 'YLT', YLT, 'YLC', YLC, 'SLX', SLX, 'SLY', SLY, 'GL12', GL12, 'NL', NL, 'NT', NT, 'A0', A0, 'PHI0', PHI0)...
                    );

                % Save the initial state of UCRT
                UCRT_initial = UCRT;

                % Run the user-defined failure criterion
                [UCRT, UCRT_MException] = abd.internal_strength.validateUserRoutine(fcnHandle, INFO, UCRT, MATERIAL_MECH, MATERIAL_FAIL, TENSORS, nPlies_points, UCRT_initial);

                % Ensure UCRT is 1xTOTAL_POINTS
                UCRT = UCRT(:).';
            end
        end

        %% GET DATA FROM OUTPUT_STRENGTH
        function [error, output] = getSettings(OUTPUT_STRENGTH)
            % Initialise output
            error = false;
            output = cell(1.0, 2.0);

            if iscell(OUTPUT_STRENGTH) == false
                % Convert to cell if necessary
                OUTPUT_STRENGTH = {OUTPUT_STRENGTH};
            end

            if cellfun(@isempty, OUTPUT_STRENGTH) == true
                % Set default values if necessary
                OUTPUT_STRENGTH = {false, 'RESERVE'};
            end

            if length(OUTPUT_STRENGTH) ~= 2.0
                % Incorrect number of arguments
                fprintf('[ERROR] The setting OUTPUT_STRENGTH requires two\narguments:\n{{[false | true] | [@<ucrt> | ''<file-name>'']}, ''<parameter>''}\n');

                % Reset the error flag and RETURN
                error = true;
                return
            end

            % Process the first argument
            argument = OUTPUT_STRENGTH{1.0};

            % Check validity of the argument
            if isempty(argument) == true
                output{1.0} = false;
            elseif (ischar(argument) == true) || (isa(argument, 'function_handle') == true)
                %{
                    The user has specified a string or a function handle.
                    If it is a file, convert it to a function handle;
                    otherwise, assume it's a function handle
                %}

                % Check if argument is a function handle
                if isa(argument, 'function_handle') == false
                    % Get function name from file path
                    [path, name, ext] = fileparts(which(argument));

                    if isempty(name) == true
                        if (isempty(path) == true) && (isempty(ext) == true) && (isempty(argument) == false) && (strcmp(argument(1.0), '@') == true)
                            %{
                                The user may have mistakenly specified a
                                function handle as a string. If so, try to
                                convert back to a valid handle
                            %}
                            argument = str2func(argument);
                        else
                            % The user routine could not be found
                            fprintf(['[ERROR] OUTPUT_STRENGTH(1) was specified as a string, but it does not match any file on the MATLAB path. The value must be a logical or a valid f',...
                                'unction handle or file name:\n{[false | true] | [@<ucrt> | ''<file-name>'']}\n']);

                            % Reset the error flag and RETURN
                            error = true;
                            return
                        end
                    else
                        % Convert to function handle
                        argument = str2func(name);
                    end
                end

                % Make sure function handle points to valid M-file
                fcnInfo = functions(argument);

                % Check for function file existence
                if (isempty(fcnInfo.file) == true) || (exist(fcnInfo.file, 'file') ~= 2.0)
                    % Incorrect argument type
                    fprintf('[ERROR] OUTPUT_STRENGTH(1) must be a valid function handle or file name:\n{@<ucrt> | ''<file-name>''}\n');

                    % Reset the error flag and RETURN
                    error = true;
                    return
                end

                if (nargin(argument) ~= 5.0) || (nargout(argument) ~= 1.0)
                    % Incorrect number of input/output arguments
                    fprintf(['[ERROR] In user routine @%s, the number of input/output arguments is incorrect: Found %.0f inputs (expected %.0f); Found %.0f outputs (expected %',...
                        '.0f)\n        The user-defined failure criterion uses the following function interface:\n        \tfunction [UCRT] = <function_handle>(INFO, UCRT, MATERIA',...
                        'L_MECH, MATERIAL_FAIL, TENSORS)\n\n        Note: In case of doubt, use the command >> abd.internal_createUcrt(''<file-name>'') to create a template user r',...
                        'outine file\n'], char(argument), nargin(argument), 5.0, nargout(argument), 1.0);

                    % Reset the error flag and RETURN
                    error = true;
                    return
                end

                % Everything is OK
                output{1.0} = argument;
            elseif (islogical(argument) == false) && (argument ~= 0.0) && (argument ~= 1.0)
                % Incorrect argument type
                fprintf('[ERROR] OUTPUT_STRENGTH(1) must be a logical:\n{false | true}\n');

                % Reset the error flag and RETURN
                error = true;
                return
            else
                % Everything is OK
                output{1.0} = argument;
            end

            % Process the second argument
            argument = OUTPUT_STRENGTH{2.0};

            % Get the parameter type
            switch lower(argument)
                case 'reserve'
                    output{2.0} = 1.0;
                case 'value'
                    output{2.0} = 2.0;
                otherwise
                    output{2.0} = 1.0;
            end
        end

        %% GET COLOUR MAP FOR FAILED SETION POINT VISUALISATION
        function [SP_COLOUR_BUFFER] = getFailedSpBuffer(MSTRS, TSAIH, HOFFMAN, TSAIW, AZZIT, MSTRN, HSNFTCRT, HSNFCCRT, HSNMTCRT, HSNMCCRT, LARPFCRT, LARMFCRT, LARKFCRT, LARSFCRT,...
                LARTFCRT, UCRT, SP_COLOUR_BUFFER, nPlies_points)
            % Collect the results of the strength analysis into a cell array
            STRENGTH_RESULTS = {MSTRS >= 1.0, TSAIH >= 1.0, HOFFMAN >= 1.0, TSAIW >= 1.0, AZZIT >= 1.0, MSTRN >= 1.0, HSNFTCRT >= 1.0, HSNFCCRT >= 1.0, HSNMTCRT >= 1.0,...
                HSNMCCRT >= 1.0, LARPFCRT >= 1.0, LARMFCRT >= 1.0, LARKFCRT >= 1.0, LARSFCRT >= 1.0, LARTFCRT >= 1.0, UCRT >= 1.0};

            %{
                Initialize failed section points buffer with logical FALSE
                values
            %}
            FAILED_SP_BUFFER = false(1.0, nPlies_points);

            %{
                For each section point, check if the value for each
                strength criterion is equal to or greater than 1
            %}
            % Looping over section points
            for i = 1:nPlies_points
                % Looping over strength criteria
                for j = 1.0:numel(STRENGTH_RESULTS)
                    % Check the current result
                    if STRENGTH_RESULTS{j}(i) == 1.0
                        % Reset the value in the failed section point buffer
                        FAILED_SP_BUFFER(i) = true;

                        % No need to check further if one value > 1
                        break;
                    end
                end
            end

            % Assign colors based on BUFFER (Red -> true; Green -> false)
            SP_COLOUR_BUFFER(FAILED_SP_BUFFER, :) = repmat([1.0, 0.0, 0.0], sum(FAILED_SP_BUFFER), 1.0);
            SP_COLOUR_BUFFER(FAILED_SP_BUFFER == false, :) = repmat([0.0, 1.0, 0.0], sum(FAILED_SP_BUFFER == false), 1.0);
        end

        %% INITIALISE FAILURE CRITERIA BUFFERS
        function [MSTRS, TSAIH, HOFFMAN, TSAIW, AZZIT, MSTRN, HSNFTCRT, HSNFCCRT, HSNMTCRT, HSNMCCRT, LARPFCRT, LARMFCRT, LARKFCRT, LARSFCRT, LARTFCRT, UCRT] = init(nPlies_points)
            % Initialise output
            MSTRS = linspace(-1.0, -1.0, nPlies_points);
            TSAIH = MSTRS;
            HOFFMAN = MSTRS;
            TSAIW = MSTRS;
            AZZIT = MSTRS;
            MSTRN = MSTRS;
            HSNFTCRT = MSTRS;
            HSNFCCRT = MSTRS;
            HSNMTCRT = MSTRS;
            HSNMCCRT = MSTRS;
            LARPFCRT = MSTRS;
            LARMFCRT = MSTRS;
            LARKFCRT = MSTRS;
            LARSFCRT = MSTRS;
            LARTFCRT = MSTRS;
            UCRT = MSTRS;
        end

        %% VALIDATE OUTPUT OF USER ROUTINE
        function [UCRT, UCRT_MException] = validateUserRoutine(fcnHandle, INFO, UCRT, MATERIAL_MECH, MATERIAL_FAIL, TENSORS, nPlies_points, UCRT_initial)
            % Initialise output
            UCRT_MException = [];
            
            % Try to run the user-defined failure criterion
            try
                % Inform the user that the routine started
                fprintf('[NOTICE] Start user routine @%s\n', char(fcnHandle));

                % Run the user routine
                UCRT = fcnHandle(INFO, UCRT, MATERIAL_MECH, MATERIAL_FAIL, TENSORS);
            catch UCRT_MException
                %{
                    Do not evaluate the user-defined failure criterion and
                    save the exception object to the workspace
                %}
                try
                    fprintf('[ERROR] Exception encountered on line %.0f in user routine file ''%s.m''\n', UCRT_MException.stack(1.0).line, char(fcnHandle));
                    fid = fopen(UCRT_MException.stack(1.0).file, 'r');

                    % Extract the line image
                    for k = 1:UCRT_MException.stack(1.0).line
                        line = fgetl(fid);
                        if ischar(line) == false
                            % Do nothing
                            break
                        end
                    end

                    % Print the line image to the command window
                    fclose(fid);
                    fprintf('        Line image: %s\n', line);
                catch
                    % Do nothing
                end

                fprintf('        MException identifier: %s\n', UCRT_MException.identifier);
                fprintf('        MException message: %s\n', UCRT_MException.message);
                fprintf('[NOTICE] The complete MATLAB Exception object will be saved in the specified output location\n');
                fprintf('[ERROR] The user-defined failure criterion has not been evaluated\n');

                % Reset the values of UCRT
                UCRT = UCRT_initial;
            end

            if (isnumeric(UCRT) == false) || (iscell(UCRT) == true)
                % Validity check
                fprintf('[ERROR] In user routine @%s, UCRT must be numeric. The user-defined failure criterion will be excluded from the output\n', char(fcnHandle));

                % Reset the values of UCRT
                UCRT = UCRT_initial;
            elseif isvector(UCRT) == false
                % Dimension check
                fprintf('[ERROR] In user routine @%s, UCRT must be one-dimensional. The user-defined failure criterion will be excludced from the output\n', char(fcnHandle));

                % Reset the values of UCRT
                UCRT = UCRT_initial;
            elseif all(UCRT == -1.0) == true
                % Validity check
                fprintf('[ERROR] In user routine @%s, UCRT values are all unset (-1). The user-defined failure criterion will be excluded from the output\n', char(fcnHandle));

                % Reset the values of UCRT
                UCRT = UCRT_initial;
            elseif numel(UCRT) ~= nPlies_points
                % NUMEL check
                fprintf('[ERROR] In user routine @%s, UCRT contains %d elements (expected %d). The user-defined failure criterion will be excludced from the output\n',...
                    char(fcnHandle), numel(UCRT), nPlies_points);

                % Reset the values of UCRT
                UCRT = UCRT_initial;
            elseif any(UCRT < 0.0) == true
                % Negative value check
                fprintf(['[WARNING] In user routine @%s, some UCRT values are negative. The expected range of output is {UCRT >= 0}. Please check the routine for programmi',...
                    'ng errors\n'], char(fcnHandle));
            elseif (any(isnan(UCRT) == true) == true) || (any(isinf(UCRT) == true) == true)
                % Validity check
                fprintf(['[WARNING] In user routine @%s, UCRT contains INF/NAN values. The expected range of output is {UCRT >= 0}. Please check the routine for programmin',...
                    'g errors\n'], char(fcnHandle));
            else
                % Inform the user that the routine ended
                fprintf('[NOTICE] End user routine @%s\n', char(fcnHandle));
            end
        end
    end
    methods(Static = true, Access = public)
        %% FAILURE CRITERION: MAXIMUM STRESS
        function [MSTRS] = getMstrs(stress, XT, XC, YT, YC, S)
            % Get stresses for each section point
            S1 = stress(1.0, :);
            S2 = stress(2.0, :);
            T12 = stress(3.0, :);

            % Tension-compression split (longitudinal)
            X(S1 >= 0.0) = XT(S1 >= 0.0);
            X(S1 < 0.0) = XC(S1 < 0.0);

            % Tension-compression split (transverse)
            Y(S2 >= 0.0) = YT(S2 >= 0.0);
            Y(S2 < 0.0) = YC(S2 < 0.0);

            % Compute individual strength terms
            MS11 = abs(S1./X);
            MS22 = abs(S2./Y);
            MS12 = abs(T12./S);

            % Get criterion from overall maximum
            MSTRS = max([MS11', MS22', MS12'], [], 2.0)';
        end

        %% FAILURE CRITERION: TSAI-HILL
        function [TSAIH] = getTsaih(parameter, stress, XT, XC, YT, YC, S)
            % Get stresses for each section point
            S1 = stress(1.0, :);
            S2 = stress(2.0, :);
            T12 = stress(3.0, :);

            % Tension-compression split (longitudinal)
            X(S1 >= 0.0) = XT(S1 >= 0.0);
            X(S1 < 0.0) = XC(S1 < 0.0);

            % Tension-compression split (transverse)
            Y(S2 >= 0.0) = YT(S2 >= 0.0);
            Y(S2 < 0.0) = YC(S2 < 0.0);

            % Compute the parameter based on user setting
            if parameter == 1.0
                % Strength reserve factor (failure index)
                TSAIH = sqrt(((S1.^2.0./X.^2.0) -...
                             ((S1.*S2)./X.^2.0) +...
                             (S2.^2.0./Y.^2.0) +...
                             (T12.^2.0./S.^2.0)));
            else
                % Criterion value
                TSAIH = ((S1.^2.0./X.^2.0) -...
                         ((S1.*S2)./X.^2.0) +...
                         (S2.^2.0./Y.^2.0) +...
                         (T12.^2.0./S.^2.0));
            end
        end

        %% FAILURE CRITERION: HOFFMAN
        function [HOFFMAN] = getHoffman(parameter, stress, XT, XC, YT, YC, S, S3, C)
            %{
                The actual HOFFMAN criterion considerrs the action of both
                the in-plane and out-of-plane shear stresses. Since Layup
                Analysis Tool only considers plane stress, the following
                simplifying assumptions are made:
                   [I]: ZT = XT; ZC = XC
                   [II]: SHR12 = SHR13 = SHR23 = S
                   [III]: S3 = S23 = S13 = 0
            %}

            ZT = XT;    ZC = XC; % Simplifying assumption [I]
            C1 = 0.5.*((1.0./(XC.*XT)) - (1.0./(YC.*YT)) - (1.0./(ZC.*ZT)));
            C2 = 0.5.*(-(1.0./(XC.*XT)) + (1.0./(YC.*YT)) - (1.0./(ZC.*ZT)));
            C3 = 0.5.*(-(1.0./(XC.*XT)) - (1.0./(YC.*YT)) + (1.0./(ZC.*ZT)));
            C4 = (1.0./XT) + (1.0./XC);
            C5 = (1.0./YT) + (1.0./YC);
            C6 = (1.0./ZT) + (1.0./ZC);

            SHR23 = S;  SHR13 = S;  SHR12 = S; % Simplifying assumption [II]
            C7 = 1.0./SHR23.^2.0;
            C8 = 1.0./SHR13.^2.0;
            C9 = 1.0./SHR12.^2.0;

            % Get stresses for each section point
            S1 = stress(1.0, :);
            S2 = stress(2.0, :);

            T12 = stress(3.0, :);
            S23 = S3;  S13 = S3;  S12 = T12; % Simplifying assumption [III]

            A = C1.*(S1 - S3).^2.0 + C2.*(S3 - S1).^2.0 + C3.*(S1 - S2).^2.0 + C7.*S23.^2.0 + C8.*S13.^2.0 + C9.*S12.^2.0;
            B = C4.*S1 + C5.*S2 + C6.*S3;

            % Compute the parameter based on user setting
            if parameter == 1.0
                % Strength reserve factor (failure index)
                HOFFMAN = abs(1.0./min([(-B + sqrt(B.^2.0 - (4.0.*A.*C)))./(2.0.*A); (-B - sqrt(B.^2.0 - (4.0.*A.*C)))./(2.0.*A)], [], 1.0));
            else
                % Criterion value
                HOFFMAN = abs(A + B);
            end
        end

        %% FAILURE CRITERION: TSAI-WU
        function [TSAIW] = getTsaiw(parameter, stress, XT, XC, YT, YC, S, C12, B12)
            % Get stresses for each section point
            S1 = stress(1.0, :);
            S2 = stress(2.0, :);
            T12 = stress(3.0, :);

            % Initialize Tsai-Wu parameters
            F1 = (1.0./XT) + (1.0./XC);
            F2 = (1.0./YT) + (1.0./YC);
            F11 = -(1.0./(XT.*XC));
            F22 = -(1.0./(YT.*YC));
            F66 = 1.0./S.^2.0;

            % Get F12 for non-zero values of B12
            B12Notzero = B12 ~= 0.0;
            F12(B12Notzero) = (1.0./(2.0.*B12(B12Notzero).^2.0)).*(1.0 - ((1.0./XT(B12Notzero)) + (1.0./XC(B12Notzero)) +...
                (1.0./YT(B12Notzero)) + (1.0./YC(B12Notzero))).*(B12(B12Notzero)) + ((1.0./(XT(B12Notzero).*XC(B12Notzero))) +...
                (1.0./(YT(B12Notzero).*YC(B12Notzero)))).*(B12(B12Notzero).^2.0));

            % Get F12 for zero values of B12
            B12Zero = B12Notzero == false;
            F12(B12Zero) = C12(B12Zero).*sqrt(F11(B12Zero).*F22(B12Zero));

            if parameter == 1.0
                % Get the quadratic coefficients
                A = (F11.*S1.*S1) + (F22.*S2.*S2) + (F66.*T12.*T12) + (2.0.*F12.*S1.*S2);
                B = (F1.*S1) + (F2.*S2);

                % Strength reserve factor (failure index)
                TSAIW = abs(1.0./min([(-B + sqrt(B.^2.0 + (4.0.*A)))./(2.0.*A);...
                    (-B - sqrt(B.^2.0 + (4.0.*A)))./(2.0.*A)], [], 1.0));
            else
                % Criterion value
                TSAIW = (F1.*S1) + (F2.*S2) + (F11.*S1.^2.0) +...
                    (F22.*S2.^2.0) + (F66.*T12.^2.0) + 2.0.*(F12.*S1.*S2);
            end

            % Reset NaN values of TSAIW to zero
            TSAIW(isnan(TSAIW) == true) = 0.0;
        end

        %% FAILURE CRITERION: AZZI-TSAI-HILL
        function [AZZIT] = getAzzit(parameter, nPlies_points, stress, XT, XC, YT, YC, S)
            % Get stresses for each section point
            S1 = stress(1.0, :);
            S2 = stress(2.0, :);
            T12 = stress(3.0, :);

            % Strength buffers
            X = zeros(1.0, nPlies_points);
            Y = zeros(1.0, nPlies_points);

            % Tension-compression split (longitudinal)
            X(S1 >= 0.0) = XT(S1 >= 0.0);
            X(S1 < 0.0) = XC(S1 < 0.0);

            % Tension-compression split (transverse)
            Y(S2 >= 0.0) = YT(S2 >= 0.0);
            Y(S2 < 0.0) = YC(S2 < 0.0);

            if parameter == 1.0
                % Strength reserve factor (failure index)
                AZZIT = sqrt((S1.^2.0./X.^2.0) -...
                             (abs((S1.*S2))./X.^2.0) +...
                             (S2.^2.0./Y.^2.0) +...
                             (T12.^2.0./S.^2.0));
            else
                % Criterion value
                AZZIT = (S1.^2.0./X.^2.0) -...
                         abs((S1.*S2)./X.^2.0) +...
                         (S2.^2.0./Y.^2.0) +...
                         (T12.^2.0./S.^2.0);
            end
        end

        %% FAILURE CRITERION: MAXIMUM STRAIN
        function [MSTRN] = getMstrn(nPlies_points, stress, E11, E22, V12, G12, XET, XEC, YET, YEC, SE)
            % Get stresses for each section point
            S1 = stress(1.0, :);
            S2 = stress(2.0, :);
            T12 = stress(3.0, :);

            % Strength buffers
            XE = zeros(1.0, nPlies_points);
            YE = zeros(1.0, nPlies_points);

            % Get the strains
            strain_11 = (S1./E11) - (V12.*(S2./E22));
            strain_22 = -V12.*(S1./E11) + (S2./E22);
            strain_12 = T12./G12;

            % Tension-compression split (longitudinal)
            XE(S1 >= 0.0) = XET(S1 >= 0.0);
            XE(S1 < 0.0) = XEC(S1 < 0.0);

            % Tension-compression split (transverse)
            YE(S2 >= 0.0) = YET(S2 >= 0.0);
            YE(S2 < 0.0) = YEC(S2 < 0.0);

            % Compute individual strength terms
            ME11 = strain_11./XE;
            ME22 = strain_22./YE;
            ME12 = abs(strain_12./SE);

            % Get criterion from overall maximum
            MSTRN = max([ME11; ME22; ME12]);
        end

        %% FAILURE CRITERION: HASHIN
        function [HSNFTCRT, HSNFCCRT, HSNMTCRT, HSNMCCRT] = getHashin(nPlies_points, stress, ALPHA, XHT, XHC, YHT, YHC, SHX, SHY)
            % Get stresses for each section point
            S1 = stress(1.0, :);
            S2 = stress(2.0, :);
            T12 = stress(3.0, :);

            % Mode I/II
            S11Pos = S1 >= 0.0;
            S11Neg = S1 < 0.0;

            % Initialise Hashin buffers
            HSNFTCRT = zeros(1.0, nPlies_points); 
            HSNFCCRT = HSNFTCRT;
            HSNMTCRT = HSNFTCRT; 
            HSNMCCRT = HSNFTCRT; 
            
            % Fibre-tension
            HSNFTCRT(S11Pos) = (S1(S11Pos)./ XHT(S11Pos)).^2.0 +...
                ALPHA(S11Pos).*(T12(S11Pos) ./ SHX(S11Pos)).^2.0;
 
            % Fibre-compression
            HSNFCCRT(S11Neg) = (S1(S11Neg) ./ XHC(S11Neg)).^2.0;
            
            % Mode III/IV
            S22Pos = S2 >= 0.0;
            S22Neg = S2 < 0.0;
            
            % Matrix-tension
            HSNMTCRT(S22Pos) = (S2(S22Pos) ./ YHT(S22Pos)).^2.0 +...
                (T12(S22Pos) ./ SHX(S22Pos)).^2.0;
            
            % Matrix-compression
            HSNMCCRT(S22Neg) = (S2(S22Neg) ./...
                (2.0*SHY(S22Neg))).^2.0 + ((YHC(S22Neg) ./...
                (2.0.*SHY(S22Neg))).^2.0 - 1.0).*(S2(S22Neg) ./...
                YHC(S22Neg)) + (T12(S22Neg) ./ SHX(S22Neg)).^2.0;
        end
    end
end
