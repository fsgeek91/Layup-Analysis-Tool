classdef internal_strength < handle
%   Perform strength calculations based on the ply stresses.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 3.0.3 Copyright Louis Vallance 2024
%   Last modified 24-Jun-2024 11:37:46 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%

    methods(Static = true, Access = public)
        %% MAIN FUNCTION FOR STRENGTH CALCULATION
        function [MSTRS, TSAIH, TSAIW, AZZIT, MSTRN, HSNFTCRT, HSNFCCRT, HSNMTCRT, HSNMCCRT, LARPFCRT, LARMFCRT, LARKFCRT, LARSFCRT, LARTFCRT, XT, XC, YT, YC, S, C12, B12, E11,...
                E22, G12, V12, XET, XEC, YET, YEC, SE, ALPHA, XHT, XHC, YHT, YHC, SHX, SHY, XLT, XLC, YLT, YLC, SLX, SLY, GL12, NL, NT, A0, PHI0, S1, S2, S3] =...
                main(noFailStress, noFailStrain, noHashin, noLaRC05, symsAvailable, XT, XC, YT, YC, S, C12, B12, E11, E22, G12, V12, XET, XEC, YET, YEC, SE, ALPHA, XHT, XHC, YHT,...
                YHC, SHX, SHY, XLT, XLC, YLT, YLC, SLX, SLY, GL12, NL, NT, A0, PHI0, stress, nPlies, nPlies_points, SECTION_POINTS, parameter)
            % Initialise the output
            S1 = [];    S2 = [];    S3 = [];

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
                if noLaRC05 == false
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
            end
            
            % Initialise output
            MSTRS = [];
            TSAIH = [];
            TSAIW = [];
            AZZIT = [];
            MSTRN = [];
            HSNFTCRT = [];
            HSNFCCRT = [];
            HSNMTCRT = [];
            HSNMCCRT = [];
            LARPFCRT = [];
            LARMFCRT = [];
            LARKFCRT = [];
            LARSFCRT = [];
            LARTFCRT = [];

            if noFailStress == false
                % Failure calculation: MSTRS
                MSTRS =...
                    ...
                    abd.internal_strength.getMstrs(stress, XT, XC, YT, YC, S);

                % Failure calculation: TSAIH
                TSAIH =...
                    ...
                    abd.internal_strength.getTsaih(parameter, stress, XT, XC, YT, YC, S);

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
                fprintf('[ERROR] The setting OUTPUT_STRENGTH requires two\narguments: {''<flag>'', ''<parameter>''}\n');

                % Reset the error flag and RETURN
                error = true;
                return
            end

            % Process the first argument
            argument = OUTPUT_STRENGTH{1.0};

            % Check validity of the argument
            if isempty(argument) == true
                output{1.0} = false;
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

        %% FAILURE CRITERION: TSAI-WU
        function [TSAIW] = getTsaiw(parameter, stress,XT, XC, YT, YC, S, C12, B12)
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
            B12Zero = ~B12Notzero;
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
