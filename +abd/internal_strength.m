classdef internal_strength < handle
%   Perform strength calculations based on the ply stresses.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 2.4 Copyright Louis Vallance 2023
%   Last modified 11-May-2023 13:34:37 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%

    methods(Static = true, Access = public)
        %% MAIN FUNCTION FOR STRENGTH CALCULATION
        function [MSTRS, TSAIH, TSAIW, AZZIT, MSTRN, HSNFTCRT, HSNFCCRT,...
                HSNMTCRT, HSNMCCRT, XT, XC, YT, YC, S, C12, B12, E11,...
                E22, G12, V12, XET, XEC, YET, YEC, SE, ALPHA, XHT, XHC,...
                YHT, YHC, SHX, SHY] =...
                main(noFailStress, noFailStrain, noHashin, XT, XC, YT,...
                YC, S, C12, B12, E11, E22, G12, V12, XET, XEC, YET, YEC,...
                SE, ALPHA, XHT, XHC, YHT, YHC, SHX, SHY, stress, nPlies,...
                nPlies_points, SECTION_POINTS)

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

            if noHashin == false
                ALPHA = abd.internal_spreadProperties(ALPHA, nPlies, SECTION_POINTS);
                XHT = abd.internal_spreadProperties(XHT, nPlies, SECTION_POINTS);
                XHC = abd.internal_spreadProperties(XHC, nPlies, SECTION_POINTS);
                YHT = abd.internal_spreadProperties(YHT, nPlies, SECTION_POINTS);
                YHC = abd.internal_spreadProperties(YHC, nPlies, SECTION_POINTS);
                SHX = abd.internal_spreadProperties(SHX, nPlies, SECTION_POINTS);
                SHY = abd.internal_spreadProperties(SHY, nPlies, SECTION_POINTS);
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

            if noFailStress == false
                % Failure calculation: MSTRS
                MSTRS =...
                    ...
                    abd.internal_strength.getMstrs(stress, XT, XC, YT,...
                    YC, S);

                % Failure calculation: TSAIH
                TSAIH =...
                    ...
                    abd.internal_strength.getTsaih(1.0, stress, XT, XC,...
                    YT, YC, S);

                % Failure calculation: TSAIW
                TSAIW =...
                    ...
                    abd.internal_strength.getTsaiw(1.0, stress, XT, XC,...
                    YT, YC, S, C12, B12);

                % Failure calculation: AZZIT
                AZZIT =...
                    ...
                    abd.internal_strength.getAzzit(1.0, nPlies_points,...
                    stress, XT, XC, YT, YC, S);
            end

            % Failure calculation: MSTRN
            if noFailStrain == false
                MSTRN = ...
                    ...
                    abd.internal_strength.getMstrn(nPlies_points,...
                    stress, E11, E22, V12, G12, XET, XEC, YET, YEC, SE);
            end

            % Failure calculation: HASHIN
            if noHashin == false
                [HSNFTCRT, HSNFCCRT, HSNMTCRT, HSNMCCRT] = ...
                    ...
                    abd.internal_strength.getHashin(nPlies_points,...
                    stress, ALPHA, XHT, XHC, YHT, YHC, SHX, SHY);
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

            % Tension-compression split
            X(S1 >= 0.0) = XT(S1 >= 0.0);
            X(S1 < 0.0) = XC(S1 < 0.0);

            Y(S2 >= 0.0) = YT(S2 >= 0.0);
            Y(S2 < 0.0) = YC(S2 < 0.0);

            MS11 = abs(S1./X);
            MS22 = abs(S2./Y);
            MS12 = abs(T12./S);

            MSTRS = max([MS11', MS22', MS12'], [], 2.0)';
        end

        %% FAILURE CRITERION: TSAI-HILL
        function [TSAIH] = getTsaih(parameter, stress, XT, XC, YT, YC, S)
            % Get stresses for each section point
            S1 = stress(1.0, :);
            S2 = stress(2.0, :);
            T12 = stress(3.0, :);

            % Tension-compression split
            X(S1 >= 0.0) = XT(S1 >= 0.0);
            X(S1 < 0.0) = XC(S1 < 0.0);

            Y(S2 >= 0.0) = YT(S2 >= 0.0);
            Y(S2 < 0.0) = YC(S2 < 0.0);

            % Compute the parameter based on user setting
            if parameter == 1.0
                % Reserve factor (failure index)
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
        function [TSAIW] = getTsaiw(parameter, stress,XT, XC, YT, YC, S,...
                C12, B12)
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

                % Reserve factor (failure index)
                TSAIW = abs(1.0./min([(-B + sqrt(B.^2.0 + (4.0.*A)))./(2.0.*A);...
                    (-B - sqrt(B.^2.0 + (4.0.*A)))./(2.0.*A)], [], 1.0));
            else
                % Criterion value
                TSAIW = (F1.*S1) + (F2.*S2) + (F11.*S1.^2.0) +...
                    (F22.*S2.^2.0) + (F66.*T12.^2.0) + 2.0.*(F12.*S1.*S2);
            end

            TSAIW(isnan(TSAIW) == true) = 0.0;
        end

        %% FAILURE CRITERION: AZZI-TSAI-HILL
        function [AZZIT] = getAzzit(parameter, nPlies_points, stress,...
                XT, XC, YT, YC, S)
            % Get stresses for each section point
            S1 = stress(1.0, :);
            S2 = stress(2.0, :);
            T12 = stress(3.0, :);

            % Strength buffers
            X = zeros(1.0, nPlies_points);
            Y = zeros(1.0, nPlies_points);

            % Tension-compression split
            X(S1 >= 0.0) = XT(S1 >= 0.0);
            X(S1 < 0.0) = XC(S1 < 0.0);

            Y(S2 >= 0.0) = YT(S2 >= 0.0);
            Y(S2 < 0.0) = YC(S2 < 0.0);

            if parameter == 1.0
                % Reserve factor (failure index)
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
        function [MSTRN] = getMstrn(nPlies_points, stress, E11, E22,...
                V12, G12, XET, XEC, YET, YEC, SE)
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

            % Tension-compression split
            XE(S1 >= 0.0) = XET(S1 >= 0.0);
            XE(S1 < 0.0) = XEC(S1 < 0.0);

            YE(S2 >= 0.0) = YET(S2 >= 0.0);
            YE(S2 < 0.0) = YEC(S2 < 0.0);

            ME11 = strain_11./XE;
            ME22 = strain_22./YE;
            ME12 = abs(strain_12./SE);

            MSTRN = max([ME11; ME22; ME12]);
        end

        %% FAILURE CRITERION: HASHIN
        function [HSNFTCRT, HSNFCCRT, HSNMTCRT, HSNMCCRT] = getHashin(...
                nPlies_points, stress, ALPHA, XHT, XHC, YHT, YHC, SHX, SHY)
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
