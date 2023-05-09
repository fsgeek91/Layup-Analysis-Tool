classdef internal_strength < handle
%   Perform strength calculations based on the ply stresses.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 2.3 Copyright Louis Vallance 2023
%   Last modified 09-May-2023 07:31:07 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%

    methods(Static = true, Access = public)
        %% MAIN FUNCTION FOR STRENGTH CALCULATION
        function [MSTRS, TSAIH, TSAIW, AZZIT, MSTRN, Xt, Xc, Yt, Yc, S,...
                C12, B12, E11, E22, G12, V12, XEt, XEc, YEt, YEc, SE] =...
                main(noFailStress, noFailStrain, Xt, Xc, Yt, Yc, S, C12,...
                B12, E11, E22, G12, V12, XEt, XEc, YEt, YEc, SE, stress,...
                nPlies, nPlies_points, SECTION_POINTS)

            % Spread material data over section points
            if noFailStress == false
                Xt = abd.internal_spreadProperties(Xt, nPlies, SECTION_POINTS);
                Xc = abd.internal_spreadProperties(Xc, nPlies, SECTION_POINTS);
                Yt = abd.internal_spreadProperties(Yt, nPlies, SECTION_POINTS);
                Yc = abd.internal_spreadProperties(Yc, nPlies, SECTION_POINTS);
                S = abd.internal_spreadProperties(S, nPlies, SECTION_POINTS);
                C12 = abd.internal_spreadProperties(C12, nPlies, SECTION_POINTS);
                B12 = abd.internal_spreadProperties(B12, nPlies, SECTION_POINTS);
            end

            if noFailStrain == false
                E11 = abd.internal_spreadProperties(E11, nPlies, SECTION_POINTS);
                E22 = abd.internal_spreadProperties(E22, nPlies, SECTION_POINTS);
                G12 = abd.internal_spreadProperties(G12, nPlies, SECTION_POINTS);
                V12 = abd.internal_spreadProperties(V12, nPlies, SECTION_POINTS);
                XEt = abd.internal_spreadProperties(XEt, nPlies, SECTION_POINTS);
                XEc = abd.internal_spreadProperties(XEc, nPlies, SECTION_POINTS);
                YEt = abd.internal_spreadProperties(YEt, nPlies, SECTION_POINTS);
                YEc = abd.internal_spreadProperties(YEc, nPlies, SECTION_POINTS);
                SE = abd.internal_spreadProperties(SE, nPlies, SECTION_POINTS);
            end
            
            % Initialise output
            MSTRS = [];
            TSAIH = [];
            TSAIW = [];
            AZZIT = [];
            MSTRN = [];

            if noFailStress == false
                % Failure calculation: MSTRS
                MSTRS =...
                    ...
                    abd.internal_strength.getMstrs(stress, Xt, Xc, Yt,...
                    Yc, S);

                % Failure calculation: TSAIH
                TSAIH =...
                    ...
                    abd.internal_strength.getTsaih(1.0, stress, Xt, Xc,...
                    Yt, Yc, S);

                % Failure calculation: TSAIW
                TSAIW =...
                    ...
                    abd.internal_strength.getTsaiw(1.0, nPlies_points,...
                    stress, Xt, Xc, Yt, Yc, S, C12, B12);

                % Failure calculation: AZZIT
                AZZIT =...
                    ...
                    abd.internal_strength.getAzzit(1.0, nPlies_points,...
                    stress, Xt, Xc, Yt, Yc, S);
            end

            % Failure calculation: MSTRN
            if noFailStrain == false
                MSTRN = ...
                    ...
                    abd.internal_strength.getMstrn(nPlies_points,...
                    stress, E11, E22, V12, G12, XEt, XEc, YEt, YEc, SE);
            end
        end
    end
    methods(Static = true, Access = public)
        %% FAILURE CRITERION: MAXIMUM STRESS
        function [MSTRS] = getMstrs(stress, Xt, Xc, Yt, Yc, S)
            % Get stresses for each section point
            S1 = stress(1.0, :);
            S2 = stress(2.0, :);
            T12 = stress(3.0, :);

            % Tension-compression split
            X(S1 >= 0.0) = Xt(S1 >= 0.0);
            X(S1 < 0.0) = Xc(S1 < 0.0);

            Y(S2 >= 0.0) = Yt(S2 >= 0.0);
            Y(S2 < 0.0) = Yc(S2 < 0.0);

            MS11 = abs(S1./X);
            MS22 = abs(S2./Y);
            MS12 = abs(T12./S);

            MSTRS = max([MS11', MS22', MS12'], [], 2.0)';
        end

        %% FAILURE CRITERION: TSAI-HILL
        function [TSAIH] = getTsaih(parameter, stress, Xt, Xc, Yt, Yc, S)
            % Get stresses for each section point
            S1 = stress(1.0, :);
            S2 = stress(2.0, :);
            T12 = stress(3.0, :);

            % Tension-compression split
            X(S1 >= 0.0) = Xt(S1 >= 0.0);
            X(S1 < 0.0) = Xc(S1 < 0.0);

            Y(S2 >= 0.0) = Yt(S2 >= 0.0);
            Y(S2 < 0.0) = Yc(S2 < 0.0);

            % Compute the parameter based on user setting
            if parameter == 1.0
                % Reserve factor (failure index)
                TSAIH = sqrt(((S1.^2.0./X.^2.0) - ((S1.*S2)./X.^2.0) +...
                    (S2.^2.0./Y.^2.0) + (T12.^2.0./S.^2.0)));
            else
                % Criterion value
                TSAIH = ((S1.^2.0./X.^2.0) - ((S1.*S2)./X.^2.0) +...
                    (S2.^2.0./Y.^2.0) + (T12.^2.0./S.^2.0));
            end
        end

        %% FAILURE CRITERION: TSAI-WU
        function [TSAIW] = getTsaiw(parameter, nPlies_points, stress,...
                Xt, Xc, Yt, Yc, S, B12, C12)
            % Get stresses for each section point
            S1 = stress(1.0, :);
            S2 = stress(2.0, :);
            T12 = stress(3.0, :);

            % Strength buffers
            X = zeros(1.0, nPlies_points);

            % Initialize Tsai-Wu parameters
            F1 = (1.0./Xt) + (1.0./Xc);
            F2 = (1.0./Yt) + (1.0./Yc);
            F11 = -(1.0./(Xt.*Xc));
            F22 = -(1.0./(Yt.*Yc));
            F66 = 1.0./S.^2.0;
            F12 = X;

            for i = 1:nPlies_points
                if B12(i) ~= 0.0
                    F12(i) = (1.0/(2.0*B12(i)^2.0))*(1.0 - ((1.0/Xt(i)) + (1.0/Xc(i)) +...
                        (1.0/Yt(i)) + (1.0/Yc(i)))*(B12(i)) + ((1.0/(Xt(i)*Xc(i))) +...
                        (1.0/(Yt(i)*Yc(i))))*(B12(i)^2.0));
                else
                    F12(i) = C12(i)*sqrt(F11(i)*F22(i));
                end
            end

            A = (F11.*S1.*S1) + (F22.*S2.*S2) + (F66.*T12.*T12) + (2.0.*F12.*S1.*S2);
            B = (F1.*S1) + (F2.*S2);
            TSAIW = -1.0.*ones(1.0, nPlies_points);

            for i = 1:nPlies_points
                Ai = A(i);
                Bi = B(i);
                Ci = -1.0;

                % Compute the parameter based on user setting
                if parameter == 1.0
                    % Reserve factor (failure index)
                    TSAIW(i) = abs(1.0./min([(-Bi + sqrt(Bi.^2.0 - (4.0.*Ai.*Ci)))./(2.0.*Ai),...
                        (-Bi - sqrt(Bi.^2.0 - (4.0.*Ai.*Ci)))./(2.0.*Ai)]));
                else
                    % Criterion value
                    TSAIW(i) = (F1.*S1) + (F2.*S2) + (F11.*S1^2.0) +...
                        (F22.*S2^2.0) + (F66.*T12^2.0) + 2.0.*(F12.*S1.*S2);
                end
            end

            TSAIW(isnan(TSAIW) == true) = 0.0;
        end

        %% FAILURE CRITERION: AZZI-TSAI-HILL
        function [AZZIT] = getAzzit(parameter, nPlies_points, stress,...
                Xt, Xc, Yt, Yc, S)
            % Get stresses for each section point
            S1 = stress(1.0, :);
            S2 = stress(2.0, :);
            T12 = stress(3.0, :);

            % Strength buffers
            X = zeros(1.0, nPlies_points);
            Y = zeros(1.0, nPlies_points);

            % Tension-compression split
            X(S1 >= 0.0) = Xt(S1 >= 0.0);
            X(S1 < 0.0) = Xc(S1 < 0.0);

            Y(S2 >= 0.0) = Yt(S2 >= 0.0);
            Y(S2 < 0.0) = Yc(S2 < 0.0);

            AZZIT = -1.0.*ones(1.0, nPlies_points);

            for i = 1:nPlies_points
                S1i = S1(i);
                S2i = S2(i);
                T12i = T12(i);

                Xi = X(i);
                Yi = Y(i);
                Si = S(i);

                % Compute the parameter based on user setting
                if parameter == 1.0
                    % Reserve factor (failure index)
                    AZZIT(i) = sqrt(max((S1i.^2.0./Xi.^2.0) -...
                        (abs((S1i.*S2i))/Xi.^2.0) + (S2i.^2.0./Yi.^2.0) +...
                        (T12i.^2.0./Si.^2.0)));
                else
                    % Criterion value
                    AZZIT(i) = (S1^2.0./Xi^2.0) - abs((S1.*S2)./Xi^2.0) +...
                        (S2^2.0./Yi^2.0) + (T12^2.0./Si^2.0);
                end
            end
        end

        %% FAILURE CRITERION: MAXIMUM STRAIN
        function [MSTRN] = getMstrn(nPlies_points, stress, E11, E22,...
                V12, G12, XEt, XEc, YEt, YEc, SE)
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
            XE(S1 >= 0.0) = XEt(S1 >= 0.0);
            XE(S1 < 0.0) = XEc(S1 < 0.0);

            YE(S2 >= 0.0) = YEt(S2 >= 0.0);
            YE(S2 < 0.0) = YEc(S2 < 0.0);

            ME11 = strain_11./XE;
            ME22 = strain_22./YE;
            ME12 = abs(strain_12./SE);

            MSTRN = max([ME11; ME22; ME12]);
        end

        %% FAILURE CRITERION: HASHIN
        function [HSNFTCRT, HSNFCCRT, HSNMTCRT, HSNMCCRT] = getHashin(...
                stress, Xht, Xhc, Yht, Yhc, Sl, St, alpha)
            % Get stresses for each section point
            S1 = stress(1.0, :);
            S2 = stress(2.0, :);
            T12 = stress(3.0, :);

            % Mode I/II
            S11Pos = S1 >= 0.0;
            S11Neg = S1 < 0.0;
            
            if any(S11Pos) == false
                HSNFTCRT = 0.0;
            else
                HSNFTCRT = max((S1(S11Pos)./ Xht).^2.0 +...
                    alpha.*(T12(S11Pos) ./ Sl).^2.0);
            end
            if any(S11Neg) == false
                HSNFCCRT = 0.0;
            else
                HSNFCCRT = max((S1(S11Neg) ./ Xhc).^2.0);
            end
            
            % Mode III/IV
            S22Pos = S2 >= 0.0;
            S22Neg = S2 < 0.0;
            
            if any(S22Pos) == false
                HSNMTCRT = 0.0;
            else
                HSNMTCRT = max((S2(S22Pos) ./ Yht).^2.0 +...
                    (T12(S22Pos) ./ Sl).^2.0);
            end
            if any(S22Neg) == false
                HSNMCCRT = 0.0;
            else
                HSNMCCRT = max((S2(S22Neg) ./...
                    (2.0*St)).^2.0 + ((Yhc ./ (2.0.*St)).^2.0 - 1.0).*...
                    (S2(S22Neg) ./ Yhc) + (T12(S22Neg) ./ Sl).^2.0);
            end
        end
    end
end
