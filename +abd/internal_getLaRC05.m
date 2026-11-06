function [LARPFCRT, LARMFCRT, LARKFCRT, LARSFCRT, LARTFCRT] = internal_getLaRC05(nPlies_points, stress, ~, S1, S2, S3, G12, Xt, Xc, Yt,Yc, Sl, St, alpha0, phi0, nl,...
    nt, ~, step, mode)
%   Perform strength calculations based on the ply stresses.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 5.1.5 Copyright Louis Vallance 2026
%   Last modified 04-Sep-2026 12:21:27 UTC
%
    
    %%

%% Initialise variables
% Basic
S11 = stress(1.0, :);
S22 = stress(2.0, :);
S33 = zeros(1.0, nPlies_points);
S12 = stress(3.0, :);
S13 = S33;
S23 = S33;

% Critical plane step size
if ischar(step) == true
    switch mode
        case 1.0 % Normal
            step = 5.0;
        case 2.0 % Stacking sequence optimisation
            step = 15.0;
        otherwise % Default case (should never be reached!)
            step = 15.0;
    end
end
    
%% Initialise St/nt/phi0 if applicable
% ST
deriveElements = St == -1.0;
St(deriveElements) = Yc(deriveElements).*cosd(alpha0(deriveElements)).*(sind(alpha0(deriveElements)) + ((cosd(alpha0(deriveElements))) ./ (tand(2.0.*(alpha0(deriveElements))))));

% NL
deriveElements = nl == -1.0;
nl(deriveElements) = (Sl(deriveElements).*cosd(2.0.*alpha0(deriveElements))) ./ (Yc(deriveElements).*(cosd(alpha0(deriveElements))).^2.0);

% NT
deriveElements = nt == -1.0;
nt(deriveElements) = -1.0 ./ (tand(2.0.*(alpha0(deriveElements))));

% PHI0 (initial fibre misalignment, in degrees)
%{
    Derive only undefined values, once per distinct material. Calibration
    uses pure longitudinal compression, not the current section-point load.
    The legacy toolbox/section-point arguments are retained in the function
    interface, but a numerical solution does not require either of them.
%}
deriveElements = phi0 == -1.0;
if any(deriveElements) == true
    [materialData, ~, materialIndexes] = unique([Xc(deriveElements)', Sl(deriveElements)', G12(deriveElements)', nl(deriveElements)'], 'rows');
    phi0_material = zeros(size(materialData, 1.0), 1.0);

    for i = 1.0:size(materialData, 1.0)
        Xc_i = materialData(i, 1.0);
        Sl_i = materialData(i, 2.0);
        G12_i = materialData(i, 3.0);
        nl_i = materialData(i, 4.0);

        %{
            Critical angle at failure under pure compression (in radians).
            Rationalise the quadratic root to avoid cancellation and 0/0
        %}
        ratio = Sl_i/Xc_i;
        discriminant = 1.0 - 4.0*(ratio + nl_i)*ratio;
        if (isfinite(discriminant) == false) || (discriminant < 0.0)
            error('abd:LaRC05:InvalidMisalignment', 'Cannot derive PHI0: the LaRC05 strengths and friction coefficient give no real critical angle.');
        end
        phic = atan(2.0*ratio/(1.0 + sqrt(discriminant)));

        %{
            Linear shear response in the initial misalignment frame: phi0 +
            Xc*sin(2*phi0)/(2*G12) = phic (all angles in radians)
        %}
        phi0_material(i) = (180.0/pi)*fzero(@(angle) angle + (Xc_i/(2.0*G12_i))*sin(2.0*angle) - phic, [0.0, phic]);
    end

    phi0(deriveElements) = phi0_material(materialIndexes);
end

%% Polymer failure
% Hydrostatic stress
Sh = (1.0/3.0).*(S11 + S22 + S33);

% k-parameter
k2 = (1.0/6.0).*((S1 - S2).^2.0 + (S2 - S3).^2.0 + (S3 - S1).^2.0);

% Failure criterion
LARPFCRT = (3.0.*(k2 - (Xt - Xc).*Sh)) ./ (Xt.*Xc);

%% Matrix failure
% Initialise critical plane search buffers
precision = floor(180.0./step) + 1.0;
precision_matrix = precision;

%{
    Under plane stress, ALPHA and 180 - ALPHA give the same matrix
    criterion. Search only the first half when the grid has reflected
    pairs; otherwise retain all samples for a user-defined step size.
    Keep the full precision for the fibre kinking/splitting search
%}
if rem(180.0, step) == 0.0
    precision_matrix = floor(90.0./step) + 1.0;
end
FI_matrix = zeros(nPlies_points, precision_matrix);

% Critical plane search
FI_matrix = par_CP_MF(precision_matrix, step, S22, S33, S23, S12, S13, St, nt, Sl, nl, Yt, FI_matrix);

% Failure criterion for the critical plane
LARMFCRT = max(FI_matrix, [], 2.0)';

%% Get the fibre criterion from S11
tens = S11 > 0.0;
kink = S11 <= -0.5.*Xc;
split = S11 >= -0.5.*Xc & S11 <= 0.0;

%% Fibre kinking/splitting failure
% Initialise critical plane search buffers
FI_split = zeros(nPlies_points, precision);
FI_kink = FI_split;

% Critical plane searching
if (any(kink) == true) || (any(split) == true)
    [FI_kink, FI_split] = CP_FKS(precision, step, S22, S33, S23, S12, S13, G12, phi0, S11, FI_kink, St, nt, Sl, nl, FI_split, Yt, kink, split);

    % Failure criterion for the critical plane (kinking)
    LARKFCRT = max(FI_kink, [], 2.0)';

    % Failure criterion for the critical plane (splitting)
    LARSFCRT = max(FI_split, [], 2.0)';
else
    % Failure criterion for the critical plane (kinking)
    LARKFCRT = zeros(1.0, nPlies_points);

    % Failure criterion for the critical plane (splitting)
    LARSFCRT = LARKFCRT;
end

%% Fibre tensile failure
% Initialise LARTFCRT buffer
LARTFCRT = zeros(1.0, nPlies_points);

% Failure criterion
LARTFCRT(tens) = S11(tens)./Xt(tens);
end

%% PARFOR version of CP search for fibre matrix failure
function [FI_matrix] = par_CP_MF(precision, step, S22, S33, S23, S12, S13, St, nt, Sl, nl, Yt, FI_matrix)
parfor alphaIndex = 1.0:precision
    % Get the current value of ALPHA
    alpha_angle = (alphaIndex*step) - step;

    % Get the tractions on the current plane
    Sn = 0.5.*(S22 + S33) + 0.5.*(S22 - S33).*cosd(2.0.*alpha_angle) + S23.*sind(2.0.*alpha_angle);
    Tt = -0.5.*(S22 - S33).*sind(2.0.*alpha_angle) + S23.*cosd(2.0.*alpha_angle);
    Tl = S12.*cosd(alpha_angle) + S13.*sind(alpha_angle);

    % McCauley normal stress
    Sn_pos = Sn;
    Sn_pos(Sn_pos < 0.0) = 0.0;

    % Failure criterion on the current plane
    FI_mi_ii = sqrt((Tt ./ (St - (nt.*Sn))).^2.0 + (Tl ./ (Sl - (nl.*Sn))).^2.0 + (Sn_pos ./ Yt).^2.0);

    % Collect output
    FI_matrix(:, alphaIndex) = FI_mi_ii;
end
end

%% CP search for fibre kinking/splitting failure
function [FI_kink, FI_split] = CP_FKS(precision, step, S22, S33, S23, S12, S13, G12, phi0, S11, FI_kink, St, nt, Sl, nl, FI_split, Yt, kink, split)
for psiIndex = 1.0:precision
    % Initialise CP buffer
    FI_kink_ii = zeros(1.0, length(S22));
    FI_split_ii = FI_kink_ii;

    % Get the current value of PSI
    psii = (psiIndex*step) - step;

    % Get the stresses on the kink band
    S2_psi = cosd(psii).^2.0.*S22 + sind(psii).^2.0.*S33 + 2.0.*sind(psii).*cosd(psii).*S23;
    Tau12_psi = S12.*cosd(psii) + S13.*sind(psii);
    Tau23_psi = sind(psii).*cosd(psii).*(S33 - S22) + (cosd(psii).^2.0 - sind(psii).^2.0).*S23;
    Tau13_psi = S13.*cosd(psii) - S12.*sind(psii);

    %{
        Get the shear strain magnitude in the initial misalignment frame.
        Pinho et al. (2006), Part II, Eqs. (7)-(10), with linear shear
        response
    %}
    Tau12_initial = 0.5.*(S2_psi - S11).*sind(2.0.*phi0) + abs(Tau12_psi).*cosd(2.0.*phi0);
    gamma0 = abs(Tau12_initial)./G12;

    %{
        Preserve the initial imperfection at zero shear. Either sign gives
        the same failure index for this plane-stress case
    %}
    shearSign = sign(Tau12_psi);
    shearSign(shearSign == 0.0) = 1.0;

    % Convert shear strain (radians) before adding to angles in degrees.
    phi = shearSign.*(phi0 + (180.0/pi).*gamma0);

    % Get the stresses on the fibre misalignment plane
    S2_m = sind(phi).^2.0.*S11 + cosd(phi).^2.0.*S2_psi - 2.0.*sind(phi).*cosd(phi).*Tau12_psi;
    Tau12_m = sind(phi).*cosd(phi).*(S2_psi - S11) + (cosd(phi).^2.0 - sind(phi).^2.0).*Tau12_psi;
    Tau23_m = Tau23_psi.*cosd(phi) - Tau13_psi.*sind(phi);

    % McCauley normal stress
    S2_m_pos = S2_m;
    S2_m_pos(S2_m_pos < 0.0) = 0.0;

    % Failure criterion on the current plane
    FI_kink_ii(kink) = sqrt((Tau23_m(kink) ./ (St(kink) - (nt(kink).*S2_m(kink)))).^2.0 + (Tau12_m(kink) ./ (Sl(kink) - (nl(kink).*S2_m(kink)))).^2.0 + ((S2_m_pos(kink)) ./ Yt(kink)).^2.0);
    FI_split_ii(split) = sqrt((Tau23_m(split) ./ (St(split) - (nt(split).*S2_m(split)))).^2.0 + (Tau12_m(split) ./ (Sl(split) - (nl(split).*S2_m(split)))).^2.0 + ((S2_m_pos(split)) ./ Yt(split)).^2.0);

    % Collect output
    FI_kink(:, psiIndex) = FI_kink_ii;
    FI_split(:, psiIndex) = FI_split_ii;
end
end
