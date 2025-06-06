function [LARPFCRT, LARMFCRT, LARKFCRT, LARSFCRT, LARTFCRT] = internal_getLaRC05(nPlies_points, stress, symsAvailable, S1, S2, S3, G12, Xt, Xc, Yt,Yc, Sl, St, alpha0, phi0, nl,...
    nt, SECTION_POINTS)
%   Perform strength calculations based on the ply stresses.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 4.0.0 Copyright Louis Vallance 2025
%   Last modified 06-Jun-2025 05:42:50 UTC
%
    
    %%

%% Initialise variables
step = 15.0;
iterate = true;
S11 = stress(1.0, :);
S22 = stress(2.0, :);
S33 = zeros(1.0, nPlies_points);
S12 = stress(3.0, :);
S13 = S33;
S23 = S33;
    
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

% PHI0 (get PHIC only)
if any(phi0 == -1.0) == true
    %{
        If the Symbolic Math Toolbox is not installed, do not try to
        calculate the value of phi0 iteratively.
    %}
    if (iterate == true) && (symsAvailable == true)
        % Initialise the intermediate variable PHIC
        phic = zeros(1.0, nPlies_points);

        % Get the indexes of undefined plies
        deriveElements = phi0 == -1.0;

        % Evaluate the intermediate variable PHIC
        phic(deriveElements) = atand((1.0 - sqrt(1.0 - 4.0.*((Sl(deriveElements)./Xc(deriveElements)) + nl(deriveElements)).*(Sl(deriveElements)./Xc(deriveElements)))) ./...
            (2.0.*((Sl(deriveElements)./Xc(deriveElements)) + nl(deriveElements))));

        % Reset NaN values
        phi0(isreal(phic) == false) = 0.0;
        phi0j = phi0;

        syms phi0sym;
    else
        phi0 = zeros(1.0, nPlies_points);
        phi0j = phi0;
    end
else
    phi0j = phi0;
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
FI_matrix = zeros(nPlies_points, precision);

% Critical plane search
FI_matrix = par_CP_MF(precision, step, S22, S33, S23, S12, S13, St, nt, Sl, nl, Yt, FI_matrix);

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
    if any(phi0 == -1.0) == true
        %{
            PHI0 is undefined, so compute the value iteratively during
            critical plane searching
        %}
        [FI_kink, FI_split] = CP_FKS(precision, step, S22, S33, S23, S12, S13, G12, nPlies_points, phi0j, phi0sym, phic, Xc, S11, FI_kink, St, nt, Sl, nl, FI_split, Yt,...
            SECTION_POINTS, kink, split);
    else
        %{
            PHI0 is already defined, so use the faster PARFOR version of
            the critical plane search algorithm
        %}
        [FI_kink, FI_split] = par_CP_FKS(precision, step, S22, S33, S23, S12, S13, G12, phi0j, S11, FI_kink, St, nt, Sl, nl, FI_split, Yt, kink, split);
    end

    % Failure crtiterion for the critical plane (kinking)
    LARKFCRT = max(FI_kink, [], 2.0)';

    % Failure crtiterion for the critical plane (splitting)
    LARSFCRT = max(FI_split, [], 2.0)';
else
    % Failure crtiterion for the critical plane (kinking)
    LARKFCRT = zeros(1.0, nPlies_points);

    % Failure crtiterion for the critical plane (splitting)
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

%% CP search for fibre kinking/splitting failure (PHI0 computation)
function [FI_kink, FI_split] = CP_FKS(precision, step, S22, S33, S23, S12, S13, G12, L, phi0j, phi0sym, phic, Xc, S11, FI_kink, St, nt, Sl, nl, FI_split, Yt, SECTION_POINTS, kink,...
    split)
% Get the plies (first section points) which are derived
loopIndexes = SECTION_POINTS.*find(phi0j(1.0:SECTION_POINTS:L) == -1.0) - (SECTION_POINTS - 1.0);

for psiIndex = 1.0:precision
    % Get the current value of PSI
    psii = (psiIndex*step) - step;

    % Initialise PHI0J on this plane
    phi0j_psi = phi0j;

    % Get the stresses on the kink band
    S2_psi = cosd(psii.*S22).^2.0 + sind(psii.*S33).^2.0 + 2.0.*sind(psii).*cosd(psii.*S23);
    Tau12_psi = S12.*cosd(psii) + S13.*sind(psii);
    Tau23_psi = -sind(psii).*cosd(psii.*S22) + sind(psii).*cosd(psii.*S33) + (cosd(psii).^2.0 - sind(psii).^2.0).*S23;
    Tau13_psi = S13.*cosd(psii) - S12.*sind(psii);

    % Get the shear strain in the initial fibre misalignment frame
    gamma0 = Tau12_psi./G12;

    % Calculate the initial fibre misalignment angle iteratively, if applicable
    phi0j_ii = zeros(1.0, length(loopIndexes));

    % Get consecutive indexes for PARFOR loop
    parforIndexes = 1:length(loopIndexes);

    % Slice the loop variables to avoid broadcasting
    phic_ii = phic(loopIndexes);
    gamma0_ii = gamma0(loopIndexes);
    Xc_ii = Xc(loopIndexes);

    % Calculate the initial fibre misalignment angle iteratively
    parfor j = parforIndexes
        warning('off', 'all')

        % Get the current variable slice
        phic_i = phic_ii(j);
        gamma0_i = gamma0_ii(j);
        Xc_i = Xc_ii(j);

        % Solve symbolically for PHI0
        phi0j_ii(j) = solve(phi0sym == phic_i - gamma0_i*sin(2.0*phi0sym*(pi/180.0))*Xc_i, phi0sym);
    end
    
    warning('on','all')

    %{
        Add the solved PHI0J values (one per ply) to the remainder of the
        section points
    %}
    for k = 1:length(loopIndexes)
        phi0j_psi(loopIndexes(k):loopIndexes(k) + (SECTION_POINTS - 1.0)) = phi0j_ii(k);
    end

    % Get the misalignment angle
    phi = phi0j_psi.*sign(Tau12_psi) + gamma0;

    % Get the stresses on the fibre misalignment plane
    S2_m = sind(phi.*S11).^2.0 + cosd(phi.*S2_psi).^2.0 - 2.0.*sind(phi).*cosd(phi.*Tau12_psi);
    Tau12_m = -sind(phi).*cosd(phi.*S11) + sind(phi).*cosd(phi.*S2_psi) + (cosd(phi).^2.0 - sind(phi).^2.0).*Tau12_psi;
    Tau23_m = Tau23_psi.*cosd(phi) - Tau13_psi.*sind(phi);

    % McCauley normal stress
    S2_m_pos = S2_m;
    S2_m_pos(S2_m_pos < 0.0) = 0.0;

    % Failure criterion on the current plane
    FI_kink(kink, psiIndex) = sqrt((Tau23_m(kink) ./ (St(kink) - (nt(kink).*S2_m(kink)))).^2.0 + (Tau12_m(kink) ./ (Sl(kink) - (nl(kink).*S2_m(kink)))).^2.0 + ((S2_m_pos(kink)) ./ Yt(kink)).^2.0);
    FI_split(split, psiIndex) = sqrt((Tau23_m(split) ./ (St(split) - (nt(split).*S2_m(split)))).^2.0 + (Tau12_m(split) ./ (Sl(split) - (nl(split).*S2_m(split)))).^2.0 + ((S2_m_pos(split)) ./ Yt(split)).^2.0);
end
end

%% CP search for fibre kinking/splitting failure (omits PHI0 computation)
function [FI_kink, FI_split] = par_CP_FKS(precision, step, S22, S33, S23, S12, S13, G12, phi0j, S11, FI_kink, St, nt, Sl, nl, FI_split, Yt, kink, split)
for psiIndex = 1.0:precision
    % Initialise CP buffer
    FI_kink_ii = zeros(1.0, length(S22));
    FI_split_ii = FI_kink_ii;

    % Get the current value of PSI
    psii = (psiIndex*step) - step;

    % Get the stresses on the kink band
    S2_psi = cosd(psii.*S22).^2.0 + sind(psii.*S33).^2.0 + 2.0.*sind(psii).*cosd(psii.*S23);
    Tau12_psi = S12.*cosd(psii) + S13.*sind(psii);
    Tau23_psi = -sind(psii).*cosd(psii.*S22) + sind(psii).*cosd(psii.*S33) + (cosd(psii).^2.0 - sind(psii).^2.0).*S23;
    Tau13_psi = S13.*cosd(psii) - S12.*sind(psii);

    % Get the shear strain in the initial fibre misalignment frame
    gamma0 = Tau12_psi./G12;

    % Get the misalignment angle
    phi = phi0j.*sign(Tau12_psi) + gamma0;

    % Get the stresses on the fibre misalignment plane
    S2_m = sind(phi.*S11).^2.0 + cosd(phi.*S2_psi).^2.0 - 2.0.*sind(phi).*cosd(phi.*Tau12_psi);
    Tau12_m = -sind(phi).*cosd(phi.*S11) + sind(phi).*cosd(phi.*S2_psi) + (cosd(phi).^2.0 - sind(phi).^2.0).*Tau12_psi;
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