function [LARPFCRT, LARMFCRT, LARKFCRT, LARSFCRT, LARTFCRT] = internal_getLaRC05(nPlies_points, stress, symsAvailable, S1, S2, S3, G12, Xt, Xc, Yt,Yc, Sl, St, alpha0, phi0, nl,...
    nt, SECTION_POINTS)
%   Perform strength calculations based on the ply stresses.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 3.0.0 Copyright Louis Vallance 2024
%   Last modified 14-Feb-2024 15:05:03 UTC
%
    
    %%

%% Initialise variables
step = 1.0;
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
FI_mi = zeros(nPlies_points, precision);

% Critical plane search
FI_mi = par_CP_MF(precision, step, S2, S3, S23, S12, S13, St, nt, Sl, nl, Yt, FI_mi);

% Failure criterion for the critical plane
LARMFCRT = max(FI_mi, [], 2.0)';

%% Fibre kinking/splitting failure
% Initialise critical plane search buffers
FI_si = zeros(nPlies_points, precision);
FI_ki = FI_si;

% Critical plane searching
if any(phi0 == -1.0) == true
    %{
        PHI0 is undefined, so compute the value iteratively during
        critical plane searching
    %}
    [FI_ki, FI_si] = CP_FKS(precision, step, S2, S3, S23, S12, S13, G12, nPlies_points, phi0j, phi0sym, phic, Xc, S1, FI_ki, St, nt, Sl, nl,FI_si, Yt, SECTION_POINTS);
else
    %{
        PHI0 is already defined, so use the faster PARFOR version of the
        critical plane search algorithm
    %}
    [FI_ki, FI_si] = par_CP_FKS(precision, step, S2, S3, S23, S12, S13, G12, phi0j, S1, FI_ki, St, nt, Sl, nl, FI_si, Yt);
end

% Failure crtiterion for the critical plane (kinking)
LARKFCRT = max(FI_ki, [], 2.0)';

% Failure crtiterion for the critical plane (splitting)
LARSFCRT = max(FI_si, [], 2.0)';

%% Fibre tensile failure
% Initialise LARTFCRT buffer
LARTFCRT = zeros(1.0, nPlies_points);

% Get positive stress points
sPos = S1 > 0.0;

% Failure criterion
LARTFCRT(sPos) = S1(sPos)./Xt(sPos);
end

%% PARFOR version of CP search for fibre matrix failure
function [FI_mi] = par_CP_MF(precision, step, S2, S3, S23, S12, S13, St, nt, Sl, nl, Yt, FI_mi)
parfor alphaIndex = 1.0:precision
    % Initialise CP buffer
    FI_mi_ii = zeros(1.0, length(S2));

    % Get the current value of ALPHA
    alpha_angle = (alphaIndex*step) - step;

    % Get the tractions on the current plane
    Sn = 0.5.*(S2 + S3) + 0.5.*(S2 - S3).*cosd(2.0.*alpha_angle) + S23.*sind(2.0.*alpha_angle);
    Tt = -0.5.*(S2 - S3).*sind(2.0.*alpha_angle) + S23.*cosd(2.0.*alpha_angle);
    Tl = S12.*cosd(alpha_angle) + S13.*sind(alpha_angle);

    % Split normal traction into positive and negative components
    SnPos = Sn > 0.0;
    SnNeg = Sn <= 0.0;

    % Failure criterion on the current plane (tensile)
    FI_mi_ii(SnPos) = ((Tt(SnPos)) ./ (St(SnPos) - (nt(SnPos).*Sn(SnPos)))).^2.0 + ((Tl(SnPos)) ./ (Sl(SnPos) - (nl(SnPos).*Sn(SnPos)))).^2.0 + ((Sn(SnPos)) ./ (Yt(SnPos))).^2.0; %#ok<PFBNS>

    % Failure criterion on the current plane (compressive)
    FI_mi_ii(SnNeg) = ((Tt(SnNeg)) ./ (St(SnNeg) - (nt(SnNeg).*Sn(SnNeg)))).^2.0 + ((Tl(SnNeg)) ./ (Sl(SnNeg) - (nl(SnNeg).*Sn(SnNeg)))).^2.0;

    % Collect output
    FI_mi(:, alphaIndex) = FI_mi_ii;
end
end

%% Regular version of CP search for fibre kinking/splitting failure
function [FI_ki, FI_si] = CP_FKS(precision, step, S2, S3, S23, S12, S13, G12, L, phi0j, phi0sym, phic, Xc, S1, FI_ki, St, nt, Sl, nl, FI_si, Yt, SECTION_POINTS)
for i = 1.0:precision
    % Get the current value of PSI
    psii = (i*step) - step;

    % Get the stresses on the kink band
    S2_psi = cosd(psii.*S2).^2.0 + sind(psii.*S3).^2.0 + 2.0*sind(psii).*cosd(psii.*S23);
    Tau12_psi = S12.*cosd(psii) + S13.*sind(psii);
    Tau23_psi = -sind(psii).*cosd(psii.*S2) + sind(psii).*cosd(psii.*S3) + (cosd(psii).^2.0 - sind(psii).^2.0).*S23;
    Tau13_psi = S13.*cosd(psii) - S12.*sind(psii);

    % Get the GAMMA0 value on the current plane
    gamma0 = Tau12_psi./G12;

    % Get the plies (first section points) which are derived
    loopIndexes = SECTION_POINTS.*find(phi0j(1.0:SECTION_POINTS:L) == -1.0) - (SECTION_POINTS - 1.0);

    % Calculate the initial fibre misalignment angle iteratively, if applicable
    phi0j_ii = zeros(1.0, length(loopIndexes));

    % Get consecutive indexes for PARFOR loop
    parforIndexes = 1:length(loopIndexes);

    % Slice the loop variables to avoid broadcasting
    phic_ii = phic(loopIndexes);
    gamma0_ii = gamma0(loopIndexes);
    Xc_ii = Xc(loopIndexes);

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
        phi0j(loopIndexes(k):loopIndexes(k) + (SECTION_POINTS - 1.0)) = phi0j_ii(k);
    end

    % Get the misalignment angle
    phi = phi0j.*sign(Tau12_psi) + gamma0;

    % Get the stresses on the fibre misalignment plane
    S2_m = sind(phi.*S1).^2.0 + cosd(phi.*S2_psi).^2.0 - 2.0.*sind(phi).*cosd(phi.*Tau12_psi);
    Tau12_m = -sind(phi).*cosd(phi.*S1) + sind(phi).*cosd(phi.*S2_psi) + (cosd(phi).^2.0 - sind(phi).^2.0).*Tau12_psi;
    Tau23_m = Tau23_psi.*cosd(phi) - Tau13_psi.*sind(phi);

    S2Pos = S2 > 0.0;
    S2Neg = S2 <= 0.0;

    % Failure criterion on the current plane (compressive)
    FI_ki(S2Neg, i) = ((Tau23_m(S2Neg)) ./ (St(S2Neg) - (nt(S2Neg).*S2_m(S2Neg)))).^2.0 + ((Tau12_m(S2Neg)) ./ (Sl(S2Neg) - (nl(S2Neg).*S2_m(S2Neg)))).^2.0;

    % Failure criterion on the current plane (tensile)
    FI_si(S2Pos, i) = ((Tau23_m(S2Pos)) ./ (St(S2Pos) - (nt(S2Pos).*S2_m(S2Pos)))).^2.0 + ((Tau12_m(S2Pos)) ./ (Sl(S2Pos) - (nl(S2Pos).*S2_m(S2Pos)))).^2.0 + ((S2_m(S2Pos)) ./ (Yt(S2Pos))).^2.0;
end
end

%% PARFOR version of CP search for fibre kinking/splitting failure
function [FI_ki, FI_si] = par_CP_FKS(precision, step, S2, S3, S23, S12, S13, G12, phi0j, S1, FI_ki, St, nt, Sl, nl, FI_si, Yt)
parfor alphaIndex = 1.0:precision
    % Initialise CP buffer
    FI_ki_ii = zeros(1.0, length(S2));
    FI_sp_ii = FI_ki_ii;

    % Get the current value of PSI
    psii = (alphaIndex*step) - step;

    % Get the stresses on the kink band
    S2_psi = cosd(psii*S2).^2.0 + sind(psii*S3).^2.0 + 2.0*sind(psii)*cosd(psii*S23);
    Tau12_psi = S12*cosd(psii) + S13*sind(psii);
    Tau23_psi = -sind(psii)*cosd(psii*S2) + sind(psii)*cosd(psii*S3) + (cosd(psii)^2.0 - sind(psii)^2.0)*S23;
    Tau13_psi = S13*cosd(psii) - S12*sind(psii);

    % Get the GAMMA0 value on the current plane
    gamma0 = Tau12_psi./G12;

    % Get the misalignment angle
    phi = phi0j.*sign(Tau12_psi) + gamma0;

    % Get the stresses on the fibre misalignment plane
    S2_m = sind(phi.*S1).^2.0 + cosd(phi.*S2_psi).^2.0 - 2.0.*sind(phi).*cosd(phi.*Tau12_psi);
    Tau12_m = -sind(phi).*cosd(phi.*S1) + sind(phi).*cosd(phi.*S2_psi) + (cosd(phi).^2.0 - sind(phi).^2.0).*Tau12_psi;
    Tau23_m = Tau23_psi.*cosd(phi) - Tau13_psi.*sind(phi);

    S2Pos = S2 > 0.0;
    S2Neg = S2 <= 0.0;

    % Failure criterion on the current plane (tensile)
    FI_ki_ii(S2Neg) = ((Tau23_m(S2Neg)) ./ (St(S2Neg) - (nt(S2Neg).*S2_m(S2Neg)))).^2.0 + ((Tau12_m(S2Neg)) ./ (Sl(S2Neg) - (nl(S2Neg).*S2_m(S2Neg)))).^2.0; %#ok<PFBNS>

    % Failure criterion on the current plane (compressive)
    FI_sp_ii(S2Pos) = ((Tau23_m(S2Pos)) ./ (St(S2Pos) - (nt(S2Pos).*S2_m(S2Pos)))).^2.0 + ((Tau12_m(S2Pos)) ./ (Sl(S2Pos) - (nl(S2Pos).*S2_m(S2Pos)))).^2.0 + ((S2_m(S2Pos)) ./ (Yt(S2Pos))).^2.0; %#ok<PFBNS>

    % Collect output
    FI_ki(:, alphaIndex) = FI_ki_ii;
    FI_si(:, alphaIndex) = FI_sp_ii;
end
end