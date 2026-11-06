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

% Use adaptive refinement only for ordinary analysis with the default step
adaptive = false;
if ischar(step) == true
    if mode == 2.0
        step = 15.0;
    else
        step = 5.0;
        adaptive = true;
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
    Reuse the calibration when successive candidates share the materials.
    The legacy toolbox/section-point arguments are retained in the function
    interface, but a numerical solution does not require either of them.
%}
deriveElements = phi0 == -1.0;
if any(deriveElements) == true
    [materialData, ~, materialIndexes] = unique([Xc(deriveElements)', Sl(deriveElements)', G12(deriveElements)', nl(deriveElements)'], 'rows');
    phi0_material = getInitialMisalignment(materialData);
    phi0(deriveElements) = phi0_material(materialIndexes);
end

%% Polymer failure
% Hydrostatic stress
Sh = (1.0/3.0).*(S11 + S22 + S33);

% k-parameter
k2 = (1.0/6.0).*((S1 - S2).^2.0 + (S2 - S3).^2.0 + (S3 - S1).^2.0);

% Failure criterion
LARPFCRT = (3.0.*(k2 - (Xt - Xc).*Sh)) ./ (Xt.*Xc);

%% Critical plane search coefficients
% Get the number of planes for each search
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

% Reuse the angular coefficients when successive candidates share a grid
angular = getAngularCoefficients(precision, step);

%% Matrix failure
% Retain only the maximum over the angular blocks
LARMFCRT = CP_MF(precision_matrix, angular, S22, S33, S23, S12, S13, St, nt, Sl, nl, Yt, adaptive);

%% Get the fibre criterion from S11
tens = S11 > 0.0;
kink = S11 <= -0.5.*Xc;
split = S11 >= -0.5.*Xc & S11 <= 0.0;

%% Fibre kinking/splitting failure
% Initialise critical plane maxima
LARKFCRT = zeros(1.0, nPlies_points);
LARSFCRT = LARKFCRT;
active = kink | split;

% Critical plane searching
if any(active) == true
    FI_fibre = zeros(1.0, nPlies_points);
    FI_fibre(active) = CP_FKS(precision, angular, S22(active), S33(active), S23(active), S12(active), S13(active),...
        G12(active), phi0(active), S11(active), St(active), nt(active), Sl(active), nl(active), Yt(active), adaptive);

    % Both modes use the same criterion and overlap at S11 = -0.5*Xc
    LARKFCRT(kink) = FI_fibre(kink);
    LARSFCRT(split) = FI_fibre(split);
end

%% Fibre tensile failure
% Initialise LARTFCRT buffer
LARTFCRT = zeros(1.0, nPlies_points);

% Failure criterion
LARTFCRT(tens) = S11(tens)./Xt(tens);
end

%% CALIBRATE THE INITIAL FIBRE MISALIGNMENT
function [phi0_material] = getInitialMisalignment(materialData)
persistent previousMaterials previousMisalignment

% Keep only the most recent material set. Each worker has its own cache
if isequal(materialData, previousMaterials) == true
    phi0_material = previousMisalignment;
    return
end

% The key contains XC, SL, G12 and NL after resolving automatic values
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

% Replace the cache only after all materials have calibrated successfully
previousMaterials = materialData;
previousMisalignment = phi0_material;
end

%% PRECOMPUTE THE ANGULAR COEFFICIENTS
function [angular] = getAngularCoefficients(precision, step)
persistent previousStep previousPrecision coefficients

% Keep only the most recent grid. Each parallel worker has its own cache
if (isempty(previousStep) == true) || (isequal(step, previousStep) == false) || (precision ~= previousPrecision)
    angle = ((1.0:precision)' .* step) - step;
    coefficients = getCoefficients(angle);
    previousStep = step;
    previousPrecision = precision;
end

angular = coefficients;
end

%% GET COEFFICIENTS FOR A GRID OR POINT-SPECIFIC TRIAL ANGLES
function [angular] = getCoefficients(angle)
angular.cos_angle = cosd(angle);
angular.sin_angle = sind(angle);
angular.cos_double = cosd(2.0.*angle);
angular.sin_double = sind(2.0.*angle);
angular.cos_squared = angular.cos_angle.^2.0;
angular.sin_squared = angular.sin_angle.^2.0;
angular.sin_cos = angular.sin_angle.*angular.cos_angle;
angular.cos_difference = angular.cos_squared - angular.sin_squared;
end

%% CP search for matrix failure
function [FI_matrix] = CP_MF(precision, angular, S22, S33, S23, S12, S13, St, nt, Sl, nl, Yt, adaptive)
% Limit each workspace to 64 planes and approximately 65536 values
blockSize = min(64.0, max(1.0, floor(65536.0/length(S22))));
Sn_mean = 0.5.*(S22 + S33);
Sn_difference = 0.5.*(S22 - S33);
FI_matrix = NaN(1.0, length(S22));

if adaptive == true
    % A denominator crossing can create a peak between coarse samples
    risk = poleRisk(min(0.0, S22), max(0.0, S22), St, nt) | poleRisk(min(0.0, S22), max(0.0, S22), Sl, nl);
    criterion = @(coefficients, indexes, points) matrixCriterion(coefficients, indexes, Sn_mean(points), Sn_difference(points),...
        S23(points), S12(points), S13(points), St(points), nt(points), Sl(points), nl(points), Yt(points));
    FI_matrix = adaptivePlanes(precision, angular, criterion, risk);
    return
end

for alphaIndex = 1.0:blockSize:precision
    % Rows represent planes; columns represent section points
    indexes = alphaIndex:min(alphaIndex + blockSize - 1.0, precision);

    FI_mi_ii = matrixCriterion(angular, indexes, Sn_mean, Sn_difference, S23, S12, S13, St, nt, Sl, nl, Yt);

    % Preserve the full search's MAX behaviour, including all-NaN inputs
    FI_matrix = max([FI_matrix; max(FI_mi_ii, [], 1.0)], [], 1.0);
end
end

%% CP search for fibre kinking/splitting failure
function [FI_fibre] = CP_FKS(precision, angular, S22, S33, S23, S12, S13, G12, phi0, S11, St, nt, Sl, nl, Yt, adaptive)
% Precompute material-only coefficients and fixed stress differences
sin_phi0 = sind(2.0.*phi0);
cos_phi0 = cosd(2.0.*phi0);
S_difference = S33 - S22;
blockSize = min(64.0, max(1.0, floor(65536.0/length(S22))));
FI_fibre = NaN(1.0, length(S22));

if adaptive == true
    %{
        Bound the misaligned normal stress using the largest possible
        rotation. This is tighter than the principal-stress bounds for
        small imperfections and remains conservative under plane stress
    %}
    phiMaximum = abs(phi0)*(pi/180.0) + (0.5.*max(abs(S11), abs(S22 - S11)).*abs(sin_phi0) + abs(S12.*cos_phi0))./abs(G12);
    sinSquared = sin(min(phiMaximum, pi/2.0)).^2.0;
    crossTerm = abs(S12).*sin(min(2.0.*phiMaximum, pi/2.0));
    lower = min(0.0, S11.*sinSquared) + min(0.0, S22) - crossTerm;
    upper = max(0.0, S11.*sinSquared) + max(0.0, S22) + crossTerm;
    variation = (0.5.*abs(S22.*sin_phi0) + abs(S12.*cos_phi0))./abs(G12);
    risk = poleRisk(lower, upper, St, nt) | poleRisk(lower, upper, Sl, nl) | (variation > 0.25) | (isfinite(phiMaximum) == false);
    criterion = @(coefficients, indexes, points) fibreCriterion(coefficients, indexes, S22(points), S33(points), S23(points),...
        S12(points), S13(points), S_difference(points), G12(points), phi0(points), sin_phi0(points), cos_phi0(points), S11(points),...
        St(points), nt(points), Sl(points), nl(points), Yt(points));
    FI_fibre = adaptivePlanes(precision, angular, criterion, risk);
    return
end

for psiIndex = 1.0:blockSize:precision
    % Rows represent planes; columns represent active section points
    indexes = psiIndex:min(psiIndex + blockSize - 1.0, precision);

    FI_fibre_ii = fibreCriterion(angular, indexes, S22, S33, S23, S12, S13, S_difference, G12, phi0, sin_phi0, cos_phi0, S11, St, nt, Sl, nl, Yt);

    % Retain only the maximum over the planes evaluated so far
    FI_fibre = max([FI_fibre; max(FI_fibre_ii, [], 1.0)], [], 1.0);
end
end

%% MATRIX FAILURE ON THE CURRENT PLANES
function [FI_matrix] = matrixCriterion(angular, indexes, Sn_mean, Sn_difference, S23, S12, S13, St, nt, Sl, nl, Yt)
Sn = Sn_mean + Sn_difference.*angular.cos_double(indexes, :) + S23.*angular.sin_double(indexes, :);
Tt = -Sn_difference.*angular.sin_double(indexes, :) + S23.*angular.cos_double(indexes, :);
Tl = S12.*angular.cos_angle(indexes, :) + S13.*angular.sin_angle(indexes, :);

% McCauley normal stress
Sn_pos = Sn;
Sn_pos(Sn_pos < 0.0) = 0.0;
FI_matrix = sqrt((Tt ./ (St - (nt.*Sn))).^2.0 + (Tl ./ (Sl - (nl.*Sn))).^2.0 + (Sn_pos ./ Yt).^2.0);
end

%% FIBRE KINKING/SPLITTING ON THE CURRENT PLANES
function [FI_fibre] = fibreCriterion(angular, indexes, S22, S33, S23, S12, S13, S_difference, G12, phi0, sin_phi0, cos_phi0, S11, St, nt, Sl, nl, Yt)
% Get the stresses on the kink band
S2_psi = angular.cos_squared(indexes, :).*S22 + angular.sin_squared(indexes, :).*S33 + 2.0.*angular.sin_cos(indexes, :).*S23;
Tau12_psi = S12.*angular.cos_angle(indexes, :) + S13.*angular.sin_angle(indexes, :);
Tau23_psi = angular.sin_cos(indexes, :).*S_difference + angular.cos_difference(indexes, :).*S23;
Tau13_psi = S13.*angular.cos_angle(indexes, :) - S12.*angular.sin_angle(indexes, :);

%{
    Get the shear strain magnitude in the initial misalignment frame.
    Pinho et al. (2006), Part II, Eqs. (7)-(10), with linear shear response
%}
Tau12_initial = 0.5.*(S2_psi - S11).*sin_phi0 + abs(Tau12_psi).*cos_phi0;
gamma0 = abs(Tau12_initial)./G12;

% Preserve the initial imperfection when the applied shear is zero
shearSign = sign(Tau12_psi);
shearSign(shearSign == 0.0) = 1.0;
phi = shearSign.*(phi0 + (180.0/pi).*gamma0);

% These coefficients depend on the current stresses and plane
sin_phi = sind(phi);
cos_phi = cosd(phi);
sin_phi_squared = sin_phi.^2.0;
cos_phi_squared = cos_phi.^2.0;
sin_cos_phi = sin_phi.*cos_phi;

% Get the stresses on the fibre misalignment plane
S2_m = sin_phi_squared.*S11 + cos_phi_squared.*S2_psi - 2.0.*sin_cos_phi.*Tau12_psi;
Tau12_m = sin_cos_phi.*(S2_psi - S11) + (cos_phi_squared - sin_phi_squared).*Tau12_psi;
Tau23_m = Tau23_psi.*cos_phi - Tau13_psi.*sin_phi;

% McCauley normal stress
S2_m_pos = S2_m;
S2_m_pos(S2_m_pos < 0.0) = 0.0;
FI_fibre = sqrt((Tau23_m ./ (St - (nt.*S2_m))).^2.0 + (Tau12_m ./ (Sl - (nl.*S2_m))).^2.0 + (S2_m_pos ./ Yt).^2.0);
end

%% CHECK WHETHER A SHEAR DENOMINATOR CAN CROSS ZERO
function [risk] = poleRisk(lower, upper, strength, friction)
first = strength - friction.*lower;
last = strength - friction.*upper;
risk = ((min(first, last) <= 0.0) & (max(first, last) >= 0.0)) | (isfinite(first) == false) | (isfinite(last) == false);
end

%% REFINE EVERY SAMPLED LOCAL PEAK
function [maximum] = adaptivePlanes(precision, angular, criterion, risk)
persistent denseAngular

%{
    This local refinement is intended for smooth, resolved peaks. It is
    not a proof of the continuous global maximum. Potential poles, rapid
    fibre rotations and nonfinite samples use a fixed dense reference
%}
maximum = NaN(1.0, length(risk));
maximumAngle = 5.0*(precision - 1.0);

for pointIndex = 1.0:1024.0:length(risk)
    % Bound the coarse workspace to at most 37 planes by 1024 points
    points = pointIndex:min(pointIndex + 1023.0, length(risk));
    values = criterion(angular, 1.0:precision, points);
    currentMaximum = max(values, [], 1.0);
    dense = risk(points) | any(isfinite(values) == false, 1.0);

    % Keep both edges of a plateau, and include peaks at either endpoint
    peaks = false(size(values));
    peaks(1.0, :) = values(1.0, :) >= values(2.0, :);
    peaks(end, :) = values(end, :) >= values(end - 1.0, :);
    middle = values(2.0:end - 1.0, :);
    before = values(1.0:end - 2.0, :);
    after = values(3.0:end, :);
    peaks(2.0:end - 1.0, :) = (middle >= before) & (middle >= after) & ((middle > before) | (middle > after));
    peaks = peaks & (dense == false) & any(values ~= values(1.0, :), 1.0);
    [angleIndexes, pointIndexes] = find(peaks);

    if isempty(angleIndexes) == false
        angleIndexes = angleIndexes';
        pointIndexes = pointIndexes';
        peakAngle = 5.0*(angleIndexes - 1.0);
        peakValue = reshape(values(sub2ind(size(values), angleIndexes, pointIndexes)), 1.0, []);
        spacing = 5.0;

        % Retain old samples while reducing the final spacing below 0.05 deg
        while spacing > 0.05
            spacing = spacing/2.0;
            trialAngle = [max(0.0, peakAngle - spacing); min(maximumAngle, peakAngle + spacing)];
            trialAngular = getCoefficients(trialAngle);
            trialValue = criterion(trialAngular, 1.0:2.0, points(pointIndexes));
            dense(pointIndexes(any(isfinite(trialValue) == false, 1.0))) = true;
            [peakValue, selected] = max([peakValue; trialValue], [], 1.0);
            peakAngle(selected == 2.0) = trialAngle(1.0, selected == 2.0);
            peakAngle(selected == 3.0) = trialAngle(2.0, selected == 3.0);
        end

        refinedMaximum = accumarray(pointIndexes', peakValue', [length(points), 1.0], @max, NaN)';
        currentMaximum = max([currentMaximum; refinedMaximum], [], 1.0);
    end

    if any(dense) == true
        % Cache this grid separately so fallback calls retain the coarse cache
        if isempty(denseAngular) == true
            denseAngular = getCoefficients(((1.0:3601.0)' .* 0.05) - 0.05);
        end
        densePoints = points(dense);
        denseMaximum = NaN(1.0, length(densePoints));
        blockSize = min(64.0, max(1.0, floor(65536.0/length(densePoints))));
        densePrecision = floor(maximumAngle/0.05) + 1.0;
        for angleIndex = 1.0:blockSize:densePrecision
            indexes = angleIndex:min(angleIndex + blockSize - 1.0, densePrecision);
            denseValue = criterion(denseAngular, indexes, densePoints);
            denseMaximum = max([denseMaximum; max(denseValue, [], 1.0)], [], 1.0);
        end
        currentMaximum(dense) = max([currentMaximum(dense); denseMaximum], [], 1.0);
    end

    maximum(points) = currentMaximum;
end
end
