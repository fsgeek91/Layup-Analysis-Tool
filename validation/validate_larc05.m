function [] = validate_larc05()
%   Regression checks for LaRC05 rotations and plane searches (#95, #102-104).
%
%   RUN THIS FUNCTION FROM THE REPOSITORY ROOT:
%   >> addpath('validation'); validate_larc05
%
%   Independent reference stresses are obtained using rotation matrices.
%   Existing parallel pool settings apply. No analysis files are written.
%
%   Layup Analysis Tool 5.1.5 Copyright Louis Vallance 2026
%
%% - DO NOT EDIT BELOW LINE
%__________________________________________________________________________
%%
% Material rows: XT, XC, YT, YC, SL, ST, G12, NL, NT, A0, PHI0.
material = [2000.0; 1500.0; 50.0; 200.0; 80.0; 70.0; 5000.0; 0.2; 0.3; 53.0; 2.0];

%% CHECK STRESS TRANSFORMATIONS AGAINST MATRIX ROTATIONS
stress = [-1000.0, -200.0, -900.0, 100.0; -40.0, 20.0, 15.0, 20.0; 30.0, -25.0, -40.0, 10.0];
materials = repmat(material, 1.0, 4.0);
materials(11.0, :) = [1.5, 2.0, 2.5, 3.0];
actual = evaluate(stress, materials, false, 1.0);
checkValues(actual(3.0:4.0, :), referenceFailure(stress, materials), 'Matrix-rotation reference');

%% CHECK ZERO STRESS AND PURE COMPRESSION WITHOUT APPLIED SHEAR
materials = repmat(material, 1.0, 2.0);
materials(11.0, 2.0) = -1.0;
checkValues(evaluate(zeros(3.0, 2.0), materials, false, 1.0), zeros(5.0, 2.0), 'Zero stress');
stress = repmat([-material(2.0); 0.0; 0.0], 1.0, 2.0);
actual = evaluate(stress, materials, false, 1.0);
assert(actual(3.0, 1.0) > 0.0, 'Zero applied shear must not remove the supplied imperfection.');
checkValues(actual(3.0, 2.0), 1.0, 'Calibrated uniaxial compressive strength');
checkValues(actual(4.0, :), zeros(1.0, 2.0), 'Pure-compression mode selection');

%% CHECK MATRIX-SEARCH SYMMETRY AGAINST THE FULL ANGULAR GRID
searchStress = [-1000.0, -200.0, 100.0, -900.0, -300.0, 100.0, 0.0, 100.0;...
    -200.0, -80.0, -40.0, 0.0, 20.0, 100.0, 0.0, 0.0; 0.0, -30.0, 40.0, -60.0, 30.0, -20.0, 0.0, 80.0];
searchMaterials = repmat(material, 1.0, size(searchStress, 2.0));
searchMaterials(1.0:7.0, 5.0:8.0) = 1.2*searchMaterials(1.0:7.0, 5.0:8.0);
searchMaterials(11.0, 5.0:8.0) = 1.5;
stepCases = {'DEFAULT', 'DEFAULT', 5.0, 15.0, 20.0, 7.0, 40.0, 2.5, 90.0};
searchModes = [1.0, 2.0, ones(1.0, 7.0)];
referenceSteps = [5.0, 15.0, 5.0, 15.0, 20.0, 7.0, 40.0, 2.5, 90.0];
baseline = evaluate(searchStress, searchMaterials, false, 1.0);
for i = 1.0:length(stepCases)
    searchActual = evaluate(searchStress, searchMaterials, false, 1.0, stepCases{i}, searchModes(i));
    checkValues(searchActual(2.0, :), referenceMatrixFailure(searchStress, searchMaterials, referenceSteps(i), 180.0), 'Full-grid matrix reference');
    checkValues(searchActual(3.0:4.0, :), referenceFailure(searchStress, searchMaterials, referenceSteps(i)), 'Full-grid kinking/splitting reference');
    checkValues(searchActual([1.0, 5.0], :), baseline([1.0, 5.0], :), 'Polymer and tensile criteria');
end

% A non-divisor must retain its upper-half samples. At 40-degree spacing,
% the transverse-compression case peaks at 120 degrees (equivalent to 60),
% which is absent from the truncated grid of 0, 40 and 80 degrees.
fullGrid = referenceMatrixFailure(searchStress, searchMaterials, 40.0, 180.0);
truncatedGrid = referenceMatrixFailure(searchStress, searchMaterials, 40.0, 90.0);
assert(fullGrid(1.0) > truncatedGrid(1.0) + 1e-6, 'The non-divisor case must detect a lost critical plane.');

%% CHECK ANGULAR BLOCKS, RUNNING MAXIMA AND CHANGES BETWEEN CALLS
blockStress = [-750.0, -751.0, -749.0, 0.0, 100.0, -900.0, 100.0, 0.0;...
    -200.0, -40.0, 20.0, 15.0, -30.0, -10.0, 0.0, 0.0; 35.0, -35.0, 25.0, 30.0, 40.0, 0.0, -80.0, 0.0];
blockMaterials = repmat(material, 1.0, size(blockStress, 2.0));
blockMaterials(7.0, :) = [5000.0, 6000.0, 7000.0, 8000.0, 5000.0, 6000.0, 7000.0, 8000.0];
blockMaterials(11.0, :) = [2.0, 1.5, 2.5, 3.0, 0.0, 1.0, 1.5, 2.0];
changedMaterials = blockMaterials;
changedMaterials(3.0:7.0, :) = 1.3*changedMaterials(3.0:7.0, :);
changedMaterials(8.0:9.0, :) = 0.7*changedMaterials(8.0:9.0, :);
changedMaterials(11.0, :) = 0.5*changedMaterials(11.0, :);

% Fine grids span several blocks and leave a partial final block. Alternate
% materials at the same spacing, then change the spacing and return to the
% first case so that cached coefficients cannot retain stale inputs.
blockSteps = [0.5, 0.5, 1.3, 5.0, 0.5];
blockResults = cell(1.0, length(blockSteps));
for i = 1.0:length(blockSteps)
    currentMaterials = blockMaterials;
    if i == 2.0
        currentMaterials = changedMaterials;
    end
    blockResults{i} = evaluate(blockStress, currentMaterials, false, 1.0, blockSteps(i), 1.0);
    checkValues(blockResults{i}(2.0, :), referenceMatrixFailure(blockStress, currentMaterials, blockSteps(i), 180.0), 'Blocked matrix reference');
    checkValues(blockResults{i}(3.0:4.0, :), referenceFailure(blockStress, currentMaterials, blockSteps(i)), 'Blocked kinking/splitting reference');
end
checkValues(blockResults{5.0}, blockResults{1.0}, 'Repeated grid and material inputs');
materialDifference = blockResults{2.0}(2.0:4.0, :) - blockResults{1.0}(2.0:4.0, :);
assert(any(abs(materialDifference(:)) > 1e-6) == true, 'Changed materials must affect the failure indices.');
checkValues(blockResults{1.0}(3.0, 1.0), blockResults{1.0}(4.0, 1.0), 'Shared kinking/splitting boundary');
assert(blockResults{1.0}(3.0, 1.0) > 0.0, 'Both criteria must remain active at S11 = -0.5*XC.');

% Check scalar output dimensions and the path with no active compression
% criterion, including pure shear and zero stress in the mixed case above.
checkValues(evaluate(blockStress(:, 1.0), blockMaterials(:, 1.0), false, 1.0, 0.5, 1.0), blockResults{1.0}(:, 1.0), 'Single section point');
tensilePoints = [5.0, 7.0];
checkValues(evaluate(blockStress(:, tensilePoints), blockMaterials(:, tensilePoints), false, 1.0, 0.5, 1.0),...
    blockResults{1.0}(:, tensilePoints), 'Tensile-only section points');

% Repetition increases the block workspace without repeating the expensive
% independent reference. Every point must keep its own running maximum.
repeatCount = 257.0;
checkValues(evaluate(repmat(blockStress, 1.0, repeatCount), repmat(blockMaterials, 1.0, repeatCount), false, 1.0, 0.5, 1.0),...
    repmat(blockResults{1.0}, 1.0, repeatCount), 'Large section-point array');

%% PRESERVE FULL-GRID REDUCTION OF NONFINITE CRITERIA
% Zero strengths produce NaN at every angle under zero stress. An active
% criterion must retain NaN, while inactive fibre kinking remains zero.
undefinedMaterial = material;
undefinedMaterial([3.0, 5.0, 6.0]) = 0.0;
undefinedResult = evaluate(zeros(3.0, 1.0), undefinedMaterial, false, 1.0, 0.5, 1.0);
assert(all(isnan(undefinedResult([2.0, 4.0]))) == true, 'All-NaN active criteria must remain NaN.');
checkValues(undefinedResult([1.0, 3.0, 5.0]), zeros(3.0, 1.0), 'Inactive nonfinite criteria');

% At 90-degree spacing, the 0/180-degree planes have 0/0 transverse shear
% terms; the 90-degree plane has a finite zero criterion. The original
% full-grid MAX ignores the NaN samples when a finite sample is available.
singularMaterial = material;
singularMaterial(6.0) = 60.0;
singularResult = evaluate([100.0; 200.0; 0.0], singularMaterial, false, 1.0, 90.0, 1.0);
checkValues(singularResult(2.0), 0.0, 'Mixed finite and NaN matrix criteria');

%% CHECK MIXED MATERIALS, SECTION POINTS AND TOOLBOX INDEPENDENCE
materials = repmat(material, 1.0, 6.0);
materials(11.0, 1.0:3.0) = -1.0;
materials(2.0, 4.0:6.0) = 1200.0;
materials(7.0, 4.0:6.0) = 6000.0;
materials(11.0, 4.0:6.0) = 1.5;
stress = [-1100.0, -200.0, -800.0, -900.0, -300.0, -700.0;...
    -40.0, 20.0, -10.0, -25.0, 15.0, 5.0; 30.0, -25.0, 0.0, -35.0, 20.0, 0.0];
actual = evaluate(stress, materials, false, 3.0);
checkValues(evaluate(stress, materials, true, 3.0), actual, 'Symbolic toolbox independence');
checkValues(evaluate(stress, materials, false, 1.0), actual, 'Section-point grouping');

% Obtain the imperfection independently by requiring the matrix-rotation
% reference to reach failure under pure longitudinal compression.
explicitMaterials = materials;
explicitMaterials(11.0, 1.0:3.0) = (180.0/pi)*fzero(@(angle) compressionFailure(angle, material) - 1.0, [0.0, 0.1]);
checkValues(evaluate(stress, explicitMaterials, false, 3.0), actual, 'Material-only calibration');
checkValues(actual(3.0:4.0, :), referenceFailure(stress, explicitMaterials), 'Mixed-material reference');
order = [6.0, 1.0, 4.0, 2.0, 5.0, 3.0];
checkValues(evaluate(stress(:, order), materials(:, order), false, 3.0), actual(:, order), 'Section-point order');

% Derive two different materials together and check the material-to-point
% mapping after changing the order of their section points.
derivedMaterials = materials;
derivedMaterials(11.0, :) = -1.0;
derivedActual = evaluate(stress, derivedMaterials, false, 3.0);
explicitMaterials(11.0, 4.0:6.0) = (180.0/pi)*fzero(@(angle) compressionFailure(angle, materials(:, 4.0)) - 1.0, [0.0, 0.1]);
checkValues(evaluate(stress, explicitMaterials, false, 3.0), derivedActual, 'Multiple material calibrations');
checkValues(evaluate(stress(:, order), derivedMaterials(:, order), false, 3.0), derivedActual(:, order), 'Derived-material mapping');

%% CHECK CALIBRATION REUSE AND CHANGES TO EACH MATERIAL INPUT
% Candidate loads and point ordering may change without changing the
% material-only imperfection. Supplied angles must still take precedence.
candidateStress = 0.7*stress(:, order);
checkValues(evaluate(candidateStress, derivedMaterials(:, order), false, 3.0),...
    evaluate(candidateStress, explicitMaterials(:, order), false, 3.0), 'Cached calibration under changed loads');
mixedMaterials = derivedMaterials;
mixedMaterials(11.0, [2.0, 5.0]) = [1.0, 3.0];
mixedExplicit = explicitMaterials;
mixedExplicit(11.0, [2.0, 5.0]) = [1.0, 3.0];
checkValues(evaluate(stress, mixedMaterials, false, 3.0), evaluate(stress, mixedExplicit, false, 3.0), 'Supplied angles with cached materials');

% Change one calibration input at a time. Obtain each expected angle from
% the independent uniaxial-failure reference rather than the cached solver.
calibrationMaterial = material;
calibrationMaterial(11.0) = -1.0;
calibrationStress = [-1000.0; -40.0; 30.0];
calibrationRows = [2.0, 5.0, 7.0, 8.0];
calibrationValues = [1300.0, 100.0, 7000.0, 0.4];
for i = 1.0:length(calibrationRows)
    evaluate(calibrationStress, calibrationMaterial, false, 1.0);
    changedMaterial = calibrationMaterial;
    changedMaterial(calibrationRows(i)) = calibrationValues(i);
    explicitMaterial = changedMaterial;
    explicitMaterial(11.0) = (180.0/pi)*fzero(@(angle) compressionFailure(angle, changedMaterial) - 1.0, [0.0, 0.1]);
    expected = evaluate(calibrationStress, explicitMaterial, false, 1.0);
    checkValues(evaluate(calibrationStress, changedMaterial, false, 1.0), expected, 'Changed calibration input');
    checkValues(evaluate(calibrationStress, changedMaterial, false, 1.0), expected, 'Repeated changed material');
end

% The cache must use the effective NL after derivation. Both ALPHA0 and YC
% can change this value even when the supplied NL remains -1.
derivedFriction = calibrationMaterial;
derivedFriction(8.0) = -1.0;
frictionCases = [53.0, 57.0, 57.0; 200.0, 200.0, 240.0];
for i = 1.0:size(frictionCases, 2.0)
    derivedFriction([10.0, 4.0]) = frictionCases(:, i);
    explicitMaterial = derivedFriction;
    explicitMaterial(8.0) = material(5.0)*cosd(2.0*derivedFriction(10.0))/(derivedFriction(4.0)*cosd(derivedFriction(10.0))^2.0);
    explicitMaterial(11.0) = (180.0/pi)*fzero(@(angle) compressionFailure(angle, explicitMaterial) - 1.0, [0.0, 0.1]);
    checkValues(evaluate(calibrationStress, derivedFriction, false, 1.0),...
        evaluate(calibrationStress, explicitMaterial, false, 1.0), 'Derived friction calibration');
end

%% CHECK INVARIANCE UNDER A CHANGE OF STRESS UNITS
scaledMaterials = materials;
scaledMaterials(1.0:7.0, :) = 1e6*scaledMaterials(1.0:7.0, :);
checkValues(evaluate(1e6*stress, scaledMaterials, false, 3.0), actual, 'Stress-unit invariance');

%% CHECK INVALID MATERIAL CALIBRATION
% A valid first material followed by an invalid one must not partially
% overwrite the previous calibration, nor turn a repeated error into a hit.
checkValues(evaluate(stress, derivedMaterials, false, 3.0), derivedActual, 'Calibration before invalid input');
invalidMaterial = repmat(calibrationMaterial, 1.0, 2.0);
invalidMaterial(2.0, 1.0) = 1000.0;
invalidMaterial(5.0, 2.0) = 1000.0;
for i = 1.0:2.0
    failedAsExpected = false;
    try
        evaluate(repmat(calibrationStress, 1.0, 2.0), invalidMaterial, false, 1.0);
    catch MException
        failedAsExpected = strcmp(MException.identifier, 'abd:LaRC05:InvalidMisalignment');
    end
    assert(failedAsExpected == true, 'Invalid calibration must report the material error.');
end
checkValues(evaluate(stress, derivedMaterials, false, 3.0), derivedActual, 'Calibration after invalid input');

fprintf('[PASS] LaRC05 rotations, plane searches, calibration reuse, zero stress, material grouping and stress units\n');
end

%% EVALUATE ALL FIVE CRITERIA
function [result] = evaluate(stress, material, symsAvailable, sectionPoints, step, mode)
if nargin < 5.0
    step = 15.0;
end
if nargin < 6.0
    mode = 1.0;
end
nPoints = size(stress, 2.0);
meanStress = 0.5.*(stress(1.0, :) + stress(2.0, :));
radius = sqrt((0.5.*(stress(1.0, :) - stress(2.0, :))).^2.0 + stress(3.0, :).^2.0);
[PF, MF, KF, SF, TF] = abd.internal_getLaRC05(nPoints, stress, symsAvailable, meanStress + radius, meanStress - radius, zeros(1.0, nPoints),...
    material(7.0, :), material(1.0, :), material(2.0, :), material(3.0, :), material(4.0, :), material(5.0, :), material(6.0, :),...
    material(10.0, :), material(11.0, :), material(8.0, :), material(9.0, :), sectionPoints, step, mode);
result = [PF; MF; KF; SF; TF];
end

%% INDEPENDENT REFERENCE USING TWO ORTHOGONAL ROTATIONS
function [result] = referenceFailure(stress, material, step)
if nargin < 3.0
    step = 15.0;
end
result = zeros(2.0, size(stress, 2.0));
for i = 1.0:size(stress, 2.0)
    tensor = [stress(1.0, i), stress(3.0, i), 0.0; stress(3.0, i), stress(2.0, i), 0.0; 0.0, 0.0, 0.0];
    initialAngle = material(11.0, i)*pi/180.0;
    initialRotation = getRotation(initialAngle);
    for psi = (0.0:step:180.0)*pi/180.0
        rotation = [1.0, 0.0, 0.0; 0.0, cos(psi), sin(psi); 0.0, -sin(psi), cos(psi)];
        kinkTensor = rotation*tensor*rotation';
        unsignedTensor = kinkTensor;
        unsignedTensor(1.0, 2.0) = abs(kinkTensor(1.0, 2.0));
        unsignedTensor(2.0, 1.0) = unsignedTensor(1.0, 2.0);
        initialTensor = initialRotation*unsignedTensor*initialRotation';
        direction = 1.0;
        if kinkTensor(1.0, 2.0) < 0.0
            direction = -1.0;
        end
        angle = direction*(initialAngle + abs(initialTensor(1.0, 2.0))/material(7.0, i));
        rotation = getRotation(angle);
        misalignedTensor = rotation*kinkTensor*rotation';
        criterion = getCriterion(misalignedTensor, material(:, i));
        if stress(1.0, i) <= -0.5*material(2.0, i)
            result(1.0, i) = max(result(1.0, i), criterion);
        end
        if (stress(1.0, i) >= -0.5*material(2.0, i)) && (stress(1.0, i) <= 0.0)
            result(2.0, i) = max(result(2.0, i), criterion);
        end
    end
end
end

%% INDEPENDENT MATRIX REFERENCE USING PLANE TRACTIONS
function [result] = referenceMatrixFailure(stress, material, step, maximumAngle)
result = zeros(1.0, size(stress, 2.0));
for i = 1.0:size(stress, 2.0)
    tensor = [stress(1.0, i), stress(3.0, i), 0.0; stress(3.0, i), stress(2.0, i), 0.0; 0.0, 0.0, 0.0];
    for alpha = (0.0:step:maximumAngle)*pi/180.0
        normal = [0.0; cos(alpha); sin(alpha)];
        transverse = [0.0; -sin(alpha); cos(alpha)];
        traction = tensor*normal;
        Sn = normal'*traction;
        Tt = transverse'*traction;
        Tl = traction(1.0);
        criterion = sqrt((Tt/(material(6.0, i) - material(9.0, i)*Sn))^2.0 +...
            (Tl/(material(5.0, i) - material(8.0, i)*Sn))^2.0 + (max(Sn, 0.0)/material(3.0, i))^2.0);
        result(i) = max(result(i), criterion);
    end
end
end

%% CALIBRATE AGAINST UNIAXIAL FAILURE USING MATRIX ROTATIONS
function [criterion] = compressionFailure(initialAngle, material)
tensor = diag([-material(2.0), 0.0, 0.0]);
rotation = getRotation(initialAngle);
initialTensor = rotation*tensor*rotation';
angle = initialAngle + abs(initialTensor(1.0, 2.0))/material(7.0);
rotation = getRotation(angle);
criterion = getCriterion(rotation*tensor*rotation', material);
end

function [rotation] = getRotation(angle)
rotation = [cos(angle), sin(angle), 0.0; -sin(angle), cos(angle), 0.0; 0.0, 0.0, 1.0];
end

function [criterion] = getCriterion(tensor, material)
criterion = sqrt((tensor(2.0, 3.0)/(material(6.0) - material(9.0)*tensor(2.0, 2.0)))^2.0 +...
    (tensor(1.0, 2.0)/(material(5.0) - material(8.0)*tensor(2.0, 2.0)))^2.0 + (max(tensor(2.0, 2.0), 0.0)/material(3.0))^2.0);
end

%% COMPARE NUMERIC OUTPUT WITH THE EXPECTED VALUES
function [] = checkValues(actual, expected, description)
tolerance = 1e-9*max(1.0, max(abs(expected(:))));
assert(isequal(size(actual), size(expected)) == true, '%s: output dimensions differ.', description);
assert(all(isfinite(actual(:))) == true, '%s: output contains non-finite values.', description);
assert(all(abs(actual(:) - expected(:)) <= tolerance) == true, '%s: incorrect output.', description);
end
