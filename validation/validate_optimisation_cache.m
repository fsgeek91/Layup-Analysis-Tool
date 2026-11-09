function [] = validate_optimisation_cache()
%   Regression checks for cached ply/angle properties (issue #106).
%
%   RUN THIS FUNCTION FROM THE REPOSITORY ROOT:
%   >> addpath('validation'); validate_optimisation_cache
%
%   Existing parallel pool settings apply. No analysis files are written.
%
%   Layup Analysis Tool 5.2.0 Copyright Louis Vallance 2026
%
%% - DO NOT EDIT BELOW LINE
%__________________________________________________________________________
%%
%% CHECK RAW LOOKUPS AGAINST THE ORIGINAL PROPERTY CALCULATIONS
for nPlies = [1.0, 3.0]
    material = getMaterial(nPlies);
    % Exercise the expansion-coefficient cleanup as well as ordinary values.
    material(5.0:8.0, 1.0) = [0.5e-12; 1.5e-12; 0.25e-12; 2.0e-12];
    [Q11, Q22, Q66, Q12] = abd.internal_getReducedQ(nPlies, material(1.0, :), material(2.0, :), material(4.0, :), material(3.0, :));
    angleLists = {[90.0, 17.0, 0.0, -43.0], 23.0};
    for nPoints = [1.0, 3.0, 5.0]
        expansion = repelem(material(5.0:8.0, :), 1.0, nPoints);
        for angleList = 1.0:length(angleLists)
            angles = angleLists{angleList};
            cache = abd.internal_optimise.getCache(angles, nPlies, nPoints, Q11, Q12, Q66, Q22,...
                expansion(1.0, :), expansion(2.0, :), expansion(3.0, :), expansion(4.0, :));
            for candidate = 1.0:numel(angles)^nPlies
                indexes = getIndexes(candidate, numel(angles), nPlies);
                theta = angles(indexes);
                actual = cell(1.0, 13.0);
                expected = actual;
                [actual{:}] = abd.internal_optimise.getCachedProperties(cache, indexes);
                expected{1.0} = repelem(theta, nPoints);
                [expected{2.0:7.0}] = abd.internal_getTransformedQ(nPlies, theta, Q11, Q12, Q66, Q22);
                [expected{8.0:13.0}] = abd.internal_getThermoHydro(expected{1.0}, expansion(1.0, :), expansion(2.0, :), expansion(3.0, :), expansion(4.0, :));
                for property = 1.0:length(expected)
                    assert(isequal(actual{property}, expected{property}) == true, 'Cached property %.0f differs from direct calculation.', property);
                end
            end
        end
    end
end

%% CHECK ALL CANDIDATE SCORES, BOTH OBJECTIVES AND ALL THREE MODES
% The first two plies share stiffness but have different expansion data.
% The single-ply and single-angle cases exercise enumeration dimensions.
cases = [3.0, 1.0, 3.0; 3.0, 3.0, 3.0; 3.0, 5.0, 3.0; 1.0, 3.0, 3.0; 3.0, 3.0, 1.0];
criteria = {'MSTRS', 'UCRT'};
for caseIndex = 1.0:size(cases, 1.0)
    nPlies = cases(caseIndex, 1.0);
    nPoints = cases(caseIndex, 2.0);
    angles = [0.0, 17.0, 90.0];
    if cases(caseIndex, 3.0) == 1.0
        angles = 17.0;
    end
    material = getMaterial(nPlies);
    nCandidates = numel(angles)^nPlies;
    sequences = zeros(nCandidates, nPlies);
    scores = zeros(2.0, nCandidates, length(criteria));
    tensors = cell(1.0, nCandidates);
    strength = repmat([1000.0; -500.0; 50.0; -150.0; 100.0], 1.0, nPlies*nPoints);
    for i = 1.0:nCandidates
        theta = angles(getIndexes(i, numel(angles), nPlies));
        sequences(i, :) = theta;
        tensors{i} = getUncachedTensors(material, theta, nPoints);
        stressCriterion = abd.internal_strength.getMstrs(tensors{i}.S_PLY_ALIGNED, strength(1.0, :), strength(2.0, :),...
            strength(3.0, :), strength(4.0, :), strength(5.0, :));
        userCriterion = max(abs(tensors{i}.E_PLY_ALIGNED), [], 1.0)/0.01;
        scores(:, i, 1.0) = [max(stressCriterion); mean(stressCriterion)];
        scores(:, i, 2.0) = [max(userCriterion); mean(userCriterion)];
    end
    for mode = 1.0:3.0
        for objective = 1.0:2.0
            for criterion = 1.0:length(criteria)
                [best, actual] = runOptimisation(material, nPoints, angles, mode, objective, criteria{criterion}, strength);
                expected = scores(objective, :, criterion);
                [bestScore, bestIndex] = min(expected);
                assert(isempty(best{5.0}) == true, 'Cached optimisation failed.');
                checkValues(actual, expected, 'Complete candidate buffer');
                checkValues(best{2.0}, bestScore, 'Winning score');
                checkValues(best{1.0}, sequences(bestIndex, :), 'Winning sequence');
                checkValues(best{6.0}.STRESS_XY, tensors{bestIndex}.S_PLY_XY, 'Winning global stress');
                checkValues(best{6.0}.STRESS_PLY, tensors{bestIndex}.S_PLY_ALIGNED, 'Winning ply stress');
                checkValues(best{6.0}.STRAIN_XY, tensors{bestIndex}.E_PLY_XY, 'Winning global strain');
                checkValues(best{6.0}.STRAIN_PLY, tensors{bestIndex}.E_PLY_ALIGNED, 'Winning ply strain');
            end
        end
    end
end

fprintf('[PASS] Optimisation cache: property lookups, section points, full candidate buffers and user tensors in all modes\n');
end

%% GET MECHANICAL AND EXPANSION PROPERTIES FOR THE VALIDATION PLIES
function [material] = getMaterial(nPlies)
material = [100000.0, 100000.0, 70000.0; 10000.0, 10000.0, 18000.0;...
    5000.0, 5000.0, 8000.0; 0.3, 0.3, 0.25; 1e-5, 2e-5, 3e-5;...
    2e-5, 4e-5, 5e-5; 0.002, 0.003, 0.004; 0.003, 0.005, 0.006];
material = material(:, 1.0:nPlies);
end

%% ENUMERATE ANGLE INDICES IN THE ORIGINAL CANDIDATE ORDER
function [indexes] = getIndexes(candidate, nAngles, nPlies)
indexes = zeros(1.0, nPlies);
candidate = candidate - 1.0;
for i = nPlies:-1.0:1.0
    indexes(i) = mod(candidate, nAngles) + 1.0;
    candidate = floor(candidate/nAngles);
end
end

%% COMPUTE ALL NINE TENSORS WITHOUT USING THE OPTIMISATION CACHE
function [tensors, expansion] = getUncachedTensors(material, theta, nPoints)
nPlies = size(material, 2.0);
thickness = [0.3, 0.5, 0.4];
[z, ~, error] = abd.internal_getThickness(nPlies, thickness(1.0:nPlies), 1e-6);
assert(error == false, 'Invalid validation thickness.');
[error, zPoints, thetaPoints, totalPoints, A11, A22, B11, B22] = abd.internal_getSectionPoints(nPoints, '', nPlies, theta, z,...
    material(5.0, :), material(6.0, :), material(7.0, :), material(8.0, :), 1e-6, []);
assert(error == false, 'Invalid validation section points.');
[Q11, Q22, Q66, Q12] = abd.internal_getReducedQ(nPlies, material(1.0, :), material(2.0, :), material(4.0, :), material(3.0, :));
[Q11t, Q12t, Q16t, Q22t, Q26t, Q66t] = abd.internal_getTransformedQ(nPlies, theta, Q11, Q12, Q66, Q22);
[axx, ayy, axy, bxx, byy, bxy] = abd.internal_getThermoHydro(thetaPoints, A11, A22, B11, B22);
expansion = [axx; ayy; axy; bxx; byy; bxy];
[ABD, ~, Qt, NxxT, NyyT, NxyT, MxxT, MyyT, MxyT, NxxM, NyyM, NxyM, MxxM, MyyM, MxyM] =...
    abd.internal_getABD(nPlies, Q11t, Q12t, Q16t, Q22t, Q26t, Q66t, z, 80.0, 0.15, axx, ayy, axy, bxx, byy, bxy, nPoints);
[tensors.E_MIDSPAN, tensors.E_PLY_XY, tensors.S_PLY_XY, tensors.E_PLY_ALIGNED, tensors.S_PLY_ALIGNED,...
    tensors.E_THERM_XY, tensors.E_HYDRO_XY, tensors.E_THERM_ALIGNED, tensors.E_HYDRO_ALIGNED] =...
    abd.internal_getTensor(ABD, 100.0, NxxT, NxxM, 30.0, NyyT, NyyM, 20.0, NxyT, NxyM, 5.0, MxxT, MxxM,...
    -2.0, MyyT, MyyM, 1.0, MxyT, MxyM, totalPoints, zPoints, thetaPoints, Qt, 80.0, 0.15, axx, ayy, axy, bxx, byy, bxy, 1e-6);
end

%% RUN THE OPTIMISER AND RETAIN ITS COMPLETE CANDIDATE BUFFER
function [best, scores] = runOptimisation(material, nPoints, angles, mode, objective, criterion, strength)
nPlies = size(material, 2.0);
thickness = [0.3, 0.5, 0.4];
[z, ~, error] = abd.internal_getThickness(nPlies, thickness(1.0:nPlies), 1e-6);
assert(error == false, 'Invalid validation thickness.');
[error, zPoints, ~, totalPoints, A11, A22, B11, B22] = abd.internal_getSectionPoints(nPoints, '', nPlies, zeros(1.0, nPlies), z,...
    material(5.0, :), material(6.0, :), material(7.0, :), material(8.0, :), 1e-6, []);
assert(error == false, 'Invalid validation section points.');
[Q11, Q22, Q66, Q12] = abd.internal_getReducedQ(nPlies, material(1.0, :), material(2.0, :), material(4.0, :), material(3.0, :));
mechanical = repelem(material, 1.0, nPoints);
unused = cell(1.0, 23.0);
principal = zeros(1.0, totalPoints);
% A chunk size of two exercises the partial final chunk.
inputs = [{{true, criterion, 1.0, objective, angles}, nPlies, totalPoints, nPoints, z, zPoints, Q11, Q22, Q66, Q12,...
    A11, A22, B11, B22, 1e-6, strength(1.0, :), strength(2.0, :), strength(3.0, :), strength(4.0, :), strength(5.0, :), [], []},...
    unused, {80.0, 0.15, 100.0, 30.0, 20.0, 5.0, -2.0, 1.0, mechanical(1.0, :), mechanical(2.0, :), mechanical(4.0, :), mechanical(3.0, :),...
    false, principal, principal, principal, nPoints, @checkUserInputs, {mode, 2.0, 'DEFAULT'}, false, 15.0}]; %#ok<NASGU>
best = [];
scores = [];
optimisationOutput = evalc('[best, scores] = abd.internal_optimise.main(inputs{:});');
if isempty(best) == true || isempty(best{5.0}) == false
    fprintf('%s', optimisationOutput);
end
end

%% CHECK ALL USER-CRITERION TENSORS AND MATERIAL ARRAYS
function [UCRT] = checkUserInputs(INFO, UCRT, MATERIAL_MECH, MATERIAL_FAIL, TENSORS) %#ok<INUSD>
material = getMaterial(INFO.PLY_COUNT);
points = 1.0:INFO.SECTION_POINTS:INFO.TOTAL_POINTS;
% Recover the candidate orientation from its rotated expansion tensor.
% All validation materials have distinct longitudinal/transverse values.
difference = material(5.0, :) - material(6.0, :);
theta = 0.5*atan2d(MATERIAL_MECH.AXY(points)./difference,...
    (MATERIAL_MECH.AXX(points) - MATERIAL_MECH.AYY(points))./difference);
[expected, expansion] = getUncachedTensors(material, theta, INFO.SECTION_POINTS);
names = fieldnames(expected);
for i = 1.0:length(names)
    checkValues(TENSORS.(names{i}), expected.(names{i}), ['User tensor ', names{i}]);
end
names = {'E11', 'E22', 'G12', 'V12'};
for i = 1.0:length(names)
    checkValues(MATERIAL_MECH.(names{i}), repelem(material(i, :), INFO.SECTION_POINTS), ['User property ', names{i}]);
end
names = {'AXX', 'AYY', 'AXY', 'BXX', 'BYY', 'BXY'};
for i = 1.0:length(names)
    checkValues(MATERIAL_MECH.(names{i}), expansion(i, :), ['User expansion ', names{i}]);
end
UCRT = max(abs(TENSORS.E_PLY_ALIGNED), [], 1.0)/0.01;
end

%% COMPARE NUMERIC OUTPUT WITH THE EXPECTED VALUES
function [] = checkValues(actual, expected, description)
tolerance = 1e-9*max(1.0, max(abs(expected(:))));
assert(isequal(size(actual), size(expected)) == true, '%s: output dimensions differ.', description);
assert(all(isfinite(actual(:))) == true, '%s: output contains non-finite values.', description);
assert(all(abs(actual(:) - expected(:)) <= tolerance) == true, '%s: incorrect output.', description);
end
