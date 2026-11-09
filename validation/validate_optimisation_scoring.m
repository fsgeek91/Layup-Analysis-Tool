function [] = validate_optimisation_scoring()
%   Regression checks for lightweight optimisation scoring (issue #107).
%
%   RUN THIS FUNCTION FROM THE REPOSITORY ROOT:
%   >> addpath('validation'); validate_optimisation_scoring
%
%   Existing parallel pool settings apply. No analysis files are written.
%
%   Layup Analysis Tool 5.2.0 Copyright Louis Vallance 2026
%
%% - DO NOT EDIT BELOW LINE
%__________________________________________________________________________
%%
%% COMPARE PLY STRESS WITH THE FULL TENSOR PATH
% Rows contain mechanical-load scale, temperature change and moisture change.
loadCases = [1.0, 0.0, 0.0; 0.0, 80.0, 0.0; 0.0, 0.0, 0.15;...
    1.0, 80.0, 0.15; 0.0, 0.0, 0.0; 1.0, -30.0, -0.05; 1e-8, 0.0, 0.0];
for nPlies = [1.0, 2.0]
    material = getMaterial(nPlies);
    theta = [17.0, -43.0];
    for nPoints = [1.0, 3.0, 5.0]
        for i = 1.0:size(loadCases, 1.0)
            inputs = getTensorInputs(material, theta(1.0:nPlies), nPoints, loadCases(i, :));
            checkStress(inputs, 'Mechanical/environmental stress recovery');
        end
    end
end

%% PRESERVE THE ORDER OF STRAIN AND STRESS CLEANUP
% A sub-threshold midspan strain is removed before stress recovery.
inputs = getPrescribedInputs([0.5e-12; 2e-12; zeros(4.0, 1.0)], 0.0, 0.0, 1e8*eye(3.0));
actual = checkStress(inputs, 'Midspan strain cleanup');
checkValues(actual, [0.0; 2e-4; 0.0], 'Midspan cleanup result');

% A small point strain produced by curvature must still generate stress.
inputs = getPrescribedInputs([zeros(3.0, 1.0); 1e-10; 0.0; 0.0], 0.005, 0.0, 1e8*eye(3.0));
actual = checkStress(inputs, 'Point strain retained during stress recovery');
checkValues(actual, [5e-5; 0.0; 0.0], 'Curvature-generated stress');

% Global components below the stress tolerance can combine above it when
% rotated. Only the final ply-aligned stress is cleaned in the light path.
inputs = getPrescribedInputs([0.75e-6; 0.75e-6; 0.75e-6; zeros(3.0, 1.0)], 0.0, 45.0, eye(3.0));
actual = checkStress(inputs, 'Stress cleanup after rotation');
checkValues(actual, [1.5e-6; 0.0; 0.0], 'Rotated stress retained');

%% COMPARE EVERY BUILT-IN CRITERION AGAINST FULL-TENSOR SCORING
material = getMaterial(2.0);
nPoints = 3.0;
failure = getFailure(nPoints);
angles = [0.0, 37.0, 90.0];
criteria = {'MSTRS', 'TSAIH', 'HOFFMAN', 'TSAIW', 'AZZIT', 'MSTRN', 'HASHIN', 'LARC05'};
nCandidates = numel(angles)^2.0;
sequences = zeros(nCandidates, 2.0);
tensors = cell(1.0, nCandidates);
for i = 1.0:nCandidates
    indexes = [floor((i - 1.0)/numel(angles)) + 1.0, mod(i - 1.0, numel(angles)) + 1.0];
    sequences(i, :) = angles(indexes);
    inputs = getTensorInputs(material, sequences(i, :), nPoints, [1.0, 80.0, 0.15]);
    tensors{i} = getFullTensors(inputs);
end
for parameter = 1.0:2.0
    for criterion = 1.0:length(criteria)
        scores = zeros(2.0, nCandidates);
        for i = 1.0:nCandidates
            values = getCriterion(criteria{criterion}, parameter, tensors{i}.STRESS_PLY, material, failure, nPoints);
            scores(:, i) = [max(values); mean(values)];
        end
        for mode = 1.0:3.0
            for objective = 1.0:2.0
                [best, actual] = runOptimisation(material, failure, nPoints, angles, mode, objective, criteria{criterion}, parameter);
                expected = scores(objective, :);
                [bestScore, bestIndex] = min(expected);
                assert(isempty(best{5.0}) == true, 'Lightweight optimisation failed.');
                checkValues(actual, expected, ['Complete candidate buffer: ', criteria{criterion}]);
                checkValues(best{2.0}, bestScore, 'Winning score');
                checkValues(best{1.0}, sequences(bestIndex, :), 'Winning sequence');
                names = fieldnames(tensors{bestIndex});
                for i = 1.0:length(names)
                    checkValues(best{6.0}.(names{i}), tensors{bestIndex}.(names{i}), ['Winning tensor: ', names{i}]);
                end
            end
        end
    end
end

fprintf('[PASS] Lightweight scoring: stress recovery, cleanup order, all built-in criteria and full winning tensors in all modes\n');
end

%% GET HETEROGENEOUS MECHANICAL AND EXPANSION PROPERTIES
function [material] = getMaterial(nPlies)
material = [100000.0, 70000.0; 10000.0, 18000.0; 5000.0, 8000.0; 0.3, 0.25;...
    1e-5, 2e-5; 2e-5, 3e-5; 0.002, 0.003; 0.003, 0.004];
material = material(:, 1.0:nPlies);
end

%% GET SECTION-POINT PROPERTIES FOR ALL BUILT-IN CRITERIA
function [failure] = getFailure(nPoints)
failure.STRESS = repmat([1000.0; -500.0; 50.0; -150.0; 100.0; -0.5; 0.0], 1.0, 2.0);
failure.STRESS(1.0:5.0, 2.0) = 0.8*failure.STRESS(1.0:5.0, 2.0);
failure.STRAIN = [0.01, 0.012; -0.006, -0.008; 0.02, 0.025; -0.01, -0.012; 0.04, 0.03];
failure.HASHIN = repmat([0.7; 1000.0; 500.0; 50.0; 150.0; 100.0; 80.0], 1.0, 2.0);
failure.HASHIN(2.0:7.0, 2.0) = 0.8*failure.HASHIN(2.0:7.0, 2.0);
failure.LARC05 = repmat([2000.0; 1500.0; 50.0; 200.0; 80.0; 70.0; 5000.0; 0.2; 0.3; 53.0; 2.0], 1.0, 2.0);
failure.LARC05(1.0:7.0, 2.0) = 1.2*failure.LARC05(1.0:7.0, 2.0);
names = fieldnames(failure);
for i = 1.0:length(names)
    failure.(names{i}) = repelem(failure.(names{i}), 1.0, nPoints);
end
end

%% ASSEMBLE THE ORIGINAL FULL-TENSOR INPUTS
function [inputs] = getTensorInputs(material, theta, nPoints, loadCase)
nPlies = size(material, 2.0);
thickness = [0.4, 0.6];
[z, ~, error] = abd.internal_getThickness(nPlies, thickness(1.0:nPlies), 1e-6);
assert(error == false, 'Invalid validation thickness.');
[error, zPoints, thetaPoints, totalPoints, A11, A22, B11, B22] = abd.internal_getSectionPoints(nPoints, '', nPlies, theta, z,...
    material(5.0, :), material(6.0, :), material(7.0, :), material(8.0, :), 1e-6, []);
assert(error == false, 'Invalid validation section points.');
[Q11, Q22, Q66, Q12] = abd.internal_getReducedQ(nPlies, material(1.0, :), material(2.0, :), material(4.0, :), material(3.0, :));
[Q11t, Q12t, Q16t, Q22t, Q26t, Q66t] = abd.internal_getTransformedQ(nPlies, theta, Q11, Q12, Q66, Q22);
[axx, ayy, axy, bxx, byy, bxy] = abd.internal_getThermoHydro(thetaPoints, A11, A22, B11, B22);
deltaT = loadCase(2.0);
deltaM = loadCase(3.0);
[ABD, ~, Qt, NxxT, NyyT, NxyT, MxxT, MyyT, MxyT, NxxM, NyyM, NxyM, MxxM, MyyM, MxyM] =...
    abd.internal_getABD(nPlies, Q11t, Q12t, Q16t, Q22t, Q26t, Q66t, z, deltaT, deltaM, axx, ayy, axy, bxx, byy, bxy, nPoints);
loadMech = loadCase(1.0)*[-150.0, 30.0, 20.0; 20.0, -10.0, 5.0];
inputs = {ABD, loadMech(1.0, 1.0), NxxT, NxxM, loadMech(1.0, 2.0), NyyT, NyyM, loadMech(1.0, 3.0), NxyT, NxyM,...
    loadMech(2.0, 1.0), MxxT, MxxM, loadMech(2.0, 2.0), MyyT, MyyM, loadMech(2.0, 3.0), MxyT, MxyM,...
    totalPoints, zPoints, thetaPoints, Qt, deltaT, deltaM, axx, ayy, axy, bxx, byy, bxy, 1e-6};
end

%% PRESCRIBE MIDSPAN STRAINS AND CURVATURES FOR CLEANUP CHECKS
function [inputs] = getPrescribedInputs(midspan, z, theta, Qt)
inputs = {eye(6.0), midspan(1.0), 0.0, 0.0, midspan(2.0), 0.0, 0.0, midspan(3.0), 0.0, 0.0,...
    midspan(4.0), 0.0, 0.0, midspan(5.0), 0.0, 0.0, midspan(6.0), 0.0, 0.0, 1.0, z, theta, Qt,...
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1e-6};
end

%% COMPARE LIGHTWEIGHT AND FULL STRESS RECOVERY
function [actual] = checkStress(inputs, description)
[~, ~, ~, ~, expected] = abd.internal_getTensor(inputs{:});
actual = abd.internal_optimise.getPlyStress(inputs{:});
checkValues(actual, expected, description);
assert(isequal(actual == 0.0, expected == 0.0) == true, '%s: zero cleanup differs.', description);
end

%% RECOVER THE COMPLETE TENSORS RETAINED FOR THE WINNER
function [tensors] = getFullTensors(inputs)
[~, tensors.STRAIN_XY, tensors.STRESS_XY, tensors.STRAIN_PLY, tensors.STRESS_PLY] = abd.internal_getTensor(inputs{:});
end

%% EVALUATE CRITERIA USING THE ORIGINAL FULL-TENSOR STRESSES
function [values] = getCriterion(criterion, parameter, stress, material, failure, nPoints)
nTotal = size(stress, 2.0);
mechanical = repelem(material, 1.0, nPoints);
F = failure.STRESS;
switch criterion
    case 'MSTRS'
        values = abd.internal_strength.getMstrs(stress, F(1.0, :), F(2.0, :), F(3.0, :), F(4.0, :), F(5.0, :));
    case 'TSAIH'
        values = abd.internal_strength.getTsaih(parameter, stress, F(1.0, :), F(2.0, :), F(3.0, :), F(4.0, :), F(5.0, :));
    case 'HOFFMAN'
        values = abd.internal_strength.getHoffman(parameter, stress, F(1.0, :), F(2.0, :), F(3.0, :), F(4.0, :), F(5.0, :), zeros(1.0, nTotal), -ones(1.0, nTotal));
    case 'TSAIW'
        values = abd.internal_strength.getTsaiw(parameter, stress, F(1.0, :), F(2.0, :), F(3.0, :), F(4.0, :), F(5.0, :), F(6.0, :), F(7.0, :));
    case 'AZZIT'
        values = abd.internal_strength.getAzzit(parameter, nTotal, stress, F(1.0, :), F(2.0, :), F(3.0, :), F(4.0, :), F(5.0, :));
    case 'MSTRN'
        F = failure.STRAIN;
        values = abd.internal_strength.getMstrn(nTotal, stress, mechanical(1.0, :), mechanical(2.0, :), mechanical(4.0, :), mechanical(3.0, :),...
            F(1.0, :), F(2.0, :), F(3.0, :), F(4.0, :), F(5.0, :));
    case 'HASHIN'
        F = failure.HASHIN;
        [FT, FC, MT, MC] = abd.internal_strength.getHashin(nTotal, stress, F(1.0, :), F(2.0, :), F(3.0, :), F(4.0, :), F(5.0, :), F(6.0, :), F(7.0, :));
        values = max([FT; FC; MT; MC], [], 1.0);
    case 'LARC05'
        F = failure.LARC05;
        meanStress = 0.5*(stress(1.0, :) + stress(2.0, :));
        radius = sqrt((0.5*(stress(1.0, :) - stress(2.0, :))).^2.0 + stress(3.0, :).^2.0);
        [PF, MF, KF, SF, TF] = abd.internal_getLaRC05(nTotal, stress, false, meanStress + radius, meanStress - radius, zeros(1.0, nTotal),...
            F(7.0, :), F(1.0, :), F(2.0, :), F(3.0, :), F(4.0, :), F(5.0, :), F(6.0, :), F(10.0, :), F(11.0, :), F(8.0, :), F(9.0, :), nPoints, 30.0, 1.0);
        values = max([PF; MF; KF; SF; TF], [], 1.0);
end
end

%% RUN THE OPTIMISER WITH ALL FAILURE PROPERTY SETS AVAILABLE
function [best, scores] = runOptimisation(material, failure, nPoints, angles, mode, objective, criterion, parameter)
nPlies = size(material, 2.0);
[z, ~, error] = abd.internal_getThickness(nPlies, [0.4, 0.6], 1e-6);
assert(error == false, 'Invalid validation thickness.');
[error, zPoints, ~, totalPoints, A11, A22, B11, B22] = abd.internal_getSectionPoints(nPoints, '', nPlies, zeros(1.0, nPlies), z,...
    material(5.0, :), material(6.0, :), material(7.0, :), material(8.0, :), 1e-6, []);
assert(error == false, 'Invalid validation section points.');
[Q11, Q22, Q66, Q12] = abd.internal_getReducedQ(nPlies, material(1.0, :), material(2.0, :), material(4.0, :), material(3.0, :));
mechanical = repelem(material, 1.0, nPoints);
properties = [reshape(num2cell(failure.STRESS, 2.0), 1.0, []), reshape(num2cell(failure.STRAIN, 2.0), 1.0, []),...
    reshape(num2cell(failure.HASHIN, 2.0), 1.0, []), reshape(num2cell(failure.LARC05, 2.0), 1.0, [])];
inputs = [{{true, criterion, parameter, objective, angles}, nPlies, totalPoints, nPoints, z, zPoints, Q11, Q22, Q66, Q12,...
    A11, A22, B11, B22, 1e-6}, properties, {80.0, 0.15, -150.0, 30.0, 20.0, 20.0, -10.0, 5.0,...
    mechanical(1.0, :), mechanical(2.0, :), mechanical(4.0, :), mechanical(3.0, :), false, [], [], [], nPoints,...
    @unexpectedCallback, {mode, 2.0, 'DEFAULT'}, false, 30.0}]; %#ok<NASGU>
best = [];
scores = [];
optimisationOutput = evalc('[best, scores] = abd.internal_optimise.main(inputs{:});');
if isempty(best) == true || isempty(best{5.0}) == false
    fprintf('%s', optimisationOutput);
end
end

%% A REGISTERED CALLBACK MUST NOT OVERRIDE THE SELECTED BUILT-IN CRITERION
function [UCRT] = unexpectedCallback(INFO, UCRT, MATERIAL_MECH, MATERIAL_FAIL, TENSORS) %#ok<INUSD>
error('A built-in scoring request must not evaluate the registered user criterion.');
end

%% COMPARE NUMERIC OUTPUT WITH THE EXPECTED VALUES
function [] = checkValues(actual, expected, description)
tolerance = 1e-10*max(1.0, max(abs(expected(:))));
assert(isequal(size(actual), size(expected)) == true, '%s: output dimensions differ.', description);
assert(all(isfinite(actual(:))) == true, '%s: output contains non-finite values.', description);
assert(all(abs(actual(:) - expected(:)) <= tolerance) == true, '%s: incorrect output.', description);
end
