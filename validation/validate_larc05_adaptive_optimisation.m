function [] = validate_larc05_adaptive_optimisation()
%   Regression checks for LaRC05 DEFAULT step selection (issue #105).
%
%   RUN THIS FUNCTION FROM THE REPOSITORY ROOT:
%   >> addpath('validation'); validate_larc05_adaptive_optimisation
%
%   Existing parallel pool settings apply. No analysis files are written.
%
%   Layup Analysis Tool 5.2.0 Copyright Louis Vallance 2026
%
%% - DO NOT EDIT BELOW LINE
%__________________________________________________________________________
%%
material = [100000.0, 70000.0; 10000.0, 18000.0; 5000.0, 8000.0; 0.3, 0.25; zeros(4.0, 2.0)];
failure = repmat([2000.0; 1500.0; 50.0; 200.0; 80.0; 70.0; 5000.0; 0.0; 0.0; 53.0; 2.0], 1.0, 2.0);
failure(1.0:7.0, 2.0) = 1.2*failure(1.0:7.0, 2.0);
failure = repelem(failure, 1.0, 3.0);
mechanical = repelem(material, 1.0, 3.0);
loadMech = [-37.7573071261914, -155.527922595307, 65.2910196661923;...
    -7.14923275767514, -29.1095619668018, -37.0620787477615];
sequences = [0.0, 0.0; 0.0, 45.0; 0.0, 90.0; 45.0, 0.0; 45.0, 45.0; 45.0, 90.0; 90.0, 0.0; 90.0, 45.0; 90.0, 90.0];
[z, ~, error] = abd.internal_getThickness(2.0, [0.4, 0.6], 1e-6);
assert(error == false, 'Invalid validation thickness.');
[error, zPoints, ~, totalPoints, A11, A22, B11, B22] = abd.internal_getSectionPoints(3.0, '', 2.0, [0.0, 0.0], z,...
    material(5.0, :), material(6.0, :), material(7.0, :), material(8.0, :), 1e-6, []);
assert(error == false, 'Invalid validation section points.');
[Q11, Q22, Q66, Q12] = abd.internal_getReducedQ(2.0, material(1.0, :), material(2.0, :), material(4.0, :), material(3.0, :));

%% COMPARE ADAPTIVE ANALYSIS AND FIXED OPTIMISATION SEARCHES
% This case changes winner between the 15-degree optimisation grid and the
% adaptive analysis. Each setting must retain its own candidate scores.
scores = zeros(2.0, size(sequences, 1.0));
denseScores = scores;
coarseScores = scores;
customScores = scores;
candidateStress = zeros(3.0, totalPoints, size(sequences, 1.0));
for i = 1.0:size(sequences, 1.0)
    [~, ~, ~, stress] = abd.internal_optimise.getOptiStressStrain(sequences(i, :), 3.0, 2.0, z, [0.0, 0.0], 1e-6,...
        Q11, Q12, Q66, Q22, A11, A22, B11, B22, 0.0, 0.0, loadMech(1.0, 1.0), loadMech(1.0, 2.0), loadMech(1.0, 3.0),...
        loadMech(2.0, 1.0), loadMech(2.0, 2.0), loadMech(2.0, 3.0), totalPoints, zPoints);
    candidateStress(:, :, i) = stress;
    values = evaluate(stress, failure, 'DEFAULT', 1.0);
    dense = evaluate(stress, failure, 0.05, 1.0);
    checkValues(values, dense, 1e-4, 'Dense angular reference');
    coarse = evaluate(stress, failure, 15.0, 1.0);
    custom = evaluate(stress, failure, 7.0, 1.0);
    checkValues(evaluate(stress, failure, 'DEFAULT', 2.0), coarse, 1e-9, 'Default optimisation step');
    checkValues(evaluate(stress, failure, 7.0, 2.0), custom, 1e-9, 'Explicit optimisation step');
    scores(:, i) = getObjectives(values);
    denseScores(:, i) = getObjectives(dense);
    coarseScores(:, i) = getObjectives(coarse);
    customScores(:, i) = getObjectives(custom);
end
[~, coarseIndex] = min(coarseScores(1.0, :));
[~, adaptiveIndex] = min(scores(1.0, :));
[~, denseIndex] = min(denseScores(1.0, :));
assert(coarseIndex ~= denseIndex, 'The case must distinguish fixed and adaptive searches.');
assert(adaptiveIndex == denseIndex, 'Adaptive analysis must recover the reference winner.');

%% CHECK ALL OPTIMISATION MODES AND BOTH OBJECTIVES
% Explicit chunk size also covers a partial final chunk (nine candidates).
unused = cell(1.0, 19.0);
principal = zeros(1.0, totalPoints);
steps = {'DEFAULT', 15.0, 7.0};
expectedScores = {coarseScores, coarseScores, customScores};
notice = 'The optimiser will use a fixed default step size of 15 degrees for LaRC05';
for mode = 1.0:3.0
    for objective = 1.0:2.0
        inputs = [{{true, 'LARC05', 2.0, objective, [0.0, 45.0, 90.0]}, 2.0, totalPoints, 3.0, z, zPoints, Q11, Q22, Q66, Q12,...
            A11, A22, B11, B22, 1e-6}, unused, {failure(1.0, :), failure(2.0, :), failure(3.0, :), failure(4.0, :), failure(5.0, :),...
            failure(6.0, :), failure(7.0, :), failure(8.0, :), failure(9.0, :), failure(10.0, :), failure(11.0, :), 0.0, 0.0,...
            loadMech(1.0, 1.0), loadMech(1.0, 2.0), loadMech(1.0, 3.0), loadMech(2.0, 1.0), loadMech(2.0, 2.0), loadMech(2.0, 3.0),...
            mechanical(1.0, :), mechanical(2.0, :), mechanical(4.0, :), mechanical(3.0, :), false, principal, principal, principal,...
            3.0, [], {mode, 2.0, 'DEFAULT'}, false, 'DEFAULT'}]; %#ok<NASGU>
        for stepIndex = 1.0:length(steps)
            inputs{end} = steps{stepIndex};
            best = [];
            actual = [];
            optimisationOutput = evalc('[best, actual] = abd.internal_optimise.main(inputs{:});');
            if isempty(best) == true || isempty(best{5.0}) == false
                fprintf('%s', optimisationOutput);
            end
            assert(isempty(best) == false && isempty(best{5.0}) == true, 'Optimisation failed.');
            assert(count(optimisationOutput, notice) == double(stepIndex == 1.0), 'Incorrect number of default step notices.');
            expected = expectedScores{stepIndex}(objective, :);
            [bestScore, bestIndex] = min(expected);
            checkValues(actual, expected, 1e-9, 'All fixed-grid candidate scores');
            checkValues(best{2.0}, bestScore, 1e-9, 'Fixed-grid winning score');
            checkValues(best{1.0}, sequences(bestIndex, :), 1e-9, 'Fixed-grid winning sequence');
            checkValues(best{6.0}.STRESS_PLY, candidateStress(:, :, bestIndex), 1e-9, 'Fixed-grid winning stress');
        end
    end
end

%% DO NOT REPORT A DEFAULT OPTIMISATION STEP WHEN OPTIMISATION IS DISABLED
inputs{1.0}{1.0} = false;
inputs{end} = 'DEFAULT';
optimisationOutput = evalc('[best, actual] = abd.internal_optimise.main(inputs{:});');
assert(isempty(actual) == true && all(cellfun(@isempty, best)) == true, 'Disabled optimisation must not evaluate candidates.');
assert(contains(optimisationOutput, notice) == false, 'Disabled optimisation must not report a default step.');

fprintf('[PASS] LaRC05 DEFAULT selection: adaptive analysis, fixed optimisation and single default step notice in all modes\n');
end

%% EVALUATE ALL FIVE CRITERIA
function [values] = evaluate(stress, material, step, mode)
nPoints = size(stress, 2.0);
meanStress = 0.5.*(stress(1.0, :) + stress(2.0, :));
radius = sqrt((0.5.*(stress(1.0, :) - stress(2.0, :))).^2.0 + stress(3.0, :).^2.0);
[PF, MF, KF, SF, TF] = abd.internal_getLaRC05(nPoints, stress, false, meanStress + radius, meanStress - radius, zeros(1.0, nPoints),...
    material(7.0, :), material(1.0, :), material(2.0, :), material(3.0, :), material(4.0, :), material(5.0, :), material(6.0, :),...
    material(10.0, :), material(11.0, :), material(8.0, :), material(9.0, :), 3.0, step, mode);
values = [PF; MF; KF; SF; TF];
end

%% GET THE MINMAX AND MINMEAN OBJECTIVES
function [scores] = getObjectives(values)
worst = max(values, [], 1.0);
scores = [max(worst); mean(worst)];
end

%% COMPARE NUMERIC OUTPUT WITH THE EXPECTED VALUES
function [] = checkValues(actual, expected, relativeTolerance, description)
tolerance = relativeTolerance*max(1.0, max(abs(expected(:))));
assert(isequal(size(actual), size(expected)) == true, '%s: output dimensions differ.', description);
assert(all(isfinite(actual(:))) == true, '%s: output contains non-finite values.', description);
assert(all(abs(actual(:) - expected(:)) <= tolerance) == true, '%s: incorrect output.', description);
end
