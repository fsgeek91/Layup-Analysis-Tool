function [] = validate_optimisation_retention()
%   Regression checks for bounded optimisation score retention (issue #111).
%
%   RUN THIS FUNCTION FROM THE REPOSITORY ROOT:
%   >> addpath('validation'); validate_optimisation_retention
%
%   Existing parallel pool settings apply. No analysis files are written.
%
%   Layup Analysis Tool 5.1.5 Copyright Louis Vallance 2026
%
%% - DO NOT EDIT BELOW LINE
%__________________________________________________________________________
%%
%% PRESERVE EXACT WINNERS WITH COMPLETE, SAMPLED AND OMITTED SCORES
angles = [0.0, 17.0, 90.0];
nCandidates = numel(angles)^4.0;
options = getOptions(nCandidates);
limits = [0.0, 1.0, 7.0, Inf, nCandidates, nCandidates + 5.0];
for objective = 1.0:2.0
    inputs = getInputs(4.0, angles, objective, []);
    reference = runOptimisation(inputs, options{2.0}, []);
    for method = 1.0:length(options)
        complete = runOptimisation(inputs, options{method}, []);
        checkWinner(complete, reference);
        checkValues(complete{2.0}, reference{2.0}, 'Default complete score buffer');
        for limit = limits
            actual = runOptimisation(inputs, options{method}, limit);
            checkWinner(actual, reference);
            checkSamples(actual{2.0}, reference{2.0}, min(limit, nCandidates));
            if limit >= nCandidates
                checkValues(actual{2.0}, reference{2.0}, 'Complete retained score buffer');
            end
        end
    end
end

%% PRESERVE SINGLE-CANDIDATE SEARCHES
inputs = getInputs(4.0, 23.0, 1.0, []);
reference = runOptimisation(inputs, options{2.0}, []);
for method = 1.0:length(options)
    for limit = [0.0, 1.0, 7.0]
        actual = runOptimisation(inputs, options{method}, limit);
        checkWinner(actual, reference);
        checkSamples(actual{2.0}, reference{2.0}, min(limit, 1.0));
    end
end

%% PRESERVE USER-CRITERION TIES AND NONFINITE SCORES
% The callbacks derive the candidate index from its expansion coefficients,
% so their results do not depend on worker scheduling or callback counters.
callbacks = {@zeroCriterion, @infiniteCriterion, @negativeInfiniteCriterion, @nanCriterion, @indexCriterion};
expectedIndexes = [1.0, 1.0, 14.0, 50.0, 1.0];
for callback = 1.0:length(callbacks)
    for objective = 1.0:2.0
        inputs = getInputs(4.0, [0.0, 45.0, 90.0], objective, callbacks{callback});
        reference = runOptimisation(inputs, options{2.0}, []);
        assert(reference{3.0} == expectedIndexes(callback), 'The user criterion must exercise its expected winner.');
        for method = 1.0:length(options)
            for limit = [0.0, 1.0, 7.0]
                actual = runOptimisation(inputs, options{method}, limit);
                checkWinner(actual, reference);
                checkSamples(actual{2.0}, reference{2.0}, limit);
            end
        end
    end
end

%% KEEP SAMPLES DETERMINISTIC WITHOUT CHANGING THE GLOBAL RANDOM STATE
inputs = getInputs(4.0, [0.0, 45.0, 90.0], 1.0, @indexCriterion);
reference = runOptimisation(inputs, options{2.0}, 7.0);
randomState = rng;
for method = 1.0:length(options)
    actual = runOptimisation(inputs, options{method}, 7.0);
    checkValues(actual{2.0}, reference{2.0}, 'Deterministic samples across optimisation modes');
end
assert(isequal(rng, randomState) == true, 'Sampling must not change the global random state.');
assert(all(diff(reference{2.0}) > 0.0) == true, 'Samples must use distinct candidates in their original order.');

%% PRESERVE THE FIRST WINNER ACROSS BOUNDED CHUNK BATCHES
% 2187 one-item chunks span more than two batches of 1024. Equal minima
% at 1025 and 2049 must select 1025, regardless of worker completion order.
for objective = 1.0:2.0
    inputs = getInputs(7.0, [0.0, 45.0, 90.0], objective, @batchCriterion);
    reference = runOptimisation(inputs, options{2.0}, []);
    assert(reference{3.0} == 1025.0, 'The batch case must select the first tied minimum.');
    for limit = [0.0, 7.0]
        actual = runOptimisation(inputs, {3.0, 1.0, 'DEFAULT'}, limit);
        checkWinner(actual, reference);
        checkSamples(actual{2.0}, reference{2.0}, limit);
    end
end

fprintf('[PASS] Score retention: exact winners, bounded ordered samples, user criteria, nonfinite scores and chunk-batch ties\n');
end

%% EXERCISE ALL SOLVERS AND SMALL, PARTIAL AND OVERSIZED CHUNKS
function [options] = getOptions(nCandidates)
options = {{1.0, 'DEFAULT', 'DEFAULT'}, {2.0, 'DEFAULT', 'DEFAULT'},...
    {3.0, 1.0, 'DEFAULT'}, {3.0, 7.0, 'DEFAULT'}, {3.0, nCandidates + 5.0, 'DEFAULT'}};
end

%% GET HETEROGENEOUS OPTIMISATION INPUTS
function [inputs] = getInputs(nPlies, angles, objective, callback)
nPoints = 3.0;
material = [linspace(100000.0, 70000.0, nPlies); linspace(10000.0, 18000.0, nPlies);...
    linspace(5000.0, 8000.0, nPlies); linspace(0.3, 0.25, nPlies);...
    linspace(1e-5, 2e-5, nPlies); linspace(2e-5, 4e-5, nPlies);...
    linspace(0.002, 0.003, nPlies); linspace(0.003, 0.004, nPlies)];
[z, ~, error] = abd.internal_getThickness(nPlies, linspace(0.2, 0.5, nPlies), 1e-6);
assert(error == false, 'Invalid validation thickness.');
[error, zPoints, ~, totalPoints, A11, A22, B11, B22] = abd.internal_getSectionPoints(nPoints, '', nPlies, zeros(1.0, nPlies), z,...
    material(5.0, :), material(6.0, :), material(7.0, :), material(8.0, :), 1e-6, []);
assert(error == false, 'Invalid validation section points.');
[Q11, Q22, Q66, Q12] = abd.internal_getReducedQ(nPlies, material(1.0, :), material(2.0, :), material(4.0, :), material(3.0, :));
mechanical = repelem(material, 1.0, nPoints);
strength = repmat([1000.0; -500.0; 50.0; -150.0; 100.0; -0.5; 0.0], 1.0, totalPoints);
properties = reshape(num2cell(strength, 2.0), 1.0, []);
unused = cell(1.0, 23.0);
criterion = 'MSTRS';
if isempty(callback) == false
    criterion = 'UCRT';
end
inputs = [{{true, criterion, 1.0, objective, angles}, nPlies, totalPoints, nPoints, z, zPoints, Q11, Q22, Q66, Q12,...
    A11, A22, B11, B22, 1e-6}, properties, unused, {80.0, 0.15, -150.0, 30.0, 20.0, 20.0, -10.0, 5.0,...
    mechanical(1.0, :), mechanical(2.0, :), mechanical(4.0, :), mechanical(3.0, :), false, [], [], [], nPoints,...
    callback, {2.0, 'DEFAULT', 'DEFAULT'}, false, 15.0}];
end

%% RUN WITH THE LEGACY INTERFACE OR AN EXPLICIT RETENTION LIMIT
function [result] = runOptimisation(inputs, options, scoreLimit)
inputs{end - 2.0} = options;
if isempty(scoreLimit) == false
    inputs{end + 1.0} = scoreLimit; %#ok<NASGU>
end
result = cell(1.0, 6.0);
optimisationOutput = evalc('[result{:}] = abd.internal_optimise.main(inputs{:});');
if isempty(result{1.0}) == true || isempty(result{1.0}{5.0}) == false
    fprintf('%s', optimisationOutput);
end
assert(isempty(result{1.0}) == false && isempty(result{1.0}{5.0}) == true, 'Optimisation failed.');
end

%% COMPARE THE EXACT WINNER AND ITS COMPLETE TENSORS
function [] = checkWinner(actual, expected)
checkValues(actual{3.0}, expected{3.0}, 'Winning candidate index');
for i = 1.0:3.0
    checkValues(actual{1.0}{i}, expected{1.0}{i}, 'Winning sequence, score or candidate count');
end
names = fieldnames(expected{1.0}{6.0});
for i = 1.0:length(names)
    checkValues(actual{1.0}{6.0}.(names{i}), expected{1.0}{6.0}.(names{i}), ['Winning tensor: ', names{i}]);
end
end

%% CHECK THAT RETAINED SCORES FORM AN ORDERED SUBSET OF THE FULL RESULTS
function [] = checkSamples(actual, complete, expectedCount)
assert(numel(actual) == expectedCount, 'Incorrect retained score count.');
if isempty(actual) == true
    return
end
assert(isrow(actual) == true, 'Retained scores must form a row vector.');
nextIndex = 1.0;
for i = 1.0:length(actual)
    matches = (complete(nextIndex:end) == actual(i)) | (isnan(complete(nextIndex:end)) & isnan(actual(i)));
    offset = find(matches, 1.0);
    assert(isempty(offset) == false, 'Retained scores must follow the original candidate order.');
    nextIndex = nextIndex + offset;
end
end

%% GET THE CANDIDATE INDEX FROM ITS THERMAL EXPANSION TENSORS
function [index] = getCandidateIndex(INFO, MATERIAL_MECH)
points = 1.0:INFO.SECTION_POINTS:INFO.TOTAL_POINTS;
theta = mod(0.5*atan2d(-MATERIAL_MECH.AXY(points), MATERIAL_MECH.AYY(points) - MATERIAL_MECH.AXX(points)), 180.0);
digits = round(theta/45.0);
index = 0.0;
for i = 1.0:INFO.PLY_COUNT
    index = 3.0*index + digits(i);
end
index = index + 1.0;
end

%% USER CRITERIA WITH CONTROLLED TIES, INFINITIES AND MISSING SCORES
function [UCRT] = zeroCriterion(INFO, UCRT, MATERIAL_MECH, MATERIAL_FAIL, TENSORS) %#ok<INUSD>
UCRT = zeros(1.0, INFO.TOTAL_POINTS);
end

function [UCRT] = infiniteCriterion(INFO, UCRT, MATERIAL_MECH, MATERIAL_FAIL, TENSORS) %#ok<INUSD>
UCRT = Inf(1.0, INFO.TOTAL_POINTS);
end

function [UCRT] = negativeInfiniteCriterion(INFO, UCRT, MATERIAL_MECH, MATERIAL_FAIL, TENSORS) %#ok<INUSD>
index = getCandidateIndex(INFO, MATERIAL_MECH);
value = index;
if index == 14.0 || index == 17.0
    value = -Inf;
end
UCRT = repmat(value, 1.0, INFO.TOTAL_POINTS);
end

function [UCRT] = nanCriterion(INFO, UCRT, MATERIAL_MECH, MATERIAL_FAIL, TENSORS) %#ok<INUSD>
index = getCandidateIndex(INFO, MATERIAL_MECH);
value = abs(index - 50.0);
if index <= 21.0
    value = NaN;
elseif index == 61.0
    value = 0.0;
end
UCRT = repmat(value, 1.0, INFO.TOTAL_POINTS);
end

function [UCRT] = indexCriterion(INFO, UCRT, MATERIAL_MECH, MATERIAL_FAIL, TENSORS) %#ok<INUSD>
UCRT = repmat(getCandidateIndex(INFO, MATERIAL_MECH), 1.0, INFO.TOTAL_POINTS);
end

function [UCRT] = batchCriterion(INFO, UCRT, MATERIAL_MECH, MATERIAL_FAIL, TENSORS) %#ok<INUSD>
index = getCandidateIndex(INFO, MATERIAL_MECH);
UCRT = repmat(double(index ~= 1025.0 && index ~= 2049.0), 1.0, INFO.TOTAL_POINTS);
end

%% COMPARE FINITE VALUES AND PRESERVE NONFINITE RESULTS
function [] = checkValues(actual, expected, description)
assert(isequal(size(actual), size(expected)) == true, '%s: output dimensions differ.', description);
special = isfinite(expected) == false;
assert(isequaln(actual(special), expected(special)) == true, '%s: nonfinite values differ.', description);
finite = expected(special == false);
if isempty(finite) == false
    tolerance = 1e-10*max(1.0, max(abs(finite(:))));
    assert(all(abs(actual(special == false) - finite) <= tolerance) == true, '%s: incorrect output.', description);
end
end
