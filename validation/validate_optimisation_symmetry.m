function [] = validate_optimisation_symmetry()
%   Regression checks for explicitly symmetric optimisation (issue #112).
%
%   RUN THIS FUNCTION FROM THE REPOSITORY ROOT:
%   >> addpath('validation'); validate_optimisation_symmetry
%
%   Existing parallel pool settings apply. No analysis files are written.
%
%   Layup Analysis Tool 5.2.0 Copyright Louis Vallance 2026
%
%% - DO NOT EDIT BELOW LINE
%__________________________________________________________________________
%%
%% COMPARE AGAINST AN INDEPENDENTLY FILTERED EXHAUSTIVE SEARCH
% Rows give physical plies, mirrored properties and section points. Fully
% specified asymmetric properties must stay at their physical ply positions.
cases = [4.0, 1.0, 3.0; 4.0, 0.0, 3.0; 2.0, 1.0, 1.0; 4.0, 1.0, 3.0];
for caseIndex = 1.0:size(cases, 1.0)
    angles = [90.0, 17.0, 0.0];
    if caseIndex == 4.0
        angles = 23.0;
    end
    nPlies = cases(caseIndex, 1.0);
    [sequences, selected] = getSymmetricCandidates(nPlies, angles);
    nCandidates = size(sequences, 1.0);
    options = getOptions(nCandidates);
    for criterion = 1.0:2.0
        for objective = 1.0:2.0
            [inputs, data] = getInputs(nPlies, cases(caseIndex, 2.0), cases(caseIndex, 3.0), angles, objective, criterion);
            complete = runOptimisation(inputs, options{2.0}, Inf, false);
            legacy = runOptimisation(inputs, options{2.0}, [], []);
            checkValues(legacy{2.0}, complete{2.0}, 'Unconstrained legacy score buffer');
            assert(complete{1.0}{3.0} == numel(angles)^nPlies, 'Explicit false must retain the full search.');

            expected = complete{2.0}(selected);
            [bestScore, bestIndex] = min(expected);
            tensors = getTensors(data, sequences(bestIndex, :));
            for method = 1.0:length(options)
                for limit = [Inf, 0.0, 4.0]
                    actual = runOptimisation(inputs, options{method}, limit, true);
                    checkValues(actual{1.0}{1.0}, sequences(bestIndex, :), 'Full mirrored winning sequence');
                    checkValues(actual{1.0}{2.0}, bestScore, 'Constrained winning score');
                    checkValues(actual{3.0}, bestIndex, 'Constrained winning index');
                    assert(actual{1.0}{3.0} == nCandidates, 'Incorrect reduced candidate count.');
                    assert(isequal(actual{1.0}{1.0}, fliplr(actual{1.0}{1.0})) == true, 'The winner must satisfy the requested symmetry.');
                    checkWinnerTensors(actual{1.0}{6.0}, tensors);
                    checkSamples(actual{2.0}, expected, min(limit, nCandidates));
                    if limit >= nCandidates
                        checkValues(actual{2.0}, expected, 'Complete constrained score order');
                    end
                    if options{method}{1.0} == 3.0
                        assert(actual{4.0} == options{method}{2.0}, 'The configured chunk size must be preserved.');
                        assert(actual{5.0} == ceil(nCandidates/actual{4.0}), 'Chunk count must use the reduced search.');
                    else
                        assert(isempty(actual{4.0}) == true && isempty(actual{5.0}) == true, 'Internal partitions must not change solver metadata.');
                    end
                end
            end
        end
    end
end

%% RETAIN THE FIRST CONSTRAINED CANDIDATE WHEN SCORES TIE
angles = [90.0, 17.0, 0.0];
options = getOptions(9.0);
for objective = 1.0:2.0
    inputs = getInputs(4.0, true, 3.0, angles, objective, 2.0);
    inputs{end - 3.0} = @zeroCriterion;
    for method = 1.0:length(options)
        for limit = [Inf, 0.0, 4.0]
            actual = runOptimisation(inputs, options{method}, limit, true);
            checkValues(actual{1.0}{1.0}, repmat(angles(1.0), 1.0, 4.0), 'First tied mirrored sequence');
            assert(actual{1.0}{2.0} == 0.0 && actual{3.0} == 1.0, 'Ties must select the first constrained candidate.');
        end
    end
end

%% REJECT AN ODD PHYSICAL PLY COUNT FOR AN EXPLICIT MIRRORED SEARCH
inputs = getInputs(3.0, false, 3.0, angles, 1.0, 1.0);
rejected = false;
try
    runOptimisation(inputs, options{2.0}, 0.0, true);
catch exception
    rejected = strcmp(exception.identifier, 'abd:Optimise:SymmetryPlyCount');
end
assert(rejected == true, 'Mirrored optimisation requires an even physical ply count.');

fprintf('[PASS] Symmetric optimisation: filtered exhaustive scores, exact winners, physical ply data, retention and stable ties\n');
end

%% ENUMERATE THE FULL SEARCH BEFORE FILTERING PALINDROMIC INDEX ROWS
function [sequences, selected] = getSymmetricCandidates(nPlies, angles)
nAngles = numel(angles);
fullIndexes = zeros(nAngles^nPlies, nPlies);
for candidate = 1.0:size(fullIndexes, 1.0)
    value = candidate - 1.0;
    for ply = nPlies:-1.0:1.0
        fullIndexes(candidate, ply) = mod(value, nAngles) + 1.0;
        value = floor(value/nAngles);
    end
end
selected = find(all(fullIndexes == fliplr(fullIndexes), 2.0));
sequences = reshape(angles(fullIndexes(selected, :)), length(selected), nPlies);
assert(length(selected) == nAngles^(0.5*nPlies), 'Incorrect exhaustive symmetry reference.');
end

%% EXERCISE ALL SOLVERS AND SMALL, PARTIAL AND OVERSIZED CHUNKS
function [options] = getOptions(nCandidates)
options = {{1.0, 'DEFAULT', 'DEFAULT'}, {2.0, 'DEFAULT', 'DEFAULT'},...
    {3.0, 1.0, 'DEFAULT'}, {3.0, 2.0, 'DEFAULT'}, {3.0, nCandidates + 5.0, 'DEFAULT'}};
end

%% BUILD PHYSICAL PLY PROPERTIES WITHOUT ASSUMING MATERIAL SYMMETRY
function [inputs, data] = getInputs(nPlies, mirrored, nPoints, angles, objective, criterion)
material = [linspace(100000.0, 70000.0, nPlies); linspace(10000.0, 18000.0, nPlies);...
    linspace(5000.0, 8000.0, nPlies); linspace(0.3, 0.25, nPlies);...
    linspace(1e-5, 2e-5, nPlies); linspace(2e-5, 4e-5, nPlies);...
    linspace(0.002, 0.003, nPlies); linspace(0.003, 0.004, nPlies)];
thickness = linspace(0.2, 0.5, nPlies);
if mirrored == true
    material = [material(:, 1.0:0.5*nPlies), fliplr(material(:, 1.0:0.5*nPlies))];
    thickness = [thickness(1.0:0.5*nPlies), fliplr(thickness(1.0:0.5*nPlies))];
end
[z, ~, error] = abd.internal_getThickness(nPlies, thickness, 1e-6);
assert(error == false, 'Invalid validation thickness.');
[error, zPoints, ~, totalPoints, A11, A22, B11, B22] = abd.internal_getSectionPoints(nPoints, '', nPlies, zeros(1.0, nPlies), z,...
    material(5.0, :), material(6.0, :), material(7.0, :), material(8.0, :), 1e-6, []);
assert(error == false, 'Invalid validation section points.');
[Q11, Q22, Q66, Q12] = abd.internal_getReducedQ(nPlies, material(1.0, :), material(2.0, :), material(4.0, :), material(3.0, :));
data = struct('MATERIAL', material, 'N_POINTS', nPoints, 'Z', z, 'Z_POINTS', zPoints, 'Q11', Q11, 'Q22', Q22, 'Q66', Q66, 'Q12', Q12,...
    'A11', A11, 'A22', A22, 'B11', B11, 'B22', B22);
mechanical = repelem(material, 1.0, nPoints);
strength = repmat([1000.0; -500.0; 50.0; -150.0; 100.0; -0.5; 0.0], 1.0, totalPoints);
properties = reshape(num2cell(strength, 2.0), 1.0, []);
unused = cell(1.0, 23.0);
failureCriterion = 'MSTRS';
callback = [];
if criterion == 2.0
    failureCriterion = 'UCRT';
    callback = @(INFO, UCRT, MATERIAL_MECH, MATERIAL_FAIL, TENSORS) checkUserInputs(INFO, UCRT, MATERIAL_MECH, MATERIAL_FAIL, TENSORS, data);
end
inputs = [{{true, failureCriterion, 1.0, objective, angles}, nPlies, totalPoints, nPoints, z, zPoints, Q11, Q22, Q66, Q12,...
    A11, A22, B11, B22, 1e-6}, properties, unused, {80.0, 0.15, -150.0, 30.0, 20.0, 20.0, -10.0, 5.0,...
    mechanical(1.0, :), mechanical(2.0, :), mechanical(4.0, :), mechanical(3.0, :), false, [], [], [], nPoints,...
    callback, {2.0, 'DEFAULT', 'DEFAULT'}, false, 15.0}];
end

%% RUN THE LEGACY OR EXPLICIT SYMMETRY INTERFACE
function [result] = runOptimisation(inputs, options, scoreLimit, symmetricPly)
inputs{end - 2.0} = options;
if isempty(symmetricPly) == false
    inputs = [inputs, {scoreLimit, symmetricPly}]; %#ok<NASGU>
end
result = cell(1.0, 6.0);
optimisationOutput = evalc('[result{:}] = abd.internal_optimise.main(inputs{:});');
if isempty(result{1.0}) == true || isempty(result{1.0}{5.0}) == false
    fprintf('%s', optimisationOutput);
end
assert(isempty(result{1.0}) == false && isempty(result{1.0}{5.0}) == true, 'Optimisation failed.');
end

%% RECOVER THE FULL TENSOR DATA WITHOUT THE OPTIMISATION CACHE
function [tensors] = getTensors(data, theta)
nPlies = size(data.MATERIAL, 2.0);
thetaPoints = repelem(theta, data.N_POINTS);
[Q11t, Q12t, Q16t, Q22t, Q26t, Q66t] = abd.internal_getTransformedQ(nPlies, theta, data.Q11, data.Q12, data.Q66, data.Q22);
[axx, ayy, axy, bxx, byy, bxy] = abd.internal_getThermoHydro(thetaPoints, data.A11, data.A22, data.B11, data.B22);
[ABD, ~, Qt, NxxT, NyyT, NxyT, MxxT, MyyT, MxyT, NxxM, NyyM, NxyM, MxxM, MyyM, MxyM] =...
    abd.internal_getABD(nPlies, Q11t, Q12t, Q16t, Q22t, Q26t, Q66t, data.Z, 80.0, 0.15, axx, ayy, axy, bxx, byy, bxy, data.N_POINTS);
[tensors.E_MIDSPAN, tensors.E_PLY_XY, tensors.S_PLY_XY, tensors.E_PLY_ALIGNED, tensors.S_PLY_ALIGNED,...
    tensors.E_THERM_XY, tensors.E_HYDRO_XY, tensors.E_THERM_ALIGNED, tensors.E_HYDRO_ALIGNED] =...
    abd.internal_getTensor(ABD, -150.0, NxxT, NxxM, 30.0, NyyT, NyyM, 20.0, NxyT, NxyM, 20.0, MxxT, MxxM,...
    -10.0, MyyT, MyyM, 5.0, MxyT, MxyM, nPlies*data.N_POINTS, data.Z_POINTS, thetaPoints, Qt, 80.0, 0.15, axx, ayy, axy, bxx, byy, bxy, 1e-6);
tensors.SYMMETRIC_ABD = abd.internal_getSymmetry(ABD, 1e-6);
end

%% CHECK THAT USER CRITERIA STILL RECEIVE ALL PHYSICAL PLIES
function [UCRT] = checkUserInputs(INFO, UCRT, MATERIAL_MECH, MATERIAL_FAIL, TENSORS, data) %#ok<INUSD>
assert(INFO.PLY_COUNT == size(data.MATERIAL, 2.0), 'The user criterion must receive the full physical ply count.');
assert(INFO.TOTAL_POINTS == INFO.PLY_COUNT*data.N_POINTS, 'The user criterion must receive every section point.');
names = {'E11', 'E22', 'G12', 'V12'};
for i = 1.0:length(names)
    checkValues(MATERIAL_MECH.(names{i}), repelem(data.MATERIAL(i, :), data.N_POINTS), ['User material: ', names{i}]);
end
points = 1.0:data.N_POINTS:INFO.TOTAL_POINTS;
difference = data.MATERIAL(5.0, :) - data.MATERIAL(6.0, :);
theta = 0.5*atan2d(MATERIAL_MECH.AXY(points)./difference, (MATERIAL_MECH.AXX(points) - MATERIAL_MECH.AYY(points))./difference);
expected = getTensors(data, theta);
names = fieldnames(TENSORS);
for i = 1.0:length(names)
    checkValues(TENSORS.(names{i}), expected.(names{i}), ['User tensor: ', names{i}]);
end
UCRT = max(abs(TENSORS.E_PLY_ALIGNED), [], 1.0)/0.01;
end

function [UCRT] = zeroCriterion(INFO, UCRT, MATERIAL_MECH, MATERIAL_FAIL, TENSORS) %#ok<INUSD>
UCRT = zeros(1.0, INFO.TOTAL_POINTS);
end

%% COMPARE THE RETAINED WINNER TENSORS WITH DIRECT RECOVERY
function [] = checkWinnerTensors(actual, expected)
names = {'STRESS_XY', 'STRESS_PLY', 'STRAIN_XY', 'STRAIN_PLY', 'SYMMETRIC_ABD'};
references = {'S_PLY_XY', 'S_PLY_ALIGNED', 'E_PLY_XY', 'E_PLY_ALIGNED', 'SYMMETRIC_ABD'};
for i = 1.0:length(names)
    checkValues(actual.(names{i}), expected.(references{i}), ['Winning tensor: ', names{i}]);
end
end

%% CHECK RETAINED SCORES AGAINST THE FILTERED CANDIDATE ORDER
function [] = checkSamples(actual, complete, expectedCount)
assert(numel(actual) == expectedCount, 'Incorrect retained score count.');
if isempty(actual) == true
    return
end
assert(isrow(actual) == true, 'Retained scores must form a row vector.');
nextIndex = 1.0;
for i = 1.0:length(actual)
    tolerance = 1e-10*max(1.0, abs(actual(i)));
    offset = find(abs(complete(nextIndex:end) - actual(i)) <= tolerance, 1.0);
    assert(isempty(offset) == false, 'Samples must follow the constrained candidate order.');
    nextIndex = nextIndex + offset;
end
end

%% COMPARE NUMERIC OUTPUT WITH THE EXPECTED VALUES
function [] = checkValues(actual, expected, description)
tolerance = 1e-10*max(1.0, max(abs(expected(:))));
assert(isequal(size(actual), size(expected)) == true, '%s: output dimensions differ.', description);
assert(all(isfinite(actual(:))) == true, '%s: output contains nonfinite values.', description);
assert(all(abs(actual(:) - expected(:)) <= tolerance) == true, '%s: incorrect output.', description);
end
