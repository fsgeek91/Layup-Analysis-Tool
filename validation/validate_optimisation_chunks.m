function [] = validate_optimisation_chunks()
%   Regression checks for chunk sizing and local score buffers (issue #110).
%
%   RUN THIS FUNCTION FROM THE REPOSITORY ROOT:
%   >> addpath('validation'); validate_optimisation_chunks
%
%   Existing parallel pool settings apply. No analysis files are written.
%
%   Layup Analysis Tool 5.1.5 Copyright Louis Vallance 2026
%
%% - DO NOT EDIT BELOW LINE
%__________________________________________________________________________
%%
%% CHECK POOL DETECTION WITHOUT CREATING A POOL
pctAvailable = (exist('gcp', 'file') ~= 0.0) && (license('test', 'Distrib_Computing_Toolbox') == true);
flags = false;
if pctAvailable == true
    flags = [false, true];
end
for pctAvail = flags
    poolBefore = [];
    if pctAvailable == true
        poolBefore = gcp('nocreate');
    end
    mode = 'N/A';
    pool = [];
    configurationOutput = evalc('[mode, pool] = abd.internal_optimise.checkParforConfig(pctAvail, mode);'); %#ok<NASGU>
    if (pctAvail == false) || (isempty(poolBefore) == true)
        assert(strcmpi(mode, 'Serial') == true && isempty(pool) == true, 'Serial pool detection failed.');
    else
        assert(strcmpi(mode, 'Parallel') == true && isequal(pool, poolBefore) == true, 'Connected pool detection failed.');
    end
    if pctAvailable == true
        assert(isequal(gcp('nocreate'), poolBefore) == true, 'Pool detection must not create or replace a pool.');
    end
end

%% COMPARE EVERY CANDIDATE, BOTH OBJECTIVES AND CHUNK BOUNDARIES
% Four heterogeneous plies give 81 candidates. The other cases exercise a
% single candidate and tied scores, which must retain the first candidate.
angleLists = {[0.0, 17.0, 90.0], 23.0, [0.0, 17.0, 90.0]};
loadScales = [1.0, 1.0, 0.0];
for caseIndex = 1.0:length(angleLists)
    nCandidates = numel(angleLists{caseIndex})^4.0;
    chunkSettings = {{3.0, 'DEFAULT', 'DEFAULT'}, {3.0, 'DEFAULT', 3.0},...
        {3.0, 1.0, 'DEFAULT'}, {3.0, 2.0, 'DEFAULT'}, {3.0, 7.0, 'DEFAULT'},...
        {3.0, nCandidates, 'DEFAULT'}, {3.0, nCandidates + 5.0, 'DEFAULT'}};
    for objective = 1.0:2.0
        inputs = getInputs(angleLists{caseIndex}, objective, loadScales(caseIndex));
        reference = runOptimisation(inputs, {2.0, 'DEFAULT', 'DEFAULT'}, pctAvailable);
        full = runOptimisation(inputs, {1.0, 'DEFAULT', 'DEFAULT'}, pctAvailable);
        checkResults(full, reference);
        if loadScales(caseIndex) == 0.0
            assert(all(reference{2.0} == 0.0) == true && reference{3.0} == 1.0, 'Zero-load ties must select the first candidate.');
        end
        for pctAvail = flags
            for setting = 1.0:length(chunkSettings)
                options = chunkSettings{setting};
                pool = [];
                if pctAvail == true
                    pool = gcp('nocreate');
                end
                expectedMode = 'Serial';
                if isempty(pool) == false
                    expectedMode = 'Parallel';
                end
                if setting <= 2.0
                    expectedSize = ceil(nCandidates/40.0);
                    if isempty(pool) == false
                        tuning = 5.0;
                        if setting == 2.0
                            tuning = 3.0;
                        end
                        expectedSize = ceil(nCandidates/(pool.NumWorkers*tuning));
                    end
                else
                    expectedSize = options{2.0};
                end
                actual = runOptimisation(inputs, options, pctAvail);
                checkResults(actual, reference);
                assert(actual{4.0} == expectedSize, 'Incorrect chunk size.');
                assert(actual{5.0} == ceil(nCandidates/expectedSize), 'Incorrect number of chunks.');
                assert(strcmpi(actual{6.0}, expectedMode) == true, 'Incorrect execution mode.');
            end
        end
    end
end

fprintf('[PASS] Chunk handling: pool detection, default/custom sizing, candidate order, partial chunks and tied winners\n');
end

%% GET HETEROGENEOUS OPTIMISATION INPUTS
function [inputs] = getInputs(angles, objective, loadScale)
nPlies = 4.0;
nPoints = 3.0;
material = [100000.0, 90000.0, 80000.0, 70000.0; 10000.0, 12000.0, 15000.0, 18000.0;...
    5000.0, 6000.0, 7000.0, 8000.0; 0.3, 0.28, 0.26, 0.25;...
    1e-5, 2e-5, 3e-5, 4e-5; 2e-5, 3e-5, 4e-5, 5e-5;...
    0.002, 0.003, 0.004, 0.005; 0.003, 0.004, 0.005, 0.006];
[z, ~, error] = abd.internal_getThickness(nPlies, [0.2, 0.3, 0.4, 0.5], 1e-6);
assert(error == false, 'Invalid validation thickness.');
[error, zPoints, ~, totalPoints, A11, A22, B11, B22] = abd.internal_getSectionPoints(nPoints, '', nPlies, zeros(1.0, nPlies), z,...
    material(5.0, :), material(6.0, :), material(7.0, :), material(8.0, :), 1e-6, []);
assert(error == false, 'Invalid validation section points.');
[Q11, Q22, Q66, Q12] = abd.internal_getReducedQ(nPlies, material(1.0, :), material(2.0, :), material(4.0, :), material(3.0, :));
mechanical = repelem(material, 1.0, nPoints);
strength = repmat([1000.0; -500.0; 50.0; -150.0; 100.0; -0.5; 0.0], 1.0, totalPoints);
properties = reshape(num2cell(strength, 2.0), 1.0, []);
unused = cell(1.0, 23.0);
loads = num2cell(loadScale*[80.0, 0.15, -150.0, 30.0, 20.0, 20.0, -10.0, 5.0]);
inputs = [{{true, 'MSTRS', 1.0, objective, angles}, nPlies, totalPoints, nPoints, z, zPoints, Q11, Q22, Q66, Q12,...
    A11, A22, B11, B22, 1e-6}, properties, unused, loads,...
    {mechanical(1.0, :), mechanical(2.0, :), mechanical(4.0, :), mechanical(3.0, :), false, [], [], [], nPoints,...
    [], {3.0, 'DEFAULT', 'DEFAULT'}, false, 15.0}];
end

%% RETAIN SCORES, WINNER AND CHUNK METADATA
function [result] = runOptimisation(inputs, options, pctAvail)
inputs{end - 2.0} = options;
inputs{end - 1.0} = pctAvail; %#ok<NASGU>
result = cell(1.0, 6.0);
optimisationOutput = evalc('[result{:}] = abd.internal_optimise.main(inputs{:});');
if isempty(result{1.0}) == true || isempty(result{1.0}{5.0}) == false
    fprintf('%s', optimisationOutput);
end
assert(isempty(result{1.0}) == false && isempty(result{1.0}{5.0}) == true, 'Optimisation failed.');
end

%% COMPARE ALL SCORES AND THE COMPLETE WINNING RESULT
function [] = checkResults(actual, expected)
checkValues(actual{2.0}, expected{2.0}, 'Complete candidate score row');
checkValues(actual{3.0}, expected{3.0}, 'Winning candidate index');
for i = 1.0:3.0
    checkValues(actual{1.0}{i}, expected{1.0}{i}, 'Winning sequence, score or candidate count');
end
names = fieldnames(expected{1.0}{6.0});
for i = 1.0:length(names)
    checkValues(actual{1.0}{6.0}.(names{i}), expected{1.0}{6.0}.(names{i}), ['Winning tensor: ', names{i}]);
end
end

%% COMPARE NUMERIC OUTPUT WITH THE EXPECTED VALUES
function [] = checkValues(actual, expected, description)
tolerance = 1e-10*max(1.0, max(abs(expected(:))));
assert(isequal(size(actual), size(expected)) == true, '%s: output dimensions differ.', description);
assert(all(isfinite(actual(:))) == true, '%s: output contains non-finite values.', description);
assert(all(abs(actual(:) - expected(:)) <= tolerance) == true, '%s: incorrect output.', description);
end
