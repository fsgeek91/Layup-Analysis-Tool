function [] = validate_optimisation_retention_public()
%   Public figure/score-retention checks for issue #111.
%
%   RUN THIS FUNCTION FROM THE REPOSITORY ROOT:
%   >> addpath('validation'); validate_optimisation_retention_public
%
%   Analysis files are written to a temporary directory. The MATLAB path is
%   restored and only figures created by these checks are closed on exit.
%
%   Layup Analysis Tool 5.2.0 Copyright Louis Vallance 2026
%
%% - DO NOT EDIT BELOW LINE
%__________________________________________________________________________
%%
originalPath = path;
pathCleanup = onCleanup(@() path(originalPath));
originalFigures = findall(groot, 'Type', 'figure');
figureCleanup = onCleanup(@() closeNewFigures(originalFigures));
outputLocation = tempname;
mkdir(outputLocation)

settings = struct('material', {{[100000.0, 10000.0, 5000.0, 0.3, 1e-5, 2e-5, 0.002, 0.003]}},...
    'fail_stress', {{[1000.0, 500.0, 50.0, 150.0, 100.0, -0.5, 0.0]}}, 'ply_thickness', 0.2, 'section_points', 3.0,...
    'load_mech', [-150.0, 30.0, 20.0; 5.0, -2.0, 1.0], 'output_strength', {{true, 'VALUE', 'DEFAULT'}},...
    'optimiser_settings', {{'MIXED-RADIX', 'DEFAULT', 'DEFAULT'}}, 'output_location', {{outputLocation, false}});

%% COMPARE NO-FIGURE, COMPLETE-CDF AND SAMPLED-CDF JOBS
% Two plies and three angles produce nine scores; six plies and five angles
% produce 15625 scores, requiring the bounded sample for the CDF.
cases = [2.0, 45.0, 9.0; 6.0, 22.5, 15625.0];
for i = 1.0:size(cases, 1.0)
    settings.stacking_sequence = linspace(0.0, 90.0, cases(i, 1.0));
    settings.output_optimised = {'MSTRS', 'VALUE', 'MINMAX', cases(i, 2.0)};
    settings.output_figure = {[], 'POINTS', 'COMPACT'};
    settings.job_name = sprintf('retention_no_figures_%.0f', i);
    reference = runAnalysis(settings);
    assert(reference.OPT.N_PERMUTATIONS == cases(i, 3.0), 'Incorrect public candidate count.');
    assert(exist(fullfile(outputLocation, settings.job_name, 'figure'), 'dir') == 0.0, 'Disabled figures must not be written.');

    settings.output_figure{1.0} = 'DEFAULT';
    settings.job_name = sprintf('retention_figures_%.0f', i);
    actual = runAnalysis(settings);
    checkWinner(actual.OPT, reference.OPT);
    expectedTitle = 'CDF plot of optimiser criterion for all stacking permutations';
    if cases(i, 3.0) > 10000.0
        expectedTitle = 'CDF plot of optimiser criterion (sample of 10000 of 15625 stacking permutations)';
    end
    checkCdf(fullfile(outputLocation, settings.job_name, 'figure'), expectedTitle);
    closeNewFigures(originalFigures)
end

%% SUPPRESS FIGURES FOR A SINGLE SECTION POINT
settings.stacking_sequence = 17.0;
settings.section_points = 1.0;
settings.output_optimised{4.0} = 45.0;
settings.output_figure{1.0} = [];
settings.job_name = 'retention_one_point_no_figures';
reference = runAnalysis(settings);
settings.output_figure{1.0} = 'DEFAULT';
settings.job_name = 'retention_one_point_figures';
actual = runAnalysis(settings);
checkWinner(actual.OPT, reference.OPT);
assert(exist(fullfile(outputLocation, settings.job_name, 'figure'), 'dir') == 0.0, 'A single section point must not create figures.');

fprintf('[PASS] Public score retention: exact winners, complete/sampled CDF labels and figure suppression\n');
fprintf('-> Validation output: %s\n', outputLocation);
end

%% RUN ANALYSIS WITHOUT PRINTING THE COMPLETE SUMMARY
function [result] = runAnalysis(settings) %#ok<INUSD>
result = [];
analysisOutput = evalc('result = abd.main(settings);');
if isstruct(result) == false
    fprintf('%s', analysisOutput);
end
assert(isstruct(result) == true, 'Analysis did not return results.');
assert(isempty(result.OPT.EXCEPTION) == true, 'Public optimisation failed.');
end

%% COMPARE THE PUBLIC WINNER AND ALL RETAINED TENSORS
function [] = checkWinner(actual, expected)
checkValues(actual.SEQUENCE, expected.SEQUENCE, 'Public winning sequence');
checkValues(actual.CRITICAL_VALUE, expected.CRITICAL_VALUE, 'Public winning score');
checkValues(actual.N_PERMUTATIONS, expected.N_PERMUTATIONS, 'Public candidate count');
names = fieldnames(expected.TENSOR);
for i = 1.0:length(names)
    actualTensor = actual.TENSOR.(names{i});
    expectedTensor = expected.TENSOR.(names{i});
    if istable(actualTensor) == true
        actualTensor = table2array(actualTensor);
        expectedTensor = table2array(expectedTensor);
    end
    checkValues(actualTensor, expectedTensor, ['Public winning tensor: ', names{i}]);
end
end

%% INSPECT THE SAVED CDF FIGURE
function [] = checkCdf(outputLocation, expectedTitle)
files = dir(fullfile(outputLocation, 'CB, *.fig'));
assert(isscalar(files) == true, 'The analysis must save exactly one optimiser CDF.');
f = openfig(fullfile(files(1.0).folder, files(1.0).name), 'invisible');
figureCleanup = onCleanup(@() close(f));
axesHandles = findall(f, 'Type', 'axes');
assert(isscalar(axesHandles) == true, 'Unexpected optimiser CDF layout.');
assert(strcmp(get(get(axesHandles, 'Title'), 'String'), expectedTitle) == true, 'The CDF title must identify complete or sampled scores.');
lines = findall(axesHandles, 'Type', 'line');
assert(isscalar(lines) == true, 'The optimiser CDF must contain one distribution.');
probability = get(lines, 'YData');
assert(min(probability) == 0.0 && max(probability) == 1.0, 'The saved CDF must span probabilities zero to one.');
end

%% CLOSE ONLY FIGURES CREATED BY THIS VALIDATION
function [] = closeNewFigures(originalFigures)
figures = findall(groot, 'Type', 'figure');
for i = 1.0:numel(figures)
    if any(figures(i) == originalFigures) == false
        close(figures(i))
    end
end
end

%% COMPARE NUMERIC OUTPUT WITH THE EXPECTED VALUES
function [] = checkValues(actual, expected, description)
tolerance = 1e-10*max(1.0, max(abs(expected(:))));
assert(isequal(size(actual), size(expected)) == true, '%s: output dimensions differ.', description);
assert(all(isfinite(actual(:))) == true, '%s: nonfinite output.', description);
assert(all(abs(actual(:) - expected(:)) <= tolerance) == true, '%s: incorrect output.', description);
end
