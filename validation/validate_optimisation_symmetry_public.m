function [] = validate_optimisation_symmetry_public()
%   Public symmetry-constrained optimisation checks for issue #112.
%
%   RUN THIS FUNCTION FROM THE REPOSITORY ROOT:
%   >> addpath('validation'); validate_optimisation_symmetry_public
%
%   Analysis files are written to a temporary directory. The MATLAB path is
%   restored and only figures created by these checks are closed on exit.
%
%   Layup Analysis Tool 5.1.5 Copyright Louis Vallance 2026
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

material = [100000.0, 70000.0; 10000.0, 18000.0; 5000.0, 8000.0; 0.3, 0.25;...
    1e-5, 2e-5; 2e-5, 3e-5; 0.002, 0.003; 0.003, 0.004];
strength = [1000.0, 800.0; 500.0, 450.0; 50.0, 60.0; 150.0, 130.0; 100.0, 90.0; -0.5, -0.5; 0.0, 0.0];
settings = struct('material', {num2cell(material, 1.0)}, 'fail_stress', {num2cell(strength, 1.0)},...
    'stacking_sequence', [17.0, -43.0], 'ply_thickness', [0.4, 0.6], 'symmetric_layup', true, 'section_points', 3.0,...
    'load_mech', [-150.0, 30.0, 20.0; 5.0, -2.0, 1.0], 'load_therm', 80.0, 'load_hydro', 0.15,...
    'output_strength', {{true, 'VALUE', 'DEFAULT'}}, 'output_figure', {{[], 'POINTS', 'COMPACT'}},...
    'optimiser_settings', {{'MIXED-RADIX', 'DEFAULT', 'DEFAULT'}}, 'output_location', {{outputLocation, false}});

%% GET INDEPENDENT REFERENCES FROM DIRECT FULL-LAMINATE ANALYSES
angles = [0.0, 45.0, 90.0];
sequences = zeros(9.0, 4.0);
directResults = cell(1.0, 9.0);
scores = zeros(2.0, 9.0);
directSettings = settings;
directSettings.symmetric_layup = false;
directSettings.material = num2cell([material, fliplr(material)], 1.0);
directSettings.fail_stress = num2cell([strength, fliplr(strength)], 1.0);
directSettings.ply_thickness = [0.4, 0.6, 0.6, 0.4];
directSettings.output_optimised = {'', 'VALUE', 'MINMAX', 45.0};
for i = 1.0:9.0
    half = angles([floor((i - 1.0)/3.0) + 1.0, mod(i - 1.0, 3.0) + 1.0]);
    sequences(i, :) = [half, fliplr(half)];
    directSettings.stacking_sequence = sequences(i, :);
    directSettings.job_name = sprintf('symmetry_direct_%.0f', i);
    directResults{i} = runAnalysis(directSettings);
    values = directResults{i}.CFAILURE.STRESS{'MSTRS', 1.0:end - 1.0};
    scores(:, i) = [max(values); mean(values)];
end

%% HONOUR THE EXPLICIT SYMMETRY SETTING IN ALL THREE SOLVERS
modes = {'FULL-MATRIX', 'MIXED-RADIX', 'CHUNKS'};
objectives = {'MINMAX', 'MINMEAN'};
for objective = 1.0:length(objectives)
    [bestScore, bestIndex] = min(scores(objective, :));
    settings.output_optimised = {'MSTRS', 'VALUE', objectives{objective}, 45.0};
    for mode = 1.0:length(modes)
        settings.optimiser_settings = {modes{mode}, 2.0, 'DEFAULT'};
        settings.job_name = sprintf('symmetry_optimised_%.0f_%.0f', objective, mode);
        actual = runAnalysis(settings);
        checkDirectWinner(actual.OPT, directResults{bestIndex}, sequences(bestIndex, :), bestScore);
        assert(actual.OPT.N_PERMUTATIONS == 9.0, 'Symmetric optimisation must enumerate the independent half.');
        checkReport(fullfile(outputLocation, settings.job_name), 9.0);
    end
end

%% ACCEPT EQUIVALENT HALF AND FULL PROPERTY DEFINITIONS
settings.optimiser_settings = {'MIXED-RADIX', 'DEFAULT', 'DEFAULT'};
settings.job_name = 'symmetry_half_properties';
reference = runAnalysis(settings);
fullSettings = settings;
fullSettings.material = directSettings.material;
fullSettings.fail_stress = directSettings.fail_stress;
fullSettings.ply_thickness = directSettings.ply_thickness;
fullSettings.job_name = 'symmetry_full_properties';
actual = runAnalysis(fullSettings);
checkWinner(actual.OPT, reference.OPT);

%% DO NOT INFER A CONSTRAINT FROM AN ALREADY MIRRORED STARTING LAYUP
fullSettings.symmetric_layup = false;
fullSettings.stacking_sequence = [17.0, -43.0, -43.0, 17.0];
fullSettings.job_name = 'symmetry_not_requested';
unconstrained = runAnalysis(fullSettings);
assert(unconstrained.OPT.N_PERMUTATIONS == 81.0, 'An unconstrained job must retain the complete search space.');
checkReport(fullfile(outputLocation, fullSettings.job_name), 81.0);

%% LABEL THE SMALL SYMMETRIC CDF AS A COMPLETE SEARCH
settings.output_figure{1.0} = 'DEFAULT';
settings.job_name = 'symmetry_complete_cdf';
actual = runAnalysis(settings);
checkWinner(actual.OPT, reference.OPT);
checkCdf(fullfile(outputLocation, settings.job_name, 'figure'),...
    'CDF plot of optimiser criterion for all stacking permutations');
closeNewFigures(originalFigures)

%% LABEL SAMPLED CDFS WITH THE CONSTRAINED CANDIDATE COUNT
% Six independent plies give twelve physical plies but only 5^6 candidates.
settings.material = {material(:, 1.0)'};
settings.fail_stress = {strength(:, 1.0)'};
settings.ply_thickness = 0.2;
settings.stacking_sequence = linspace(0.0, 90.0, 6.0);
settings.output_optimised = {'MSTRS', 'VALUE', 'MINMAX', 22.5};
settings.output_figure{1.0} = [];
settings.job_name = 'symmetry_large_no_figures';
reference = runAnalysis(settings);
assert(reference.OPT.N_PERMUTATIONS == 15625.0, 'Incorrect large constrained candidate count.');
assert(isequal(reference.OPT.SEQUENCE, fliplr(reference.OPT.SEQUENCE)) == true, 'The winner must retain mirrored orientations.');
settings.output_figure{1.0} = 'DEFAULT';
settings.job_name = 'symmetry_sampled_cdf';
actual = runAnalysis(settings);
checkWinner(actual.OPT, reference.OPT);
checkReport(fullfile(outputLocation, settings.job_name), 15625.0);
checkCdf(fullfile(outputLocation, settings.job_name, 'figure'),...
    'CDF plot of optimiser criterion (sample of 10000 of 15625 stacking permutations)');

fprintf('[PASS] Public symmetry: direct candidate reference, half/full properties, unconstrained jobs and constrained report/CDF counts\n');
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
if isfield(result, 'OPT') == true
    assert(isempty(result.OPT.EXCEPTION) == true, 'Public optimisation failed.');
end
end

%% COMPARE AN OPTIMISED WINNER WITH DIRECT STRESS/STRAIN RECOVERY
function [] = checkDirectWinner(actual, expected, sequence, score)
checkValues(actual.SEQUENCE, sequence, 'Direct-reference winning sequence');
checkValues(actual.CRITICAL_VALUE, score, 'Direct-reference winning score');
checkValues(table2array(actual.TENSOR.STRESS_XY), table2array(expected.STRESS.TENSOR_XY), 'Direct-reference global stress');
checkValues(table2array(actual.TENSOR.STRESS_PLY), table2array(expected.STRESS.TENSOR_PLY), 'Direct-reference ply stress');
checkValues(table2array(actual.TENSOR.STRAIN_XY), table2array(expected.STRAIN.TENSOR_XY), 'Direct-reference global strain');
checkValues(table2array(actual.TENSOR.STRAIN_PLY), table2array(expected.STRAIN.TENSOR_PLY), 'Direct-reference ply strain');
end

%% COMPARE TWO PUBLIC OPTIMISATION RESULTS
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

%% CHECK THE NUMBER OF CANDIDATES RECORDED IN THE TEXT REPORT
function [] = checkReport(outputLocation, expectedCount)
report = fileread(fullfile(outputLocation, 'summary.log'));
match = regexp(report, 'Checked ([0-9,]+) stacking sequence permutations', 'tokens', 'once');
assert(isempty(match) == false, 'The optimisation candidate count is missing from the report.');
assert(str2double(strrep(match{1.0}, ',', '')) == expectedCount, 'The report must use the constrained candidate count.');
end

%% INSPECT THE SAVED CDF TITLE
function [] = checkCdf(outputLocation, expectedTitle)
files = dir(fullfile(outputLocation, 'CB, *.fig'));
assert(isscalar(files) == true, 'The analysis must save exactly one optimiser CDF.');
f = openfig(fullfile(files(1.0).folder, files(1.0).name), 'invisible');
figureCleanup = onCleanup(@() close(f));
axesHandles = findall(f, 'Type', 'axes');
assert(isscalar(axesHandles) == true, 'Unexpected optimiser CDF layout.');
assert(strcmp(get(get(axesHandles, 'Title'), 'String'), expectedTitle) == true, 'The CDF must use the constrained candidate count.');
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
