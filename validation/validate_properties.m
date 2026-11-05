function [] = validate_properties()
%   Regression checks for section-point material properties (issue #99).
%
%   RUN THIS FUNCTION FROM THE REPOSITORY ROOT:
%   >> addpath('validation'); validate_properties
%
%   Analysis files are written to a temporary directory. The MATLAB path is
%   restored when the checks finish. Existing parallel pool settings apply.
%
%   Layup Analysis Tool 5.1.5 Copyright Louis Vallance 2026
%
%% - DO NOT EDIT BELOW LINE
%__________________________________________________________________________
%%
originalPath = path;
pathCleanup = onCleanup(@() path(originalPath));
outputLocation = tempname;
mkdir(outputLocation)
[material, limits] = getMaterials();
settings = struct('material', {num2cell(material, 1.0)}, 'fail_strain', {num2cell(limits, 1.0)},...
    'fail_stress', {{[1000.0, 500.0, 50.0, 150.0, 100.0, -0.5, 0.0]}},...
    'stacking_sequence', [0.0, 35.0, 90.0], 'ply_thickness', [0.3, 0.5, 0.4],...
    'load_mech', [100.0, 30.0, 20.0; 5.0, -2.0, 1.0], 'output_location', {{outputLocation, false}});

%% CHECK HETEROGENEOUS PLIES AND DIFFERENT SECTION-POINT COUNTS
for nPoints = [3.0, 1.0, 2.0, 5.0]
    settings.section_points = nPoints;
    settings.fail_strain = num2cell(limits, 1.0);
    settings.output_strength = {true, 'RESERVE', 'DEFAULT'};
    settings.job_name = sprintf('properties_%.0f_points', nPoints);
    result = runAnalysis(settings);
    expected = getExpected(table2array(result.STRESS.TENSOR_PLY), material, limits, nPoints);
    checkValues(result.CFAILURE.STRAIN{'MSTRN', 1.0:end - 1.0}, expected, 'Heterogeneous maximum-strain output');

    % Check all four mechanical properties in the user routine, with and
    % without FAIL_STRAIN. Its output also confirms that it ran successfully.
    settings.output_strength = {@checkUserProperties, 'RESERVE', 'DEFAULT'};
    settings.job_name = sprintf('user_properties_%.0f_points', nPoints);
    result = runAnalysis(settings);
    assert(isfield(result.CFAILURE, 'UCRT') == true, 'User criterion did not return results.');
    checkValues(result.CFAILURE.UCRT{1.0, 1.0:end - 1.0}, expected, 'User criterion with strain failure');
    settings.fail_strain = {[]};
    settings.job_name = sprintf('user_properties_no_strain_%.0f_points', nPoints);
    result = runAnalysis(settings);
    assert(isfield(result.CFAILURE, 'UCRT') == true, 'User criterion did not return results.');
    checkValues(result.CFAILURE.UCRT{1.0, 1.0:end - 1.0}, expected, 'User criterion without strain failure');
end

%% RETAIN THE HOMOGENEOUS CASE
settings.material = {material(:, 1.0)};
settings.fail_strain = {limits(:, 1.0)};
settings.output_strength = {true, 'RESERVE', 'DEFAULT'};
settings.section_points = 3.0;
settings.job_name = 'homogeneous_properties';
result = runAnalysis(settings);
expected = getExpected(table2array(result.STRESS.TENSOR_PLY), repmat(material(:, 1.0), 1.0, 3.0),...
    repmat(limits(:, 1.0), 1.0, 3.0), settings.section_points);
checkValues(result.CFAILURE.STRAIN{'MSTRN', 1.0:end - 1.0}, expected, 'Homogeneous maximum-strain output');

%% CHECK ALL OPTIMISATION MODES AGAINST INDEPENDENT CANDIDATE SCORES
settings.material = num2cell(material, 1.0);
settings.fail_strain = num2cell(limits, 1.0);
sequences = [0.0, 0.0, 0.0; 0.0, 0.0, 90.0; 0.0, 90.0, 0.0; 0.0, 90.0, 90.0;...
    90.0, 0.0, 0.0; 90.0, 0.0, 90.0; 90.0, 90.0, 0.0; 90.0, 90.0, 90.0];
scores = zeros(1.0, size(sequences, 1.0));
for i = 1.0:size(sequences, 1.0)
    settings.stacking_sequence = sequences(i, :);
    settings.job_name = sprintf('property_candidate_%.0f', i);
    result = runAnalysis(settings);
    expected = getExpected(table2array(result.STRESS.TENSOR_PLY), material, limits, settings.section_points);
    checkValues(result.CFAILURE.STRAIN{'MSTRN', 1.0:end - 1.0}, expected, 'Candidate maximum-strain output');
    scores(i) = max(expected);
end

[bestScore, bestIndex] = min(scores);
modes = {'FULL-MATRIX', 'MIXED-RADIX', 'CHUNKS'};
criteria = {'MSTRN', 'UCRT'};
settings.output_strength = {@checkUserProperties, 'RESERVE', 'DEFAULT'};
for i = 1.0:length(modes)
    settings.optimiser_settings = {modes{i}, 2.0, 'DEFAULT'};
    for j = 1.0:length(criteria)
        settings.output_optimised = {criteria{j}, 'RESERVE', 'MINMAX', 90.0};
        settings.job_name = sprintf('property_optimisation_%.0f_%.0f', i, j);
        result = runAnalysis(settings);
        assert(isempty(result.OPT.EXCEPTION) == true, 'Optimisation failed.');
        checkValues(result.OPT.CRITICAL_VALUE, bestScore, 'Optimised failure index');
        checkValues(result.OPT.SEQUENCE, sequences(bestIndex, :), 'Optimised stacking sequence');
    end
end

fprintf('[PASS] Section-point properties: heterogeneous plies, user criteria and all optimisation modes\n');
fprintf('-> Validation output: %s\n', outputLocation);
end

%% GET DISTINCT MECHANICAL PROPERTIES AND STRAIN LIMITS FOR THREE PLIES
function [material, limits] = getMaterials()
material = [100000.0, 70000.0, 150000.0; 10000.0, 18000.0, 12000.0;...
    5000.0, 8000.0, 6000.0; 0.3, 0.25, 0.2; 1e-5, 2e-5, 3e-5;...
    2e-5, 3e-5, 4e-5; 0.002, 0.003, 0.004; 0.003, 0.004, 0.005];
limits = [0.01, 0.012, 0.009; 0.006, 0.008, 0.005;...
    0.02, 0.025, 0.015; 0.01, 0.012, 0.008; 0.04, 0.03, 0.05];
end

%% CHECK THE PROPERTIES RECEIVED BY A USER-DEFINED CRITERION
function [UCRT] = checkUserProperties(INFO, UCRT, MATERIAL_MECH, MATERIAL_FAIL, TENSORS) %#ok<INUSD>
[material, limits] = getMaterials();
checkValues(INFO.PLY_COUNT, 3.0, 'User criterion ply count');
checkValues(INFO.TOTAL_POINTS, 3.0*INFO.SECTION_POINTS, 'User criterion point count');
names = {'E11', 'E22', 'G12', 'V12'};
for i = 1.0:length(names)
    checkValues(MATERIAL_MECH.(names{i}), repelem(material(i, :), INFO.SECTION_POINTS), names{i});
end
if isempty(MATERIAL_FAIL.STRAIN.XET) == false
    signedLimits = limits;
    signedLimits([2.0, 4.0], :) = -signedLimits([2.0, 4.0], :);
    names = {'XET', 'XEC', 'YET', 'YEC', 'SE'};
    for i = 1.0:length(names)
        checkValues(MATERIAL_FAIL.STRAIN.(names{i}), repelem(signedLimits(i, :), INFO.SECTION_POINTS), names{i});
    end
end
UCRT = getExpected(TENSORS.S_PLY_ALIGNED, material, limits, INFO.SECTION_POINTS);
end

%% EVALUATE EACH PLY USING ITS OWN COMPLIANCE AND STRAIN LIMITS
function [expected] = getExpected(stress, material, limits, nPoints)
expected = zeros(1.0, size(stress, 2.0));
for i = 1.0:size(material, 2.0)
    points = (i - 1.0)*nPoints + (1.0:nPoints);
    E11 = material(1.0, i);
    E22 = material(2.0, i);
    G12 = material(3.0, i);
    V12 = material(4.0, i);
    compliance = [1.0/E11, -V12/E11, 0.0; -V12/E11, 1.0/E22, 0.0; 0.0, 0.0, 1.0/G12];
    strain = compliance*stress(:, points);
    expected(points) = max([strain(1.0, :)/limits(1.0, i); -strain(1.0, :)/limits(2.0, i);...
        strain(2.0, :)/limits(3.0, i); -strain(2.0, :)/limits(4.0, i); abs(strain(3.0, :))/limits(5.0, i)], [], 1.0);
end
end

%% RUN ANALYSIS WITHOUT PRINTING THE COMPLETE SUMMARY
function [result] = runAnalysis(settings) %#ok<INUSD>
result = [];
analysisOutput = evalc('result = abd.main(settings);');
if isstruct(result) == false
    fprintf('%s', analysisOutput);
end
assert(isstruct(result) == true, 'Analysis did not return results.');
end

%% COMPARE NUMERIC OUTPUT WITH THE EXPECTED VALUES
function [] = checkValues(actual, expected, description)
tolerance = 1e-10*max(1.0, max(abs(expected(:))));
assert(isequal(size(actual), size(expected)) == true, '%s: output dimensions differ.', description);
assert(all(isfinite(actual(:))) == true, '%s: output contains non-finite values.', description);
assert(all(abs(actual(:) - expected(:)) <= tolerance) == true, '%s: incorrect output.', description);
end
