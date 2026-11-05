function [] = validate_environmental()
%   Regression checks for thermal and moisture resultants (issue #93).
%
%   RUN THIS FUNCTION FROM THE REPOSITORY ROOT:
%   >> addpath('validation'); validate_environmental
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

%% CHECK ENVIRONMENTAL FORCES AND MOMENTS BY DIRECT INTEGRATION
% Use coupled stiffness matrices and different expansion coefficients in
% each ply so that all six components of each resultant are exercised.
z = [-1.0, -0.2, 1.0];
Q = cat(3.0, [100000.0, 3000.0, 1000.0; 3000.0, 10000.0, 500.0; 1000.0, 500.0, 5000.0],...
    [80000.0, 2000.0, -800.0; 2000.0, 12000.0, -400.0; -800.0, -400.0, 6000.0]);
alpha = [1e-5, 2e-5; 2e-5, 3e-5; -5e-6, 8e-6];
beta = [0.002, 0.003; 0.003, 0.004; -0.001, 0.0005];
loadCases = [100.0, 0.0; 0.0, 0.1; 100.0, 0.1; -50.0, -0.02; 0.0, 0.0];
expectedABD = zeros(6.0, 6.0);
expectedThermal = zeros(6.0, 1.0);
expectedMoisture = expectedThermal;

for i = 1.0:2.0
    % Integrate a constant ply stress and its first moment analytically.
    weights = [z(i + 1.0) - z(i), 0.5*(z(i + 1.0)^2.0 - z(i)^2.0)];
    cubicWeight = (z(i + 1.0)^3.0 - z(i)^3.0)/3.0;
    expectedABD = expectedABD + [weights(1.0)*Q(:, :, i), weights(2.0)*Q(:, :, i);...
        weights(2.0)*Q(:, :, i), cubicWeight*Q(:, :, i)];
    expectedThermal = expectedThermal + reshape((Q(:, :, i)*alpha(:, i))*weights, [], 1.0);
    expectedMoisture = expectedMoisture + reshape((Q(:, :, i)*beta(:, i))*weights, [], 1.0);
end

% Supply the global coefficients directly to isolate resultant assembly
% from the separate moisture-coefficient issue (#94).
alphaPoints = repelem(alpha, 1.0, 2.0);
betaPoints = repelem(beta, 1.0, 2.0);
for i = 1.0:size(loadCases, 1.0)
    deltaT = loadCases(i, 1.0);
    deltaM = loadCases(i, 2.0);
    [ABD, ~, ~, NxxT, NyyT, NxyT, MxxT, MyyT, MxyT, NxxM, NyyM, NxyM, MxxM, MyyM, MxyM, error] =...
        abd.internal_getABD(2.0, reshape(Q(1.0, 1.0, :), 1.0, []), reshape(Q(1.0, 2.0, :), 1.0, []),...
        reshape(Q(1.0, 3.0, :), 1.0, []), reshape(Q(2.0, 2.0, :), 1.0, []), reshape(Q(2.0, 3.0, :), 1.0, []),...
        reshape(Q(3.0, 3.0, :), 1.0, []), z, deltaT, deltaM, alphaPoints(1.0, :), alphaPoints(2.0, :), alphaPoints(3.0, :),...
        betaPoints(1.0, :), betaPoints(2.0, :), betaPoints(3.0, :), 2.0);
    assert(error == false, 'Environmental resultant assembly failed.');
    checkValues(ABD, expectedABD, 'ABD stiffness');
    checkValues([NxxT; NyyT; NxyT; MxxT; MyyT; MxyT], deltaT*expectedThermal, 'Thermal resultants');
    checkValues([NxxM; NyyM; NxyM; MxxM; MyyM; MxyM], deltaM*expectedMoisture, 'Moisture resultants');
end

%% CHECK FREE EXPANSION THROUGH THE PUBLIC SETTINGS INTERFACE
settings = struct('material', {{[100000.0, 10000.0, 5000.0, 0.3, 1e-5, 2e-5, 0.002, 0.003]}},...
    'stacking_sequence', 0.0, 'ply_thickness', 1.0, 'section_points', 2.0, 'output_location', {{outputLocation, false}});
% Off-axis moisture cases have nonzero BXY, avoiding the independent #94
% defect. Check global strains because aligned strains have a separate bug.
freeCases = [0.0, 100.0, 0.0; 45.0, 0.0, 0.1; 45.0, 100.0, 0.1];
for i = 1.0:size(freeCases, 1.0)
    theta = freeCases(i, 1.0);
    settings.job_name = sprintf('free_expansion_%.0f', i);
    settings.stacking_sequence = theta;
    settings.load_therm = freeCases(i, 2.0);
    settings.load_hydro = freeCases(i, 3.0);
    expansion = settings.load_therm*[1e-5; 2e-5] + settings.load_hydro*[0.002; 0.003];
    expectedStrain = [cosd(theta)^2.0, sind(theta)^2.0; sind(theta)^2.0, cosd(theta)^2.0;...
        2.0*sind(theta)*cosd(theta), -2.0*sind(theta)*cosd(theta)]*expansion;
    result = runAnalysis(settings);
    checkValues(table2array(result.E_MIDSPAN), [expectedStrain; zeros(3.0, 1.0)], 'Free expansion');
    checkValues(table2array(result.STRESS.TENSOR_XY), zeros(3.0, 2.0), 'Free-expansion stress');
end

%% CHECK THE MECHANICAL-ONLY CASE
settings.job_name = 'mechanical_only';
settings.stacking_sequence = 0.0;
settings.load_therm = 0.0;
settings.load_hydro = 0.0;
settings.load_mech = [10.0, 5.0, 2.0; 0.0, 0.0, 0.0];
result = runAnalysis(settings);
checkValues(table2array(result.STRESS.TENSOR_XY), repmat([10.0; 5.0; 2.0], 1.0, 2.0), 'Mechanical-only stress');
checkValues(table2array(result.E_MIDSPAN), [(10.0 - 0.3*5.0)/100000.0; 5.0/10000.0 - 0.3*10.0/100000.0;...
    2.0/5000.0; 0.0; 0.0; 0.0], 'Mechanical-only strain');

%% CHECK ALL THREE OPTIMISATION PATHS AGAINST DIRECT ANALYSIS
% Use thermal loading here so that the separate moisture-coefficient defect
% cannot affect the 0/90-degree candidate set. MSTRS avoids other criteria.
settings.material = {[100000.0, 10000.0, 5000.0, 0.3, 1e-5, 2e-5, 0.002, 0.003],...
    [100000.0, 10000.0, 5000.0, 0.3, 2e-5, 4e-5, 0.002, 0.003]};
settings.fail_stress = {[1000.0, 500.0, 50.0, 150.0, 100.0, -0.5, 0.0]};
settings.ply_thickness = 0.5;
settings.load_therm = 100.0;
settings.load_mech = [10.0, 5.0, 2.0; 1.0, -0.5, 0.2];
settings.output_strength = {true, 'RESERVE', 'DEFAULT'};
sequences = [0.0, 0.0; 0.0, 90.0; 90.0, 0.0; 90.0, 90.0];
scores = zeros(1.0, 4.0);
directResults = cell(1.0, 4.0);
for i = 1.0:4.0
    settings.job_name = sprintf('direct_candidate_%.0f', i);
    settings.stacking_sequence = sequences(i, :);
    directResults{i} = runAnalysis(settings);
    criterion = directResults{i}.CFAILURE.STRESS{'MSTRS', 1.0:end - 1.0};
    scores(i) = max(criterion);
end

[bestScore, bestIndex] = min(scores);
modes = {'FULL-MATRIX', 'MIXED-RADIX', 'CHUNKS'};
settings.output_optimised = {'MSTRS', 'RESERVE', 'MINMAX', 90.0};
for i = 1.0:length(modes)
    settings.job_name = sprintf('optimisation_%.0f', i);
    settings.optimiser_settings = {modes{i}, 2.0, 'DEFAULT'};
    result = runAnalysis(settings);
    assert(isempty(result.OPT.EXCEPTION) == true, 'Optimisation failed.');
    checkValues(result.OPT.CRITICAL_VALUE, bestScore, 'Optimised failure index');
    checkValues(result.OPT.SEQUENCE, sequences(bestIndex, :), 'Optimised stacking sequence');
    checkValues(table2array(result.OPT.TENSOR.STRESS_XY), table2array(directResults{bestIndex}.STRESS.TENSOR_XY), 'Optimised stress');
end

fprintf('[PASS] Environmental resultants, free expansion, mechanical loading and all optimisation modes\n');
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
end

%% COMPARE NUMERIC OUTPUT WITH THE EXPECTED VALUES
function [] = checkValues(actual, expected, description)
tolerance = 1e-9*max(1.0, max(abs(expected(:))));
assert(isequal(size(actual), size(expected)) == true, '%s: output dimensions differ.', description);
assert(all(isfinite(actual(:))) == true, '%s: output contains non-finite values.', description);
assert(all(abs(actual(:) - expected(:)) <= tolerance) == true, '%s: incorrect output.', description);
end
