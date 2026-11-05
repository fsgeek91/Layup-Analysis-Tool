function [] = validate_strain_transform()
%   Regression checks for ply-aligned engineering strains (issue #100).
%
%   RUN THIS FUNCTION FROM THE REPOSITORY ROOT:
%   >> addpath('validation'); validate_strain_transform
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

%% CHECK THE REPORTED NORMAL-STRAIN CASE AND PURE ENGINEERING SHEAR
result = evaluate([0.001; zeros(5.0, 1.0)], 45.0, 0.0, zeros(3.0, 1.0), zeros(3.0, 1.0));
checkValues(result.E_PLY_ALIGNED, [0.0005; 0.0005; -0.001], 'Normal strain at 45 degrees');
result = evaluate([0.0; 0.0; 0.001; zeros(3.0, 1.0)], [45.0, -45.0], [0.0, 0.0], zeros(3.0, 2.0), zeros(3.0, 2.0));
checkValues(result.E_PLY_ALIGNED, [0.0005, -0.0005; -0.0005, 0.0005; 0.0, 0.0], 'Engineering shear at +/-45 degrees');

%% CHECK MIXED STRAINS, CURVATURES AND ENVIRONMENTAL STRAINS
theta = [0.0, 90.0, 45.0, -45.0, 23.0, -67.0, 180.0];
z = [-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6];
midspan = [0.001; -0.0004; 0.0008; 0.0002; -0.0003; 0.0005];
thermal = [0.0001; 0.0003; -0.0002]*[1.0, 0.8, 1.2, 0.6, 1.4, 0.9, 1.1];
moisture = [-0.0002; 0.0001; 0.0003]*[0.5, 1.0, 0.8, 1.2, 0.7, 1.1, 0.9];
result = evaluate(midspan, theta, z, thermal, moisture);
globalStrain = midspan(1.0:3.0) + midspan(4.0:6.0)*z;
checkValues(result.E_PLY_XY, globalStrain, 'Global strains and curvatures');
checkValues(result.E_PLY_ALIGNED, rotateTensor(globalStrain, theta, 2.0), 'Mixed engineering-strain rotation');
checkValues(result.E_THERM_ALIGNED, rotateTensor(thermal, theta, 2.0), 'Thermal-strain rotation');
checkValues(result.E_HYDRO_ALIGNED, rotateTensor(moisture, theta, 2.0), 'Moisture-strain rotation');
expectedStress = [100000.0, 2000.0, 3000.0; 2000.0, 10000.0, 500.0; 3000.0, 500.0, 5000.0]*(globalStrain - thermal - moisture);
checkValues(result.S_PLY_XY, expectedStress, 'Global stress');
checkValues(result.S_PLY_ALIGNED, rotateTensor(expectedStress, theta, 1.0), 'Stress rotation');
checkValues(sum(result.S_PLY_ALIGNED.*result.E_PLY_ALIGNED, 1.0), sum(expectedStress.*globalStrain, 1.0), 'Stress-strain work invariance');

%% CHECK FREE EXPANSION THROUGH THE PUBLIC SETTINGS INTERFACE
material = [100000.0, 10000.0, 5000.0, 0.3, 1e-5, 2e-5, 0.002, 0.003];
settings = struct('material', {{material}}, 'ply_thickness', 1.0, 'section_points', 3.0,...
    'output_location', {{outputLocation, false}});
cases = [45.0, 80.0, 0.0; -45.0, 0.0, 0.15; 30.0, 80.0, 0.15; 0.0, 80.0, 0.15; 90.0, 80.0, 0.15];
for i = 1.0:size(cases, 1.0)
    settings.stacking_sequence = cases(i, 1.0);
    settings.load_therm = cases(i, 2.0);
    settings.load_hydro = cases(i, 3.0);
    settings.job_name = sprintf('aligned_free_expansion_%.0f', i);
    result = runAnalysis(settings);
    thermal = repmat(settings.load_therm*[material(5.0:6.0), 0.0]', 1.0, 3.0);
    moisture = repmat(settings.load_hydro*[material(7.0:8.0), 0.0]', 1.0, 3.0);
    checkValues(table2array(result.STRAIN.THERM_PLY), thermal, 'Ply thermal expansion');
    checkValues(table2array(result.STRAIN.HYDRO_PLY), moisture, 'Ply moisture expansion');
    checkValues(table2array(result.STRAIN.TENSOR_PLY), thermal + moisture, 'Ply free expansion');
    checkValues(table2array(result.STRESS.TENSOR_PLY), zeros(3.0, 3.0), 'Stress-free expansion');
end

%% CHECK USER-CRITERION INPUTS WITH MECHANICAL AND ENVIRONMENTAL LOADS
settings.material = {material, [70000.0, 18000.0, 8000.0, 0.25, 2e-5, 3e-5, 0.003, 0.004]};
settings.ply_thickness = [0.4, 0.6];
settings.stacking_sequence = [25.0, -40.0];
settings.load_mech = [100.0, 30.0, 20.0; 5.0, -2.0, 1.0];
settings.fail_strain = {[0.01, 0.006, 0.02, 0.01, 0.04]};
settings.output_strength = {@checkUserStrains, 'RESERVE', 'DEFAULT'};
settings.job_name = 'aligned_user_strains';
result = runAnalysis(settings);
checkPublicStrains(result, settings.stacking_sequence, settings.section_points);

%% CHECK ALL OPTIMISATION MODES WITH A STRAIN-BASED USER CRITERION
sequences = [0.0, 0.0; 0.0, 45.0; 0.0, 90.0; 45.0, 0.0; 45.0, 45.0; 45.0, 90.0; 90.0, 0.0; 90.0, 45.0; 90.0, 90.0];
scores = zeros(1.0, size(sequences, 1.0));
for i = 1.0:size(sequences, 1.0)
    settings.stacking_sequence = sequences(i, :);
    settings.job_name = sprintf('aligned_candidate_%.0f', i);
    result = runAnalysis(settings);
    scores(i) = checkPublicStrains(result, settings.stacking_sequence, settings.section_points);
end
[bestScore, bestIndex] = min(scores);
settings.output_optimised = {'UCRT', 'RESERVE', 'MINMAX', 45.0};
modes = {'FULL-MATRIX', 'MIXED-RADIX', 'CHUNKS'};
for i = 1.0:length(modes)
    settings.optimiser_settings = {modes{i}, 2.0, 'DEFAULT'};
    settings.job_name = sprintf('aligned_optimisation_%.0f', i);
    result = runAnalysis(settings);
    assert(isempty(result.OPT.EXCEPTION) == true, 'Optimisation failed.');
    checkValues(result.OPT.CRITICAL_VALUE, bestScore, 'Optimised strain criterion');
    checkValues(result.OPT.SEQUENCE, sequences(bestIndex, :), 'Optimised stacking sequence');
    expected = rotateTensor(table2array(result.OPT.TENSOR.STRAIN_XY), repelem(result.OPT.SEQUENCE, settings.section_points), 2.0);
    checkValues(table2array(result.OPT.TENSOR.STRAIN_PLY), expected, 'Optimised ply strains');
end

fprintf('[PASS] Engineering-strain transformations: mechanical, thermal, moisture, user criteria and optimisation\n');
fprintf('-> Validation output: %s\n', outputLocation);
end

%% PRESCRIBE GLOBAL STRAINS AND CURVATURES DIRECTLY
function [result] = evaluate(midspan, theta, z, thermal, moisture)
nPoints = numel(theta);
Qt = repmat([100000.0, 2000.0, 3000.0; 2000.0, 10000.0, 500.0; 3000.0, 500.0, 5000.0], 1.0, 1.0, nPoints);
[result.E_MIDSPAN, result.E_PLY_XY, result.S_PLY_XY, result.E_PLY_ALIGNED, result.S_PLY_ALIGNED,...
    result.E_THERM_XY, result.E_HYDRO_XY, result.E_THERM_ALIGNED, result.E_HYDRO_ALIGNED] =...
    abd.internal_getTensor(eye(6.0), midspan(1.0), 0.0, 0.0, midspan(2.0), 0.0, 0.0, midspan(3.0), 0.0, 0.0,...
    midspan(4.0), 0.0, 0.0, midspan(5.0), 0.0, 0.0, midspan(6.0), 0.0, 0.0, nPoints, z, theta, Qt, 1.0, 1.0,...
    thermal(1.0, :), thermal(2.0, :), thermal(3.0, :), moisture(1.0, :), moisture(2.0, :), moisture(3.0, :), 1e-6);
end

%% ROTATE SYMMETRIC SECOND-ORDER TENSORS INDEPENDENTLY
function [aligned] = rotateTensor(values, theta, shearFactor)
aligned = zeros(size(values));
for i = 1.0:numel(theta)
    rotation = [cosd(theta(i)), sind(theta(i)); -sind(theta(i)), cosd(theta(i))];
    tensor = [values(1.0, i), values(3.0, i)/shearFactor; values(3.0, i)/shearFactor, values(2.0, i)];
    tensor = rotation*tensor*rotation';
    aligned(:, i) = [tensor(1.0, 1.0); tensor(2.0, 2.0); shearFactor*tensor(1.0, 2.0)];
end
end

%% CHECK PUBLIC STRAINS AND THE USER-CRITERION OUTPUT
function [score] = checkPublicStrains(result, theta, nPoints)
expected = rotateTensor(table2array(result.STRAIN.TENSOR_XY), repelem(theta, nPoints), 2.0);
checkValues(table2array(result.STRAIN.TENSOR_PLY), expected, 'Public ply strain output');
assert(isfield(result.CFAILURE, 'UCRT') == true, 'User criterion did not return results.');
criterion = max(abs(expected), [], 1.0)/0.01;
checkValues(result.CFAILURE.UCRT{1.0, 1.0:end - 1.0}, criterion, 'User strain criterion');
score = max(criterion);
end

%% VERIFY LOCAL CONSTITUTIVE CONSISTENCY INSIDE THE USER CRITERION
function [UCRT] = checkUserStrains(INFO, UCRT, MATERIAL_MECH, MATERIAL_FAIL, TENSORS) %#ok<INUSD>
stress = TENSORS.S_PLY_ALIGNED;
elastic = [(stress(1.0, :) - MATERIAL_MECH.V12.*stress(2.0, :))./MATERIAL_MECH.E11;...
    stress(2.0, :)./MATERIAL_MECH.E22 - MATERIAL_MECH.V12.*stress(1.0, :)./MATERIAL_MECH.E11;...
    stress(3.0, :)./MATERIAL_MECH.G12];
checkValues(TENSORS.E_PLY_ALIGNED - TENSORS.E_THERM_ALIGNED - TENSORS.E_HYDRO_ALIGNED, elastic, 'Local elastic strain');
checkValues(TENSORS.E_THERM_ALIGNED, repelem(80.0*[1e-5, 2e-5; 2e-5, 3e-5; 0.0, 0.0], 1.0, INFO.SECTION_POINTS), 'User thermal strains');
checkValues(TENSORS.E_HYDRO_ALIGNED, repelem(0.15*[0.002, 0.003; 0.003, 0.004; 0.0, 0.0], 1.0, INFO.SECTION_POINTS), 'User moisture strains');
UCRT = max(abs(TENSORS.E_PLY_ALIGNED), [], 1.0)/0.01;
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
