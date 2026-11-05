function [] = validate_mstrn()
%   Regression checks for maximum-strain failure (issue #98).
%
%   RUN THIS FUNCTION FROM THE REPOSITORY ROOT:
%   >> addpath('validation'); validate_mstrn
%
%   The public-interface check writes analysis files to a temporary folder.
%   The MATLAB path is restored when the checks finish.
%
%   Layup Analysis Tool 5.1.5 Copyright Louis Vallance 2026
%
%% - DO NOT EDIT BELOW LINE
%__________________________________________________________________________
%%
originalPath = path;
pathCleanup = onCleanup(@() path(originalPath));

% Material rows: E11, E22, V12, G12, XET, XEC, YET, YEC, SE.
% Internal compression limits are negative.
material = [100000.0; 10000.0; 0.3; 5000.0; 0.01; -0.0001; 0.1; -0.1; 0.1];

%% CHECK THE REPORTED POISSON-COUPLING FAILURE
checkValues(evaluate([0.0; 100.0; 0.0], material), 3.0, 'Poisson-induced compressive failure');
checkValues(evaluate([1000.0; 100.0; 0.0], material), 0.97, 'Longitudinal compliance');

%% CHECK BOTH DIRECTIONS OF STRESS/STRAIN SIGN REVERSAL
stress = [-10.0, 10.0, 1000.0, -1000.0; -100.0, 100.0, 10.0, -10.0; zeros(1.0, 4.0)];
materials = repmat(material, 1.0, 4.0);
materials(5.0:6.0, 1.0) = [0.0001; -0.01];
materials(5.0:6.0, 3.0:4.0) = repmat([1.0; -1.0], 1.0, 2.0);
materials(8.0, 3.0) = -0.0005;
materials(7.0, 4.0) = 0.0005;
checkValues(evaluate(stress, materials), [2.0, 2.0, 4.0, 4.0], 'Strain-sign allowable selection');

%% CHECK EACH STRAIN LIMIT, HALF LOADS AND ZERO STRESS
material(5.0:9.0) = [0.01; -0.005; 0.02; -0.01; 0.04];
strain = [material(5.0), material(6.0), zeros(1.0, 5.0);...
    0.0, 0.0, material(7.0), material(8.0), 0.0, 0.0, 0.0;...
    zeros(1.0, 4.0), material(9.0), -material(9.0), 0.0];
materials = repmat(material, 1.0, 7.0);
stress = getStress(strain, materials);
checkValues(evaluate(stress, materials), [ones(1.0, 6.0), 0.0], 'Strain strength limits');
checkValues(evaluate(0.5*stress, materials), [0.5*ones(1.0, 6.0), 0.0], 'Half-load strain indices');

%% CHECK HETEROGENEOUS SECTION-POINT PROPERTIES
materials = repmat(material, 1.0, 4.0);
materials(1.0:4.0, :) = [100000.0, 80000.0, 150000.0, 60000.0;...
    10000.0, 20000.0, 12000.0, 15000.0; 0.3, 0.25, 0.28, 0.2; 5000.0, 7000.0, 4500.0, 6000.0];
materials(5.0:9.0, :) = materials(5.0:9.0, :).*repmat([1.0, 0.5, 1.5, 2.0], 5.0, 1.0);
strain = [0.001, -0.002, 0.0004, -0.001; -0.003, 0.001, -0.0008, 0.002; 0.005, -0.003, 0.0, 0.001];
stress = getStress(strain, materials);
expected = max([strain(1.0, :)./materials(5.0, :); strain(1.0, :)./materials(6.0, :);...
    strain(2.0, :)./materials(7.0, :); strain(2.0, :)./materials(8.0, :); abs(strain(3.0, :)./materials(9.0, :))], [], 1.0);
checkValues(evaluate(stress, materials), expected, 'Prescribed strains and heterogeneous properties');
for i = 1.0:size(stress, 2.0)
    checkValues(evaluate(stress(:, i), materials(:, i)), expected(i), 'Scalar/vector agreement');
end

%% CHECK THE PUBLIC SETTINGS INTERFACE
% One material isolates this correction from the separate property-spread
% issue (#99). Positive input limits exercise the normal sign conversion.
outputLocation = tempname;
mkdir(outputLocation)
settings = struct('job_name', 'maximum_strain_poisson',...
    'material', {{[100000.0, 10000.0, 5000.0, 0.3, 1e-5, 2e-5, 0.002, 0.003]}},...
    'fail_strain', {{[0.01, 0.0001, 0.1, 0.1, 0.1]}}, 'stacking_sequence', 0.0, 'ply_thickness', 1.0,...
    'section_points', 2.0, 'load_mech', [0.0, 100.0, 0.0; 0.0, 0.0, 0.0],...
    'output_strength', {{true, 'RESERVE', 'DEFAULT'}}, 'output_location', {{outputLocation, false}});
result = runAnalysis(settings);
checkValues(result.CFAILURE.STRAIN{'MSTRN', 1.0:end - 1.0}, [3.0, 3.0], 'Public maximum-strain output');

fprintf('[PASS] Maximum strain: compliance, strain signs, limits, section points and public output\n');
fprintf('-> Validation output: %s\n', outputLocation);
end

%% EVALUATE THE MAXIMUM-STRAIN CRITERION
function [result] = evaluate(stress, material)
result = abd.internal_strength.getMstrn(size(stress, 2.0), stress, material(1.0, :), material(2.0, :), material(3.0, :),...
    material(4.0, :), material(5.0, :), material(6.0, :), material(7.0, :), material(8.0, :), material(9.0, :));
end

%% CONVERT PRESCRIBED STRAINS TO STRESSES USING PLY STIFFNESS
function [stress] = getStress(strain, material)
nPoints = size(strain, 2.0);
[~, ~, ~, ~, Q] = abd.internal_getReducedQ(nPoints, material(1.0, :), material(2.0, :), material(3.0, :), material(4.0, :));
stress = zeros(3.0, nPoints);
for i = 1.0:nPoints
    stress(:, i) = Q(:, :, i)*strain(:, i);
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
