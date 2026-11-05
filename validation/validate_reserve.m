function [] = validate_reserve()
%   Regression checks for forward-load failure indices (issue #96).
%
%   RUN THIS FUNCTION FROM THE REPOSITORY ROOT:
%   >> addpath('validation'); validate_reserve
%
%   No analysis files are written.
%
%   Layup Analysis Tool 5.1.5 Copyright Louis Vallance 2026
%
%% - DO NOT EDIT BELOW LINE
%__________________________________________________________________________
%%
% Material rows: XT, XC, YT, YC, S, C12, B12. Compression is negative
% because these checks call the internal strength functions directly.
material = [1000.0; -500.0; 50.0; -150.0; 100.0; -0.5; 0.0];

%% CHECK THE REPORTED TSAI-WU TENSION/COMPRESSION LIMITS
stress = [1000.0, -500.0, 500.0, -250.0, 0.0, 0.0, 0.0;...
    0.0, 0.0, 0.0, 0.0, 50.0, -150.0, 0.0; zeros(1.0, 7.0)];
materials = repmat(material, 1.0, 7.0);
checkValues(evaluate(1.0, 1.0, stress, materials), [1.0, 1.0, 0.5, 0.5, 1.0, 1.0, 0.0], 'Tsai-Wu strength limits');

%% CHECK HOFFMAN WITHOUT DEPENDING ON THE SEPARATE EQUATION ISSUE (#97)
% Equal biaxial stresses make the disputed squared stress differences
% identical. With these strengths, failure occurs at 500 +/- 500*sqrt(3).
hoffmanMaterial = material;
hoffmanMaterial(3.0:4.0) = material(1.0:2.0);
limits = [500.0 + 500.0*sqrt(3.0), 500.0 - 500.0*sqrt(3.0)];
biaxial = [limits, 0.5*limits, 0.0];
stress = [biaxial; biaxial; zeros(1.0, 5.0)];
checkValues(evaluate(2.0, 1.0, stress, repmat(hoffmanMaterial, 1.0, 5.0)), [1.0, 1.0, 0.5, 0.5, 0.0], 'Hoffman forward roots');

%% CHECK ZERO STRESS, SHEAR AND EQUAL TENSILE/COMPRESSIVE STRENGTHS
stress = [zeros(2.0, 5.0); 0.0, 100.0, -100.0, 50.0, -50.0];
for criterion = 1.0:2.0
    checkValues(evaluate(criterion, 1.0, stress, repmat(material, 1.0, 5.0)), [0.0, 1.0, 1.0, 0.5, 0.5], 'Zero stress and shear');
end
symmetricMaterial = material;
symmetricMaterial(2.0) = -symmetricMaterial(1.0);
stress = [300.0, -300.0; zeros(2.0, 2.0)];
checkValues(evaluate(1.0, 1.0, stress, repmat(symmetricMaterial, 1.0, 2.0)), [0.3, 0.3], 'Equal tension/compression strengths');

%% CHECK VECTOR INPUTS AND PROPORTIONAL LOADING
stress = [250.0, -200.0, 120.0, -300.0; 15.0, -20.0, -30.0, 10.0; 20.0, -15.0, 35.0, -25.0];
materials = repmat(material, 1.0, 4.0);
materials(1.0:5.0, :) = materials(1.0:5.0, :).*repmat([1.0, 1.5, 0.8, 2.0], 5.0, 1.0);
for criterion = 1.0:2.0
    actual = evaluate(criterion, 1.0, stress, materials);
    assert(all(actual > 0.0), 'The selected load cases must have positive failure indices.');
    for i = 1.0:size(stress, 2.0)
        checkValues(actual(i), evaluate(criterion, 1.0, stress(:, i), materials(:, i)), 'Section-point indexing');
    end
    checkValues(evaluate(criterion, 1.0, 2.5*stress, materials), 2.5*actual, 'Proportional scaling');
    failureStress = stress./repmat(actual, 3.0, 1.0);
    checkValues(evaluate(criterion, 2.0, failureStress, materials), ones(1.0, 4.0), 'Forward failure surface');
end

%% CHECK LINEAR, TANGENT AND NONCONVEX FAILURE POLYNOMIALS
% Each pair defines A*lambda^2 + B*lambda - 1 = 0. For A=-2 and
% B=3, both roots are positive (1/2 and 1); the first failure is at 1/2.
A = [0.0, 0.0, 0.0, -2.0, -1.0, -2.0, -2.0];
B = [2.0, -2.0, 0.0, 3.0, 2.0, 1.0, -3.0];
checkValues(abd.internal_strength.getReserve(A, B), [2.0, 0.0, 0.0, 2.0, 1.0, 0.0, 0.0], 'Degenerate and positive-root selection');

% A valid Tsai-Wu coupling of -1 can also produce a linear polynomial.
linearMaterial = [2.0; -1.0; 2.0; -1.0; 1.0; -1.0; 0.0];
stress = [-1.0, 1.0; -1.0, 1.0; 0.0, 0.0];
checkValues(evaluate(1.0, 1.0, stress, repmat(linearMaterial, 1.0, 2.0)), [1.0, 0.0], 'Tsai-Wu linear limit');

%% CHECK CANCELLATION AND INVALID COEFFICIENTS
% Use relative error here: an absolute tolerance would accept zero and
% miss cancellation of the small positive inverse root.
actual = abd.internal_strength.getReserve([1.0, 1e-200], [-1e20, -1.0]);
expected = [1e-20, 1e-200];
assert(all(abs(actual./expected - 1.0) < 1e-12), 'Small positive roots must survive cancellation.');
assert(all(isnan(abd.internal_strength.getReserve([NaN, Inf, 1.0], [1.0, 1.0, NaN]))), 'Invalid coefficients must remain invalid.');

fprintf('[PASS] Tsai-Wu and Hoffman forward roots, zero stress, shear, scaling and degenerate polynomials\n');
end

%% EVALUATE THE SELECTED CRITERION
function [result] = evaluate(criterion, parameter, stress, material)
if criterion == 1.0
    result = abd.internal_strength.getTsaiw(parameter, stress, material(1.0, :), material(2.0, :), material(3.0, :),...
        material(4.0, :), material(5.0, :), material(6.0, :), material(7.0, :));
else
    nPoints = size(stress, 2.0);
    result = abd.internal_strength.getHoffman(parameter, stress, material(1.0, :), material(2.0, :), material(3.0, :),...
        material(4.0, :), material(5.0, :), zeros(1.0, nPoints), -ones(1.0, nPoints));
end
end

%% COMPARE NUMERIC OUTPUT WITH THE EXPECTED VALUES
function [] = checkValues(actual, expected, description)
tolerance = 1e-10*max(1.0, max(abs(expected(:))));
assert(isequal(size(actual), size(expected)) == true, '%s: output dimensions differ.', description);
assert(all(isfinite(actual(:))) == true, '%s: output contains non-finite values.', description);
assert(all(abs(actual(:) - expected(:)) <= tolerance) == true, '%s: incorrect output.', description);
end
