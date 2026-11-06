function [] = validate_abd_inverse()
%   Regression checks for optional ABD inversion (issue #108).
%
%   RUN THIS FUNCTION FROM THE REPOSITORY ROOT:
%   >> addpath('validation'); validate_abd_inverse
%
%   No analysis files are written.
%
%   Layup Analysis Tool 5.1.5 Copyright Louis Vallance 2026
%
%% - DO NOT EDIT BELOW LINE
%__________________________________________________________________________
%%
%% PRESERVE ABD ASSEMBLY WITH AND WITHOUT THE INVERSE
loadCases = [0.0, 0.0; 80.0, 0.0; 0.0, 0.15; 80.0, 0.15; -30.0, -0.05];
for nPoints = [1.0, 3.0]
    for i = 1.0:size(loadCases, 1.0)
        inputs = getInputs(nPoints, loadCases(i, :));
        original = cell(1.0, 16.0);
        enabled = original;
        disabled = original;
        [original{:}] = abd.internal_getABD(inputs{:});
        [enabled{:}] = abd.internal_getABD(inputs{:}, true);
        [disabled{:}] = abd.internal_getABD(inputs{:}, false);

        assert(original{16.0} == false, 'ABD assembly failed.');
        assert(isequal(original, enabled) == true, 'Explicit inversion differs from the default.');
        assert(isequal(size(original{2.0}), [6.0, 6.0]) == true, 'The default inverse is missing.');
        assert(norm(original{1.0}*original{2.0} - eye(6.0), inf) <= 1e-10, 'The default inverse is incorrect.');
        assert(isempty(disabled{2.0}) == true, 'Inversion was not disabled.');
        remaining = [1.0, 3.0:16.0];
        assert(isequal(original(remaining), disabled(remaining)) == true,...
            'Skipping the inverse changed stiffness, section-point matrices, resultants or the error flag.');
    end
end

%% PRESERVE SINGULAR MATRIX DETECTION WHEN INVERSION IS DISABLED
inputs = getInputs(3.0, [80.0, 0.15]);
for i = 2.0:7.0
    inputs{i} = zeros(1.0, 2.0);
end
original = cell(1.0, 16.0);
enabled = original;
disabled = original;
originalMessage = evalc('[original{:}] = abd.internal_getABD(inputs{:});');
enabledMessage = evalc('[enabled{:}] = abd.internal_getABD(inputs{:}, true);');
disabledMessage = evalc('[disabled{:}] = abd.internal_getABD(inputs{:}, false);');
assert(original{16.0} == true, 'A singular ABD matrix was not detected.');
assert(isempty(original{2.0}) == true, 'A singular ABD matrix must not return an inverse.');
assert(contains(originalMessage, '[ERROR] ABD matrix is singular to working precision') == true,...
    'The singular-matrix diagnostic is missing.');
assert(isequal(original, enabled, disabled) == true, 'Inversion selection changed singular-matrix outputs.');
assert(isequal(originalMessage, enabledMessage, disabledMessage) == true, 'Inversion selection changed the singular-matrix diagnostic.');

fprintf('[PASS] ABD inverse: default compatibility, optional inversion, environmental resultants and singular-matrix handling\n');
end

%% GET A HETEROGENEOUS OFF-AXIS LAMINATE WITH ENVIRONMENTAL LOADING
function [inputs] = getInputs(nPoints, loadCase)
nPlies = 2.0;
theta = [17.0, -43.0];
z = [-0.5, -0.1, 0.5];
[Q11, Q22, Q66, Q12] = abd.internal_getReducedQ(nPlies, [100000.0, 70000.0], [10000.0, 18000.0], [0.3, 0.25], [5000.0, 8000.0]);
[Q11t, Q12t, Q16t, Q22t, Q26t, Q66t] = abd.internal_getTransformedQ(nPlies, theta, Q11, Q12, Q66, Q22);
expansion = repelem([1e-5, 2e-5; 2e-5, 3e-5; 0.002, 0.003; 0.003, 0.004], 1.0, nPoints);
[axx, ayy, axy, bxx, byy, bxy] = abd.internal_getThermoHydro(repelem(theta, nPoints),...
    expansion(1.0, :), expansion(2.0, :), expansion(3.0, :), expansion(4.0, :));
inputs = {nPlies, Q11t, Q12t, Q16t, Q22t, Q26t, Q66t, z, loadCase(1.0), loadCase(2.0), axx, ayy, axy, bxx, byy, bxy, nPoints};
end
