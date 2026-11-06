function [] = validate_larc05_optimisation()
%   Regression checks for LaRC05 optimisation scores and steps (#101, #105).
%
%   RUN THIS FUNCTION FROM THE REPOSITORY ROOT:
%   >> addpath('validation'); validate_larc05_optimisation
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
material = [100000.0, 70000.0; 10000.0, 18000.0; 5000.0, 8000.0; 0.3, 0.25;...
    1e-5, 2e-5; 2e-5, 3e-5; 0.002, 0.003; 0.003, 0.004];
settings = struct('material', {num2cell(material, 1.0)}, 'ply_thickness', [0.4, 0.6], 'section_points', 3.0,...
    'load_therm', 20.0, 'load_hydro', 0.01, 'output_strength', {{true, 'VALUE', 15.0}}, 'output_location', {{outputLocation, false}});
sequences = [0.0, 0.0; 0.0, 45.0; 0.0, 90.0; 45.0, 0.0; 45.0, 45.0; 45.0, 90.0; 90.0, 0.0; 90.0, 45.0; 90.0, 90.0];
modes = {'FULL-MATRIX', 'MIXED-RADIX', 'CHUNKS'};
objectives = {'MINMAX', 'MINMEAN'};

% The first case makes polymer failure govern; the second makes matrix
% failure exceed one. Explicit plane spacing matches direct and optimised
% analysis. DEFAULT step selection is checked separately.
failureMaterials = {[50.0; 50.0; 100000.0; 100000.0; 100000.0; 100000.0; 5000.0; 0.0; 0.0; 53.0; 2.0],...
    [1e6; 1e6; 25.0; 200.0; 40.0; 50.0; 5000.0; 0.0; 0.0; 53.0; 2.0]};
loads = {[-150.0, 30.0, 20.0; 20.0, -10.0, 5.0], [100.0, 150.0, 30.0; 5.0, -2.0, 1.0]};
for caseIndex = 1.0:length(failureMaterials)
    failure = repmat(failureMaterials{caseIndex}, 1.0, 2.0);
    failure(1.0:7.0, 2.0) = 1.2*failure(1.0:7.0, 2.0);
    settings.larc05 = num2cell(failure, 1.0);
    settings.load_mech = loads{caseIndex};
    settings.output_optimised = {'', 'VALUE', 'MINMAX', 45.0};
    scores = zeros(2.0, size(sequences, 1.0));
    directResults = cell(1.0, size(sequences, 1.0));
    polymerGoverns = false;
    clippingChangesScore = false;
    averagingModesChangesScore = false;

    %% GET ALL CANDIDATE SCORES FROM DIRECT ANALYSIS
    for i = 1.0:size(sequences, 1.0)
        settings.stacking_sequence = sequences(i, :);
        settings.job_name = sprintf('larc05_candidate_%.0f_%.0f', caseIndex, i);
        result = runAnalysis(settings);
        directResults{i} = result;
        values = result.CFAILURE.LARC05{:, 1.0:end - 1.0};
        worst = max(values, [], 1.0);
        scores(:, i) = [max(worst); mean(worst)];

        % Check the principal-stress contribution independently using the
        % plane-stress deviatoric invariant, without an eigensystem.
        stress = table2array(result.STRESS.TENSOR_PLY);
        Xt = repelem(failure(1.0, :), settings.section_points);
        Xc = repelem(failure(2.0, :), settings.section_points);
        polymer = (stress(1.0, :).^2.0 - stress(1.0, :).*stress(2.0, :) + stress(2.0, :).^2.0 +...
            3.0*stress(3.0, :).^2.0 - (Xt - Xc).*(stress(1.0, :) + stress(2.0, :)))./(Xt.*Xc);
        checkValues(values(1.0, :), polymer, 'Candidate polymer failure');
        polymerGoverns = polymerGoverns || any(values(1.0, :) > max(values(2.0:end, :), [], 1.0));
        clipped = values;
        clipped(2.0, :) = min(clipped(2.0, :), 1.0);
        clippingChangesScore = clippingChangesScore || any(abs(worst - max(clipped, [], 1.0)) > 1e-6);
        averagingModesChangesScore = averagingModesChangesScore || abs(mean(worst) - mean(values(:))) > 1e-6;
    end
    if caseIndex == 1.0
        assert(polymerGoverns == true, 'The first case must exercise governing polymer failure.');
    else
        assert(clippingChangesScore == true, 'The second case must exercise uncapped matrix failure.');
    end
    assert(averagingModesChangesScore == true, 'The cases must distinguish worst-mode and all-mode averaging.');

    %% CHECK EVERY SCORE IN ALL THREE OPTIMISATION MODES
    for mode = 1.0:length(modes)
        for objective = 1.0:length(objectives)
            [bestScore, bestIndex] = min(scores(objective, :));
            % Supply principal stresses from two different initial layups.
            % Neither may affect any candidate's score or the winner.
            for initialIndex = [1.0, size(sequences, 1.0)]
                stress = table2array(directResults{initialIndex}.STRESS.TENSOR_PLY);
                [best, actual] = runOptimisation(settings, material, failure, mode, objective, stress);
                assert(isempty(best{5.0}) == true, 'Optimisation failed.');
                checkValues(actual, scores(objective, :), 'All candidate scores');
                checkValues(best{2.0}, bestScore, 'Optimised LaRC05 score');
                checkValues(best{1.0}, sequences(bestIndex, :), 'Optimised stacking sequence');
                checkValues(best{6.0}.STRESS_PLY, table2array(directResults{bestIndex}.STRESS.TENSOR_PLY), 'Winning candidate stress');
            end
        end
    end

    %% CHECK THE PUBLIC OPTIMISATION INTERFACE
    settings.stacking_sequence = [25.0, -40.0];
    settings.optimiser_settings = {'MIXED-RADIX', 2.0, 'DEFAULT'};
    for objective = 1.0:length(objectives)
        settings.output_optimised = {'LARC05', 'VALUE', objectives{objective}, 45.0};
        settings.job_name = sprintf('larc05_public_%.0f_%.0f', caseIndex, objective);
        result = runAnalysis(settings);
        [bestScore, bestIndex] = min(scores(objective, :));
        assert(isempty(result.OPT.EXCEPTION) == true, 'Public optimisation failed.');
        checkValues(result.OPT.CRITICAL_VALUE, bestScore, 'Public optimised LaRC05 score');
        checkValues(result.OPT.SEQUENCE, sequences(bestIndex, :), 'Public optimised stacking sequence');
    end
end

%% CHECK DEFAULT AND EXPLICIT STEP ROUTING THROUGH THE PUBLIC INTERFACE
notice = 'The optimiser will use a fixed default step size of 15 degrees for LaRC05';
settings.output_strength = {true, 'RESERVE', 15.0};
settings.output_optimised = {'LARC05', 'VALUE', 'MINMAX', 45.0};
settings.job_name = 'larc05_public_fixed_default';
[fixedResult, analysisOutput] = runAnalysis(settings);
assert(contains(analysisOutput, notice) == false, 'An explicit step must not report the default step notice.');
settings.output_strength{3.0} = 'DEFAULT';
settings.job_name = 'larc05_public_default_optimised';
[defaultResult, analysisOutput] = runAnalysis(settings);
assert(count(analysisOutput, notice) == 1.0, 'DEFAULT optimisation must report its step exactly once.');
checkValues(defaultResult.CFAILURE.LARC05{:, 1.0:end - 1.0}, fixedResult.CFAILURE.LARC05{:, 1.0:end - 1.0}, 'DEFAULT initial-layup fixed step');
checkValues(defaultResult.OPT.CRITICAL_VALUE, fixedResult.OPT.CRITICAL_VALUE, 'DEFAULT public optimisation score');
checkValues(defaultResult.OPT.SEQUENCE, fixedResult.OPT.SEQUENCE, 'DEFAULT public optimisation sequence');

settings.output_optimised{1.0} = '';
settings.job_name = 'larc05_public_default_analysis';
[adaptiveResult, analysisOutput] = runAnalysis(settings);
assert(contains(analysisOutput, notice) == false, 'Ordinary analysis must not report an optimisation step.');
stress = table2array(adaptiveResult.STRESS.TENSOR_PLY);
checkValues(adaptiveResult.CFAILURE.LARC05{:, 1.0:end - 1.0}, evaluate(stress, failure, settings.section_points, 'DEFAULT'), 'Public adaptive analysis');

settings.output_strength{3.0} = 7.0;
settings.job_name = 'larc05_public_custom_analysis';
[customResult, analysisOutput] = runAnalysis(settings);
assert(contains(analysisOutput, notice) == false, 'An explicit step must not report the default step notice.');
checkValues(customResult.CFAILURE.LARC05{:, 1.0:end - 1.0}, evaluate(stress, failure, settings.section_points, 7.0), 'Public custom fixed step');
[customBest, ~] = runOptimisation(settings, material, failure, 2.0, 1.0, stress);
settings.output_optimised{1.0} = 'LARC05';
settings.job_name = 'larc05_public_custom_optimised';
[customOptimised, analysisOutput] = runAnalysis(settings);
assert(contains(analysisOutput, notice) == false, 'Custom optimisation must not report the default step notice.');
checkValues(customOptimised.CFAILURE.LARC05{:, 1.0:end - 1.0}, customResult.CFAILURE.LARC05{:, 1.0:end - 1.0}, 'Custom initial-layup fixed step');
checkValues(customOptimised.OPT.CRITICAL_VALUE, customBest{2.0}, 'Custom public optimisation score');
checkValues(customOptimised.OPT.SEQUENCE, customBest{1.0}, 'Custom public optimisation sequence');

fprintf('[PASS] LaRC05 optimisation: candidate principal stresses, uncapped matrix failure and both objectives in all modes\n');
fprintf('[PASS] LaRC05 public DEFAULT/custom step selection and command-window notices\n');
fprintf('-> Validation output: %s\n', outputLocation);
end

%% RUN THE OPTIMISER AND RETURN ITS COMPLETE CANDIDATE BUFFER
function [best, scores] = runOptimisation(settings, material, failure, mode, objective, initialStress)
nPlies = size(material, 2.0);
nPoints = settings.section_points;
[z, ~, error] = abd.internal_getThickness(nPlies, settings.ply_thickness, 1e-6);
assert(error == false, 'Invalid validation thickness.');
[error, zPoints, ~, totalPoints, A11, A22, B11, B22] = abd.internal_getSectionPoints(nPoints, '', nPlies, zeros(1.0, nPlies), z,...
    material(5.0, :), material(6.0, :), material(7.0, :), material(8.0, :), 1e-6, []);
assert(error == false, 'Invalid validation section points.');
[Q11, Q22, Q66, Q12] = abd.internal_getReducedQ(nPlies, material(1.0, :), material(2.0, :), material(4.0, :), material(3.0, :));
mechanical = repelem(material, 1.0, nPoints);
failure = repelem(failure, 1.0, nPoints);
principal = zeros(3.0, totalPoints);
for i = 1.0:totalPoints
    principal(:, i) = eig([initialStress(1.0, i), initialStress(3.0, i), 0.0;...
        initialStress(3.0, i), initialStress(2.0, i), 0.0; 0.0, 0.0, 0.0]);
end
loadMech = settings.load_mech;
best = [];
scores = [];
% Stress, strain and Hashin properties are unused in this LaRC05 test.
unused = cell(1.0, 19.0);
% Explicit chunk size also covers a partial final chunk (nine candidates).
inputs = [{{true, 'LARC05', 2.0, objective, [0.0, 45.0, 90.0]}, nPlies, totalPoints, nPoints, z, zPoints, Q11, Q22, Q66, Q12,...
    A11, A22, B11, B22, 1e-6}, unused, {failure(1.0, :), failure(2.0, :), failure(3.0, :), failure(4.0, :), failure(5.0, :),...
    failure(6.0, :), failure(7.0, :), failure(8.0, :), failure(9.0, :), failure(10.0, :), failure(11.0, :), settings.load_therm, settings.load_hydro,...
    loadMech(1.0, 1.0), loadMech(1.0, 2.0), loadMech(1.0, 3.0), loadMech(2.0, 1.0), loadMech(2.0, 2.0), loadMech(2.0, 3.0),...
    mechanical(1.0, :), mechanical(2.0, :), mechanical(4.0, :), mechanical(3.0, :), false, principal(1.0, :), principal(2.0, :), principal(3.0, :),...
    nPoints, [], {mode, 2.0, 'DEFAULT'}, false, settings.output_strength{3.0}}]; %#ok<NASGU>
optimisationOutput = evalc('[best, scores] = abd.internal_optimise.main(inputs{:});');
if isempty(best) == true || isempty(best{5.0}) == false
    fprintf('%s', optimisationOutput);
end
end

%% RUN ANALYSIS WITHOUT PRINTING THE COMPLETE SUMMARY
function [result, analysisOutput] = runAnalysis(settings) %#ok<INUSD>
result = [];
analysisOutput = evalc('result = abd.main(settings);');
if isstruct(result) == false
    fprintf('%s', analysisOutput);
end
assert(isstruct(result) == true, 'Analysis did not return results.');
end

%% EVALUATE THE INITIAL LAYUP USING THE REQUESTED PLANE SEARCH
function [values] = evaluate(stress, material, sectionPoints, step)
material = repelem(material, 1.0, sectionPoints);
nPoints = size(stress, 2.0);
meanStress = 0.5.*(stress(1.0, :) + stress(2.0, :));
radius = sqrt((0.5.*(stress(1.0, :) - stress(2.0, :))).^2.0 + stress(3.0, :).^2.0);
[PF, MF, KF, SF, TF] = abd.internal_getLaRC05(nPoints, stress, false, meanStress + radius, meanStress - radius, zeros(1.0, nPoints),...
    material(7.0, :), material(1.0, :), material(2.0, :), material(3.0, :), material(4.0, :), material(5.0, :), material(6.0, :),...
    material(10.0, :), material(11.0, :), material(8.0, :), material(9.0, :), sectionPoints, step, 1.0);
values = [PF; MF; KF; SF; TF];
end

%% COMPARE NUMERIC OUTPUT WITH THE EXPECTED VALUES
function [] = checkValues(actual, expected, description)
tolerance = 1e-9*max(1.0, max(abs(expected(:))));
assert(isequal(size(actual), size(expected)) == true, '%s: output dimensions differ.', description);
assert(all(isfinite(actual(:))) == true, '%s: output contains non-finite values.', description);
assert(all(abs(actual(:) - expected(:)) <= tolerance) == true, '%s: incorrect output.', description);
end
