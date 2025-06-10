function [varargout] = internal_getCriticalPly(DATA, symmetricAbd, plyBuffer, nPlies)
%   Get the worst ply from the static failure assessment.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 4.2.0 Copyright Louis Vallance 2025
%   Last modified 10-Jun-2025 08:28:19 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%

%% Adjust failure values
%{
    Identify failure values that are very close to each other (witin a
    toleance) and make them equal. This ensures that plies with the same
    failure value are not ignored due to small rounding errors
%}
% Get the number of columns to process individually
cols = width(DATA);

% Initialise VARARGOUT
varargout = cell(1.0, 3.0*cols + 1.0);

% Get indexes for worst plies and symmetric failure flags
dataIndexes = 1.0:3.0:3.0*cols;
valueIndexes = 2.0:3.0:3.0*cols;
symmetryIndexes = 3.0:3.0:3.0*cols;
sfailratioIndex = 3.0*cols + 1.0;

% Get list of all section points
nPoints = 1.0:length(plyBuffer);

% Initialise buffer for SFAILRATIO
SFAILRATIO = zeros(1.0, cols);

% Buffer to hold fail ply information for each criterion
FAILED_PLY_OUTER = false(nPlies, cols);

for i = 1.0:cols
    % Get the current column from DATA
    DATA_i = DATA(:, i);

    % Get the unique values within tolerance
    [C, ~, IC] = uniquetol(DATA_i, 1e-6);

    % Replace elements with adjusted values
    DATA_i = C(IC);

    % Compute the single worst value over all section points for each ply
    DATA_ply = zeros(nPlies, 1.0);

    % Buffer to record failed plies
    FAILED_PLY_INNER = false(1.0, nPlies);

    for p = 1.0:nPlies
        % Get maximum value in current ply based on output section points
        currentOutputPoints = nPoints(plyBuffer == p);

        % Get the values over all section points for the current ply
        DATA_i_all = DATA_i(currentOutputPoints);

        if isempty(currentOutputPoints) == true
            %{
                There is no data at the requested location. Since the
                calculation considers ALL section points, this condition
                should never be met!
            %}
            DATA_ply(p) = -1.0;
        else
            DATA_ply(p) = max(DATA_i_all);
        end

        if all(DATA_i_all >= 1.0)
            %{
                All section points in the ply have failed. Update the
                failed ply inner buffer
            %}
            FAILED_PLY_INNER(p) = true;
        end
    end

    % Update the value of SFAILRATIO
    SFAILRATIO(i) = length(FAILED_PLY_INNER(FAILED_PLY_INNER == true))/nPlies;

    % Update the failed ply outer buffer
    FAILED_PLY_OUTER(:, i) = FAILED_PLY_INNER;

    % Update VARARGOUT
    varargout{dataIndexes(i)} = find(DATA_ply == max(DATA_ply), 1.0);
    varargout{valueIndexes(i)} = DATA_ply;

    %% Check if the critical ply has a symmetric failure
    if symmetricAbd == true
        %{
            This check is only performed if the ABD matrix was previously
            confirmed to be symmetric
        %}
        if (all(find(flip(DATA_i) == max(DATA_i)) == find(DATA_i == max(DATA_i))) == true)
            %{
                After flipping the results for the current ply, the indexes
                of the worst plies are identical. This indicates that the
                corresponding plies on the other side of the symmetry plane
                have also failed.
            %}
            varargout{symmetryIndexes(i)} = 'YES';
        else
            % The ply failure is not symmetric
            varargout{symmetryIndexes(i)} = 'NO';
        end
    else
        % The symmetric failure flag does not apply in this case
        varargout{symmetryIndexes(i)} = 'N/A';
    end
end

% Output SFAILRATIO
varargout{sfailratioIndex} = SFAILRATIO;

% Output the failed ply outer buffer
varargout{sfailratioIndex + 1.0} = max(FAILED_PLY_OUTER, [], 2.0)';
end