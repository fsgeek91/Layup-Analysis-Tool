function [varargout] = internal_getCriticalPly(DATA, symmetricAbd,...
    outputPoints, plyBuffer, nPlies)
%   Get the worst ply from the static failure assessment.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 2.3 Copyright Louis Vallance 2023
%   Last modified 09-May-2023 07:31:07 UTC
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
[rows, cols] = size(DATA);

% Initialise VARARGOUT
varargout = cell(1.0, 3.0*cols);

% Get indexes for worst plies and symmetric failure flags
dataIndexes = 1.0:3.0:3.0*cols;
valueIndexes = 2.0:3.0:3.0*cols;
symmetryIndexes = 3.0:3.0:3.0*cols;

for i = 1.0:cols
    % Get the current column from DATA
    DATA_i = DATA(:, i);

    % Get the unique values within tolerance
    [C, IA, IC] = uniquetol(DATA_i, 1e-6);

    % Replace elements with adjusted values
    DATA_i = C(IC);

    % Compute the single worst value over all section points for each ply
    DATA_ply = zeros(nPlies, 1.0);
    for p = 1.0:nPlies
        currentOutputPoints = outputPoints(plyBuffer == p);

        if isempty(currentOutputPoints) == true
            %{
                There is no data at the requested location
            %}
            DATA_ply(p) = -1.0;
        else
            DATA_ply(p) = max(DATA_i(currentOutputPoints));
        end
    end

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
end