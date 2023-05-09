function [error, OUTPUT_PLY_POINTS, plyBuffer, OUTPUT_ENVELOPE,...
    ENVELOPE_MODE, outputApproximate] = internal_getOutputPoints(...
    OUTPUT_PLY, z, z_points, nPlies, nPlies_points, plyBuffer,...
    SECTION_POINTS, tolerance)
%   Get list of section points for stress/strain output.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 2.2 Copyright Louis Vallance 2023
%   Last modified 09-May-2023 07:31:07 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%
% Initialise output
error = false;
OUTPUT_PLY_POINTS = [];
OUTPUT_ENVELOPE =  false;
ENVELOPE_MODE = 1.0;
outputApproximate = false;

%% If cell, convert to CHAR/NUM
if (iscell(OUTPUT_PLY) == 1.0) && (isempty(OUTPUT_PLY) == 0.0)
    % Extract element from CELL (if applicable)
    OUTPUT_PLY = OUTPUT_PLY{1.0};
end

%% If definition is empty, use the default value and RETURN
if isempty(OUTPUT_PLY) == 1.0
    OUTPUT_PLY = 'DEFAULT';
end

if isnumeric(OUTPUT_PLY) == 1.0
    %% Process numeric definition
    invalidCondition = OUTPUT_PLY <= 0.0 | mod(OUTPUT_PLY, 1.0) ~= 0.0 | OUTPUT_PLY > length(z_points);
    OUTPUT_PLY(OUTPUT_PLY <= 0.0 | mod(OUTPUT_PLY, 1.0) ~= 0.0 | OUTPUT_PLY > length(z_points)) = [];

    if (any(invalidCondition) == true) && (isempty(OUTPUT_PLY) == false)
        fprintf(['[ABD WARNING] Invalid section point numbers found in',...
            ' OUTPUT_PLY have been removed\n'])
    end

    OUTPUT_PLY_POINTS = OUTPUT_PLY;
elseif ischar(OUTPUT_PLY) == 1.0
    %% Process string definition
    % Remove white space and convert to lower case
    OUTPUT_PLY(ismember(OUTPUT_PLY, ' ')) = [];
    OUTPUT_PLY(ismember(OUTPUT_PLY, '-')) = [];
    OUTPUT_PLY = lower(OUTPUT_PLY);

    % Check for incompatible output combinations
    switch SECTION_POINTS
        case 1.0
            %{
                There is only one section point per ply, so output for TOP
                and BOTTOM faces is not supported
            %}
            if (strcmpi(OUTPUT_PLY, 'top') == true) ||...
                    (strcmpi(OUTPUT_PLY, 'bottom') == true) ||...
                    (strcmpi(OUTPUT_PLY, 'default') == true)
                fprintf(['[ABD ERROR] At least two section points are ',...
                    'required for output to locations TOP and BOTTOM\n'])
                error = true;
                return
            end
        case 2.0
            %{
                Output for MIDDLE faces is not supported when two section
                points are requested
            %}
            if (strcmpi(OUTPUT_PLY, 'middle') == true)
                fprintf(['[ABD ERROR] Output to location MIDDLE is not',...
                    ' available when SECTION_POINTS = 2.0\n'])
                error = true;
                return
            end
        otherwise
    end

    switch OUTPUT_PLY
        case 'default' % Top and bottom faces
            OUTPUT_PLY_POINTS = find(ismember(z_points, z) == true);
        case 'top' % Top face only
            OUTPUT_PLY_POINTS = nPlies_points/nPlies:nPlies_points/nPlies:nPlies_points;
        case 'middle' % Middle face only
            % Get z-values at midspan of each ply
            meanValues = (z(1.0:end - 1.0) + z(2.0:end))/2.0;

            % Reset very small values to zero
            meanValues(abs(meanValues) < tolerance) = 0.0;
            
            % Get section points which lie at midspan locations
            OUTPUT_PLY_POINTS = find(ismember(z_points, meanValues) == true);

            if isempty(OUTPUT_PLY_POINTS) == true
                %{
                    There are no available section points for location
                    MIDDLE, so accept the section points at the ply average
                    z-value for each ply
                %}
                % Initialise section point loop index
                spIndex = 1.0;

                for i = 1.0:nPlies
                    % Get all section point location for the current ply
                    z_points_ply = z_points(spIndex:spIndex + (SECTION_POINTS - 1.0));

                    
                    % Assign the section point closest to the middle
                    [~, indexOfMin] = min(abs(z_points_ply - mean(z_points_ply)));
                    OUTPUT_PLY_POINTS(i) = indexOfMin + (spIndex - 1.0);

                    % Update the section point loop index
                    spIndex = spIndex + SECTION_POINTS;
                end

                % Warn user later that output is at approximate locations
                outputApproximate = true;
            end
        case 'bottom' % Bottom face only
            OUTPUT_PLY_POINTS = 1.0:nPlies_points/nPlies:nPlies_points;
        case 'all' % All section points
            OUTPUT_PLY_POINTS = 1.0:nPlies_points;
        case 'envelopeabsmax' % Numerically largest (+ve/-ve) values over all section points in each ply
            %{
                The selected section points depend on the computed
                stresses. All section points will be requested and a flag
                set to perform the evnelope calculation during
                post-processing
            %}
            OUTPUT_PLY_POINTS = 1.0:nPlies_points;
            OUTPUT_ENVELOPE = true;
            ENVELOPE_MODE = 1.0;
        case 'envelopemax' % Largest (+ve) values over all section points in each ply
            %{
                The selected section points depend on the computed
                stresses. All section points will be requested and a flag
                set to perform the evnelope calculation during
                post-processing
            %}
            OUTPUT_PLY_POINTS = 1.0:nPlies_points;
            OUTPUT_ENVELOPE = true;
            ENVELOPE_MODE = 2.0;
        case 'envelopemin' % Largest (-ve) values over all section points in each ply
            %{
                The selected section points depend on the computed
                stresses. All section points will be requested and a flag
                set to perform the evnelope calculation during
                post-processing
            %}
            OUTPUT_PLY_POINTS = 1.0:nPlies_points;
            OUTPUT_ENVELOPE = true;
            ENVELOPE_MODE = 3.0;
        otherwise
            fprintf('[ABD ERROR] Invalid parameter in OUTPUT_PLY: ''%s''\n', OUTPUT_PLY)
            error = true;
            return
    end
else
    fprintf('[ABD ERROR] Invalid value of OUTPUT_PLY\n')
    error = true;
    return
end

if isempty(OUTPUT_PLY_POINTS) == true
    %{
        The user-selected section point output has resulted in no result
        locations, so exit with an error
    %}
    fprintf(['[ABD ERROR] There are no locations for result output. Ei',...
        'ther change the\n            section points for output with O',...
        'UTPUT_PLY, or increase the\n            number of section poi',...
        'nts for the calculation with\n            SECTION_POINTS\n'])
    error = true;
    return
end

% Update the ply buffer
plyBuffer = plyBuffer(OUTPUT_PLY_POINTS);