function [error, OUTPUT_PLY_POINTS, plyBuffer, OUTPUT_ENVELOPE, ENVELOPE_MODE, outputApproximate, plyBuffer_sfailratio] =...
    internal_getOutputPoints(OUTPUT_PLY, z, z_points, nPlies, nPlies_points, plyBuffer, SECTION_POINTS, tolerance, enableTensor)
%   Get list of section points for stress/strain output.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 4.1.0 Copyright Louis Vallance 2025
%   Last modified 06-Jun-2025 11:07:25 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%
% Initialise output
error = false;
OUTPUT_PLY_POINTS = [];
OUTPUT_ENVELOPE =  false;
ENVELOPE_MODE = 1.0;
outputApproximate = [];

%{
    For strength output, the ply is considered to have failed when all of
    the section points in the layer have failed. Therefore, a version of
    the ply buffer containing all section points is required
%}
plyBuffer_sfailratio = [];

%% If cell, convert to CHAR/NUM
if (iscell(OUTPUT_PLY) == true) && (isempty(OUTPUT_PLY) == false)
    % Extract element from CELL (if applicable)
    OUTPUT_PLY = OUTPUT_PLY{1.0};
end

%% If definition is empty, use the default value and RETURN
if isempty(OUTPUT_PLY) == true
    OUTPUT_PLY = 'DEFAULT';
end

if isnumeric(OUTPUT_PLY) == true
    %% Process numeric definition
    invalidCondition = OUTPUT_PLY <= 0.0 | mod(OUTPUT_PLY, 1.0) ~= 0.0 | OUTPUT_PLY > length(z_points);
    OUTPUT_PLY(OUTPUT_PLY <= 0.0 | mod(OUTPUT_PLY, 1.0) ~= 0.0 | OUTPUT_PLY > length(z_points)) = [];

    % Number of section points must be a positive integer
    if (any(invalidCondition) == true) && (isempty(OUTPUT_PLY) == false)
        fprintf('[WARNING] Invalid section point numbers found in OUTPUT_PLY have been removed\n')
    end

    % Update the section point list
    OUTPUT_PLY_POINTS = OUTPUT_PLY;
elseif ischar(OUTPUT_PLY) == true
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
            if (strcmpi(OUTPUT_PLY, 'top') == true) || (strcmpi(OUTPUT_PLY, 'bottom') == true) || (strcmpi(OUTPUT_PLY, 'default') == true)
                if enableTensor == true
                    fprintf(['[ERROR] At least two section points are required for output to locations TOP\nand BOTTOM\n-> The available options for single section point output ar',...
                        'e:\n   MIDDLE, ALL, ENVELOPEABSMAX, ENVELOPEMAX, ENVELOPEMIN, [SECTION-POINT-LIST]\n'])

                    % Reset the error flag and RETURN
                    error = true;
                else
                    OUTPUT_PLY_POINTS = [];
                    plyBuffer = [];
                end
                
                return
            end
        case 2.0
            %{
                Output for MIDDLE faces is not supported when two section
                points are requested
            %}
            if (strcmpi(OUTPUT_PLY, 'middle') == true)
                if enableTensor == true
                    fprintf(['[ERROR] Output to location MIDDLE is not available when SECTION_POINTS = 2.0\n-> The available options for two section point output are:\n   DEFAULT,',...
                        ' TOP, BOTTOM, ALL, ENVELOPEABSMAX, ENVELOPEMAX, ENVELOPEMIN,\n   [SECTION-POINT-LIST]\n'])

                    % Reset the error flag and RETURN
                    error = true;
                else
                    OUTPUT_PLY_POINTS = [];
                    plyBuffer = [];
                end
                
                return
            end
        otherwise
            % Do not check other values of SECTION_POINTS
    end

    switch OUTPUT_PLY
        case 'default' % Top and bottom faces
            % Get z-points at ply boundaries
            OUTPUT_PLY_POINTS = find(ismember(z_points, z) == true);
        case 'top' % Top face only
            % Get z-points at ply boundaries (ignore bottom faces)
            OUTPUT_PLY_POINTS = nPlies_points/nPlies:nPlies_points/nPlies:nPlies_points;
        case 'middle' % Middle face only
            % Get z-values at midspan of each ply
            meanValues = (z(1.0:end - 1.0) + z(2.0:end))/2.0;

            % Reset very small values to zero
            meanValues(abs(meanValues) < tolerance) = 0.0;
            z_points(abs(z_points) < tolerance) = 0.0;
            
            % Get section points which lie at midspan locations
            OUTPUT_PLY_POINTS = find(any(abs(z_points' - meanValues) < tolerance, 2.0) == true)';

            if isempty(OUTPUT_PLY_POINTS) == true
                %{
                    There are no available section points for location
                    MIDDLE, so accept the section points at the ply average
                    z-value for each ply
                %}
                % Initialise section point loop index/buffer
                spIndex = 1.0;
                deviation_buffer = zeros(1.0, nPlies);

                for i = 1.0:nPlies
                    % Get all section point locations for the current ply
                    z_points_ply = z_points(spIndex:spIndex + (SECTION_POINTS - 1.0));
                    
                    % Assign the section point closest to the middle
                    [~, indexOfMin] = min(abs(z_points_ply - mean(z_points_ply)));
                    OUTPUT_PLY_POINTS(i) = indexOfMin + (spIndex - 1.0);

                    % Update the section point loop index
                    spIndex = spIndex + SECTION_POINTS;

                    %{
                        For the current ply, get the z-coordinates of one
                        of the two nearest middle points
                    %}
                    middle_z_point = z_points_ply(0.5*SECTION_POINTS);

                    % Get the z-coordinate of the mid plane of the first ply
                    middle_z_coord = meanValues(i);

                    % Compute the deviation for the current ply
                    if middle_z_point/middle_z_coord > 1.05
                        deviation_buffer(i) = abs(middle_z_point - middle_z_coord);
                    end
                end

                %{
                    If the maximum computed deviation exceeds 5%, set a
                    notice/warning and print this to the output file later
                %}
                if any(deviation_buffer > 0.0) == true
                    % Set the warning for the output file
                    outputApproximate = sprintf(['\nWarning: The maximum deviation of the midspan from its nearest section\npoints exceeds 5%%%% (%gmm) of the ply thickness.\n-> P',...
                        'recise output at location MIDDLE is unavailable\n-> For output that lies exactly on the layup midspan, specify an odd number\n   of section points\n-> Res',...
                        'ults will be written at approximate locations\n'], max(deviation_buffer));
                else
                    % Set the notice for the output file
                    outputApproximate = sprintf(['\nNote: Precise output at location MIDDLE is unavailable\n-> For output that lies exactly on the layup midspan, specify an odd nu',...
                        'mber\n   of section points\n-> Results will be written at approximate locations\n']);
                end
            end
        case 'bottom' % Bottom face only
            % Get z-points at ply boundaries (ignore top faces)
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
            % An invalid parameter was specified
            if enableTensor == true
                fprintf('[ERROR] Invalid parameter in OUTPUT_PLY: ''%s''\n', OUTPUT_PLY)

                % Reset the error flag and RETURN
                error = true;
            else
                OUTPUT_PLY_POINTS = [];
                plyBuffer = [];
            end
            
            return
    end
else
    % An invalid parameter was specified
    if enableTensor == true
        fprintf('[ERROR] Invalid value of OUTPUT_PLY\n')

        % Reset the error flag and RETURN
        error = true;
    else
        OUTPUT_PLY_POINTS = [];
        plyBuffer = [];
    end
    
    return
end

if isempty(OUTPUT_PLY_POINTS) == true
    %{
        The user-selected section point output has resulted in no result
        locations, so exit with an error
    %}
    fprintf(['[ERROR] There are no locations for result output.\nEither change the section points for output with OUTPUT_PLY, or increase\nthe number of section points for the cal',...
        'culation with SECTION_POINTS\n'])

    % Reset the error flag and RETURN
    error = true;
    return
end

% Update the ply buffer
plyBuffer_sfailratio = plyBuffer;
plyBuffer = plyBuffer(OUTPUT_PLY_POINTS);