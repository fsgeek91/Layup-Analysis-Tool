function [error, z_points, theta, nPlies_points, A11, A22, B11, B22, plyBuffer, thickness, SECTION_POINTS, OUTPUT_PLY] =...
    internal_getSectionPoints(SECTION_POINTS, TAG, nPlies, theta, z, A11, A22, B11, B22, tolerance, OUTPUT_PLY)
%   Get list of section points from user definition.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 4.2.3 Copyright Louis Vallance 2025
%   Last modified 23-Jun-2025 14:28:39 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%
% Initialise output
error = false;
z_points = [];
nPlies_points = nPlies;
plyBuffer = [];
thickness = [];

if isempty(SECTION_POINTS) == true
    % The number of section points was not specified, so use minimum (2)
    SECTION_POINTS = 2.0;
end

% Get the number of section points from the "DEFAULT" parameter
if (ischar(SECTION_POINTS) == true) && (strcmpi(SECTION_POINTS, 'DEFAULT') == true)
    % Process the value of OUTPUT_PLY
    if (iscell(OUTPUT_PLY) == true) && (isempty(OUTPUT_PLY) == false)
        % Extract element from CELL (if applicable)
        OUTPUT_PLY = OUTPUT_PLY{1.0};
    end

    % Set the number of section points based on the definition of OUTPUT_PLY
    if isempty(OUTPUT_PLY) == true
        % OUTPUT_PLY is empty, so use default value of 2
        SECTION_POINTS = 2.0;
    else
        if ischar(OUTPUT_PLY) == true
            % OUTPUT_PLY is a parameter
            switch lower(OUTPUT_PLY)
                case 'middle'
                    % Ensure that there is a section point on the midspan
                    SECTION_POINTS = 3.0;
                otherwise
                    % Only 2 section points are required
                    SECTION_POINTS = 2.0;
            end
        elseif isnumeric(OUTPUT_PLY) == true
            %{
                OUTPUT_PLY is a section point list. Inform the user that
                they need to specify the number of section points directly
            %}
            str = strtrim(sprintf('%d, ', OUTPUT_PLY));
            fprintf('[ERROR] Invalid value of OUTPUT_PLY ([%s])\n-> For user-defined section point lists, the number of section points must\n   be specified with SECTION_POINTS\n',...
                str(1.0:end - 1.0));
            error = true;
            return
        else
            % Only 2 section points are required
            SECTION_POINTS = 2.0;
        end
    end
elseif (SECTION_POINTS == 1.0) && (strcmpi(OUTPUT_PLY, 'DEFAULT') == true)
    %{
        The user specified a single section point and the output location
        is set to DEFAULT, so reset the value of OUTPUT_PLY to MIDDLE in
        order to satisfy the request
    %}
    OUTPUT_PLY = 'MIDDLE';
end

if (SECTION_POINTS <= 0.0) || (mod(SECTION_POINTS, 1.0) ~= 0.0) || (isnan(SECTION_POINTS) == true) || (isinf(SECTION_POINTS) == true)
    % Number of section points must be a positive integer
    fprintf('[ERROR] Invalid value of %s. The number of section points\nmust be a positive integer\n', TAG);
    error = true;
    return
else
    % Reset output
    z_points = zeros(1.0, nPlies*SECTION_POINTS);
    theta_i = z_points;
    A11_i = z_points;
    A22_i = z_points;
    B11_i = z_points;
    B22_i = z_points;
    plyBuffer = z_points;
    index = 1.0;

    switch SECTION_POINTS
        case 1.0
            % One section point per ply: Use mid-span as ply location
            for i = 1:nPlies
                z_points(i) = mean([z(i), z(i + 1.0)]);
            end

            % Record the plies containing the current section points
            plyBuffer = 1.0:nPlies;
        otherwise
            % N section points per ply: Evenly distribute points over layup
            for i = 1.0:nPlies
                %{
                    Create the evenly distributed section point list for
                    the current ply
                %}
                z_points(index:index + (SECTION_POINTS - 1.0)) =...
                    linspace(z(i), z(i + 1.0), SECTION_POINTS);

                % Update buffers to include section points
                theta_i(index:index + (SECTION_POINTS - 1.0)) = theta(i);
                A11_i(index:index + (SECTION_POINTS - 1.0)) = A11(i);
                A22_i(index:index + (SECTION_POINTS - 1.0)) = A22(i);
                B11_i(index:index + (SECTION_POINTS - 1.0)) = B11(i);
                B22_i(index:index + (SECTION_POINTS - 1.0)) = B22(i);

                % Record the plies containing the current section points
                plyBuffer(index:index + (SECTION_POINTS - 1.0)) = i;

                % Update the loop index
                index = index + SECTION_POINTS;
            end

            % Re-assign buffers to original variables
            theta = theta_i;
            A11 = A11_i;
            A22 = A22_i;
            B11 = B11_i;
            B22 = B22_i;
            nPlies_points = length(theta);
    end
end

% Reset very small values of z_points to zero
z_points(abs(z_points) < tolerance) = 0.0;

if nPlies_points == 1.0
    % Indicate that there is only one section point in total
    thickness = 0.5;
else
    % Get the normalised thickness list (for tensor output)
    thickness = z_points/max(z);
end