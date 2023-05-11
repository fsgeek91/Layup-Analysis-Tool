function [error, z_points, theta, nPlies_points, A11, A22, B11, B22,...
    plyBuffer, thickness]...
    = internal_getSectionPoints(DEFINITION, TAG, nPlies, theta, z, A11,...
    A22, B11, B22, tolerance)
%   Get list of section points from user definition.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 2.4 Copyright Louis Vallance 2023
%   Last modified 11-May-2023 13:34:37 UTC
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

if (DEFINITION <= 0.0) || (mod(DEFINITION, 1.0) ~= 0.0) ||...
        (isnan(DEFINITION) == true) || (isinf(DEFINITION) == true)
    fprintf(['[ABD ERROR] Invalid value of %s. The number of section p',...
        'oints\nmust be a positive integer\n'], TAG);
    error = true;
    return
else
    z_points = zeros(1.0, nPlies*DEFINITION);
    theta_i = z_points;
    A11_i = z_points;
    A22_i = z_points;
    B11_i = z_points;
    B22_i = z_points;
    plyBuffer = z_points;
    index = 1.0;

    switch DEFINITION
        case 1.0
            for i = 1:nPlies
                z_points(i) = mean([z(i), z(i + 1.0)]);
            end

            % Record the plies containing the current section points
            plyBuffer = 1.0:nPlies;
        otherwise
            for i = 1.0:nPlies
                z_points(index:index + (DEFINITION - 1.0)) =...
                    linspace(z(i), z(i + 1.0), DEFINITION);

                % Update buffers to include section points
                theta_i(index:index + (DEFINITION - 1.0)) = theta(i);
                A11_i(index:index + (DEFINITION - 1.0)) = A11(i);
                A22_i(index:index + (DEFINITION - 1.0)) = A22(i);
                B11_i(index:index + (DEFINITION - 1.0)) = B11(i);
                B22_i(index:index + (DEFINITION - 1.0)) = B22(i);

                % Record the plies containing the current section points
                plyBuffer(index:index + (DEFINITION - 1.0)) = i;

                index = index + DEFINITION;
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
    thickness = z_points/max(z_points);
end