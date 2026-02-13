function [E_midspan, E_ply_xy, S_ply_xy, E_ply_aligned, S_ply_aligned, E_therm_xy, E_moist_xy, E_therm_aligned, E_moist_aligned] =...
    internal_getTensor(ABD, Nxx, NxxT, NxxM, Nyy, NyyT, NyyM, Nxy, NxyT, NxyM, Mxx, MxxT, MxxM, Myy, MyyT, MyyM, Mxy, MxyT, MxyM, nPlies_points, z, theta_points, Qt, deltaT,...
    deltaM, axx, ayy, axy, bxx, byy, bxy, tolerance)
%   Get tensor quantities from ABD matrix.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 5.1.1 Copyright Louis Vallance 2026
%   Last modified 13-Feb-2026 11:01:37 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%

%% GET INDUCED MIDSPAN STRAIN
% Induced midspan strains and curvatures for the laminate
%{
    E_midspan = [epsilon_xx_0;
                  epsilon_yy_0;
                  epsilon_xy_0;
                  kappa_xx;
                  kappa_yy;
                  kappa_xy]
%}
E_midspan = ABD\[(Nxx + NxxT + NxxM);...
                  (Nyy + NyyT + NyyM);...
                  (Nxy + NxyT + NxyM);...
                  (Mxx + MxxT + MxxM);...
                  (Myy + MyyT + MyyM);...
                  (Mxy + MxyT + MxyM)];

% Reset very small values to zero
E_midspan(abs(E_midspan) < 1e-12) = 0.0;

%% GET X-Y STRAIN FOR EACH PLY
% Initialise stress/strain buffers for mechanical quantities
E_ply_xy = zeros(3.0, nPlies_points);
S_ply_xy = E_ply_xy;
E_ply_aligned = E_ply_xy;
S_ply_aligned = E_ply_xy;

% Initialise strain buffers for thermal/moisture quantities
E_therm_xy = E_ply_xy;
E_moist_xy = E_ply_xy;
E_therm_aligned = E_ply_xy;
E_moist_aligned = E_ply_xy;

% Engineering/tensor conversion
%A = 2.0;
A = 1.0;
%B = 0.5;
B = 1.0;

for i = 1:nPlies_points
    z_points = z(i);
    theta_point = theta_points(i);

    % Ply strains in X-Y coordinate system (evaluate at z)
    %{
        E_ply_xy = [epsilon_xx;
                    epsilon_yy;
                    gamma_xy]

        Note: gamma_xy = 2.0.*epsilon_xy
    %}
    E_ply_xy(:, i) = [E_midspan(1.0);
                      E_midspan(2.0);
                      A*E_midspan(3.0)] +...
                      ...
                      z_points*[E_midspan(4.0);
                             E_midspan(5.0);
                             E_midspan(6.0)];

    % Thermal/moiture strain components (XY)
    E_therm_xy(:, i) = deltaT*[axx(i); ayy(i); axy(i)];
    E_moist_xy(:, i) = deltaM*[bxx(i); byy(i); bxy(i)];

    % Ply stresses in X-Y coordinate system (evaluate at ply number)
    S_ply_xy(:, i) = Qt(:, :, i)*...
        [(E_ply_xy(1.0, i) - deltaT*axx(i) - deltaM*bxx(i));...
         (E_ply_xy(2.0, i) - deltaT*ayy(i) - deltaM*byy(i));...
         (B*E_ply_xy(3.0, i) - deltaT*axy(i) - deltaM*bxy(i))];

    % Transformation matrix for current ply
    T = [cosd(theta_point).^2.0, sind(theta_point).^2.0, 2.0.*cosd(theta_point).*sind(theta_point);...
        sind(theta_point).^2.0, cosd(theta_point).^2.0, -2.0.*cosd(theta_point).*sind(theta_point);...
        -cosd(theta_point).*sind(theta_point), cosd(theta_point).*sind(theta_point), cosd(theta_point).^2.0 - sind(theta_point).^2.0];

    % Compute the fibre/transverse strains for the current ply
    E_ply_aligned(:, i) = T*E_ply_xy(:, i);

    % Compute the fibre/transverse stresses for the current ply
    S_ply_aligned(:, i) = T*S_ply_xy(:, i);

    % Compute thermal/moisture strain for the current ply
    E_therm_aligned(:, i) = T*E_therm_xy(:, i);
    E_moist_aligned(:, i) = T*E_moist_xy(:, i);
end

% Remove near-zero stress values
S_ply_xy(abs(S_ply_xy) < tolerance) = 0.0;
S_ply_aligned(abs(S_ply_aligned) < tolerance) = 0.0;
E_ply_xy(abs(E_ply_xy) < 1e-12) = 0.0;
E_ply_aligned(abs(E_ply_aligned) < 1e-12) = 0.0;