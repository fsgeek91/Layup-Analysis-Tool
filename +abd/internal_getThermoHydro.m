function [axx, ayy, axy, bxx, byy, bxy] = internal_getThermoHydro(theta_points, A11_points, A22_points, B11_points, B22_points)
%   Get list of section points from user definition.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 3.0.5 Copyright Louis Vallance 2025
%   Last modified 11-Apr-2025 10:21:25 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%
% Compute the effective thermal constants for each ply
axx = A11_points.*cosd(theta_points).^2.0 + A22_points.*sind(theta_points).^2.0;
ayy = A11_points.*sind(theta_points).^2.0 + A22_points.*cosd(theta_points).^2.0;
axy = 2.0.*cosd(theta_points).*sind(theta_points).*(A11_points - A22_points);

% Compute the effective moisture constants for each ply
bxx = B11_points.*cosd(theta_points).^2.0 + B22_points.*sind(theta_points).^2.0;
byy = B11_points.*sind(theta_points).^2.0 + B22_points.*cosd(theta_points).^2.0;
bxy = 2.0.*cosd(theta_points).*sind(theta_points).*(B11_points - B22_points);

% Reset very small values to zero
axx(abs(axx) < 1e-12) = 0.0;
ayy(abs(ayy) < 1e-12) = 0.0;
axy(abs(axy) < 1e-12) = 0.0;
bxx(abs(bxx) < 1e-12) = 0.0;
byy(abs(bxy) < 1e-12) = 0.0;
bxy(abs(bxy) < 1e-12) = 0.0;