function [axx, ayy, axy, bxx, byy, bxy] =...
    internal_getThermoHydro(theta_points, A11_points, A22_points,...
    B11_points, B22_points)
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
axx = A11_points.*cosd(theta_points).^2.0 + A22_points.*sind(theta_points).^2.0;
ayy = A11_points.*sind(theta_points).^2.0 + A22_points.*cosd(theta_points).^2.0;
axy = 2.0.*cosd(theta_points).*sind(theta_points).*(A11_points - A22_points);

bxx = B11_points.*cosd(theta_points).^2.0 + B22_points.*sind(theta_points).^2.0;
byy = B11_points.*sind(theta_points).^2.0 + B22_points.*cosd(theta_points).^2.0;
bxy = 2.0.*cosd(theta_points).*sind(theta_points).*(B11_points - B22_points);