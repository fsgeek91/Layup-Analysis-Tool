function [z, t] = internal_getThickness(nPlies, t_ply, tolerance)
%   Get z-coordinates from the layup definition.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 4.0.1 Copyright Louis Vallance 2025
%   Last modified 06-Jun-2025 11:07:25 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%
% Buffer for ply boundary locations
z = zeros(1.0, nPlies + 1.0);

% Total layup thickness
t = [];

% Check for invalid thickness values
if any(t_ply <= 0.0) == true
    fprintf('[ERROR] Zero or negative ply thickness values are not\nallowed\n');
    return
end

% Get the total ply thickness
t = sum(t_ply);

% Compute the first z-coordinate (ply half-thickness)
z(1.0) = -t/2.0;

% Get the remaining z-coordinates
for i = 2.0:nPlies + 1.0
    z(i) = z(i - 1.0) + t_ply(i-1);
end

% Reset near-zero values to zero
z(abs(z) < tolerance) = 0.0;