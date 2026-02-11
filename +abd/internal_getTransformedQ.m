function [Q11t, Q12t, Q16t, Q22t, Q26t, Q66t, Qt] = internal_getTransformedQ(N, theta, Q11, Q12, Q66, Q22)
%   Get the reduced Q-matrix.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 5.0.0 Copyright Louis Vallance 2026
%   Last modified 11-Feb-2026 08:06:52 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%
%% Get transformed stiffness matrices (ply coordinates)
Q11t = Q11.*cosd(theta).^4.0 + 2.0.*(Q12 + 2.0.*Q66).*cosd(theta).^2.0.*sind(theta).^2.0 + Q22.*sind(theta).^4.0;
Q12t = Q12.*(cosd(theta).^4.0 + sind(theta).^4.0) + (Q11 + Q22 - 4.0.*Q66).*cosd(theta).^2.0.*sind(theta).^2.0;
Q16t = (Q11 - Q12 - 2.0.*Q66).*cosd(theta).^3.0.*sind(theta) - (Q22 - Q12 - 2.0.*Q66).*cosd(theta).*sind(theta).^3.0;
Q22t = Q11.*sind(theta).^4.0 + 2.0.*(Q12 + 2.0.*Q66).*cosd(theta).^2.0.*sind(theta).^2.0 + Q22.*cosd(theta).^4.0;
Q26t = (Q11 - Q12 - 2.0.*Q66).*cosd(theta).*sind(theta).^3.0 - (Q22 - Q12 - 2.0.*Q66).*cosd(theta).^3.0.*sind(theta);
Q66t = (Q11 + Q22 - 2.0.*Q12 - 2.0.*Q66).*cosd(theta).^2.0.*sind(theta).^2.0 + Q66.*(cosd(theta).^4.0 + sind(theta).^4.0);

%% Concatenate stiffness matrix into N-D array
Qt = abd.internal_getAssembledNDArray(N, Q11t, Q22t, Q66t, Q12t, Q16t, Q26t);