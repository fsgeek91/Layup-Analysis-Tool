function [Q11t, Q12t, Q16t, Q22t, Q26t, Q66t, Qt] = internal_getTransformedQ(theta, Q11, Q12, Q66, Q22)
%   Get the reduced Q-matrix.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 3.0.3 Copyright Louis Vallance 2024
%   Last modified 24-Jun-2024 11:37:46 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%
Q11t = Q11.*cosd(theta).^4.0 + 2.0.*(Q12 + 2.0.*Q66).*cosd(theta).^2.0.*sind(theta).^2.0 + Q22.*sind(theta).^4.0;
Q12t = Q12.*(cosd(theta).^4.0 + sind(theta).^4.0) + (Q11 + Q22 - 4.0.*Q66).*cosd(theta).^2.0.*sind(theta).^2.0;
Q16t = (Q11 - Q12 - 2.0.*Q66).*cosd(theta).^3.0.*sind(theta) - (Q22 - Q12 - 2.0.*Q66).*cosd(theta).*sind(theta).^3.0;
Q22t = Q11.*sind(theta).^4.0 + 2.0.*(Q12 + 2.0.*Q66).*cosd(theta).^2.0.*sind(theta).^2.0 + Q22.*cosd(theta).^4.0;
Q26t = (Q11 - Q12 - 2.0.*Q66).*cosd(theta).*sind(theta).^3.0 - (Q22 - Q12 - 2.0.*Q66).*cosd(theta).^3.0.*sind(theta);
Q66t = (Q11 + Q22 - 2.0.*Q12 - 2.0.*Q66).*cosd(theta).^2.0.*sind(theta).^2.0 + Q66.*(cosd(theta).^4.0 + sind(theta).^4.0);

% Transformed Q-matrix
N = length(theta);
Qt = zeros(3.0, 3.0, N);

% Concatenate into 4-D array
parfor i = 1.0:N
    Qt(:, :, i) = [Q11t(i), Q12t(i), Q16t(i); Q12t(i), Q22t(i), Q26t(i); Q16t(i), Q26t(i), Q66t(i)];
end