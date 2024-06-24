function [Q11, Q22, Q66, Q12, Qij] = internal_getReducedQ(E11, E22, V12, G12)
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
Q11 = (E11.^2.0)./(E11 - (V12.^2.0.*E22));
Q22 = (E11.*E22)./(E11 - (V12.^2.0.*E22));
Q66 = G12;
Q12 = (V12.*E11.*E22)./(E11 - (V12.^2.0.*E22));

% Q-matrix
N = length(E11);
Qij = zeros(3.0, 3.0, N);

% Concatenate into 4-D array
parfor i = 1.0:N
    Qij(:, :, i) = [Q11(i), Q12(i), 0.0; Q12(i), Q22(i), 0.0; 0.0, 0.0, Q66(i)];
end