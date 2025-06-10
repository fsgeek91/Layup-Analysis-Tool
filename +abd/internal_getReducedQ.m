function [Q11, Q22, Q66, Q12, Qij] = internal_getReducedQ(N, E11, E22, V12, G12)
%   Get the reduced Q-matrix.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 4.2.0 Copyright Louis Vallance 2025
%   Last modified 10-Jun-2025 08:28:19 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%
%% Get stiffness matrices (global coordinates)
Q11 = (E11.^2.0)./(E11 - (V12.^2.0.*E22));
Q22 = (E11.*E22)./(E11 - (V12.^2.0.*E22));
Q66 = G12;
Q12 = (V12.*E11.*E22)./(E11 - (V12.^2.0.*E22));

%% Concatenate stiffness matrix into N-D array
Z = zeros(1.0, N);
Qij = abd.internal_getAssembledNDArray(N, Q11, Q22, Q66, Q12, Z, Z);