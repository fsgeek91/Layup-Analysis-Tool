function [Q11, Q22, Q66, Q12] = internal_getReducedQ(E11, E22, V12, G12)
%   Get the reduced Q-matrix.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 3.0.0 Copyright Louis Vallance 2024
%   Last modified 14-Feb-2024 15:05:03 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%
Q11 = (E11.^2.0)./(E11 - (V12.^2.0.*E22));
Q22 = (E11.*E22)./(E11 - (V12.^2.0.*E22));
Q66 = G12;
Q12 = (V12.*E11.*E22)./(E11 - (V12.^2.0.*E22));

% Q-matrix
%Qij = [Q11, Q12, 0.0; Q12, Q22, 0.0; 0.0, 0.0, Q66];