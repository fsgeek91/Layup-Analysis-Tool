function [Q] = internal_getAssembledNDArray(N, N1, N2, N3, S12, S16, S26)
%   Concatenate N 3x3 arrays into an N-dimensional array.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 3.0.5 Copyright Louis Vallance 2025
%   Last modified 11-Apr-2025 10:21:25 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%
%% Concatenate stiffness matrix into N-D array
% Assemble matrix terms
normals = [N1; N2; N3]';
shears = [S12; S16; S26]';

% Assign diagonal elements to each slice
diagonals = repmat(eye([3.0, 3.0]), 1.0, 1.0, N);
diagonals(diagonals > 0.0) = normals';

% Assign off-diagonal elements to each slice
nonDiagonals = repmat(tril(ones(3.0), -1.0), 1.0, 1.0, N);
nonDiagonals(nonDiagonals > 0.0) = shears';

% Assemble stiffness matrix
Q = diagonals + nonDiagonals + permute(nonDiagonals, [2.0, 1.0, 3.0]);