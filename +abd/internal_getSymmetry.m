function [symmetricAbd] = internal_getSymmetry(ABD, tolerance)
%   Gather variables from user inputs.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 2.7.3 Copyright Louis Vallance 2024
%   Last modified 12-Feb-2024 14:08:48 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%
% Initialise output
symmetricAbd = true;

% Loop over ABD elements to check for symmetry
for i = 4.0:6.0
    for j = 1.0:3.0
        if (abs(ABD(i, j)) > tolerance) || (abs(ABD(j, i)) > tolerance)
            symmetricAbd = false;
            break
        end
    end
end