function [symmetricAbd] = internal_getSymmetry(ABD, tolerance)
%   Gather variables from user inputs.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 3.0.5 Copyright Louis Vallance 2025
%   Last modified 11-Apr-2025 10:21:25 UTC
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