function [symmetricAbd] = internal_getSymmetry(ABD, tolerance)
%   Gather variables from user inputs.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 2.3 Copyright Louis Vallance 2023
%   Last modified 09-May-2023 07:31:07 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%
symmetricAbd = true;

for i = 4.0:6.0
    for j = 1.0:3.0
        if (abs(ABD(i, j)) > tolerance) || (abs(ABD(j, i)) > tolerance)
            symmetricAbd = false;
            break
        end
    end
end