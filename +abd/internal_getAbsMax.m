function [VARIABLE_MAX] = internal_getAbsMax(VARIABLE, MODE)
%   Get the absolute maximum value of an array.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 5.1.1 Copyright Louis Vallance 2026
%   Last modified 13-Feb-2026 11:01:37 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%
switch MODE
    case 1.0 % ABS MAX
        VARIABLE_ABS = abs(VARIABLE);
        [~, POSITION] = max(VARIABLE_ABS);
        VARIABLE_MAX = VARIABLE(POSITION);
    case 2.0 % MAX
        VARIABLE_MAX = max(VARIABLE);
    case 3.0 % MIN
        VARIABLE_MAX = min(VARIABLE);
    otherwise
end