function [PARAM] = internal_correctSign(PARAM, MODE)
%   Correct the sign of strength properties.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 3.0.1 Copyright Louis Vallance 2024
%   Last modified 23-Feb-2024 13:20:04 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%

switch MODE
    case 1.0
        % Make negative properties positive
        badValues = PARAM < 0.0;
        PARAM(badValues) = abs(PARAM(badValues));
    case -1.0
        % Make positive properties negative
        badValues = PARAM > 0.0;
        PARAM(badValues) = -PARAM(badValues);
    case 2.0
        % Same as CASE(1.0), ignore -1.0 as this is treated as default

        % Make negative properties positive
        badValues = (PARAM < 0.0) & (PARAM ~= -1.0);
        PARAM(badValues) = abs(PARAM(badValues));
    otherwise
end