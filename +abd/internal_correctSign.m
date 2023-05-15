function [PARAM] = internal_correctSign(PARAM, MODE)
%   Correct the sign of strength properties.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 2.5 Copyright Louis Vallance 2023
%   Last modified 15-May-2023 07:15:38 UTC
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
    otherwise
end