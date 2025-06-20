function [EXT, EYT, GXYT, NUXYT, NUYXT, EXB, EYB, GXYB, NUXYB, NUYXB] =...
    internal_getModuli(t, ABD_INV)
%   Get the equivalent extensional and bending moduli.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 4.2.2 Copyright Louis Vallance 2025
%   Last modified 20-Jun-2025 07:44:10 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%
% Extensional moduli
EXT = 1.0/(t*ABD_INV(1.0, 1.0));
EYT = 1.0/(t*ABD_INV(2.0, 2.0));
GXYT = 1.0/(t*ABD_INV(3.0, 3.0));
NUXYT = -1.0*ABD_INV(2.0, 1.0)/ABD_INV(1.0, 1.0);
NUYXT = -1.0*ABD_INV(2.0, 1.0)/ABD_INV(2.0, 2.0);

% Bending moduli
EXB = 12.0/(t^3.0*ABD_INV(4.0, 4.0));
EYB = 12.0/(t^3.0*ABD_INV(5.0, 5.0));
GXYB = 12.0/(t^3.0*ABD_INV(6.0, 6.0));
NUXYB = -1.0*ABD_INV(5.0, 4.0)/ABD_INV(4.0, 4.0);
NUYXB = -1.0*ABD_INV(5.0, 4.0)/ABD_INV(5.0, 5.0);