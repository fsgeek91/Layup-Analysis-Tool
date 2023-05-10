%   Validation case for Joyce. See page 22 of "Macro-mechanics of
%   lamina.pdf".
%
%   RUN THIS SCRIPT.
%
%   See also abd.main.
%
%   Layup Analysis Tool 2.4 Copyright Louis Vallance 2023
%   Last modified 10-May-2023 10:16:13 UTC

%% 1: MATERIAL DATA
% MATERIAL  Mechanical material properties
%{
    Note: When SYMMETRIC_LAYUP = true, material definitions are only
    required on one side of the symmetry plane.
%}
%{
TIP:

    Specify material properties for n plies as a 1xn cell array:
    Modulus of Elasticity (longitudinal/transverse);
    In-plane Shear Modulus; In-plane Poisson's Ratio;
    Coefficient of Thermal Expansion (longitudinal/transverse);
    Coefficient of Hydroscopic Expansion (longitudinal/transverse).

    MATERIAL = {[E11, E22, G12, V12, A11, A22, B11, B22](1),
                ...,
                [E11, E22, G12, V12, A11, A22, B11, B22](n)}

    Note: MATERIAL(1) = Bottom; MATERIAL(n) = Top.

	Units:
    Stress - [N/mm2]
    Thermal expansion - [1/degC]
    Hydroscopic expansion - [1/mm]
%}
MATERIAL = [181e3, 10.3e3, 7179.0, 0.28, 1e-5, 1e-5, 2e-3, 2e-3];

% FAIL_STRESS  Strength properties for stress-based failure criteria
%{
TIP:

    Specify strength properties for n plies as a 1xn cell array:
    Tensile/compressive stress limit (longitudinal);
    Tensile/compressive stress limit (transverse);
    Shear strength in the XY-plane; Cross-product coefficient;
    Biaxial stress limit.

    FAIL_STRESS = {[XT, XC, YT, YC, S, C, B](1),
                   ...,
                   [XT, XC, YT, YC, S, C, B](n)}

    Note: FAIL_STRESS(1) = Bottom; FAIL_STRESS(n) = Top.

    Note: If B = 0, the coupling term is computed from C.

    Units:
    Stress - [N/mm2]
%}
FAIL_STRESS = [];

% FAIL_STRAIN  Strength properties for strain-based failure criteria
%{
TIP:

    Specify strength properties for n plies as a 1xn cell array:
    Tensile/compressive strain limit (longitudinal);
    Tensile/compressive strain limit (transverse);
    Shear strain limit in the XY-plane.

    FAIL_STRAIN = {[XET, XEC, YET, YEC, SE](1),
                   ...,
                   [XET, XEC, YET, YEC, SE](n)}

    Note: FAIL_STRAIN(1) = Bottom; FAIL_STRAIN(n) = Top.

    Units:
    Strain - [mm/mm]
%}
FAIL_STRAIN = [];

% HASHIN  Strength properties for Hashin damage initiation criteria
%{
TIP:

    Specify strength properties for n plies as a 1xn cell array:
    Shear influence parameter;
    Lamina tensile/compressive strength (longitudinal);
    Lamina tensile/compressive strength (transverse);
    Lamina in-plane/transverse shear strength;

    HASHIN = {[ALPHA, XHT, XHC, YHT, YHC, SHX, SHY](1),
              ...,
              [ALPHA, XHT, XHC, YHT, YHC, SHX, SHY](n)}

    Note: HASHIN(1) = Bottom; HASHIN(n) = Top.

    Units:
    Stress - [N/mm2]
%}
HASHIN = [];

%% 2: LAYUP PROPERTIES
% STACKING_SEQUENCE  Layup stacking sequence (bottom-up) [degrees]
STACKING_SEQUENCE = 0.0;

% PLY_THICKNESS  Ply thickness list [mm]
%{
    [t1,... , tn]: Specify thickness values per ply
    t: Specify constant ply thickness

    Note: (t1) = Bottom; (tn) = Top.

	Note: When SYMMETRIC_LAYUP = true, ply thickness values are only
    required on one side of the symmetry plane.
%}
PLY_THICKNESS = 1.0;

% SYMMETRIC_LAYUP  Make the calculated section symmetric
SYMMETRIC_LAYUP = false;

% SECTION_POINTS  Number of stress/strain calculation points per ply
%{
    Note: Section points are evenly distributed over the layup. The total
    number of section points is SECTION_POINTS*length(STACKING_SEQUENCE).
%}
SECTION_POINTS = 2.0;

%% 3: LOAD MATRIX
% Mechanical load (forces) [N]
NXX = 2.0;
NYY = -3.0;
NXY = 4.0;

% Mechanical load (moments) [N.mm]
MXX = 0.0;
MYY = 0.0;
MXY = 0.0;

% Thermal/hydroscopic load
DELTA_T = 0.0; % [degC]
DELTA_M = 0.0; % [%/100 moisture weight content change]

%% 4: OUTPUT DEFINITION
%{
    '<location>': Default (top and bottom); Top; Middle; Bottom; All;
    EnvelopeAbsMax; EnvelopeMax; EnvelopeMin
    [SP1,..., SPn]: Section point list (1 = Bottom; n = Top)
%}
OUTPUT_PLY = 'DEFAULT';

% OUTPUT_FIGURE  MATLAB figure of stress and strain
%{
    '<mode>': [] (do not create figures); Default (no smoothing); Smooth
    (smooth output at ply boundaries)
%}
OUTPUT_FIGURE = [];

% OUTPUT_STRENGTH  Evaluate static failure criteria at each ply
OUTPUT_STRENGTH = false;

% OUTPUT_OPTIMISED  Compute the optimised stacking sequence
%{
    Note: The stacking optimisation properties are taken from the number of
    plies and section points in the layup definition.

    First argument (failure/damage initiation criterion):
    '<criterion>': Mstrs (Maximum stress); Tsaih (Tsai-Hill);
    Tsaiw (Tsai-Wu); Azzit (Azzi-Tsai-Hill); Mstrn (Maximum strain);
    Hashin;

    Second argument (failure parameter):
    '<param>': Reserve (reserve factor); Value (criterion value)

    Note: The reserve factor applies only to Tsai-Hill, Tsai-Wu and
    Azzi-Tsai-Hill failure criteria.

    Third argument (objective function):
    '<fun>': MinMax (minimise the maximum criterion value); MinMean
    (minimise the average criterion value)

    Fourth argument:
    theta: Angular step size [degrees]
%}
OUTPUT_OPTIMISED = {'', 'RESERVE', 'MINMAX', 10.0};

% OUTPUT_LOCATION  Results output location
%{
    'DEFAULT': Default output location
    '<location>': User-specified output location
%}
OUTPUT_LOCATION = 'DEFAULT';

%% - DO NOT EDIT BELOW LINE
%__________________________________________________________________________
%%

% Submit the layup for analysis!
[ABD1, ABD2, E_MIDPLANE, E_PLY, S_PLY, EQ_MODULI, CFAILURE, OPT_SEQ] =...
    ...
    abd.main({MATERIAL, FAIL_STRESS, FAIL_STRAIN, HASHIN},...
    {STACKING_SEQUENCE, PLY_THICKNESS, SYMMETRIC_LAYUP, SECTION_POINTS},...
    {OUTPUT_PLY, OUTPUT_FIGURE, OUTPUT_STRENGTH, OUTPUT_OPTIMISED,...
    OUTPUT_LOCATION}, [NXX, NYY, NXY, MXX, MYY, MXY], [DELTA_T, DELTA_M]);

clear MATERIAL FAIL_STRESS FAIL_STRAIN HASHIN STACKING_SEQUENCE...
    PLY_THICKNESS SYMMETRIC_LAYUP SECTION_POINTS NXX NYY NXY MXX MYY MXY...
    DELTA_T DELTA_M OUTPUT_PLY OUTPUT_FIGURE OUTPUT_STRENGTH...
    OUTPUT_OPTIMISED OUTPUT_LOCATION