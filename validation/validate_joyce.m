%   Validation case for Joyce. See page 22 of "Macro-mechanics of
%   lamina.pdf".
%
%   RUN THIS SCRIPT.
%
%   Units:
%   Force [N]
%   Length [mm]
%   Angle [degrees]
%   Moment [N.mm]
%   Stress [N/mm2]
%   Temperature [degC]
%   Thermal expansion [1/degC]
%   Moisture [%/100 moisture weight content change]
%   Moisture expansion [1/mm]
%
%   See also abd.main.
%
%   Layup Analysis Tool 4.0.1 Copyright Louis Vallance 2025
%   Last modified 06-Jun-2025 11:07:25 UTC

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
%}
HASHIN = [];

% LARC05  Strength properties for LaRC05 damage initiation criteria
%{
TIP:

    Specify strength properties for n plies as a 1xn cell array:
    Tensile/compressive strength (longitudinal);
    Tensile/compressive strength (transverse);
    In-plane/transverse shear strength;
    Shear modulus in the 12-plane;
    Longitudinal/transverse shear friction coefficient;
    Fracture plane angle for pure compression;
    Misalignment angle at failure for pure compression.

    LARC05 = {[XLT, XLC, YLT, YLC*, SLX, SLY*, GL12, NL*, NT*, A0*, PHI0*](1),
              ...,
              [XLT, XLC, YLT, YLC*, SLX, SLY*, GL12, NL*, NT*, A0*, PHI0*](n)}

    Note: LARC05(1) = Bottom; LARC05(n) = Top.
    
    Note: Parameters marked with an asterisk (*) are derived if a value of
    -1 is specified.
%}
LARC05 = [];

%% 2: LAYUP PROPERTIES
% STACKING_SEQUENCE  Layup stacking sequence (bottom-up)
STACKING_SEQUENCE = 0.0;

% PLY_THICKNESS  Ply thickness list
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

% SECTION_POINTS  Number of stress/strain section points per ply
%{
    'Default': Program controlled
    sp: User-defined

    Note: The layup section is integrated once before the stress analysis;
    section points are thus treated as sample points.

    Note: Section points are evenly distributed over the layup. The total
    number of section points is SECTION_POINTS*length(STACKING_SEQUENCE).
%}
SECTION_POINTS = 'DEFAULT';

%% 3: LOAD MATRIX
% Mechanical load (forces)
NXX = 2.0;
NYY = -3.0;
NXY = 4.0;

% Mechanical load (moments)
MXX = 0.0;
MYY = 0.0;
MXY = 0.0;

% Thermal/hydroscopic load
DELTA_T = 0.0;
DELTA_M = 0.0;

%% 4: OUTPUT DEFINITION
% OUTPUT_PLY  Section points for stress/strain output
%{
    '<location>': Default (top and bottom); Top; Middle; Bottom; All;
    EnvelopeAbsMax; EnvelopeMax; EnvelopeMin
    [SP1,..., SPn]: Section point list (1 = Bottom; n = Top)
%}
OUTPUT_PLY = 'DEFAULT';

% OUTPUT_FIGURE  MATLAB figure of stress and strain
%{
    First argument (MATLAB figures):
    '<mode>': [] (do not create figures); Default (no smoothing); Smooth
    (smooth output at ply boundaries)
    
    Second argument (section points):
    '<mode>': [] (do not visualise section points); Points
    (plot all section points)

    Third argument (MATLAB figure layout):
    '<layout>': Split (one figure per tensor component); Compact (single
    figure for all tensor components)

    Note: MATLAB figures are only generated when at least one load matrix
    component is specified with NXX/NYY/NXY or MXX/MYY/MXY.
%}
OUTPUT_FIGURE = {'DEFAULT', 'POINTS', 'SPLIT'};

% OUTPUT_STRENGTH  Evaluate static failure/damage initiation criteria
%{
    First argument (strength assessment):
    status: false (do not evaluate); true (evaluate based on available
    material data)

    Second argument (failure parameter):
    '<param>': Reserve (strength reserve factor); Value (criterion value)

    Note: The setting of the failure parameter only applies to the
    Tsai-Hill, Tsai-Wu and Azzi-Tsai-Hill failure criteria.
%}
OUTPUT_STRENGTH = {false, 'RESERVE'};

% OUTPUT_OPTIMISED  Compute the optimised stacking sequence
%{
    Note: The stacking optimisation properties are taken from the number of
    plies and section points in the layup definition. The optimisation
    requires a load matrix definition (see Section 3) and the results of a
    strength evaluation using OUTPUT_STRENGTH = {true, <param>}.

    First argument (failure/damage initiation criterion):
    '<criterion>': Mstrs (Maximum stress); Tsaih (Tsai-Hill);
    Tsaiw (Tsai-Wu); Azzit (Azzi-Tsai-Hill); Mstrn (Maximum strain);
    Hashin; LaRC05

    Second argument (failure parameter):
    '<param>': Reserve (strength reserve factor); Value (criterion value)

    Third argument (objective function):
    '<fun>': MinMax (minimise the maximum criterion value); MinMean
    (minimise the average criterion value)

    Fourth argument (precision):
    theta: Angular step size
%}
OUTPUT_OPTIMISED = {'', 'RESERVE', 'MINMAX', 5.0};

% OPTIMISER_SETTINGS  Optimiser solver settings (advanced)
%{
    First argument (solver):
    '<param>': Full matrix; Mixed-radix; Chunks
    
    Second argument (chunk size):
    'DEFAULT': Program controlled
    n: User-defined
    
    Third argument (tuning constant):
    'DEFAULT': Program-controlled
    k: User-defined
%}
OPTIMISER_SETTINGS = {'MIXED-RADIX', 'DEFAULT', 'DEFAULT'};

% OUTPUT_LOCATION  Results output location
%{
    First argument (output location):
    'DEFAULT': Default output location
    '<location>': User-specified output location
    
    Second argument (open results file);
    true: Open results file after analysis
    false: Do not open results file after analysis
%}
OUTPUT_LOCATION = {'DEFAULT', true};

%% - DO NOT EDIT BELOW LINE
%__________________________________________________________________________
%%

% Submit the layup for analysis!
[ABD1, ABD2, Q, E_MIDPLANE, E_PLY, S_PLY, EQ_MODULI, CFAILURE, OPT_SEQ] =...
    ...
    abd.main({MATERIAL, FAIL_STRESS, FAIL_STRAIN, HASHIN, LARC05}, {STACKING_SEQUENCE, PLY_THICKNESS, SYMMETRIC_LAYUP, SECTION_POINTS}, {OUTPUT_PLY, OUTPUT_FIGURE,...
    OUTPUT_STRENGTH, OUTPUT_OPTIMISED, OPTIMISER_SETTINGS, OUTPUT_LOCATION}, [NXX, NYY, NXY, MXX, MYY, MXY], [DELTA_T, DELTA_M]);

clear MATERIAL FAIL_STRESS FAIL_STRAIN HASHIN LARC05 STACKING_SEQUENCE PLY_THICKNESS SYMMETRIC_LAYUP SECTION_POINTS NXX NYY NXY MXX MYY MXY DELTA_T DELTA_M OUTPUT_PLY...
    OUTPUT_FIGURE OUTPUT_STRENGTH OUTPUT_OPTIMISED OPTIMISER_SETTINGS OUTPUT_LOCATION