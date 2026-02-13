function [] = user_definitions(varargin)
%USER_DEFINITIONS    Helper script for Layup Analysis Tool.
%   Fill out this script with your layup definitions and analysis settings.
%   Read the tips above each option for usage hints.
%
%   Select "Run" from the Run section of the Editor, or hit the F5 key to
%   start the analysis.
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
%   See also abd.main, examples.
%
%   Layup Analysis Tool 5.1.1 Copyright Louis Vallance 2026
%   Last modified 13-Feb-2026 11:01:37 UTC

%__________________________________________________________________________
%% 1: JOB

JOB_NAME = mfilename;   JOB_DESCRIPTION = '';

%__________________________________________________________________________
%% 2: MATERIAL DATA
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
MATERIAL = [2e5, 7e4, 5e3, 0.3, 1e-5, 1e-5, 2e-3, 2e-3];

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
    Lamina in-plane/transverse shear strength.

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

%__________________________________________________________________________
%% 3: LAYUP PROPERTIES
% STACKING_SEQUENCE  Layup stacking sequence (bottom-up)
STACKING_SEQUENCE = [0.0, 45.0, 90.0];

% PLY_THICKNESS  Ply thickness list
%{
    [t1,... , tn]: Specify thickness values per ply
    t: Specify constant ply thickness

    Note: (t1) = Bottom; (tn) = Top.

    Note: When SYMMETRIC_LAYUP = true, ply thickness values are only
    required on one side of the symmetry plane.
%}
PLY_THICKNESS = 0.1;

% SYMMETRIC_LAYUP  Make the calculated section symmetric
SYMMETRIC_LAYUP = false;

% SECTION_POINTS  Number of stress/strain section points per ply
%{
    'Default': Program controlled
    sp: User-defined number of section points

    Note: The layup section is integrated once before the stress analysis;
    section points are thus treated as sample points.

    Note: Section points are evenly distributed over the layup. The total
    number of section points is SECTION_POINTS*length(STACKING_SEQUENCE).
%}
SECTION_POINTS = 'DEFAULT';

%__________________________________________________________________________
%% 4: LOAD MATRIX
% Mechanical load (forces)
NXX = 0.0;
NYY = 0.0;
NXY = 0.0;

% Mechanical load (moments)
MXX = 0.0;
MYY = 0.0;
MXY = 0.0;

% Thermal/hydroscopic load
DELTA_T = 0.0;
DELTA_M = 0.0;

%__________________________________________________________________________
%% 5: OUTPUT
% OUTPUT_PLY  Section points for stress/strain output
%{
    '<location>': Default (program controlled); Top; Middle (midspan/single
    section point); Bottom; All; EnvelopeAbsMax; EnvelopeMax; EnvelopeMin
    [SP1,..., SPn]: User-defined section point list
    (SP1 = Bottom; SPn = Top)
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
    '<layout>': Split (single MATLAB figure, separate plot for each tensor
    component); Compact (single MATLAB figure, one plot for all tensor
    components); DETACHED (separate MATLAB figure for each tensor
    component)

    Note: MATLAB figures are only generated when at least one load matrix
    component is specified with NXX/NYY/NXY or MXX/MYY/MXY.
%}
OUTPUT_FIGURE = {'DEFAULT', 'POINTS', 'SPLIT'};

% OUTPUT_STRENGTH  Evaluate static failure/damage initiation criteria
%{
    First argument (strength assessment):
    status: false (do not evaluate); true (evaluate based on available
    material data); @<ucrt> (function handle of user routine containing 
    user-defined failure criterion)
    
    Note: When the strength assessment is a user-defined failure criterion,
    the user criterion is evaluated in addition to all previously evaluated
    critera. Run the following command to generate a template user routine
    file: >> abd.createUcrt('<criterion-name>').

    Second argument (failure parameter):
    '<param>': Reserve (strength reserve factor); Value (criterion value)

    Note: The setting of the failure parameter only applies to the
    Tsai-Hill, 2D Hoffman, Tsai-Wu and Azzi-Tsai-Hill failure criteria.
%}
OUTPUT_STRENGTH = {false, 'RESERVE'};

% OUTPUT_OPTIMISED  Compute the optimised stacking sequence
%{
    Note: The stacking optimisation properties are taken from the number of
    plies and section points in the layup definition. The optimisation
    requires a load matrix definition (see Section 3) and the results of a
    strength evaluation using OUTPUT_STRENGTH = {true, <param>}.

    First argument (failure/damage initiation criterion):
    '<criterion>': Mstrs (Maximum stress); Tsaih (Tsai-Hill); Hoffman;
    Tsaiw (Tsai-Wu); Azzit (Azzi-Tsai-Hill); Mstrn (Maximum strain);
    Hashin; LaRC05; Ucrt (User-defined failure criterion)

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
    s: User-defined chunk size
    
    Third argument (tuning constant):
    'DEFAULT': Program-controlled
    k: User-defined constant (typically in the range 2-10)
%}
OPTIMISER_SETTINGS = {'MIXED-RADIX', 'DEFAULT', 'DEFAULT'};

% OUTPUT_LOCATION  Results output location
%{
    First argument (output location):
    'DEFAULT': Default output location
    <path>: User-defined results location
    
    Second argument (open results file);
    true: Open summary file after analysis
    false: Do not open summary file after analysis
%}
OUTPUT_LOCATION = {'DEFAULT', true};

%% - DO NOT EDIT BELOW LINE
%__________________________________________________________________________
%%

% Submit the layup for analysis!
[~] = abd.main(struct(...
    ...
    'job_name', JOB_NAME, 'job_description', JOB_DESCRIPTION,...                                                                                       %% 1: JOB
    ...
    'material', {MATERIAL}, 'fail_stress', {FAIL_STRESS}, 'fail_strain', {FAIL_STRAIN}, 'hashin', {HASHIN}, 'larc05', {LARC05},...                     %% 2: MATERIAL DATA
    ...
    'stacking_sequence', STACKING_SEQUENCE, 'ply_thickness', PLY_THICKNESS, 'symmetric_layup', SYMMETRIC_LAYUP, 'section_points', SECTION_POINTS,...   %% 3: LAYUP PROPERTIES
    ...
    'load_mech', [NXX, NYY, NXY; MXX, MYY, MXY], 'load_therm', DELTA_T, 'load_hydro', DELTA_M,...                                                      %% 4: LOAD MATRIX
    ...
    'output_ply', OUTPUT_PLY, 'output_figure', {OUTPUT_FIGURE}, 'output_strength', {OUTPUT_STRENGTH}, 'output_optimised', {OUTPUT_OPTIMISED},...       %% 5: OUTPUT
    'optimiser_settings', {OPTIMISER_SETTINGS}, 'output_location', {OUTPUT_LOCATION})...
    ...
    );