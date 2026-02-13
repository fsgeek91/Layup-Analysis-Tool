function [] = validate_joyce(varargin)
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
%   Layup Analysis Tool 5.1.1 Copyright Louis Vallance 2026
%   Last modified 13-Feb-2026 11:01:37 UTC

%__________________________________________________________________________
%% 1: JOB

JOB_NAME = mfilename;   JOB_DESCRIPTION = '';

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
MATERIAL = [181e3, 10.3e3, 7179.0, 0.28, 1e-5, 1e-5, 2e-3, 2e-3];

%% 3: LAYUP PROPERTIES
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
    sp: User-defined number of section points

    Note: The layup section is integrated once before the stress analysis;
    section points are thus treated as sample points.

    Note: Section points are evenly distributed over the layup. The total
    number of section points is SECTION_POINTS*length(STACKING_SEQUENCE).
%}
SECTION_POINTS = 'DEFAULT';

%% 4: LOAD MATRIX
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
    'job_name', JOB_NAME, 'job_description', JOB_DESCRIPTION,...           %% 1: JOB
    ...
    'material', {MATERIAL}, 'fail_stress', {[]}, 'fail_strain', {[]},...   %% 2: MATERIAL DATA
    'hashin', {[]}, 'larc05', {[]},...
    ...
    'stacking_sequence', STACKING_SEQUENCE, 'ply_thickness',...            %% 3: LAYUP PROPERTIES
    PLY_THICKNESS, 'symmetric_layup', SYMMETRIC_LAYUP, 'section_points',...
    SECTION_POINTS,...
    ...
    'load_mech', [NXX, NYY, NXY; MXX, MYY, MXY], 'load_therm', DELTA_T,... %% 4: LOAD MATRIX
    'load_hydro', DELTA_M,...
    ...
    'output_ply', OUTPUT_PLY, 'output_figure', {OUTPUT_FIGURE},...         %% 5: OUTPUT
    'output_strength', {{false, 'RESERVE'}}, 'output_optimised',...
    {{'', 'RESERVE', 'MINMAX', 5.0}},...
    'optimiser_settings', {{'MIXED-RADIX', 'DEFAULT', 'DEFAULT'}},...
    'output_location', {OUTPUT_LOCATION})...
    ...
    );