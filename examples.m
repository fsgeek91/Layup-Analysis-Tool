%EXAMPLES    Tutorial examples for Layup Analysis Tool.
%   Highlight the function call, then select the "Evaluate Selection in
%   Command Window" context option, or hit the F9 key to run an example.
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
%   Theory References:
%   Python code (GitHub): https://tinyurl.com/48phdvv8
%   Coefficient of Hydroscopic Expansion: https://tinyurl.com/2p9b4mvp
%   Definition of Laminate Plane Element Behaviour:
%   https://tinyurl.com/suwrhfcs
%
%   Resources:
%   Interactive Composite Laminate Calculator: https://tinyurl.com/2nww8yft
%   Finding Stiffness Matrices A, B and D (eFunda):
%   https://tinyurl.com/4tvscjk7
%
%   See also abd.main, user_definitions.
%
%   Layup Analysis Tool 4.3.0 Copyright Louis Vallance 2025
%   Last modified 17-Jun-2025 08:14:14 UTC
%
%==========================================================================
%__________________________________________________________________________
%   USE CASE I - Stiffness calculation only (ABD matrix output):
%
%   ABD calculation of a 0.1mm thick UD [45(1)/90(1)] layup; return ABD
%   matrix and its inverse:
[~] = abd.main(struct(...
    'job_name', 'Use Case Ia',...                                           % Job name
    'job_description', sprintf(['ABD calculation of a 0.1mm thick UD [',... % Job description
    '45(1)/90(1)] layup;\nreturn ABD matrix and its inverse']),...
    'material', {[2e5, 7e4, 5e3, 0.3, 1e-5, 1e-5, 2e-3, 2e-3]},...          % Mechanical material properties
    'stacking_sequence', [45, 90],...                                       % Stacking sequence
    'ply_thickness', 0.1,...                                                % Ply thickness values                                   
    'output_location', {{'DEFAULT', false}}));                              % Analysis output location
%
%   ABD calculation of a variable thickness UD [45(1)/90(1)]-Symmetric
%   layup, different materials; return ABD matrices and the equivalent
%   bending/tension moduli:
[~] = abd.main(struct(...
    'job_name', 'Use Case Ib',...                                           % Job name
    'job_description', sprintf(['ABD calculation of a variable thickne',... % Job description
    'ss UD\n[45(1)/90(1)]-Symmetric layup, different materials; return',...
    ' ABD matrices and\nthe equivalent bending/tension moduli']),...
    'material', {{[2e5, 7e4, 5e3, 0.3, 1e-5, 1e-5, 2e-3, 2e-3],...          % Mechanical material properties
    [1.5e5, 5e4, 3.5e3, 0.29, 1.2e-5, 1.2e-5, 1.5e-3, 1.5e-3]}},...
    'stacking_sequence', [45, 90],...                                       % Stacking sequence
    'ply_thickness', [0.1, 0.15],...                                        % Ply thickness values
    'symmetric_layup', true,...                                             % Symmetric layup  
    'output_location', {{'DEFAULT', false}}));                              % Analysis output location
%__________________________________________________________________________
%   USE CASE II - Stress analysis:
%
%   Stress analysis of a 0.1mm thick UD [45(1)/90(1)] layup in pure bending
%   (x-axis); output for 3 section points per ply, output stresses at ply
%   midspans, return midspan strains and curvatures and stress/strain
%   tensors: 
[~] = abd.main(struct(...
    'job_name', 'Use Case II',...                                           % Job name
    'job_description', sprintf(['Stress analysis of a 0.1mm thick UD [',... % Job description
    '45(1)/90(1)] layup in\npure bending (x-axis); output for 3 sectio',...
    'n points per ply, output stresses\nat ply midspans, return midspa',...
    'n strains and curvatures and stress/strain\ntensors']),...
    'material', {{[2e5, 7e4, 5e3, 0.3, 1e-5, 1e-5, 2e-3, 2e-3]}},...        % Mechanical material properties
    'stacking_sequence', [45, 90],...                                       % Stacking sequence
    'ply_thickness', 0.1,...                                                % Ply thickness values 
    'section_points', 3,...                                                 % Number of section points
    'load_mech', [0, 0, 0; 100, 0, 0],...                                   % Load matrix (mechanical)
    'output_ply', 'MIDDLE',...                                              % Ply output location
    'output_figure', {{'DEFAULT', [], 'SPLIT'}},...                         % MATLAB figures
    'output_location', {{'DEFAULT', false}}));                              % Analysis output location

% Load matrix
%__________________________________________________________________________
%   USE CASE III - Strength analysis:
%
%   Strength analysis of a variable thickness UD [45(1)/90(2)]-Symmetric
%   layup in combined bending (x-axis) and tension (y-direciton); output
%   for 2 section points per ply, output at bottom faces, disable MATLAB
%   figures, return failure measure components:
[~] = abd.main(struct(...
    'job_name', 'Use Case III',...                                          % Job name
    'job_description', sprintf(['Strength analysis of a variable thick',... % Job description
    'ness UD\n[45(1)/90(2)]-Symmetric layup in combined bending (x-axi',...
    's) and tension\n(y-direciton); output for 2 section points per pl',...
    'y, output at bottom faces,\ndisable MATLAB figures, return failur',...
    'e measure components']),...
    'material', {{[2e5, 7e4, 5e3, 0.3, 1e-5, 1e-5, 2e-3, 2e-3]}},...        % Mechanical material properties
    'fail_stress', {{[400, 300, 400, 300, 150, 1, 0]}},...                  % Fail stress properties
    'stacking_sequence', [45, 90, 90],...                                   % Stacking sequence
    'ply_thickness', [0.15, 0.2, 0.2],...                                   % Ply thickness values 
    'symmetric_layup', true,...                                             % Symmetric layup  
    'section_points', 2,...                                                 % Number of section points
    'load_mech', [0, 150, 0; -100, 0, 0],...                                % Load matrix (mechanical)
    'output_ply', 'BOTTOM',...                                              % Ply output location
    'output_figure', {{'DEFAULT', 'POINTS', 'SPLIT'}},...                   % MATLAB figures
    'output_strength', {{true, 'RESERVE'}},...                              % Strength calculation
    'output_location', {{'DEFAULT', false}}));                              % Analysis output location
%__________________________________________________________________________
%   USE CASE IV - Stacking sequence optimisation:
%
%   Stacking sequence optimisation of the layup example from USE CASE III
%   using the Tsai-Wu failure criterion (reserve); return results of the
%   stacking sequence optimisation:
[~] = abd.main(struct(...
    'job_name', 'Use Case IV',...                                           % Job name
    'job_description', sprintf(['Stacking sequence optimisation of the',... % Job description
    ' layup example from\nUSE CASE III using the Tsai-Wu failure crite',...
    'rion (reserve); return results\nof the stacking sequence optimisa',...
    'tion']),...
    'material', {{[2e5, 7e4, 5e3, 0.3, 1e-5, 1e-5, 2e-3, 2e-3]}},...        % Mechanical material properties
    'fail_stress', {{[400, 300, 400, 300, 150, 1, 0]}},...                  % Fail stress properties
    'stacking_sequence', [45, 90, 90],...                                   % Stacking sequence
    'ply_thickness', [0.15, 0.2, 0.2],...                                   % Ply thickness values 
    'symmetric_layup', true,...                                             % Symmetric layup  
    'section_points', 2,...                                                 % Number of section points
    'load_mech', [0, 150, 0; -100, 0, 0],...                                % Load matrix (mechanical)
    'output_ply', 'BOTTOM',...                                              % Ply output location
    'output_figure', {{'DEFAULT', 'POINTS', 'SPLIT'}},...                   % MATLAB figures
    'output_strength', {{true, 'RESERVE'}},...                              % Strength calculation
    'output_optimised', {{'TSAIW', 'RESERVE', 'MINMAX', 10.0}},...          % Stacking sequence optimisation
    'output_location', {{'DEFAULT', false}}));                              % Analysis output location