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
%   Layup Analysis Tool 4.2.0 Copyright Louis Vallance 2025
%   Last modified 10-Jun-2025 08:28:19 UTC
%
%==========================================================================
%__________________________________________________________________________
%   USE CASE I - Stiffness calculation only (ABD matrix output):
%
%   ABD calculation of a 0.1mm thick UD [45(1)/90(1)] layup; return ABD
%   matrix and its inverse:
[~] = abd.main(struct(...
    'jobname', 'Use Case Ia',...                                           % Job name
    'jobdescription', sprintf(['ABD calculation of a 0.1mm thick UD [4',...% Job description
    '5(1)/90(1)] layup;\nreturn ABD matrix and its inverse']),...
    'material', {[2e5, 7e4, 5e3, 0.3, 1e-5, 1e-5, 2e-3, 2e-3]},...         % Mechanical material properties
    'stackingsequence', [45, 90],...                                       % Stacking sequence
    'plythickness', 0.1,...                                                % Ply thickness values                                   
    'outputlocation', {{'DEFAULT', false}}));                              % Analysis output location
%
%   ABD calculation of a variable thickness UD [45(1)/90(1)]-Symmetric
%   layup, different materials; return ABD matrices and the equivalent
%   bending/tension moduli:
[~] = abd.main(struct(...
    'jobname', 'Use Case Ib',...                                           % Job name
    'jobdescription', sprintf(['ABD calculation of a variable thicknes',...% Job description
    's UD\n[45(1)/90(1)]-Symmetric layup, different materials; return ',...
    'ABD matrices and\nthe equivalent bending/tension moduli']),...
    'material', {{[2e5, 7e4, 5e3, 0.3, 1e-5, 1e-5, 2e-3, 2e-3],...         % Mechanical material properties
    [1.5e5, 5e4, 3.5e3, 0.29, 1.2e-5, 1.2e-5, 1.5e-3, 1.5e-3]}},...
    'stackingsequence', [45, 90],...                                       % Stacking sequence
    'plythickness', [0.1, 0.15],...                                        % Ply thickness values
    'symmetriclayup', true,...                                             % Symmetric layup  
    'outputlocation', {{'DEFAULT', false}}));                              % Analysis output location
%__________________________________________________________________________
%   USE CASE II - Stress analysis:
%
%   Stress analysis of a 0.1mm thick UD [45(1)/90(1)] layup in pure bending
%   (x-axis); output for 3 section points per ply, output stresses at ply
%   midspans, return midspan strains and curvatures and stress/strain
%   tensors: 
[~] = abd.main(struct(...
    'jobname', 'Use Case II',...                                           % Job name
    'jobdescription', sprintf(['Stress analysis of a 0.1mm thick UD [4',...% Job description
    '5(1)/90(1)] layup in\npure bending (x-axis); output for 3 section',...
    ' points per ply, output stresses\nat ply midspans, return midspan',...
    ' strains and curvatures and stress/strain\ntensors']),...
    'material', {{[2e5, 7e4, 5e3, 0.3, 1e-5, 1e-5, 2e-3, 2e-3]}},...       % Mechanical material properties
    'stackingsequence', [45, 90],...                                       % Stacking sequence
    'plythickness', 0.1,...                                                % Ply thickness values 
    'sectionpoints', 3,...                                                 % Number of section points
    'loadmech', [0, 0, 0; 100, 0, 0],...                                   % Load matrix (mechanical)
    'outputply', 'MIDDLE',...                                              % Ply output location
    'outputfigure', {{'DEFAULT', [], 'SPLIT'}},...                         % MATLAB figures
    'outputlocation', {{'DEFAULT', false}}));                              % Analysis output location

% Load matrix
%__________________________________________________________________________
%   USE CASE III - Strength analysis:
%
%   Strength analysis of a variable thickness UD [45(1)/90(2)]-Symmetric
%   layup in combined bending (x-axis) and tension (y-direciton); output
%   for 2 section points per ply, output at bottom faces, disable MATLAB
%   figures, return failure measure components:
[~] = abd.main(struct(...
    'jobname', 'Use Case III',...                                          % Job name
    'jobdescription', sprintf(['Strength analysis of a variable thickn',...% Job description
    'ess UD\n[45(1)/90(2)]-Symmetric layup in combined bending (x-axis',...
    ') and tension\n(y-direciton); output for 2 section points per ply',...
    ', output at bottom faces,\ndisable MATLAB figures, return failure',...
    ' measure components']),...
    'material', {{[2e5, 7e4, 5e3, 0.3, 1e-5, 1e-5, 2e-3, 2e-3]}},...       % Mechanical material properties
    'failstress', {{[400, 300, 400, 300, 150, 1, 0]}},...                  % Fail stress properties
    'stackingsequence', [45, 90, 90],...                                   % Stacking sequence
    'plythickness', [0.15, 0.2, 0.2],...                                   % Ply thickness values 
    'symmetriclayup', true,...                                             % Symmetric layup  
    'sectionpoints', 2,...                                                 % Number of section points
    'loadmech', [0, 150, 0; -100, 0, 0],...                                % Load matrix (mechanical)
    'outputply', 'BOTTOM',...                                              % Ply output location
    'outputstrength', {{true, 'RESERVE'}},...                              % Strength calculation
    'outputlocation', {{'DEFAULT', false}}));                              % Analysis output location
%__________________________________________________________________________
%   USE CASE IV - Stacking sequence optimisation:
%
%   Stacking sequence optimisation of the layup example from USE CASE III
%   using the Tsai-Wu failure criterion (reserve); return results of the
%   stacking sequence optimisation:
[~] = abd.main(struct(...
    'jobname', 'Use Case IV',...                                           % Job name
    'jobdescription', sprintf(['Stacking sequence optimisation of the ',...% Job description
    'layup example from\nUSE CASE III using the Tsai-Wu failure criter',...
    'ion (reserve); return results\nof the stacking sequence optimisat',...
    'ion']),...
    'material', {{[2e5, 7e4, 5e3, 0.3, 1e-5, 1e-5, 2e-3, 2e-3]}},...       % Mechanical material properties
    'failstress', {{[400, 300, 400, 300, 150, 1, 0]}},...                  % Fail stress properties
    'stackingsequence', [45, 90, 90],...                                   % Stacking sequence
    'plythickness', [0.15, 0.2, 0.2],...                                   % Ply thickness values 
    'symmetriclayup', true,...                                             % Symmetric layup  
    'sectionpoints', 2,...                                                 % Number of section points
    'loadmech', [0, 150, 0; -100, 0, 0],...                                % Load matrix (mechanical)
    'outputply', 'BOTTOM',...                                              % Ply output location
    'outputfigure', {{'DEFAULT', [], 'SPLIT'}},...                         % MATLAB figures
    'outputstrength', {{true, 'RESERVE'}},...                              % Strength calculation
    'outputoptimised', {{'TSAIW', 'RESERVE', 'MINMAX', 10.0}},...          % Stacking sequence optimisation
    'outputlocation', {{'DEFAULT', false}}));                              % Analysis output location