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
%   Layup Analysis Tool 2.7.3 Copyright Louis Vallance 2024
%   Last modified 12-Feb-2024 14:08:48 UTC
%
%==========================================================================
%__________________________________________________________________________
%   USE CASE I - Stiffness calculation only (ABD matrix output):
%
%   ABD calculation of a 0.1mm thick UD [45(1)/90(1)] layup; return ABD
%   matrix and its inverse:
[ABD, ABD_INV] =...                                                        % Output requests
    abd.main(...                                                           % Function name
    {[2e5, 7e4, 5e3, 0.3, 1e-5, 1e-5, 2e-3, 2e-3], [], [], []},...         % Material properties
    {[45, 90], 0.1, false, []},...                                         % Stacking sequence, ply thickness, symmetry flag
    {[], [], [], [], 'DEFAULT'});                                          % Output location
%
%   ABD calculation of a variable thickness UD [45(1)/90(1)]-Symmetric
%   layup, different materials; return ABD matrices and the equivalent
%   bending/tension moduli:
[ABD, ABD_INV, ~, ~, ~, EMTB] =...
    abd.main(...
    {{[2e5, 7e4, 5e3, 0.3, 1e-5, 1e-5, 2e-3, 2e-3], [1.5e5, 5e4, 3.5e3, 0.29, 1.2e-5, 1.2e-5, 1.5e-3, 1.5e-3]}, [], [], []},...
    {[45, 90], [0.1, 0.15], true, []},...
    {[], [], [], [], 'DEFAULT'});
%__________________________________________________________________________
%   USE CASE II - Stress analysis:
%
%   Stress analysis of a 0.1mm thick UD [45(1)/90(1)] layup in pure bending
%   (x-axis); output for 3 section points per ply, output stresses at ply
%   midspans, return midplane strains and curvatures and stress/strain
%   tensors:
[~, ~, EI, EP, SP] =...                                                    % Output requests
    abd.main(...                                                           % Function name
    {[2e5, 7e4, 5e3, 0.3, 1e-5, 1e-5, 2e-3, 2e-3], [], [], []},...         % Material properties
    {[45, 90], 0.1, false, 3},...                                          % Stacking sequence, ply thickness, symmetry flag, number of section points
    {'MIDDLE', {'DEFAULT', 'SPLIT'}, [], [], 'DEFAULT'},...                % Output location, MATLAB figure appearance, output location
    [0, 0, 0, 100, 0, 0]);                                                 % Load matrix
%__________________________________________________________________________
%   USE CASE III - Strength analysis:
%
%   Strength analysis of a variable thickness UD [45(1)/90(2)]-Symmetric
%   layup in combined bending (x-axis) and tension (y-direciton); output
%   for 2 section points per ply, output at bottom faces, disable MATLAB
%   figures, return failure measure components:
[~, ~, ~, ~, ~, ~, CFAILURE] =...
    abd.main(...
    {[2e5, 7e4, 5e3, 0.3, 1e-5, 1e-5, 2e-3, 2e-3], [400, 300, 400, 300, 150, 1, 0], [], []},...
    {[45, 90, 90], [0.15, 0.2, 0.2], true, 2},...
    {'BOTTOM', {[], 'SPLIT'}, {true, 'RESERVE'}, [], 'DEFAULT'},...
    [0, 150, 0, -100, 0, 0]);
%__________________________________________________________________________
%   USE CASE IV - Layup optimisation:
%
%   Stacking sequence optimisation of the layup example from USE CASE III
%   using the Tsai-Wu failure criterion (reserve); return results of the
%   stacking sequence optimisation:
[~, ~, ~, ~, ~, ~, ~, OPT] =...
    abd.main({[2e5, 7e4, 5e3, 0.3, 1e-5, 1e-5, 2e-3, 2e-3], [400, 300, 400, 300, 150, 1, 0], [], []},...
    {[45, 90, 90], [0.15, 0.2, 0.2], true, 2},...
    {'BOTTOM', {'DEFAULT', 'SPLIT'}, {true, 'RESERVE'}, {'TSAIW', 'RESERVE', 'MINMAX', 10.0}, 'DEFAULT'},...
    [0, 150, 0, -100, 0, 0]);