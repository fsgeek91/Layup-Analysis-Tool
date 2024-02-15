function [enableTensor, printTensor, materialDataMechanical, materialDataFailStress, materialDataFailStrain, materialDataHashin, materialDataLaRC05, theta, t_ply, symmetricPly,...
    SECTION_POINTS, OUTPUT_PLY, OUTPUT_FIGURE, OUTPUT_STRENGTH, OUTPUT_OPTIMISED, OUTPUT_LOCATION, Nxx, Nyy, Nxy, Mxx, Myy, Mxy, deltaT, deltaM, error] =...
    internal_initialise(nargin, USER_INPUTS)
%   Gather variables from user inputs.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 3.0.0 Copyright Louis Vallance 2024
%   Last modified 14-Feb-2024 15:05:03 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%
% Initialise output
materialDataMechanical = [];
materialDataFailStress = [];
materialDataFailStrain = [];
materialDataHashin = [];
materialDataLaRC05 = [];
theta = [];
t_ply = [];
symmetricPly = [];
SECTION_POINTS = [];
OUTPUT_PLY = [];
OUTPUT_FIGURE = [];
OUTPUT_STRENGTH = [];
OUTPUT_OPTIMISED = [];
OUTPUT_LOCATION = [];
Nxx = [];
Nyy = [];
Nxy = [];
Mxx = [];
Myy = [];
Mxy = [];
deltaT = [];
deltaM = [];
error = false;

% Flag to enable tensor output
enableTensor = true;
printTensor = 1.0;

switch nargin
    case 3.0 % ABD only
        % Material data
        materialData = USER_INPUTS{1.0};
        materialDataMechanical = materialData(1.0);
        materialDataFailStress = materialData(2.0);
        materialDataFailStrain = materialData(3.0);
        materialDataHashin = materialData(4.0);
        materialDataLaRC05 = materialData(5.0);

        % Layup data
        [theta, t_ply, symmetricPly, SECTION_POINTS] = deal(USER_INPUTS{2.0}{1.0}, USER_INPUTS{2.0}{2.0}, USER_INPUTS{2.0}{3.0}, USER_INPUTS{2.0}{4.0});

        % Output data
        [OUTPUT_PLY, OUTPUT_FIGURE, OUTPUT_STRENGTH, OUTPUT_OPTIMISED, OUTPUT_LOCATION] = deal(USER_INPUTS{3.0}{1.0}, USER_INPUTS{3.0}{2.0}, USER_INPUTS{3.0}{3.0},...
            USER_INPUTS{3.0}{4.0}, USER_INPUTS{3.0}{5.0});

        % Disable tensor output
        enableTensor = false;
    case 4.0 % ABD + Load matrix
        % Material data
        materialData = USER_INPUTS{1.0};
        materialDataMechanical = materialData(1.0);
        materialDataFailStress = materialData(2.0);
        materialDataFailStrain = materialData(3.0);
        materialDataHashin = materialData(4.0);
        materialDataLaRC05 = materialData(5.0);

        % Layup data
        [theta, t_ply, symmetricPly, SECTION_POINTS] = deal(USER_INPUTS{2.0}{1.0}, USER_INPUTS{2.0}{2.0}, USER_INPUTS{2.0}{3.0}, USER_INPUTS{2.0}{4.0});

        % Output data
        [OUTPUT_PLY, OUTPUT_FIGURE, OUTPUT_STRENGTH, OUTPUT_OPTIMISED, OUTPUT_LOCATION] = deal(USER_INPUTS{3.0}{1.0}, USER_INPUTS{3.0}{2.0}, USER_INPUTS{3.0}{3.0},...
            USER_INPUTS{3.0}{4.0}, USER_INPUTS{3.0}{5.0});

        % Load matrix data
        [Nxx, Nyy, Nxy, Mxx, Myy, Mxy] = deal(USER_INPUTS{4.0}(1.0), USER_INPUTS{4.0}(2.0), USER_INPUTS{4.0}(3.0), USER_INPUTS{4.0}(4.0), USER_INPUTS{4.0}(5.0),...
            USER_INPUTS{4.0}(6.0));

        % Thermo/hydro load data
        [deltaT, deltaM] = deal(0.0, 0.0);
    case 5.0 % ABD + Load matrix + Thermo/hydro
        % Material data
        materialData = USER_INPUTS{1.0};
        materialDataMechanical = materialData(1.0);
        materialDataFailStress = materialData(2.0);
        materialDataFailStrain = materialData(3.0);
        materialDataHashin = materialData(4.0);
        materialDataLaRC05 = materialData(5.0);

        % Layup data
        [theta, t_ply, symmetricPly, SECTION_POINTS] = deal(USER_INPUTS{2.0}{1.0}, USER_INPUTS{2.0}{2.0}, USER_INPUTS{2.0}{3.0}, USER_INPUTS{2.0}{4.0});

        % Output data
        [OUTPUT_PLY, OUTPUT_FIGURE, OUTPUT_STRENGTH, OUTPUT_OPTIMISED, OUTPUT_LOCATION] = deal(USER_INPUTS{3.0}{1.0}, USER_INPUTS{3.0}{2.0}, USER_INPUTS{3.0}{3.0},...
            USER_INPUTS{3.0}{4.0}, USER_INPUTS{3.0}{5.0});

        % Load matrix data
        [Nxx, Nyy, Nxy, Mxx, Myy, Mxy] = deal(USER_INPUTS{4.0}(1.0), USER_INPUTS{4.0}(2.0), USER_INPUTS{4.0}(3.0), USER_INPUTS{4.0}(4.0), USER_INPUTS{4.0}(5.0),...
            USER_INPUTS{4.0}(6.0));

        % Thermo/hydro load data
        [deltaT, deltaM] = deal(USER_INPUTS{5.0}(1.0), USER_INPUTS{5.0}(2.0));
    otherwise
        % Report error to the user
        fprintf('[ERROR] An invalid number of arguments was specified\n');

        if nargin == 0.0
            % The user probably ran main.abd by mistake
            fprintf('-> To perform a layup analysis, run the input file ''user_definitions.m '' directly');
        end

        % Set the error flag and RETURN
        error = true;
        return
end

%% Process OUTPUT_OPTIMISED
if iscell(OUTPUT_OPTIMISED) == false
    OUTPUT_OPTIMISED = {OUTPUT_OPTIMISED};
end

if cellfun(@isempty, OUTPUT_OPTIMISED) == true
    % Set default values if necessary
    OUTPUT_OPTIMISED = {'', 'RESERVE', 'MINMAX', 10.0};
end

%% Process output location
if strcmpi(OUTPUT_LOCATION, 'default') == true
    % Save results under OUTPUT folder inside user's PWD
    OUTPUT_LOCATION = [pwd, filesep, 'output'];
elseif strcmpi(OUTPUT_LOCATION, 'qft') == true
    % Save results under QFT folder structure inside user's PWD
    OUTPUT_LOCATION = [pwd, filesep, 'Project', filesep, 'output'];
end
end