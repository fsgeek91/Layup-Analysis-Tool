function [enableTensor, printTensor, materialDataMechanical,...
    materialDataFailStress, materialDataFailStrain, materialDataHashin,...
    theta, t_ply, symmetricPly, SECTION_POINTS, OUTPUT_PLY,...
    OUTPUT_FIGURE, OUTPUT_STRENGTH, OUTPUT_OPTIMISED, OUTPUT_LOCATION,...
    Nxx, Nyy, Nxy, Mxx, Myy, Mxy, deltaT, deltaM, error] =...
    internal_initialise(nargin, USER_INPUTS)
%   Gather variables from user inputs.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 2.6 Copyright Louis Vallance 2023
%   Last modified 16-Aug-2023 05:51:23 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%
% Initialise output
materialDataMechanical = [];
materialDataFailStress = [];
materialDataFailStrain = [];
materialDataHashin = [];
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

        % Layup data
        [theta, t_ply, symmetricPly, SECTION_POINTS] =...
            deal(USER_INPUTS{2.0}{1.0}, USER_INPUTS{2.0}{2.0},...
            USER_INPUTS{2.0}{3.0}, USER_INPUTS{2.0}{4.0});

        % Output data
        [OUTPUT_PLY, OUTPUT_FIGURE, OUTPUT_STRENGTH, OUTPUT_OPTIMISED,...
            OUTPUT_LOCATION] = deal(USER_INPUTS{3.0}{1.0},...
            USER_INPUTS{3.0}{2.0}, USER_INPUTS{3.0}{3.0},...
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

        % Layup data
        [theta, t_ply, symmetricPly, SECTION_POINTS] =...
            deal(USER_INPUTS{2.0}{1.0}, USER_INPUTS{2.0}{2.0},...
            USER_INPUTS{2.0}{3.0}, USER_INPUTS{2.0}{4.0});

        % Output data
        [OUTPUT_PLY, OUTPUT_FIGURE, OUTPUT_STRENGTH, OUTPUT_OPTIMISED,...
            OUTPUT_LOCATION] = deal(USER_INPUTS{3.0}{1.0},...
            USER_INPUTS{3.0}{2.0}, USER_INPUTS{3.0}{3.0},...
            USER_INPUTS{3.0}{4.0}, USER_INPUTS{3.0}{5.0});

        % Load matrix data
        [Nxx, Nyy, Nxy, Mxx, Myy, Mxy] = deal(USER_INPUTS{4.0}(1.0),...
            USER_INPUTS{4.0}(2.0), USER_INPUTS{4.0}(3.0),...
            USER_INPUTS{4.0}(4.0), USER_INPUTS{4.0}(5.0), USER_INPUTS{4.0}(6.0));
    case 5.0 % ABD + Load matrix + Thermo/hydro
        % Material data
        materialData = USER_INPUTS{1.0};
        materialDataMechanical = materialData(1.0);
        materialDataFailStress = materialData(2.0);
        materialDataFailStrain = materialData(3.0);
        materialDataHashin = materialData(4.0);

        % Layup data
        [theta, t_ply, symmetricPly, SECTION_POINTS] =...
            deal(USER_INPUTS{2.0}{1.0}, USER_INPUTS{2.0}{2.0},...
            USER_INPUTS{2.0}{3.0}, USER_INPUTS{2.0}{4.0});

        % Output data
        [OUTPUT_PLY, OUTPUT_FIGURE, OUTPUT_STRENGTH, OUTPUT_OPTIMISED,...
            OUTPUT_LOCATION] = deal(USER_INPUTS{3.0}{1.0},...
            USER_INPUTS{3.0}{2.0}, USER_INPUTS{3.0}{3.0},...
            USER_INPUTS{3.0}{4.0}, USER_INPUTS{3.0}{5.0});

        % Load matrix data
        [Nxx, Nyy, Nxy, Mxx, Myy, Mxy] = deal(USER_INPUTS{4.0}(1.0),...
            USER_INPUTS{4.0}(2.0), USER_INPUTS{4.0}(3.0),...
            USER_INPUTS{4.0}(4.0), USER_INPUTS{4.0}(5.0), USER_INPUTS{4.0}(6.0));

        % Thermo/hydro load data
        [deltaT, deltaM] = deal(USER_INPUTS{5.0}(1.0), USER_INPUTS{5.0}(2.0));
    otherwise
        % NARGIN is invalid, so RETURN
        fprintf(['[LAYUP-ANALYSIS-TOOL ERROR] An invalid number of arg',...
            'uments was specified\n']);
        error = true;
        return
end

%% Process output location
if strcmpi(OUTPUT_LOCATION, 'default') == true
    % Save results under OUTPUT folder inside user's PWD
    OUTPUT_LOCATION = [pwd, '\output'];
elseif strcmpi(OUTPUT_LOCATION, 'qft') == true
    % Save results under QFT folder structure inside user's PWD
    OUTPUT_LOCATION = [pwd, '\Project\output'];
end
end