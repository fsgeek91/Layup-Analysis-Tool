function [enableTensor, printTensor, materialDataMechanical,...
    materialDataFailStress, materialDataFailStrain, theta, t_ply,...
    symmetricPly, SECTION_POINTS, OUTPUT_PLY, OUTPUT_FIGURE,...
    OUTPUT_STRENGTH, OUTPUT_OPTIMISED, OUTPUT_LOCATION, Nxx, Nyy, Nxy,...
    Mxx, Myy, Mxy, deltaT, deltaM] = internal_initialise(nargin,...
    USER_INPUTS)
%   Gather variables from user inputs.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 2.2 Copyright Louis Vallance 2023
%   Last modified 09-May-2023 07:31:07 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%
% Flag to enable tensor output
enableTensor = true;
printTensor = 1.0;

switch nargin
    case 3.0
        materialData = USER_INPUTS{1.0};
        materialDataMechanical = materialData(1.0);
        materialDataFailStress = materialData(2.0);
        materialDataFailStrain = materialData(3.0);

        [theta, t_ply, symmetricPly, SECTION_POINTS] =...
            deal(USER_INPUTS{2.0}{1.0}, USER_INPUTS{2.0}{2.0},...
            USER_INPUTS{2.0}{3.0}, USER_INPUTS{2.0}{4.0});

        [OUTPUT_PLY, OUTPUT_FIGURE, OUTPUT_STRENGTH, OUTPUT_OPTIMISED,...
            OUTPUT_LOCATION] = deal(USER_INPUTS{3.0}{1.0},...
            USER_INPUTS{3.0}{2.0}, USER_INPUTS{3.0}{3.0},...
            USER_INPUTS{3.0}{4.0}, USER_INPUTS{3.0}{5.0});

        % Disable tensor output
        enableTensor = false;
    case 4.0
        materialData = USER_INPUTS{1.0};
        materialDataMechanical = materialData(1.0);
        materialDataFailStress = materialData(2.0);
        materialDataFailStrain = materialData(3.0);

        [theta, t_ply, symmetricPly, SECTION_POINTS] =...
            deal(USER_INPUTS{2.0}{1.0}, USER_INPUTS{2.0}{2.0},...
            USER_INPUTS{2.0}{3.0}, USER_INPUTS{2.0}{4.0});

        [OUTPUT_PLY, OUTPUT_FIGURE, OUTPUT_STRENGTH, OUTPUT_OPTIMISED,...
            OUTPUT_LOCATION] = deal(USER_INPUTS{3.0}{1.0},...
            USER_INPUTS{3.0}{2.0}, USER_INPUTS{3.0}{3.0},...
            USER_INPUTS{3.0}{4.0}, USER_INPUTS{3.0}{5.0});

        [Nxx, Nyy, Nxy, Mxx, Myy, Mxy] = deal(USER_INPUTS{4.0}(1.0),...
            USER_INPUTS{4.0}(2.0), USER_INPUTS{4.0}(3.0),...
            USER_INPUTS{4.0}(4.0), USER_INPUTS{4.0}(5.0), USER_INPUTS{4.0}(6.0));
    case 5.0
        materialData = USER_INPUTS{1.0};
        materialDataMechanical = materialData(1.0);
        materialDataFailStress = materialData(2.0);
        materialDataFailStrain = materialData(3.0);

        [theta, t_ply, symmetricPly, SECTION_POINTS] =...
            deal(USER_INPUTS{2.0}{1.0}, USER_INPUTS{2.0}{2.0},...
            USER_INPUTS{2.0}{3.0}, USER_INPUTS{2.0}{4.0});

        [OUTPUT_PLY, OUTPUT_FIGURE, OUTPUT_STRENGTH, OUTPUT_OPTIMISED,...
            OUTPUT_LOCATION] = deal(USER_INPUTS{3.0}{1.0},...
            USER_INPUTS{3.0}{2.0}, USER_INPUTS{3.0}{3.0},...
            USER_INPUTS{3.0}{4.0}, USER_INPUTS{3.0}{5.0});

        [Nxx, Nyy, Nxy, Mxx, Myy, Mxy] = deal(USER_INPUTS{4.0}(1.0),...
            USER_INPUTS{4.0}(2.0), USER_INPUTS{4.0}(3.0),...
            USER_INPUTS{4.0}(4.0), USER_INPUTS{4.0}(5.0), USER_INPUTS{4.0}(6.0));

        [deltaT, deltaM] = deal(USER_INPUTS{5.0}(1.0), USER_INPUTS{5.0}(2.0));
    otherwise
        fprintf(['[ABD ERROR] An invalid number of arguments was speci',...
            'fied\n']);
        return
end

%% Process output location
if strcmpi(OUTPUT_LOCATION, 'default') == true
    OUTPUT_LOCATION = [pwd, '\output'];
elseif strcmpi(OUTPUT_LOCATION, 'qft') == true
    OUTPUT_LOCATION = [pwd, '\Project\output'];
end
end