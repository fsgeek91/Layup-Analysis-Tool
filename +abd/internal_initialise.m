function [enableTensor, printTensor, materialDataMechanical, materialDataFailStress, materialDataFailStrain, materialDataHashin, materialDataLaRC05, theta, t_ply, symmetricPly,...
    SECTION_POINTS, OUTPUT_PLY, OUTPUT_FIGURE, OUTPUT_STRENGTH, OUTPUT_OPTIMISED, OPTIMISER_SETTINGS, OUTPUT_LOCATION, Nxx, Nyy, Nxy, Mxx, Myy, Mxy, deltaT, deltaM, error] =...
    internal_initialise(nargin, USER_INPUTS)
%   Gather variables from user inputs.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 4.2.0 Copyright Louis Vallance 2025
%   Last modified 10-Jun-2025 08:28:19 UTC
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
OPTIMISER_SETTINGS = [];
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
        [OUTPUT_PLY, OUTPUT_FIGURE, OUTPUT_STRENGTH, OUTPUT_OPTIMISED, OPTIMISER_SETTINGS, OUTPUT_LOCATION] = deal(USER_INPUTS{3.0}{1.0}, USER_INPUTS{3.0}{2.0},...
            USER_INPUTS{3.0}{3.0}, USER_INPUTS{3.0}{4.0}, USER_INPUTS{3.0}{5.0}, USER_INPUTS{3.0}{6.0});

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
        [OUTPUT_PLY, OUTPUT_FIGURE, OUTPUT_STRENGTH, OUTPUT_OPTIMISED, OPTIMISER_SETTINGS, OUTPUT_LOCATION] = deal(USER_INPUTS{3.0}{1.0}, USER_INPUTS{3.0}{2.0},...
            USER_INPUTS{3.0}{3.0}, USER_INPUTS{3.0}{4.0}, USER_INPUTS{3.0}{5.0}, USER_INPUTS{3.0}{6.0});

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
        [OUTPUT_PLY, OUTPUT_FIGURE, OUTPUT_STRENGTH, OUTPUT_OPTIMISED, OPTIMISER_SETTINGS, OUTPUT_LOCATION] = deal(USER_INPUTS{3.0}{1.0}, USER_INPUTS{3.0}{2.0},...
            USER_INPUTS{3.0}{3.0}, USER_INPUTS{3.0}{4.0}, USER_INPUTS{3.0}{5.0}, USER_INPUTS{3.0}{6.0});

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

if all(cellfun(@isempty, OUTPUT_OPTIMISED)) == true
    % Set default values if necessary
    OUTPUT_OPTIMISED = {'', 'RESERVE', 'MINMAX', 10.0};
end

%% Process OPTIMISER_SETTINGS
if iscell(OPTIMISER_SETTINGS) == false
    OPTIMISER_SETTINGS = {OPTIMISER_SETTINGS};
end

if all((cellfun(@isempty, OPTIMISER_SETTINGS) == true)) || (length(OPTIMISER_SETTINGS) < 3.0)
    % Set default values if necessary
    OPTIMISER_SETTINGS = {'MIXED-RADIX', 'DEFAULT', 'DEFAULT'};
end

% Process the first argument
argument = OPTIMISER_SETTINGS{1.0};

if ischar(argument) == false
    OPTIMISER_SETTINGS{1.0} = 2.0;
else
    argument = lower(argument);
    argument(ismember(argument, ' ')) = [];
    argument(ismember(argument, '-')) = [];

    switch argument
        case 'fullmatrix'
            OPTIMISER_SETTINGS{1.0} = 1.0;
        case 'mixedradix'
            OPTIMISER_SETTINGS{1.0} = 2.0;
        case 'chunks'
            OPTIMISER_SETTINGS{1.0} = 3.0;
        otherwise
            OPTIMISER_SETTINGS{1.0} = 2.0;
    end
end


%% Process output location
if iscell(OUTPUT_LOCATION) == false
    OUTPUT_LOCATION = {OUTPUT_LOCATION};
end

if length(OUTPUT_LOCATION) < 2.0
    OUTPUT_LOCATION{2.0} = false;
end

% Process the first argument
argument = OUTPUT_LOCATION{1.0};

if ischar(argument) == false
    OUTPUT_LOCATION{1.0} = 'DEFAULT';
else
    if strcmpi(argument, 'default') == true
        % Save results under OUTPUT folder inside user's PWD
        OUTPUT_LOCATION{1.0} = [pwd, filesep, 'output'];
    elseif strcmpi(argument, 'qft') == true
        % Save results under QFT folder structure inside user's PWD
        OUTPUT_LOCATION{1.0} = [pwd, filesep, 'Project', filesep, 'output'];
    end
end

% Process the second argument
argument = OUTPUT_LOCATION{2.0};

if islogical(argument) == false
    OUTPUT_LOCATION{2.0} = false;
end
end