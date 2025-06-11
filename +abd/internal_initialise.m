function [enableTensor, printTensor, materialDataMechanical, materialDataFailStress, materialDataFailStrain, materialDataHashin, materialDataLaRC05, theta, t_ply, symmetricPly,...
    SECTION_POINTS, OUTPUT_PLY, OUTPUT_FIGURE, OUTPUT_STRENGTH, OUTPUT_OPTIMISED, OPTIMISER_SETTINGS, OUTPUT_LOCATION, Nxx, Nyy, Nxy, Mxx, Myy, Mxy, deltaT, deltaM, jobName,...
    jobDescription, settings, error] =...
    internal_initialise(settings)
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
jobName = [];
jobDescription = [];
error = false;

requiredFields = {'jobname', 'jobdescription', 'material', 'failstress', 'failstrain', 'hashin', 'larc05', 'stackingsequence', 'plythickness', 'symmetriclayup', 'sectionpoints',...
    'loadmech', 'loadtherm', 'loadmoist', 'outputply', 'outputfigure', 'outputstrength', 'outputoptimised', 'optimisersettings', 'outputlocation'};

for i = 1:length(requiredFields)
    % Get the current field name
    currentField = requiredFields{i};

    % Check that the current field exists
    if isfield(settings, currentField) == false
        % Update the user
        fprintf('[ERROR] Undefined option ''%s''\n', upper(currentField));
        error = true;
        return
    end
end

% Flag to enable tensor output
enableTensor = true;
printTensor = true;

% Extract inputs from FLAGS structure
jobName = settings.jobname;
jobDescription = settings.jobdescription;
materialDataMechanical = settings.material;
materialDataFailStress = settings.failstress;
materialDataFailStrain = settings.failstrain;
materialDataHashin = settings.hashin;
materialDataLaRC05 = settings.larc05;
theta = settings.stackingsequence;
t_ply = settings.plythickness;
symmetricPly = settings.symmetriclayup;
SECTION_POINTS = settings.sectionpoints;
loadmech = settings.loadmech;
if isempty(loadmech) == false
    Nxx = loadmech(1.0, 1.0);
    Nyy = loadmech(1.0, 2.0);
    Nxy = loadmech(1.0, 3.0);
    Mxx = loadmech(2.0, 1.0);
    Myy = loadmech(2.0, 2.0);
    Mxy = loadmech(2.0, 3.0);
end
deltaT = settings.loadtherm;
if isempty(deltaT) == true
    deltaT = 0.0;
end
deltaM = settings.loadmoist;
if isempty(deltaM) == true
    deltaM = 0.0;
end
OUTPUT_PLY = settings.outputply;
OUTPUT_FIGURE = settings.outputfigure;
OUTPUT_STRENGTH = settings.outputstrength;
OUTPUT_OPTIMISED = settings.outputoptimised;
OPTIMISER_SETTINGS = settings.optimisersettings;
OUTPUT_LOCATION = settings.outputlocation;

% Disable tensor output if applicable
if all(all(loadmech == 0.0)) == true && all(all([deltaT, deltaM] == 0.0)) == true
    enableTensor = false;
end

%% Process OUTPUT_FIGURE
if iscell(OUTPUT_FIGURE) == false
    OUTPUT_FIGURE = {OUTPUT_FIGURE};
end

if all(cellfun(@isempty, OUTPUT_FIGURE)) == true
    % Set default values if necessary
    OUTPUT_FIGURE = {'', 'POINTS', 'SPLIT'};
end

%% Process OUTPUT_STRENGTH
if iscell(OUTPUT_STRENGTH) == false
    OUTPUT_STRENGTH = {OUTPUT_STRENGTH};
end

if all(cellfun(@isempty, OUTPUT_STRENGTH)) == true
    % Set default values if necessary
    OUTPUT_STRENGTH = {false, 'RESERVE'};
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

if all(cellfun(@isempty, OUTPUT_LOCATION)) == true
    % Set default values if necessary
    OUTPUT_LOCATION = {'DEFAULT', true};
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

%% Get the job ID from the input structure
% Append the job id to the job settings structure
job_id = simon.hash.DataHash(settings, 'SHA-512', 'array');
settings.jobid = job_id;

% Append the job date to the job settings structure
jobDate = char(datetime('now'));
jobDate(ismember(jobDate, ':')) = [];
settings.jobdate = deal(jobDate);

if exist([OUTPUT_LOCATION{1.0}, filesep, jobName], 'dir') == 7.0
    % Get the previous job ID (if applicable)
    [job_id_previous, job_date_previous] = abd.internal_getPreviousJobID(OUTPUT_LOCATION{1.0}, jobName);

    % Set the default check string
    checkString = 'default';
    
    if strcmp(job_id, job_id_previous) == true
        %{
            Job settings have not been modified since the last submission.
            Ask the user if it's OK to overwrite the previous results
        %}
        [response, tf] = uigetpref('latprefdialogues', 'checkOverwrite_resultsNotModified', 'Layup Analysis Tool', sprintf(['The settings for job ''%s'' have n',...
            'ot been modified since\nthe last submission.\n\nOK to overwrite?'], jobName), ["OK", "Cancel"], "DefaultButton", "Cancel");
    else
        %{
            Job settings have been modified since the last submission. Ask
            the user if they want to keep the previous results or overwrite
            the previous results
        %}
        [response, tf] = uigetpref('latprefdialogues', 'checkOverwrite_resultsAlreadyExist', 'Layup Analysis Tool', sprintf(['An output directory already exist',...
            's for job ''%s''.\n\nOK to overwrite?'], jobName), ["Overwrite previous", "Keep previous", "Cancel"], "DefaultButton", "Cancel");
    end
    
    if tf == false
        %{
            The dialogue was never displayed due to the user preference, so
            overwrite if the job settings are the same, and keep previous
            results if the job settings are different
        %}
        if strcmp(response, 'overwrite previous') == true
            response = 'keep previous';
        end

        %{
            The analysis will not be aborted, so match the response to the
            check string
        %}
        checkString = response;
    elseif (strcmp(response, 'cancel') == false) && (isempty(response) == false)
        %{
            The user selected an option that does not abort the analysis,
            so match the response to the check string
        %}
        checkString = response;
    elseif isempty(response) == true
        % There was no user response, so the analysis will be aborted
        response = 'cancel';
    end

    if (strcmp(response(1.0), 'k') == true) && (strcmp(job_id, job_id_previous) == false)
        % Temporarily disable all warnings
        warning('off','all')

        try
            % Old job directory name
            oldJobDir = [pwd, filesep, 'output', filesep, jobName];

            % New job directory name
            newJobDir = [pwd, filesep, 'output', filesep, sprintf('%s_%s', jobName, job_date_previous)];

            % Rename the old job directory to the new name
            movefile(oldJobDir, newJobDir)

            % Add the old job directory back to the MATLAB path
            addpath(genpath(newJobDir))
        catch

            % The old job directory could not be renamed
            response = questdlg('The old job directory could not be renamed. Results will be overwritten. OK to continue?', 'Quick Fatigue Tool', 'Yes', 'No', 'Yes');

            if (strcmp(response, 'Yes') == true) || (strcmpi(response, 'Y') == true)
                checkString = response;
            end
        end

        % Re-enable all warnings
        warning('on','all')
    end

    if strcmpi(checkString, response) == false
        % The user aborted the analysis
        fprintf('[ERROR] Job ''%s'' was aborted by the user\n', jobName);
        error = 1.0;
        return
    end
end
end