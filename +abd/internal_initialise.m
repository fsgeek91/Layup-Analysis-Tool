function [enableTensor, printTensor, material, failstress, failstrain, hashin, larc05, stackingsequence, plythickness, symmetriclayup, sectionpoints, outputply, outputfigure,...
    outputstrength, outputoptimised, optimisersettings, outputlocation, Nxx, Nyy, Nxy, Mxx, Myy, Mxy, loadtherm, loadmoist, jobname, jobdescription, settings, error] =...
    internal_initialise(settings)
%   Gather variables from user inputs.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 4.2.0 Copyright Louis Vallance 2025
%   Last modified 10-Jun-2025 08:28:19 UTC
%
%#ok<*NODEF>

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%
%% Initialise output
stack = dbstack;
jobname = stack(end).name;
jobdescription = [];
material = [];
failstress = [];
failstrain = [];
hashin = [];
larc05 = [];
stackingsequence = [];
plythickness = [];
symmetriclayup = [];
sectionpoints = [];
Nxx = 0.0;
Nyy = 0.0;
Nxy = 0.0;
Mxx = 0.0;
Myy = 0.0;
Mxy = 0.0;
loadtherm = [];
loadmoist = [];
outputply = [];
outputfigure = [];
outputstrength = [];
outputoptimised = [];
optimisersettings = [];
outputlocation = [];
enableTensor = true;
printTensor = true;
error = false;

% Set the default settings structure
defaultSettings = struct('jobname', jobname, 'jobdescription', [], 'material', {{[]}}, 'failstress', {{[]}}, 'failstrain', {{[]}}, 'hashin', {{[]}}, 'larc05', {{[]}},...
    'stackingsequence', [0.0, 45.0, 90.0], 'plythickness', 0.1, 'symmetriclayup', false, 'sectionpoints', 'DEFAULT', 'loadmech', [0.0, 0.0, 0.0; 0.0, 0.0, 0.0], 'loadtherm', 0.0,...
    'loadmoist', 0.0, 'outputply', 'DEFAULT', 'outputfigure', {{[], 'POINTS', 'SPLIT'}}, 'outputstrength', {{false, 'RESERVE'}}, 'outputoptimised',...
    {{'', 'RESERVE', 'MINMAX', 5.0}}, 'optimisersettings', {{'MIXED-RADIX', 'DEFAULT', 'DEFAULT'}}, 'outputlocation', {{'DEFAULT', true}});

% Get the field name list
settingsNames =  fieldnames(defaultSettings);

% Set fields which are compulsory
compulsorySettings = [false, false, true, false, false, false, false, true, true, false, false, false, false, false, false, false, false, false, false, false];

% Check in setting in turn
for i = 1:length(settingsNames)
    % Get the current field name
    currentField = settingsNames{i};

    % Check that the current field exists
    if isfield(settings, currentField) == false
        % The current field was not specified by the user, so give it a
        if compulsorySettings(i) == true
            % The current setting is compulsory, so exit with an error
            fprintf('[ERROR] Missing option ''%s''\n-> A value for this option is required\n', upper(currentField));

            error = 1.0;
            return
        else
            % The current setting is optional, so set a default value
            dummy = evalc (sprintf('%s = defaultSettings.(settingsNames{i})', currentField)); %#ok<NASGU>

            % Update the settings structure
            settings.(currentField) = defaultSettings.(settingsNames{i});
        end
    else
        dummy = evalc(sprintf('%s = settings.(settingsNames{i})', currentField)); %#ok<NASGU>
    end
end

%% Process NXX/NYY/NXY/MXX/MYY/MXY
% Extract inputs from SETTINGS structure
if isempty(loadmech) == false
    Nxx = loadmech(1.0, 1.0);
    Nyy = loadmech(1.0, 2.0);
    Nxy = loadmech(1.0, 3.0);
    Mxx = loadmech(2.0, 1.0);
    Myy = loadmech(2.0, 2.0);
    Mxy = loadmech(2.0, 3.0);
end

% Disable tensor output if applicable
if (all(all(loadmech == 0.0)) == true) && (all(all([loadtherm, loadmoist] == 0.0)) == true)
    enableTensor = false;
end

%% Process OUTPUT_FIGURE
if iscell(outputfigure) == false
    outputfigure = {outputfigure};
end

if all(cellfun(@isempty, outputfigure)) == true
    % Set default values if necessary
    outputfigure = {'', 'POINTS', 'SPLIT'};
end

%% Process OUTPUT_STRENGTH
if iscell(outputstrength) == false
    outputstrength = {outputstrength};
end

if all(cellfun(@isempty, outputstrength)) == true
    % Set default values if necessary
    outputstrength = {false, 'RESERVE'};
end

%% Process OUTPUT_OPTIMISED
if iscell(outputoptimised) == false
    outputoptimised = {outputoptimised};
end

if all(cellfun(@isempty, outputoptimised)) == true
    % Set default values if necessary
    outputoptimised = {'', 'RESERVE', 'MINMAX', 10.0};
end

%% Process OPTIMISER_SETTINGS
if iscell(optimisersettings) == false
    optimisersettings = {optimisersettings};
end

if all((cellfun(@isempty, optimisersettings) == true)) || (length(optimisersettings) < 3.0)
    % Set default values if necessary
    optimisersettings = {'MIXED-RADIX', 'DEFAULT', 'DEFAULT'};
end

% Process the first argument
argument = optimisersettings{1.0};

if ischar(argument) == false
    optimisersettings{1.0} = 2.0;
else
    argument = lower(argument);
    argument(ismember(argument, ' ')) = [];
    argument(ismember(argument, '-')) = [];

    switch argument
        case 'fullmatrix'
            optimisersettings{1.0} = 1.0;
        case 'mixedradix'
            optimisersettings{1.0} = 2.0;
        case 'chunks'
            optimisersettings{1.0} = 3.0;
        otherwise
            optimisersettings{1.0} = 2.0;
    end
end

%% Process output location
if iscell(outputlocation) == false
    outputlocation = {outputlocation};
end

if all(cellfun(@isempty, outputlocation)) == true
    % Set default values if necessary
    outputlocation = {'DEFAULT', true};
end

if length(outputlocation) < 2.0
    outputlocation{2.0} = false;
end

% Process the first argument
argument = outputlocation{1.0};

if ischar(argument) == false
    outputlocation{1.0} = 'DEFAULT';
else
    if strcmpi(argument, 'default') == true
        % Save results under OUTPUT folder inside user's PWD
        outputlocation{1.0} = [pwd, filesep, 'output'];
    elseif strcmpi(argument, 'qft') == true
        % Save results under QFT folder structure inside user's PWD
        outputlocation{1.0} = [pwd, filesep, 'Project', filesep, 'output'];
    end
end

% Process the second argument
argument = outputlocation{2.0};

if islogical(argument) == false
    outputlocation{2.0} = false;
end

%% Get the job ID from the input structure
% Append the job id to the job settings structure
job_id = simon.hash.DataHash(settings, 'SHA-512', 'array');
settings.jobid = job_id;

% Append the job date to the job settings structure
jobDate = char(datetime('now'));
jobDate(ismember(jobDate, ':')) = [];
settings.jobdate = deal(jobDate);

% Get the complete path to the output directory
outputlocationFull = [outputlocation{1.0}, filesep, jobname];

if exist(outputlocationFull, 'dir') == 7.0
    % Get the previous job ID (if applicable)
    [job_id_previous, job_date_previous] = abd.internal_getPreviousJobID(outputlocation{1.0}, jobname);

    % Set the default check string
    checkString = 'default';
    
    if strcmp(job_id, job_id_previous) == true
        %{
            Job settings have not been modified since the last submission.
            Ask the user if it's OK to overwrite the previous results
        %}
        [response, tf] = uigetpref('latprefdialogues', 'checkOverwrite_resultsNotModified', 'Layup Analysis Tool', sprintf(['The settings for job ''%s'' have n',...
            'ot been modified since\nthe last submission.\n\nOK to overwrite?'], jobname), ["OK", "Cancel"], "DefaultButton", "Cancel");
    else
        %{
            Job settings have been modified since the last submission. Ask
            the user if they want to keep the previous results or overwrite
            the previous results
        %}
        [response, tf] = uigetpref('latprefdialogues', 'checkOverwrite_resultsAlreadyExist', 'Layup Analysis Tool', sprintf(['An output directory already exist',...
            's for job ''%s''.\n\nOK to overwrite?'], jobname), ["Overwrite previous", "Keep previous", "Cancel"], "DefaultButton", "Cancel");
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

    % Temporarily disable all warnings
    warning('off','all')

    if (strcmp(response(1.0), 'k') == true) && (strcmp(job_id, job_id_previous) == false)
        try
            % Old job directory name
            oldJobDir = [pwd, filesep, 'output', filesep, jobname];

            % New job directory name
            newJobDir = [pwd, filesep, 'output', filesep, sprintf('%s_%s', jobname, job_date_previous)];

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
    end

    if strcmpi(checkString, response) == true
        % Remove the contents of the existing directory
        try
            if exist(outputlocationFull, 'dir') == 7.0
                rmdir(outputlocationFull, 's');
            end
        catch 
            % Do nothing
        end
    else
        % Re-enable all warnings
        warning('on','all')

        % The user aborted the analysis
        fprintf('[ERROR] Job ''%s'' was aborted by the user\n', jobname);
        error = 1.0;
        return
    end

    % Re-enable all warnings
    warning('on','all')
end
end