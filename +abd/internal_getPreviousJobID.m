function [job_id_previous, job_date_previous] = internal_getPreviousJobID(folder, jobName)
%   Get the ID of the previous job.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 5.1.2 Copyright Louis Vallance 2026
%   Last modified 16-Feb-2026 12:06:38 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%
% Initialize output
job_id_previous = '';
job_date_previous = '';

% Search for <JOB_NAME>.MAT inside the JOBNAME directory
previousJobSettings = dir(fullfile([folder, filesep, jobName, filesep, 'settings.mat']));

% If the previous job settings exist, get the job ID
if isempty(previousJobSettings) == false
    try
        % Previous output folder
        folder_previous = previousJobSettings.folder;

        % Reconstruct the full file path
        f = [folder_previous, filesep, 'settings.mat'];

        % Load the job settings data
        x = load(f);

        % Extract the job ID
        job_id_previous = x.settings.jobid;

        % Extract the date of the previous job
        job_date_previous = x.settings.jobdate;
    catch MException
        % Something went wrong while extracting the data
        
        % Save the MATLAB exception object
        abd.internal_saveMExceptionObj(MException)
    end
end
end