function [job_id_previous, job_date_previous] = internal_getPreviousJobID(folder, jobName)
%   Get the ID of the previous job.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 4.3.0 Copyright Louis Vallance 2025
%   Last modified 17-Jun-2025 08:14:14 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%
% Initialize output
job_id_previous = '';
job_date_previous = '';

% Search for <JOB_NAME>.MAT inside the JOBNAME directory
previousJobSettings = dir(fullfile([folder, filesep, jobName, filesep, jobName '.mat']));

% If the previous job settings exist, get the job ID
if isempty(previousJobSettings) == false
    try
        % Previous output folder
        folder_previous = previousJobSettings.folder;

        % Reconstruct the full file path
        f = [folder_previous, filesep, jobName, '.mat'];

        % Load the job settings data
        x = load(f);

        % Extract the job ID
        job_id_previous = x.settings.jobid;

        % Extract the date of the previous job
        job_date_previous = x.settings.jobdate;
    catch
        % Something went wrong while extracting the data
    end
end
end