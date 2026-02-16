function [] = internal_saveMExceptionObj(MExcpetion)
% Get the user output folder
outputDir = fullfile(pwd, 'output');  % Change if needed
if exist(outputDir, 'dir') == false
    mkdir(outputDir);
end

% Format the time stamp
timestamp = strrep(char(datetime), ':', '_');

% Get the function name
stack = dbstack;
if numel(stack) > 1.0
    funcName = stack(2.0).name;
else
    funcName = mfilename;
end

% Get the exception identifier (shortened)
if isempty(MExcpetion.identifier) == false
    idPart = strrep(MExcpetion.identifier, ':', '_');
else
    idPart = 'NoIdentifier';
end

% Shorten long identifiers
idPart = idPart(1.0:min(end, 40.0));

% Build the file name
fileName = sprintf('Exception_%s_%s_%s.mat', funcName, idPart, timestamp);

% Remove illegal filename characters
fileName = regexprep(fileName, '[^\w\-\.]', '_');

fullPath = fullfile(outputDir, fileName);

% Save the exception object to the output folder
save(fullPath, 'MExcpetion');

% Inform the user
fprintf('[ERROR] Layup Analysis Tool encountered an exception\n-> The exception has been saved to:\n%s\n-> Please contact the developer for assistance\n', fullPath);