function [] = internal_resetDialogues()
%   Function to reset user preference dialogues.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 4.2.0 Copyright Louis Vallance 2025
%   Last modified 10-Jun-2025 08:28:19 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%
%% First delete all QFT user preferences
if isempty(getpref) == false
    % Get the field names of existing user preferences
    names = fieldnames(getpref);

    % Filter field names by QFT user preferences only
    groupsToRemove = names(strcmp(names, 'latprefdialogues'));

    % Remove the user preference
    if isempty(groupsToRemove) == false
        rmpref(char(groupsToRemove))
    end
end
end