function [BUFFER] = internal_updateCrtBuffer(BUFFER, KEYS, VALUES)
%   Update the overall worst criterion buffer.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 5.1.5 Copyright Louis Vallance 2026
%   Last modified 04-Sep-2026 12:21:27 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%

% Update the buffer iteratively
for i = 1:length(KEYS)
    BUFFER(KEYS{i}) = VALUES(i);
end
end