function [BUFFER] = internal_updateCrtBuffer(BUFFER, KEYS, VALUES)
%   Update the overall worst criterion buffer.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 5.1.3 Copyright Louis Vallance 2026
%   Last modified 17-Feb-2026 06:40:45 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%

% Update the buffer iteratively
for i = 1:length(KEYS)
    BUFFER(KEYS{i}) = VALUES(i);
end
end