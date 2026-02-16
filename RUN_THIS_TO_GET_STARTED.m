% Print the Layup Analysis Tool help contents to the MATLAB Command Window
help lat

% Open the user definitions template file
try 
    edit user_definitions.m
catch MException
    % Save the MATLAB exception object
    abd.internal_saveMExceptionObj(MException)
end