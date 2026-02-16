help lat

try 
    edit user_definitions.m
catch MException
    % Save the MATLAB exception object
    abd.internal_saveMExceptionObj(MException)
end