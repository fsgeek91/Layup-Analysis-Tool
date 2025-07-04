function [t_ply, theta, nPlies, error] = internal_mirror(symmetricPly, t_ply, theta)
%   Mirror the ply stacking list (if applicable).
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 4.2.3 Copyright Louis Vallance 2025
%   Last modified 23-Jun-2025 14:28:39 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%
% Initialise output
error = false;
nPlies = [];

% Apply symmetry
if symmetricPly == false
    if (length(t_ply) > 1.0) && (length(t_ply) ~= length(theta))
        %{
            For layups defined without symmetry, the number of thickness
            values must equal the number of plies if more than one
            thickness value was specified
        %}
        fprintf('[ERROR] The number of ply thickness values does not match the number of\nplies in the layup definition\n');

        % Reset the error flag and RETURN
        error = true;
        return
    elseif isscalar(t_ply) == true
        % Propagate the thickness list to the number of plies
        t_ply = linspace(t_ply, t_ply, length(theta));
    end
elseif (symmetricPly == true) && (length(t_ply) > 1.0)
    %{
        For layups defined with symmetry and more than one thickness value
        is specified, the number of thickness values must either be twice
        the number of plies (before mirroring), or equal to the total
        number of plies after accounting for symmetry
    %}
    if length(t_ply) == 2.0*length(theta)
        % There is no need to mirror the ply thickness list
    elseif length(t_ply) ~= length(theta) 
        % The number of specified thickness values is invalid
        fprintf('[ERROR] The number of ply thickness values does not match the number of\nplies in the layup definition\n');

        % Reset the error flag and RETURN
        error = true;
        return
    else
        % Mirror the ply thickness values
        t_ply = [t_ply, flip(t_ply)];
    end

    % Mirror the layup stacking list
    theta = [theta, flip(theta)];
elseif symmetricPly == true
    % Mirror the layup stacking list and the ply thickness list
    theta = [theta, flip(theta)];
    t_ply = linspace(t_ply, t_ply, length(theta));
end

% Number of plies and the number of material definitions
nPlies = length(theta);