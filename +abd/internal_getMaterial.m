function [error, noStrength, varargout] = internal_getMaterial(data, nPlies, symmetricPly, mode, tag)
%   Get material data for each ply.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 3.0.1 Copyright Louis Vallance 2024
%   Last modified 23-Feb-2024 13:20:04 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%

% Initialise output
error = false;
noStrength = false;

% Process MODE
switch mode
    case 1.0 % MECHANICAL
        varargout{1.0} = [];
        varargout{2.0} = [];
        varargout{3.0} = [];
        varargout{4.0} = [];
        varargout{5.0} = [];
        varargout{6.0} = [];
        varargout{7.0} = [];
        varargout{8.0} = [];
    case 2.0 % FAIL STRESS OR HASHIN
        varargout{1.0} = [];
        varargout{2.0} = [];
        varargout{3.0} = [];
        varargout{4.0} = [];
        varargout{5.0} = [];
        varargout{6.0} = [];
        varargout{7.0} = [];
    case 3.0  % FAIL STRAIN
        varargout{1.0} = [];
        varargout{2.0} = [];
        varargout{3.0} = [];
        varargout{4.0} = [];
        varargout{5.0} = [];
    case 4.0 % LARC05
        varargout{1.0} = [];
        varargout{2.0} = [];
        varargout{3.0} = [];
        varargout{4.0} = [];
        varargout{5.0} = [];
        varargout{6.0} = [];
        varargout{7.0} = [];
        varargout{8.0} = [];
        varargout{9.0} = [];
        varargout{10.0} = [];
        varargout{11.0} = [];
    otherwise
end

% Extract nested cell (if applicable)
if iscell(data) == false
    data = {data};
elseif iscell(data{1.0}) == true
    data = data{1.0};
end

% Get the number of materials
nMaterials = length(data);

% Check for a consistent definition
if symmetricPly == false
    if (nMaterials > 1.0) && (nMaterials ~= nPlies)
        %{
            For layups defined without symmetry, the number of materials
            must equal the number of plies if more than one material was
            specified
        %}
        fprintf('[ERROR] The number of material definitions does not\nmatch the number of plies\n');

        % Reset the error flag and RETURN
        error = true;
        return
    elseif nMaterials == 1.0
        % Propagate the material list to the number of plies
        data = repmat(data, [1.0, nPlies]);
    end
elseif (symmetricPly == true) && (nMaterials > 1.0)
    %{
        For layups defined with symmetry and more than one material is
        specified, the number of materials must be equal to the number of
        plies
    %}
    if nMaterials == 0.5*nPlies
        % Mirror the material list
        data = [data, flip(data)];
    elseif nMaterials ~= nPlies
        % The number of specified materials is invalid
        fprintf('[ERROR] The number of materials does not match the\nnumber of plies in the layup definition\n');

        % Reset the error flag and RETURN
        error = true;
        return
    end
elseif symmetricPly == true
    % Propagate the material list to the number of plies
    data = repmat(data, [1.0, nPlies]);
end

if (all(cellfun(@isempty, data)) == true) && (mode ~= 1.0)
    %{
        Even if OUTPUT_STRENGTH = true, it is not compulsory to define
        FAIL_STRESS and FAIL_STRAIN properties. RETURN now and check that
        at least FAIL_STRESS or FAIL_STRAIN properties are specified
    %}

    % Reset the flag and RETURN
    noStrength = true;
    return
elseif any(cellfun(@isempty, data)) == true
    %{
        At least one ply is missing material properties, so RETURN with an
        error
    %}
    fprintf('[ERROR] One or more plies are missing material\nproperties\n');

    % Reset the error flag and RETURN
    error = true;
    return
end

switch mode
    case 1.0
        % Initialise material property buffers
        D1 = zeros(1.0, nPlies); D2 = D1; D3 = D1; D4 = D1; D5 = D1;
        D6 = D1; D7 = D1; D8 = D1;

        % Populate material property buffers
        for i = 1:nPlies
            % Get material properties for the current ply
            currentMaterial = data{i};

            % Property count check
            if length(currentMaterial) ~= 8.0
                fprintf('[ERROR] Incorrect number of properties specified in %s\n', tag);

                % Reset the error flag and RETURN
                error = true;
                return
            end

            % Assign values for the material property buffers
            [D1(i), D2(i), D3(i), D4(i), D5(i), D6(i), D7(i), D8(i)] = deal(currentMaterial(1.0), currentMaterial(2.0), currentMaterial(3.0), currentMaterial(4.0),...
                currentMaterial(5.0), currentMaterial(6.0), currentMaterial(7.0), currentMaterial(8.0));
        end

        % Assign values to VARARGOUT
        varargout{1.0} = D1;
        varargout{2.0} = D2;
        varargout{3.0} = D3;
        varargout{4.0} = D4;
        varargout{5.0} = D5;
        varargout{6.0} = D6;
        varargout{7.0} = D7;
        varargout{8.0} = D8;
    case 2.0
        % Initialise material property buffers
        D1 = zeros(1.0, nPlies); D2 = D1; D3 = D1; D4 = D1; D5 = D1;
        D6 = D1; D7 = D1;

        % Populate material property buffers
        for i = 1:nPlies
            % Get material properties for the current ply
            currentMaterial = data{i};

            % Property count check
            if length(currentMaterial) ~= 7.0
                fprintf('[ERROR] Incorrect number of properties specified in %s\n', tag);

                % Reset the error flag and RETURN
                error = true;
                return
            end

            % Assign values for the material property buffers
            [D1(i), D2(i), D3(i), D4(i), D5(i), D6(i), D7(i)] = deal(currentMaterial(1.0), currentMaterial(2.0), currentMaterial(3.0), currentMaterial(4.0), currentMaterial(5.0),...
                currentMaterial(6.0), currentMaterial(7.0));
        end

        % Assign values to VARARGOUT
        varargout{1.0} = D1;
        varargout{2.0} = D2;
        varargout{3.0} = D3;
        varargout{4.0} = D4;
        varargout{5.0} = D5;
        varargout{6.0} = D6;
        varargout{7.0} = D7;
    case 3.0
        % Initialise material property buffers
        D1 = zeros(1.0, nPlies); D2 = D1; D3 = D1; D4 = D1; D5 = D1;

        % Populate material property buffers
        for i = 1:nPlies
            % Get material properties for the current ply
            currentMaterial = data{i};

            % Property count check
            if length(currentMaterial) ~= 5.0
                fprintf('[ERROR] Incorrect number of properties specified in %s\n', tag);

                % Reset the error flag and RETURN
                error = true;
                return
            end

            % Assign values for the material property buffers
            [D1(i), D2(i), D3(i), D4(i), D5(i)] = deal(currentMaterial(1.0), currentMaterial(2.0), currentMaterial(3.0), currentMaterial(4.0), currentMaterial(5.0));
        end

        % Assign values to VARARGOUT
        varargout{1.0} = D1;
        varargout{2.0} = D2;
        varargout{3.0} = D3;
        varargout{4.0} = D4;
        varargout{5.0} = D5;
    case 4.0
        % Initialise material property buffers
        D1 = zeros(1.0, nPlies); D2 = D1; D3 = D1; D4 = D1; D5 = D1;
        D6 = D1; D7 = D1; D8 = D1; D9 = D1; D10 = D1; D11 = D1;

        % Populate material property buffers
        for i = 1:nPlies
            % Get material properties for the current ply
            currentMaterial = data{i};

            % Property count check
            if length(currentMaterial) ~= 11.0
                fprintf('[ERROR] Incorrect number of properties specified in %s\n', tag);

                % Reset the error flag and RETURN
                error = true;
                return
            end

            % Assign values for the material property buffers
            [D1(i), D2(i), D3(i), D4(i), D5(i), D6(i), D7(i), D8(i), D9(i), D10(i), D11(i)] = deal(currentMaterial(1.0), currentMaterial(2.0), currentMaterial(3.0),...
                currentMaterial(4.0), currentMaterial(5.0), currentMaterial(6.0), currentMaterial(7.0), currentMaterial(8.0), currentMaterial(9.0), currentMaterial(10.0),...
                currentMaterial(11.0));
        end

        % Replace undefined values of D4 with D3
        D4(D4 == -1.0) = D3(D4 == -1.0);

        % Replace undefined values of D10 with default value
        D10(D10 == -1.0) = 53.0;

        % Assign values to VARARGOUT
        varargout{1.0} = D1;
        varargout{2.0} = D2;
        varargout{3.0} = D3;
        varargout{4.0} = D4;
        varargout{5.0} = D5;
        varargout{6.0} = D6;
        varargout{7.0} = D7;
        varargout{8.0} = D8;
        varargout{9.0} = D9;
        varargout{10.0} = D10;
        varargout{11.0} = D11;
    otherwise
end
end