function [error, noStrength, varargout] = internal_getMaterial(data, nPlies, symmetricPly, mode, tag)
%   Get material data for each ply.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 5.1.0 Copyright Louis Vallance 2026
%   Last modified 12-Feb-2026 12:33:07 UTC
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
        fprintf('[ERROR] The number of material definitions does not match the number of\nplies\n');

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
        fprintf('[ERROR] The number of materials does not match the number of plies in\nthe layup definition\n');

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
    fprintf('[ERROR] One or more plies are missing material properties\n');

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
                fprintf('[ERROR] Incorrect number of properties specified in %s\n-> Expected 8, got %.0f\n', tag, length(currentMaterial));

                % Reset the error flag and RETURN
                error = true;
            end

            % Validity check
            if any(currentMaterial <= 0.0)
                fprintf('[ERROR] In %s, mechanical properties must be positive\n', tag);

                % Reset the error flag and RETURN
                error = true;
            end

            % If error, RETURN
            if error == true
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

            % Argument number check
            if length(currentMaterial) ~= 7.0
                fprintf('[ERROR] Incorrect number of properties specified in %s\n-> Expected 7, got %.0f\n', tag, length(currentMaterial));

                % Reset the error flag and RETURN
                error = true;
            end

            % Validity check
            switch lower(tag)
                case 'fail_stress'
                    if (length(currentMaterial) > 4.0) && (any(currentMaterial(1.0:5.0) <= 0.0) == true)
                        fprintf('[ERROR] In %s, strength properties must be positive\n', tag);

                        % Reset the error flag and RETURN
                        error = true;
                    elseif (length(currentMaterial) > 5.0) && ((currentMaterial(6.0) < -1.0) || (currentMaterial(6.0) > 1.0))
                        fprintf('[ERROR] In %s, stress coupling term must be in the range {-1 <= C <= 1}\n', tag);

                        % Reset the error flag and RETURN
                        error = true;
                    elseif (length(currentMaterial) > 6.0) && (currentMaterial(7.0) < 0.0)
                        fprintf('[ERROR] In %s, biaxiality ratio (B) cannot be negative\n', tag);

                        % Reset the error flag and RETURN
                        error = true;
                    end
                case 'hashin'
                    if (length(currentMaterial) > 6.0) && (any(currentMaterial(2.0:7.0) <= 0.0) == true)
                        fprintf('[ERROR] In %s, strength properties must be positive\n', tag);

                        % Reset the error flag and RETURN
                        error = true;
                    elseif (isempty(currentMaterial) == false) && ((currentMaterial(1.0) < 0.0) || (currentMaterial(1.0) > 1.0))
                        fprintf('[ERROR] In %s, coupling term must be in the range {0 <= ALPHA <= 1}\n', tag);

                        % Reset the error flag and RETURN
                        error = true;
                    end
            end

            % If error, RETURN
            if error == true
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

            % Argument number check
            if length(currentMaterial) ~= 5.0
                fprintf('[ERROR] Incorrect number of properties specified in %s\n-> Expeted 5, got %.0f\n', tag, length(currentMaterial));

                % Reset the error flag and RETURN
                error = true;
            end

            % Validity check
            if any(currentMaterial <= 0.0)
                fprintf('[ERROR] In %s, strength properties must be positive\n', tag);

                % Reset the error flag and RETURN
                error = true;
            end

            % If error, RETURN
            if error == true
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

            % Argument number check
            if length(currentMaterial) ~= 11.0
                fprintf('[ERROR] Incorrect number of properties specified in %s\n-> Expected 11, got %.0f\n', tag, length(currentMaterial));

                % Reset the error flag and RETURN
                error = true;
            end

            % Validity check
            if ((length(currentMaterial) > 6.0) && (any(currentMaterial([1.0, 2.0, 3.0, 5.0, 7.0]) <= 0.0))) ||...
                    ((length(currentMaterial) > 5.0) && (any(all([currentMaterial([4.0, 6.0]) <= 0.0; currentMaterial([4.0, 6.0]) ~= -1.0])) == true))
                fprintf('[ERROR] In %s, strength properties must be positive\n', tag);

                % Reset the error flag and RETURN
                error = true;
            elseif (length(currentMaterial) > 8.0) && (any(all([any([(currentMaterial(8.0:9.0) < 0.0); (currentMaterial(8.0:9.0) > 1.0)]); currentMaterial(8.0:9.0) ~= -1.0])) == true)
                fprintf('[ERROR] In %s, longitudinal/transverse shear friction coefficient\nmust be in the range {0 <= NL/T <= 1}\n', tag);

                % Reset the error flag and RETURN
                error = true;
            elseif (length(currentMaterial) > 9.0) && ((currentMaterial(10.0) < 0.0 || currentMaterial(10.0) > 180.0) && (currentMaterial(10.0) ~= -1.0))
                fprintf('[ERROR] In %s, fracture plane angle for pure compression must be in\nthe range {0 <= A0 <= 180}\n', tag);

                % Reset the error flag and RETURN
                error = true;
            elseif (length(currentMaterial) > 10.0) && ((currentMaterial(11.0) <= 0.0) && (currentMaterial(11.0) ~= -1.0))
                fprintf('[ERROR] In %s, misalignment angle at failure for pure compression\n(PHI0) must be positive\n', tag);

                % Reset the error flag and RETURN
                error = true;
            end

            % If error, RETURN
            if error == true
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