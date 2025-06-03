function [] = internal_printTensor(fid, OUTPUT_ENVELOPE, ENVELOPE_MODE, S_ply_xy, S_ply_aligned, E_ply_xy, E_ply_aligned, E_therm_xy, E_moist_xy, E_therm_aligned, E_moist_aligned,...
    nPlies, outputPoints, plyBuffer, symmetricAbd, outputApproximate, thickness, header, OUTPUT_PLY)
%   Print stress/strain tensor information to the output file.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 3.0.7 Copyright Louis Vallance 2025
%   Last modified 03-Jun-2025 10:08:33 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%
%% Get the envelope mode in case it's needed later
switch ENVELOPE_MODE
    case 1.0
        envelopeString = 'ENVELOPE ABSMAX';
    case 2.0
        envelopeString = 'ENVELOPE MAX';
    case 3.0
        envelopeString = 'ENVELOPE MIN';
end

%% Set the location string for each section points
if ischar(OUTPUT_PLY) == true
    switch lower(OUTPUT_PLY)
        case 'default'
            % Top and bottom faces
            locationString = {'BOTTOM', 'TOP'};
        case 'top'
            % Top face only
            locationString = {'TOP'};
        case 'middle'
            % Midspan only
            if isempty(outputApproximate) == false
                locationString = {'MIDDLE (APPROX.)'};
            else
                locationString = {'MIDDLE (PRECISE)'};
            end
        case 'bottom'
            % Bottom face only
            locationString = {'BOTTOM'};
        otherwise
            % Unknown
            locationString = repmat({'UNKNOWN'}, 1.0, length(outputPoints));
    end
else
    % User-specified section point list
    locationString = repmat({'USER'}, 1.0, length(outputPoints));
end

%% Inform use if output is incomplete
if length(unique(plyBuffer)) ~= nPlies
    % Inform the user of incomplete output
    fprintf(fid, '\nNote: Result output is not available at all plies\n');
elseif isempty(outputApproximate) == false
    fprintf(fid, outputApproximate);
end

%% Print stress tensor data
% Print the table header
fprintf(fid, '\n%s\n', header);

% Print the stress tensor table for each ply location
fprintf(fid, '\nStress at user-selected locations:\n');

if OUTPUT_ENVELOPE == true
    % Print the header
    fprintf(fid, 'PLY    SECTION POINT                 Sxx           Syy           Sxy           S11           S22           S12           \n');

    % Print the numerically largest stress in the current ply
    S_ply_xy_envelope_11 = abd.internal_getAbsMax(S_ply_xy(1.0, :), ENVELOPE_MODE);
    S_ply_xy_envelope_22 = abd.internal_getAbsMax(S_ply_xy(2.0, :), ENVELOPE_MODE);
    S_ply_xy_envelope_12 = abd.internal_getAbsMax(S_ply_xy(3.0, :), ENVELOPE_MODE);
    S_ply_aligned_envelope_11 = abd.internal_getAbsMax(S_ply_aligned(1.0, :), ENVELOPE_MODE);
    S_ply_aligned_envelope_22 = abd.internal_getAbsMax(S_ply_aligned(2.0, :), ENVELOPE_MODE);
    S_ply_aligned_envelope_12 = abd.internal_getAbsMax(S_ply_aligned(3.0, :), ENVELOPE_MODE);

    fprintf(fid, '%-7s%-30s%-14g%-14g%-14g%-14g%-14g%-14g\n', 'ALL', envelopeString, S_ply_xy_envelope_11, S_ply_xy_envelope_22, S_ply_xy_envelope_12, S_ply_aligned_envelope_11,...
        S_ply_aligned_envelope_22, S_ply_aligned_envelope_12);
else
    % Print the header
    fprintf(fid, 'PLY    SECTION POINT                 LOCATION            Sxx           Syy           Sxy           S11           S22           S12           \n');

    for i = 1.0:nPlies
        % Extract the stresses at the section points for the current ply
        currentOutputPoints = outputPoints(plyBuffer == i);

        if isempty(currentOutputPoints) == true
            fprintf(fid, '%-7.0f%s\n', i, 'NO RESULTS');
        else
            % Get the stresses for the current section points
            S_ply_xy_i = S_ply_xy(:, currentOutputPoints);
            S_ply_aligned_i = S_ply_aligned(:, currentOutputPoints);

            for sp = 1.0:length(currentOutputPoints)
                % Get the thickness fraction of the current section point
                thicknessFraction = sprintf('%.0f (fraction = %.5g)', currentOutputPoints(sp), thickness(currentOutputPoints(sp)));

                fprintf(fid, '%-7.0f%-30s%-20s%-14g%-14g%-14g%-14g%-14g%-14g\n', i, thicknessFraction, locationString{sp}, S_ply_xy_i(1.0, sp), S_ply_xy_i(2.0, sp), S_ply_xy_i(3.0, sp),...
                    S_ply_aligned_i(1.0, sp), S_ply_aligned_i(2.0, sp), S_ply_aligned_i(3.0, sp));
            end
        end

        % Draw symmetry plane (if applicable)
        if (symmetricAbd == true) && (i == nPlies/2.0)
            fprintf(fid, '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - SYM\n');
        end
    end
end

%% Print strain tensor data
fprintf(fid, '\nStrain at user-selected locations:\n');

if OUTPUT_ENVELOPE == true
    % Print the header
    fprintf(fid, 'PLY    SECTION POINT                 Exx           Eyy           Exy           E11           E22           E12           \n');

    % Print the numerically largest stress over all plies
    E_ply_xy_envelope_11 = abd.internal_getAbsMax(E_ply_xy(1.0, :), ENVELOPE_MODE);
    E_ply_xy_envelope_22 = abd.internal_getAbsMax(E_ply_xy(2.0, :), ENVELOPE_MODE);
    E_ply_xy_envelope_12 = abd.internal_getAbsMax(E_ply_xy(3.0, :), ENVELOPE_MODE);
    E_ply_aligned_envelope_11 = abd.internal_getAbsMax(E_ply_aligned(1.0, :), ENVELOPE_MODE);
    E_ply_aligned_envelope_22 = abd.internal_getAbsMax(E_ply_aligned(2.0, :), ENVELOPE_MODE);
    E_ply_aligned_envelope_12 = abd.internal_getAbsMax(E_ply_aligned(3.0, :), ENVELOPE_MODE);

    fprintf(fid, '%-7s%-30s%-14g%-14g%-14g%-14g%-14g%-14g\n', 'ALL', envelopeString, E_ply_xy_envelope_11, E_ply_xy_envelope_22, E_ply_xy_envelope_12, E_ply_aligned_envelope_11,...
        E_ply_aligned_envelope_22, E_ply_aligned_envelope_12);
else
    % Print the header
    fprintf(fid, 'PLY    SECTION POINT                 LOCATION            Exx           Eyy           Exy           E11           E22           E12           \n');

    for i = 1.0:nPlies
        % Extract the strains at the section points for the current ply
        currentOutputPoints = outputPoints(plyBuffer == i);

        if isempty(currentOutputPoints) == true
            fprintf(fid, '%-7.0f%s\n', i, 'NO RESULTS');
        else
            % Get the strains for the current section points
            E_ply_xy_i = E_ply_xy(:, currentOutputPoints);
            E_ply_aligned_i = E_ply_aligned(:, currentOutputPoints);

            for sp = 1.0:length(currentOutputPoints)
                thicknessFraction = sprintf('%.0f (fraction = %.5g)', currentOutputPoints(sp), thickness(currentOutputPoints(sp)));

                fprintf(fid, '%-7.0f%-30s%-20s%-14g%-14g%-14g%-14g%-14g%-14g\n', i, thicknessFraction, locationString{sp}, E_ply_xy_i(1.0, sp), E_ply_xy_i(2.0, sp),...
                    E_ply_xy_i(3.0, sp), E_ply_aligned_i(1.0, sp), E_ply_aligned_i(2.0, sp), E_ply_aligned_i(3.0, sp));
            end
        end

        % Draw symmetry plane (if applicable)
        if (symmetricAbd == true) && (i == nPlies/2.0)
            fprintf(fid, '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - SYM\n');
        end
    end
end

%% Print the stress-free thermal strain tensor data
if isempty(E_therm_xy) == false && (any(any(E_therm_xy)) == true || any(any(E_therm_aligned)) == true)
    fprintf(fid, '\nStress-free strain due to thermal process at user-selected locations:\n');
    fprintf(fid, 'PLY    SECTION POINT                 E_Therm_xx    E_Therm_yy    E_Therm_xy    E_Therm_11    E_Therm_22    E_Therm_12    \n');

    if OUTPUT_ENVELOPE == true
        % Print the numerically largest stress over all plies
        E_therm_xy_envelope_11 = abd.internal_getAbsMax(E_therm_xy(1.0, :), ENVELOPE_MODE);
        E_therm_xy_envelope_22 = abd.internal_getAbsMax(E_therm_xy(2.0, :), ENVELOPE_MODE);
        E_therm_xy_envelope_12 = abd.internal_getAbsMax(E_therm_xy(3.0, :), ENVELOPE_MODE);
        E_therm_aligned_envelope_11 = abd.internal_getAbsMax(E_therm_aligned(1.0, :), ENVELOPE_MODE);
        E_therm_aligned_envelope_22 = abd.internal_getAbsMax(E_therm_aligned(2.0, :), ENVELOPE_MODE);
        E_therm_aligned_envelope_12 = abd.internal_getAbsMax(E_therm_aligned(3.0, :), ENVELOPE_MODE);

        fprintf(fid, '%-7s%-30s%-14g%-14g%-14g%-14g%-14g%-14g\n', 'ALL', envelopeString, E_therm_xy_envelope_11, E_therm_xy_envelope_22, E_therm_xy_envelope_12,...
            E_therm_aligned_envelope_11, E_therm_aligned_envelope_22, E_therm_aligned_envelope_12);
    else
        for i = 1.0:nPlies
            %{
                Extract the strains for the current section point for the
                current ply
            %}
            currentOutputPoints = outputPoints(plyBuffer == i);

            if isempty(currentOutputPoints) == true
                fprintf(fid, '%-7.0f%s\n', i, 'NO RESULTS');
            else
                % Get the strains for the current section points
                E_therm_xy_i = E_therm_xy(:, currentOutputPoints);
                E_therm_aligned_i = E_therm_aligned(:, currentOutputPoints);

                for sp = 1.0:length(currentOutputPoints)
                    thicknessFraction = sprintf('%.0f (fraction = %.5g)', currentOutputPoints(sp), thickness(currentOutputPoints(sp)));

                    fprintf(fid, '%-7.0f%-30s%-14g%-14g%-14g%-14g%-14g%-14g\n', i, thicknessFraction, E_therm_xy_i(1.0, sp), E_therm_xy_i(2.0, sp), E_therm_xy_i(3.0, sp),...
                        E_therm_aligned_i(1.0, sp), E_therm_aligned_i(2.0, sp), E_therm_aligned_i(3.0, sp));
                end
            end

            % Draw symmetry plane (if applicable)
            if (symmetricAbd == true) && (i == nPlies/2.0)
                fprintf(fid, '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - SYM\n');
            end
        end
    end
end

%% Print the stress-free moisture strain tensor data
if isempty(E_moist_xy) == false && (any(any(E_moist_xy)) == true || any(any(E_moist_aligned)) == true)
    fprintf(fid, '\nStress-free strain due to moisture process at user-selected locations:\n');
    fprintf(fid, 'PLY    SECTION POINT                 E_Moist_xx    E_Moist_yy    E_Moist_xy    E_Moist_11    E_Moist_22    E_Moist_12    \n');

    if OUTPUT_ENVELOPE == true
        % Print the numerically largest stress over all plies
        E_moist_xy_envelope_11 = abd.internal_getAbsMax(E_moist_xy(1.0, :), ENVELOPE_MODE);
        E_moist_xy_envelope_22 = abd.internal_getAbsMax(E_moist_xy(2.0, :), ENVELOPE_MODE);
        E_moist_xy_envelope_12 = abd.internal_getAbsMax(E_moist_xy(3.0, :), ENVELOPE_MODE);
        E_moist_aligned_envelope_11 = abd.internal_getAbsMax(E_moist_aligned(1.0, :), ENVELOPE_MODE);
        E_moist_aligned_envelope_22 = abd.internal_getAbsMax(E_moist_aligned(2.0, :), ENVELOPE_MODE);
        E_moist_aligned_envelope_12 = abd.internal_getAbsMax(E_moist_aligned(3.0, :), ENVELOPE_MODE);

        fprintf(fid, '%-7s%-30s%-14g%-14g%-14g%-14g%-14g%-14g\n', 'ALL', envelopeString, E_moist_xy_envelope_11, E_moist_xy_envelope_22, E_moist_xy_envelope_12,...
            E_moist_aligned_envelope_11, E_moist_aligned_envelope_22, E_moist_aligned_envelope_12);
    else
        for i = 1.0:nPlies
            %{
                Extract the strains for the current section point for the
                current ply
            %}
            currentOutputPoints = outputPoints(plyBuffer == i);

            if isempty(currentOutputPoints) == true
                fprintf(fid, '%-7.0f%s\n', i, 'NO RESULTS');
            else
                % Get the strains for the current section points
                E_moist_xy_i = E_moist_xy(:, currentOutputPoints);
                E_moist_aligned_i = E_moist_aligned(:, currentOutputPoints);

                for sp = 1.0:length(currentOutputPoints)
                    thicknessFraction = sprintf('%.0f (fraction = %.5g)', currentOutputPoints(sp), thickness(currentOutputPoints(sp)));

                    fprintf(fid, '%-7.0f%-30s%-14g%-14g%-14g%-14g%-14g%-14g\n', i, thicknessFraction, E_moist_xy_i(1.0, sp), E_moist_xy_i(2.0, sp), E_moist_xy_i(3.0, sp),...
                        E_moist_aligned_i(1.0, sp), E_moist_aligned_i(2.0, sp), E_moist_aligned_i(3.0, sp));
                end
            end

            % Draw symmetry plane (if applicable)
            if (symmetricAbd == true) && (i == nPlies/2.0)
                fprintf(fid, '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - SYM\n');
            end
        end
    end
end