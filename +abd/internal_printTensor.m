function [] = internal_printTensor(fid, OUTPUT_ENVELOPE,...
    ENVELOPE_MODE, S_ply_xy, S_ply_aligned, E_ply_xy, E_ply_aligned,...
    nPlies, outputPoints, plyBuffer, symmetricAbd, outputApproximate,...
    thickness, header)
%   Print stress/strain tensor information to the output file.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 2.2 Copyright Louis Vallance 2023
%   Last modified 09-May-2023 07:31:07 UTC
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

%% Inform use if output is incomplete
if length(unique(plyBuffer)) ~= nPlies
    % Inform the user of incomplete output
    fprintf(fid, ['\nNote: Result output is not available at all plies',...
        '\n']);
elseif outputApproximate == true
    % Inform the user of approximate output
    fprintf(fid, ['\nNote: Output location MIDDLE is unavailble. Resul',...
        't output is at approximate\nlocations\n']);
end

%% Print stress tensor data
% Print the table header
fprintf(fid, '\n%s\n', header);

% Print the stress tensor table for each ply location
fprintf(fid, '\nStress tensors at user-selected locations:\n');
fprintf(fid, ['PLY  SECTION POINT                 Sxx           Syy   ',...
    '        Sxy           S11           S22           S12           \n']);

if OUTPUT_ENVELOPE == true
    % Print the numerically largest stress in the current ply
    S_ply_xy_envelope_11 = abd.internal_getAbsMax(S_ply_xy(1.0, :), ENVELOPE_MODE);
    S_ply_xy_envelope_22 = abd.internal_getAbsMax(S_ply_xy(2.0, :), ENVELOPE_MODE);
    S_ply_xy_envelope_12 = abd.internal_getAbsMax(S_ply_xy(3.0, :), ENVELOPE_MODE);
    S_ply_aligned_envelope_11 = abd.internal_getAbsMax(S_ply_aligned(1.0, :), ENVELOPE_MODE);
    S_ply_aligned_envelope_22 = abd.internal_getAbsMax(S_ply_aligned(2.0, :), ENVELOPE_MODE);
    S_ply_aligned_envelope_12 = abd.internal_getAbsMax(S_ply_aligned(3.0, :), ENVELOPE_MODE);

    fprintf(fid, '%-5s%-30s%-14g%-14g%-14g%-14g%-14g%-14g\n',...
        'ALL', envelopeString, S_ply_xy_envelope_11,...
        S_ply_xy_envelope_22, S_ply_xy_envelope_12,...
        S_ply_aligned_envelope_11, S_ply_aligned_envelope_22,...
        S_ply_aligned_envelope_12);
else
    for i = 1.0:nPlies
        %{
                Extract the stresses for the current section point for the
                current ply
        %}
        currentOutputPoints = outputPoints(plyBuffer == i);

        if isempty(currentOutputPoints) == true
            fprintf(fid, '%-5.0f%s\n', i, 'NO RESULTS');
        else
            % Get the stresses for the current section points
            S_ply_xy_i = S_ply_xy(:, currentOutputPoints);
            S_ply_aligned_i = S_ply_aligned(:, currentOutputPoints);

            for sp = 1.0:length(currentOutputPoints)
                thicknessFraction = sprintf('%.0f (fraction = %.5g)',...
                    currentOutputPoints(sp),...
                    thickness(currentOutputPoints(sp)));

                fprintf(fid, '%-5.0f%-30s%-14g%-14g%-14g%-14g%-14g%-14g\n',...
                    i, thicknessFraction, S_ply_xy_i(1.0, sp),...
                    S_ply_xy_i(2.0, sp), S_ply_xy_i(3.0, sp),...
                    S_ply_aligned_i(1.0, sp),...
                    S_ply_aligned_i(2.0, sp),...
                    S_ply_aligned_i(3.0, sp));
            end
        end

        % Draw symmetry plane (if applicable)
        if (symmetricAbd == true) && (i == nPlies/2.0)
            fprintf(fid, ['- - - - - - - - - - - - - - - - - - - - - -',...
                ' - - - - - - - - - - - - - - - - - - - - - - - - - - ',...
                '- - - - - - - - - - - - SYM\n']);
        end
    end
end

%% Print strain tensor data
fprintf(fid, '\nStrain tensors at user-selected locations:\n');
fprintf(fid, ['PLY  SECTION POINT                 Exx           Eyy   ',...
    '        Exy           E11           E22           E12           \n']);

if OUTPUT_ENVELOPE == true
    % Print the numerically largest stress over all plies
    E_ply_xy_envelope_11 = abd.internal_getAbsMax(E_ply_xy(1.0, :), ENVELOPE_MODE);
    E_ply_xy_envelope_22 = abd.internal_getAbsMax(E_ply_xy(2.0, :), ENVELOPE_MODE);
    E_ply_xy_envelope_12 = abd.internal_getAbsMax(E_ply_xy(3.0, :), ENVELOPE_MODE);
    E_ply_aligned_envelope_11 = abd.internal_getAbsMax(E_ply_aligned(1.0, :), ENVELOPE_MODE);
    E_ply_aligned_envelope_22 = abd.internal_getAbsMax(E_ply_aligned(2.0, :), ENVELOPE_MODE);
    E_ply_aligned_envelope_12 = abd.internal_getAbsMax(E_ply_aligned(3.0, :), ENVELOPE_MODE);

    fprintf(fid, '%-5s%-30s%-14g%-14g%-14g%-14g%-14g%-14g\n',...
        'ALL', envelopeString, E_ply_xy_envelope_11,...
        E_ply_xy_envelope_22, E_ply_xy_envelope_12,...
        E_ply_aligned_envelope_11, E_ply_aligned_envelope_22,...
        E_ply_aligned_envelope_12);
else
    for i = 1.0:nPlies
        %{
                Extract the strains for the current section point for the
                current ply
        %}
        currentOutputPoints = outputPoints(plyBuffer == i);

        if isempty(currentOutputPoints) == true
            fprintf(fid, '%-5.0f%s\n', i, 'NO RESULTS');
        else
            % Get the strains for the current section points
            E_ply_xy_i = E_ply_xy(:, currentOutputPoints);
            E_ply_aligned_i = E_ply_aligned(:, currentOutputPoints);

            for sp = 1.0:length(currentOutputPoints)
                thicknessFraction = sprintf('%.0f (fraction = %.5g)',...
                    currentOutputPoints(sp),...
                    thickness(currentOutputPoints(sp)));

                fprintf(fid, '%-5.0f%-30s%-14g%-14g%-14g%-14g%-14g%-14g\n',...
                    i, thicknessFraction, E_ply_xy_i(1.0, sp),...
                    E_ply_xy_i(2.0, sp), E_ply_xy_i(3.0, sp),...
                    E_ply_aligned_i(1.0, sp),...
                    E_ply_aligned_i(2.0, sp),...
                    E_ply_aligned_i(3.0, sp));
            end
        end

        % Draw symmetry plane (if applicable)
        if (symmetricAbd == true) && (i == nPlies/2.0)
            fprintf(fid, ['- - - - - - - - - - - - - - - - - - - - - -',...
                ' - - - - - - - - - - - - - - - - - - - - - - - - - - ',...
                '- - - - - - - - - - - - SYM\n']);
        end
    end
end