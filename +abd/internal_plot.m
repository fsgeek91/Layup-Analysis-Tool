function [] = internal_plot(OUTPUT_FIGURE, outputLocation, nPlies,...
    E_ply_xy, S_ply_xy, E_ply_aligned, S_ply_aligned, z, z_points)
%   Plot a MATLAB figure of the ply strain.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 2.3 Copyright Louis Vallance 2023
%   Last modified 09-May-2023 07:31:07 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%

% Plot parameters
fontX = 12.0;
fontY = 12.0;
fontTitle = 14.0;
fontTicks = 12.0;
lineWidth = 1.5;

% Get the plot domain
z_points_norm = z_points/max(z_points);
z_plies_norm = z/max(z);
increment = floor(length(unique(z_points))/nPlies);

%% Apply smoothing (if applicable)
if strcmpi(OUTPUT_FIGURE, 'smooth') == true
    E_ply_xy = smoothdata(E_ply_xy, 2.0, 'loess');
    E_ply_aligned = smoothdata(E_ply_aligned, 2.0, 'loess');
    S_ply_xy = smoothdata(S_ply_xy, 2.0, 'loess');
    S_ply_aligned = smoothdata(S_ply_aligned, 2.0, 'loess');
end

%% EP, Ply strains in X-Y coordinates
% Create the figure
f1 = figure('visible', 'off');
figureTitle = 'XY ply strains for all section points';

% Plot the strains
p1 = plot(E_ply_xy(1.0, :), z_points_norm, 'LineWidth', lineWidth);
hold on
p2 = plot(E_ply_xy(2.0, :), z_points_norm, 'LineWidth', lineWidth);
p3 = plot(E_ply_xy(3.0, :), z_points_norm, 'LineWidth', lineWidth);

% Plot the ply boundaries
if nPlies < 51.0
    step = 1.0;

    for i = 1:nPlies + 1.0
        line([min(min(E_ply_xy)), max(max(E_ply_xy))], [z_plies_norm(i), z_plies_norm(i)],...
            'Color', 'green', 'LineStyle', '--')

        step = step + increment;
    end
end

% Other options
legend([p1, p2, p3], '\epsilon_x_x', '\epsilon_y_y', '\gamma_x_y')
grid minor
xlabel('Strain [mm/mm]', 'FontSize', fontX);
ylabel('Thickness fraction [mm/mm]', 'FontSize', fontY);
title(figureTitle, 'FontSize', fontTitle)
set(gca, 'FontSize', fontTicks)

try
    axis tight
catch
    % Don't tighten the axis
end

% Set the figure path
dir = [outputLocation, '\EP, ', figureTitle];
saveas(f1, dir, 'fig')

% Make the figure visible
abd.internal_makeVisible([dir, '.fig'], abd.internal_getMATLABVersion)

%% EP, Ply strains in ply coordinates
% Create the figure
f2 = figure('visible', 'off');
figureTitle = 'Aligned ply strains for all section points';

% Plot the strains
p1 = plot(E_ply_aligned(1.0, :), z_points_norm, 'LineWidth', lineWidth);
hold on
p2 = plot(E_ply_aligned(2.0, :), z_points_norm, 'LineWidth', lineWidth);
p3 = plot(E_ply_aligned(3.0, :), z_points_norm, 'LineWidth', lineWidth);

% Plot the ply boundaries
if nPlies < 51.0
    step = 1.0;

    for i = 1:nPlies + 1.0
        line([min(min(E_ply_aligned)), max(max(E_ply_aligned))], [z_plies_norm(i), z_plies_norm(i)],...
            'Color', 'green', 'LineStyle', '--')

        step = step + increment;
    end
end

% Other options
legend([p1, p2, p3], '\epsilon_f_i_b_r_e', '\epsilon_t_r_a_n_s_v_e_r_s_e', '\gamma_p_l_y')
grid minor
xlabel('Strain [mm/mm]', 'FontSize', fontX);
ylabel('Thickness fraction [mm/mm]', 'FontSize', fontY);
title(figureTitle, 'FontSize', fontTitle)
set(gca, 'FontSize', fontTicks)

try
    axis tight
catch
    % Don't tighten the axis
end

% Set the figure path
dir = [outputLocation, '\EP, ', figureTitle];
saveas(f2, dir, 'fig')

% Make the figure visible
abd.internal_makeVisible([dir, '.fig'], abd.internal_getMATLABVersion)

%% SP, Ply stresses in X-Y coordinates
% Create the figure
f3 = figure('visible', 'off');
figureTitle = 'XY ply stresses for all section points';

% Plot the stresses
p1 = plot(S_ply_xy(1.0, :), z_points_norm, 'LineWidth', lineWidth);
hold on
p2 = plot(S_ply_xy(2.0, :), z_points_norm, 'LineWidth', lineWidth);
p3 = plot(S_ply_xy(3.0, :), z_points_norm, 'LineWidth', lineWidth);

% Plot the ply boundaries
if nPlies < 51.0
    step = 1.0;

    for i = 1:nPlies + 1.0
        line([min(min(S_ply_xy)), max(max(S_ply_xy))], [z_plies_norm(i), z_plies_norm(i)],...
            'Color', 'green', 'LineStyle', '--')

        step = step + increment;
    end
end

% Other options
legend([p1, p2, p3], '\sigma_x_x', '\sigma_y_y', '\tau_x_y')
grid minor
xlabel('Stress [N/mm2]', 'FontSize', fontX)
ylabel('Thickness fraction [mm/mm]', 'FontSize', fontY)
title(figureTitle, 'FontSize', fontTitle)
set(gca, 'FontSize', fontTicks)

try
    axis tight
catch
    % Don't tighten the axis
end

% Set the figure path
dir = [outputLocation, '\SP, ', figureTitle];
saveas(f3, dir, 'fig')

% Make the figure visible
abd.internal_makeVisible([dir, '.fig'], abd.internal_getMATLABVersion)

%% SP, Ply stresses in ply coordinates
% Create the figure
f4 = figure('visible', 'off');
figureTitle = 'Aligned ply stresses for all section points';

% Plot the stresses
p1 = plot(S_ply_aligned(1.0, :), z_points_norm, 'LineWidth', lineWidth);
hold on
p2 = plot(S_ply_aligned(2.0, :), z_points_norm, 'LineWidth', lineWidth);
p3 = plot(S_ply_aligned(3.0, :), z_points_norm, 'LineWidth', lineWidth);

% Plot the ply boundaries
if nPlies < 51.0
    step = 1.0;

    for i = 1:nPlies + 1.0
        line([min(min(S_ply_xy)), max(max(S_ply_xy))], [z_plies_norm(i), z_plies_norm(i)],...
            'Color', 'green', 'LineStyle', '--')

        step = step + increment;
    end
end

% Other options
legend([p1, p2, p3], '\sigma_f_i_b_r_e', '\sigma_t_r_a_n_s_v_e_r_s_e', '\tau_p_l_y')
grid minor
xlabel('Stress [N/mm2]', 'FontSize', fontX)
ylabel('Thickness fraction [mm/mm]', 'FontSize', fontY)
title(figureTitle, 'FontSize', fontTitle)
set(gca, 'FontSize', fontTicks)

try
    axis tight
catch
    % Don't tighten the axis
end

% Set the figure path
dir = [outputLocation, '\SP, ', figureTitle];
saveas(f4, dir, 'fig')

% Make the figure visible
abd.internal_makeVisible([dir, '.fig'], abd.internal_getMATLABVersion)
end