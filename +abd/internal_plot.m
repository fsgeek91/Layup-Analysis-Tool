classdef internal_plot < handle
%   Plot MATLAB figures of output variables.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 5.0.0 Copyright Louis Vallance 2026
%   Last modified 11-Feb-2026 08:06:52 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%
    methods(Static = true, Access = public)
        %% MAIN FUNCTION FOR PLOTTING
        function [] = main(OUTPUT_FIGURE, SP_VIZ, PLOT_STYLE, outputLocation, nPlies, E_ply_xy, S_ply_xy, E_ply_aligned, S_ply_aligned, z, z_points, CRITERION_BUFFER,...
                OUTPUT_OPTIMISED, SP_COLOUR_BUFFER)

            % Plot parameters
            fontX = 14.0;
            fontY = 14.0;
            fontTitle = 16.0;
            fontTicks = 12.0;
            lineWidth = 1.5;

            % Set the figure title string buffer
            switch PLOT_STYLE
                case 'detached'
                    figureTitles_strainXY = {'XY ply strain in XX-direction for all section points',...
                        'XY ply strain in YY-direction for all section points',...
                        'XY ply strain in XY-direction for all section points'};
                    figureTitles_strainPly = {'Aligned ply strain in 11-direction for all section points',...
                        'Aligned ply strain in 22-direction for all section points',...
                        'Aligned ply strain in 12-direction for all section points'};

                    figureTitles_stressXY = {'XY ply stress in XX-direction for all section points',...
                        'XY ply stress in YY-direction for all section points',...
                        'XY ply stress in XY-direction for all section points'};
                    figureTitles_stressPly = {'Aligned ply stress in 11-direction for all section points',...
                        'Aligned ply stress in 22-direction for all section points',...
                        'Aligned ply stress in 12-direction for all section points'};
                otherwise
                    figureTitles_strainXY = {'XY ply strains for all section points'};
                    figureTitles_strainPly = {'Aligned ply strains for all section points'};

                    figureTitles_stressXY = {'XY ply stresses for all section points'};
                    figureTitles_stressPly = {'Aligned ply stresses for all section points'};
            end

            % Get the plot domain
            z_points_norm = z_points/max(z);
            z_plies_norm = z/max(z);

            %% Apply smoothing (if applicable)
            if strcmpi(OUTPUT_FIGURE, 'smooth') == true
                E_ply_xy = smoothdata(E_ply_xy, 2.0, 'loess');
                E_ply_aligned = smoothdata(E_ply_aligned, 2.0, 'loess');
                S_ply_xy = smoothdata(S_ply_xy, 2.0, 'loess');
                S_ply_aligned = smoothdata(S_ply_aligned, 2.0, 'loess');
            end

            %% EP, Ply strains in X-Y coordinates
            abd.internal_plot.now(...
                figureTitles_strainXY, {'\epsilon_x_x', '\epsilon_y_y', '\gamma_x_y'}, {'XY ply strain in XX-direction', 'XY ply strain in YY-direction',...
                'XY ply strain in XY-direction'}, PLOT_STYLE, E_ply_xy, z_points_norm, lineWidth, nPlies, z_plies_norm, fontTitle, fontTicks, fontX, fontY, 'Strain [mm/mm]',...
                outputLocation, [filesep, 'EP, '], SP_COLOUR_BUFFER, SP_VIZ)

            %% EP, Ply strains in ply coordinates
            abd.internal_plot.now(...
                figureTitles_strainPly, {'\epsilon_f_i_b_r_e', '\epsilon_t_r_a_n_s_v_e_r_s_e', '\gamma_p_l_y'}, {'Aligned ply strain in 11-direction',...
                'Aligned ply strain in 22-direction', 'Aligned ply strain in 12-direction'}, PLOT_STYLE, E_ply_aligned, z_points_norm, lineWidth, nPlies, z_plies_norm, fontTitle,...
                fontTicks, fontX, fontY, 'Strain [mm/mm]', outputLocation, [filesep, 'EP, '], SP_COLOUR_BUFFER, SP_VIZ)

            %% SP, Ply stresses in X-Y coordinates
            abd.internal_plot.now(...
                figureTitles_stressXY, {'\sigma_x_x', '\sigma_y_y', '\tau_x_y'}, {'XY ply stress in XX-direction', 'XY ply stress in YY-direction',...
                'XY ply stress in XY-direction'}, PLOT_STYLE, S_ply_xy, z_points_norm, lineWidth, nPlies, z_plies_norm, fontTitle, fontTicks, fontX, fontY, 'Stress [N/mm2]',...
                outputLocation, [filesep, 'SP, '], SP_COLOUR_BUFFER, SP_VIZ)

            %% SP, Ply stresses in ply coordinates
            abd.internal_plot.now(...
                figureTitles_stressPly, {'\sigma_f_i_b_r_e', '\sigma_t_r_a_n_s_v_e_r_s_e', '\tau_p_l_y'}, {'Aligned ply stress in 11-direction',...
                'Aligned ply stress in 22-direction', 'Aligned ply stress in 12-direction'}, PLOT_STYLE, S_ply_aligned, z_points_norm, lineWidth, nPlies, z_plies_norm, fontTitle,...
                fontTicks, fontX, fontY, 'Stress [N/mm2]', outputLocation, [filesep, 'SP, '], SP_COLOUR_BUFFER, SP_VIZ)

            %% CB, Optimised criterion buffer
            if isempty(CRITERION_BUFFER) == false
                % Create the figure
                f = abd.internal_plot.createFigure();

                % Set the figure title
                figureTitle = sprintf('CDF plot of optimiser criterion for all stacking permutations');

                % Plot the criterion
                h = cdfplot(CRITERION_BUFFER);
                set(h, 'LineWidth', 2.0)
                hold on

                % Activate the grid
                grid minor

                % Set the x-label string
                if (strcmpi(OUTPUT_OPTIMISED{2.0}, 'tsaih') == true || strcmpi(OUTPUT_OPTIMISED{2.0}, 'tsaiw') == true || strcmpi(OUTPUT_OPTIMISED{2.0}, 'azzit') == true)
                    if OUTPUT_OPTIMISED{3.0} == 1.0
                        xLabelString = sprintf('%s reserve factor', upper(OUTPUT_OPTIMISED{2.0}));
                    else
                        xLabelString = sprintf('%s value', upper(OUTPUT_OPTIMISED{2.0}));
                    end
                else
                    xLabelString = sprintf('%s value', upper(OUTPUT_OPTIMISED{2.0}));
                end

                % Set axis labels
                xlabel(xLabelString, 'FontSize', fontY)
                ylabel('Cumulative probability', 'FontSize', fontX)

                % Set the figure title
                title(figureTitle, 'FontSize', fontTitle)
                set(gca, 'FontSize', fontTicks)

                % Try to tighten the axes
                try
                    axis tight
                catch
                    % Don't tighten the axis
                end
                
                % Save the MATLAB figure to a file
                abd.internal_plot.save(outputLocation, [filesep, 'CB, '], figureTitle, f)
            end
        end

        %% CREATE A MATLAB FIGURE OF THE SELECTED PLOT VARIABLE
        function [] = now(figureTitle, legendStrings, plotTitle, PLOT_STYLE, VARIABLE, RANGE, lineWidth, nPlies, z_plies_norm, fontTitle, fontTicks, fontX, fontY, xlabelString,...
                outputLocation, leadString, SP_COLOUR_BUFFER, SP_VIZ)
            % Configure the plot layout
            switch lower(PLOT_STYLE)
                case 'compact'
                    % Single figure
                    N = ones(1.0, 3.0);
                    P = 1.0;

                    % Create the figure
                    f = abd.internal_plot.createFigure();
                case 'split'
                    % Tiled figures
                    N = 1.0:3.0;
                    P = 3.0;

                    % Create the figure
                    f = abd.internal_plot.createFigure();
                case 'detached'
                    % Individual figures
                    N = ones(1.0, 3.0);
                    P = 1.0;
                otherwise
                    % This case should never be reached!

                    % Tiled figures
                    N = 1.0:3.0;
                    P = 3.0;

                    % Create the figure
                    f = abd.internal_plot.createFigure();
            end

            % Initialise figure handle buffer
            H = zeros(1.0, 3.0);

            % Plot each tensor component in turn
            for plotNumber = 1.0:3.0
                % Create a new figure each time (detached mode only)
                if strcmpi(PLOT_STYLE, 'detached') == true
                    % Create the figure
                    f = abd.internal_plot.createFigure();
                end

                % Set the current plot space
                subplot(1.0, P, N(plotNumber))

                % Get the current domain
                DOMAIN = VARIABLE(plotNumber, :);

                % Plot the current variable
                H(plotNumber) = plot(DOMAIN, RANGE, 'LineWidth', lineWidth);
                hold on

                % Plot the ply boundaries
                if (strcmpi(PLOT_STYLE, 'split') == true) || (strcmpi(PLOT_STYLE, 'detached') == true)
                    abd.internal_plot.boundaries(nPlies, DOMAIN, z_plies_norm)
                end

                % Plot the section points
                if isempty(SP_VIZ) == false
                    if (strcmpi(PLOT_STYLE, 'compact') == true) && (plotNumber == 1.0)
                        scatter(linspace(0.5*(min(min(VARIABLE, [], 2.0)) + max(max(VARIABLE, [], 2.0))), 0.5*(min(min(VARIABLE, [], 2.0)) + max(max(VARIABLE, [], 2.0))),...
                            length(RANGE)), RANGE, 18.0, SP_COLOUR_BUFFER)
                    elseif (strcmpi(PLOT_STYLE, 'split') == true) || (strcmpi(PLOT_STYLE, 'detached') == true)
                        scatter(linspace(0.5*(min(DOMAIN) + max(DOMAIN)), 0.5*(min(DOMAIN) + max(DOMAIN)), length(RANGE)), RANGE, 18.0, SP_COLOUR_BUFFER)
                    end
                end

                % Set the legend and figure title
                if strcmpi(PLOT_STYLE, 'detached') == true
                    % Set the legend
                    legend(H(plotNumber), legendStrings{plotNumber});

                    % Get the title of the current figure
                    figureTitleCurrent = figureTitle{plotNumber};

                    % Set the figure title
                    title(figureTitleCurrent, 'FontSize', fontTitle)

                    % Set font tick size for all plots
                    set(gca, 'FontSize', fontTicks)
                elseif (P == 1.0) && (plotNumber == 3.0)
                    % Set the legend
                    L = legend(H, legendStrings);

                    % Do not allow further modifications to the legend
                    L.AutoUpdate = 'off';

                    % Set the figure title
                    title(figureTitle{1.0}, 'FontSize', fontTitle)

                    % Set font tick size for all plots
                    set(gca, 'FontSize', fontTicks)
                elseif P == 3.0
                    % Set the figure title
                    title(plotTitle{plotNumber}, 'FontSize', fontTitle)

                    % Set font tick size for the current plot
                    set(gca, 'FontSize', fontTicks)
                end

                % Activate the grid
                grid minor

                % Set axis labels
                xlabel(xlabelString, 'FontSize', fontX);
                ylabel('Thickness fraction [mm/mm]', 'FontSize', fontY);

                % Try to tighten the axes
                try
                    axis tight
                catch
                    % Don't tighten the axis
                end

                % Save the MATLAB figure to a file now (detached mode only)
                if strcmpi(PLOT_STYLE, 'detached') == true
                    abd.internal_plot.save(outputLocation, leadString, figureTitleCurrent, f)
                end
            end

            % Plot the ply boundaries
            if strcmpi(PLOT_STYLE, 'compact') == true
                abd.internal_plot.boundaries(nPlies, [min(min(VARIABLE, [], 2.0)), max(max(VARIABLE, [], 2.0))], z_plies_norm)
            end

            % Save the MATLAB figure to a file (excludes detached mode)
            if strcmpi(PLOT_STYLE, 'detached') == false
                abd.internal_plot.save(outputLocation, leadString, figureTitle{1.0}, f)
            end
        end

        %% PLOT THE PLY BOUNDARIES
        function [] = boundaries(nPlies, DOMAIN, RANGE)
            %{
                Add horizontal dashed lines to the plot indicating the
                boundary between each ply in the layup
            %}

            %{
                If the value of each element in DOMAIN is identical, adjust
                the endpoints of DOMAIN to force rendering of the ply
                boundaries
            %}
            if any(abs(diff(DOMAIN)) > 1e-17) == false
                % Get arbitrary element from DOMAIN
                DOMAIN_ELEMENT = DOMAIN(1.0);

                % Adjust endpoints by a small fraction of DOMAIN_ELEMENT
                if DOMAIN_ELEMENT == 0.0
                    DOMAIN(1.0) = -1e-12;
                    DOMAIN(end) = +1e-12;
                else
                    DOMAIN(1.0) = DOMAIN_ELEMENT - 1e-12*DOMAIN_ELEMENT;
                    DOMAIN(end) = DOMAIN_ELEMENT + 1e-12*DOMAIN_ELEMENT;
                end
            end

            %{
                Do not plot ply boundaries if there are greater than 50
                plies in the layup
            %}
            if nPlies < 51.0
                for i = 1:nPlies + 1.0
                    line([min(min(DOMAIN)), max(max(DOMAIN))], [RANGE(i), RANGE(i)], 'Color', [0.5, 0.5, 0.5], 'LineStyle', '--')
                end
            end
        end

        %% SAVE THE MATLAB FIGURE TO A FILE
        function [] = save(outputLocation, leadString, figureTitle, f)
            % Save the figure
            abd.internal_plot.saveFigure([outputLocation, leadString, figureTitle], f, 'fig');
        end

        %% GET DATA FROM OUTPUT_FIGURE
        function [error, output] = getSettings(OUTPUT_FIGURE)
            % Initialise output
            error = false;
            output = cell(1.0, 3.0);

            if iscell(OUTPUT_FIGURE) == false
                % Convert to cell if necessary
                OUTPUT_FIGURE = {OUTPUT_FIGURE};
            end

            if cellfun(@isempty, OUTPUT_FIGURE) == true
                % Set default values if necessary
                OUTPUT_FIGURE = {'DEFAULT', 'POINTS', 'SPLIT'};
            end

            if length(OUTPUT_FIGURE) ~= 3.0
                % Incorrect number of arguments
                fprintf('[ERROR] The setting OUTPUT_FIGURE requires three arguments:\n{''<mode>'', ''<section-visualisation>'', ''<layout>''}\n');

                % Reset the error flag and RETURN
                error = true;
                return
            end

            % Process the first argument
            argument = OUTPUT_FIGURE{1.0};

            if ((isempty(argument) == false) && (ischar(argument) == false)) || ((isempty(argument) == false) && (strcmpi(argument, 'default') == false) &&...
                    (strcmpi(argument, 'smooth') == false))
                % Incorrect variable type
                fprintf('[ERROR] The setting OUTPUT_FIGURE(1) must be one of the following:\n{''DEFAULT'' | ''SMOOTH''}\n');

                % Reset the error flag and RETURN
                error = true;
                return
            else
                % Everything is OK
                output{1.0} = argument;
            end

            % Process the second argument
            argument = OUTPUT_FIGURE{2.0};

            if ((isempty(argument) == false) && (ischar(argument) == false)) || ((isempty(argument) == false) && (strcmpi(argument, 'points') == false))
                % Incorrect variable type
                fprintf('[ERROR] The setting OUTPUT_FIGURE(2) must be one of the following:\n{[] | ''POINTS''}\n');

                % Reset the error flag and RETURN
                error = true;
                return
            else
                % Everything is OK
                output{2.0} = argument;
            end

            % Process the third argument
            argument = OUTPUT_FIGURE{3.0};

            if ((isempty(argument) == false) && (ischar(argument) == false)) || ((strcmpi(argument, 'compact') == false) && (strcmpi(argument, 'split') == false) &&...
                    (strcmpi(argument, 'detached') == false))
                % Incorrect variable type
                fprintf('[ERROR] The setting OUTPUT_FIGURE(3) must be one of the following:\n{''COMPACT'' | ''SPLIT'' | ''DETACHED''}\n');

                % Reset the error flag and RETURN
                error = true;
                return
            else
                % Everything is OK
                output{3.0} = argument;
            end
        end

        %% CREATE A MATLAB FIGURE
        function [f] = createFigure()
            % Get the figure visibility
            try
                % Try to get the user preferences object (Quick Fatigue Tool)
                [p, ~, ~] = qpref.direct;

                % Get the figure visibility
                figureVisibility = p.qftpref_output.figureVisibility;
            catch
                %{
                    There is no QFT user preference object, so use default
                    value of 'off'
                %}
                figureVisibility = 'off';
            end

            % Create the MATLAB figure window
            f = figure('visible', figureVisibility);

            % Set the figure visibility
            f.CreateFcn = 'set(gcbo, ''Visible'', ''on'')';
        end

        %% SAVE THE MATLAB FIGURE
        function [f] = saveFigure(dir, f, figureFormat)
            % Save the figure file
            saveas(f, dir, figureFormat)
        end
    end
end
