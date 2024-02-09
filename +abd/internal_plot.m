classdef internal_plot < handle
%   Plot MATLAB figures of output variables.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 2.7 Copyright Louis Vallance 2024
%   Last modified 09-Feb-2024 09:10:19 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%
    methods(Static = true, Access = public)
        %% MAIN FUNCTION FOR PLOTTING
        function [] = main(OUTPUT_FIGURE, PLOT_STYLE, outputLocation, nPlies, E_ply_xy, S_ply_xy, E_ply_aligned, S_ply_aligned, z, z_points, CRITERION_BUFFER, OUTPUT_OPTIMISED)

            % Plot parameters
            fontX = 14.0;
            fontY = 14.0;
            fontTitle = 16.0;
            fontTicks = 12.0;
            lineWidth = 1.5;

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
                'XY ply strains for all section points', {'\epsilon_x_x', '\epsilon_y_y', '\gamma_x_y'}, {'XY ply strain in XX-direction', 'XY ply strain in YY-direction',...
                'XY ply strain in XY-direction'}, PLOT_STYLE, E_ply_xy, z_points_norm, lineWidth, nPlies, z_plies_norm, fontTitle, fontTicks, fontX, fontY, 'Strain [mm/mm]',...
                outputLocation, '\EP, ')

            %% EP, Ply strains in ply coordinates
            abd.internal_plot.now(...
                'Aligned ply strains for all section points', {'\epsilon_f_i_b_r_e', '\epsilon_t_r_a_n_s_v_e_r_s_e', '\gamma_p_l_y'}, {'Aligned ply strain in 11-direction',...
                'Aligned ply strain in 22-direction', 'Aligned ply strain in 12-direction'}, PLOT_STYLE, E_ply_aligned, z_points_norm, lineWidth, nPlies, z_plies_norm, fontTitle,...
                fontTicks, fontX, fontY, 'Strain [mm/mm]', outputLocation, '\EP, ')

            %% SP, Ply stresses in X-Y coordinates
            abd.internal_plot.now(...
                'XY ply stresses for all section points', {'\sigma_x_x', '\sigma_y_y', '\tau_x_y'}, {'XY ply stress in XX-direction', 'XY ply stress in YY-direction',...
                'XY ply stress in XY-direction'}, PLOT_STYLE, S_ply_xy, z_points_norm, lineWidth, nPlies, z_plies_norm, fontTitle, fontTicks, fontX, fontY, 'Stress [N/mm2]',...
                outputLocation, '\SP, ')

            %% SP, Ply stresses in ply coordinates
            abd.internal_plot.now(...
                'Aligned ply stresses for all section points', {'\sigma_f_i_b_r_e', '\sigma_t_r_a_n_s_v_e_r_s_e', '\tau_p_l_y'}, {'Aligned ply stress in 11-direction',...
                'Aligned ply stress in 22-direction', 'Aligned ply stress in 12-direction'}, PLOT_STYLE, S_ply_aligned, z_points_norm, lineWidth, nPlies, z_plies_norm, fontTitle,...
                fontTicks, fontX, fontY, 'Stress [N/mm2]', outputLocation, '\SP, ')

            %% CB, Optimised criterion buffer
            if isempty(CRITERION_BUFFER) == false
                % Create the figure
                f = figure('visible', 'off');

                % Set the figure title
                figureTitle = sprintf('Optimiser criterion for all stacking permutations');

                % Plot the criterion
                histogram(CRITERION_BUFFER);
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
                ylabel('Number of occurrences', 'FontSize', fontX)

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
                abd.internal_plot.save(outputLocation, '\CB, ', figureTitle, f)
            end
        end

        %% CREATE A MATLAB FIGURE OF THE SELECTED PLOT VARIABLE
        function [] = now(figureTitle, legendStrings, plotTitle, PLOT_STYLE, VARIABLE, RANGE, lineWidth, nPlies, z_plies_norm, fontTitle, fontTicks, fontX, fontY, xlabelString,...
                outputLocation, leadString)
            % Create the figure
            f = figure('visible', 'off');

            % Initialise figure handle buffer
            H = zeros(1.0, 3.0);

            % Configure the plot layout
            if strcmpi(PLOT_STYLE, 'compact') == true
                % Single figure
                N = ones(1.0, 3.0);
                P = 1.0;
            else
                % Tiled figures
                N = 1.0:3.0;
                P = 3.0;
            end

            % Plot each tensor component in turn
            for plotNumber = 1.0:3.0
                % Set the current plot space
                subplot(1.0, P, N(plotNumber))

                % Get the current domain
                DOMAIN = VARIABLE(plotNumber, :);

                % Plot the current variable
                H(plotNumber) = plot(DOMAIN, RANGE, 'LineWidth', lineWidth);
                hold on

                % Plot the ply boundaries
                abd.internal_plot.boundaries(nPlies, DOMAIN, z_plies_norm)

                % Set the legend and figure title
                if (P == 1.0) && (plotNumber == 3.0)
                    % Set the legend
                    legend(H, legendStrings)

                    % Set the figure title
                    title(figureTitle, 'FontSize', fontTitle)

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
            end

            % Save the MATLAB figure to a file
            abd.internal_plot.save(outputLocation, leadString, figureTitle, f)
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
                DOMAIN(1.0) = DOMAIN_ELEMENT - 1e-12*DOMAIN_ELEMENT;
                DOMAIN(end) = DOMAIN_ELEMENT + 1e-12*DOMAIN_ELEMENT;
            end

            %{
                Do not plot ply boundaries if there are greater than 50
                plies in the layup
            %}
            if nPlies < 51.0
                for i = 1:nPlies + 1.0
                    line([min(min(DOMAIN)), max(max(DOMAIN))], [RANGE(i), RANGE(i)], 'Color', 'green', 'LineStyle', '--')
                end
            end
        end

        %% SAVE THE MATLAB FIGURE TO A FILE
        function [] = save(outputLocation, leadString, figureTitle, f)
            % Set the figure path
            dir = [outputLocation, leadString, figureTitle];
            saveas(f, dir, 'fig')

            % Make the figure visible
            try
                abd.internal_makeVisible([dir, '.fig'], abd.internal_getMATLABVersion)
            catch
                %{
                    The contents of the figure are probably too large.
                    Accept the conequences and move on
                %}
            end
        end

        %% GET DATA FROM OUTPUT_FIGURE
        function [error, output] = getSettings(OUTPUT_FIGURE)
            % Initialise output
            error = false;
            output = cell(1.0, 2.0);

            if iscell(OUTPUT_FIGURE) == false
                % Convert to cell if necessary
                OUTPUT_FIGURE = {OUTPUT_FIGURE};
            end

            if (all(cellfun(@isempty, OUTPUT_FIGURE)) == true) || (length(OUTPUT_FIGURE) ~= 2.0)
                % Incorrect number of arguments
                fprintf('[ERROR] The setting OUTPUT_FIGURE requires two\narguments: {''<mode>'', ''<layout>''}\n');

                % Reset the error flag and RETURN
                error = true;
                return
            end

            % Process the first argument
            argument = OUTPUT_FIGURE{1.0};

            if (isempty(argument) == false) && (ischar(argument) == false)
                % Incorrect variable type
                fprintf('[ERROR] The setting OUTPUT_FIGURE(1) must be a string:\n{''DEFAULT'' | ''SMOOTH''}\n');

                % Reset the error flag and RETURN
                error = true;
                return
            elseif (isempty(argument) == false) &&...
                    (strcmpi(argument, 'default') == false) &&...
                    (strcmpi(argument, 'smooth') == false)
                % Incorrect mode tag
                fprintf('[ERROR] The setting OUTPUT_FIGURE(1) must be one of\nthe following: {''DEFAULT'' | ''SMOOTH''}\n');

                % Reset the error flag and RETURN
                error = true;
                return
            else
                % Everything is OK
                output{1.0} = argument;
            end

            % Process the second argument
            argument = OUTPUT_FIGURE{2.0};

            if ischar(argument) == false
                % Incorrect variable type
                fprintf('[ERROR] The setting OUTPUT_FIGURE(2) must be a\nstring: {''COMPACT'' | ''SPLIT''}\n');

                % Reset the error flag and RETURN
                error = true;
                return
            elseif (strcmpi(argument, 'compact') == false) && (strcmpi(argument, 'split') == false)
                % Incorrect mode tag
                fprintf('[ERROR] The setting OUTPUT_FIGURE(2) must be one of\nthe following: {''COMPACT'' | ''SPLIT''}\n');

                % Reset the error flag and RETURN
                error = true;
                return
            else
                % Everything is OK
                output{2.0} = argument;
            end
        end
    end
end
