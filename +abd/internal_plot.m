classdef internal_plot < handle
%   Plot a MATLAB figure of the ply strain.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 2.6 Copyright Louis Vallance 2023
%   Last modified 17-May-2023 07:40:13 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%
    methods(Static = true, Access = public)
        %% MAIN FUNCTION FOR PLOTTING
        function [] = main(OUTPUT_FIGURE, PLOT_STYLE, outputLocation,...
                nPlies, E_ply_xy, S_ply_xy, E_ply_aligned,...
                S_ply_aligned, z, z_points, CRITERION_BUFFER,...
                OUTPUT_OPTIMISED)

            % Plot parameters
            fontX = 14.0;
            fontY = 14.0;
            fontTitle = 16.0;
            fontTicks = 12.0;
            lineWidth = 1.5;

            % Get the plot domain
            z_points_norm = z_points/max(z);
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
            abd.internal_plot.now(...
                'XY ply strains for all section points',...
                {'\epsilon_x_x', '\epsilon_y_y', '\gamma_x_y'},...
                {'XY ply strain in XX-direction',...
                'XY ply strain in YY-direction',...
                'XY ply strain in XY-direction'},...
                PLOT_STYLE, E_ply_xy, z_points_norm, lineWidth, nPlies,...
                z_plies_norm, increment, fontTitle, fontTicks, fontX,...
                fontY, outputLocation, '\EP, ')

            %% EP, Ply strains in ply coordinates
            abd.internal_plot.now(...
                'Aligned ply strains for all section points',...
                {'\epsilon_f_i_b_r_e', '\epsilon_t_r_a_n_s_v_e_r_s_e',...
                '\gamma_p_l_y'}, {'Aligned ply strain in 11-direction',...
                'Aligned ply strain in 22-direction',...
                'Aligned ply strain in 12-direction'},...
                PLOT_STYLE, E_ply_aligned, z_points_norm, lineWidth,...
                nPlies, z_plies_norm, increment, fontTitle, fontTicks,...
                fontX, fontY, outputLocation, '\EP, ')

            %% SP, Ply stresses in X-Y coordinates
            abd.internal_plot.now(...
                'XY ply stresses for all section points',...
                {'\sigma_x_x', '\sigma_y_y', '\tau_x_y'},...
                {'XY ply stress in XX-direction',...
                'XY ply stress in YY-direction',...
                'XY ply stress in XY-direction'},...
                PLOT_STYLE, S_ply_xy, z_points_norm, lineWidth, nPlies,...
                z_plies_norm, increment, fontTitle, fontTicks, fontX,...
                fontY, outputLocation, '\SP, ')

            %% SP, Ply stresses in ply coordinates
            abd.internal_plot.now(...
                'Aligned ply stresses for all section points',...
                {'\sigma_f_i_b_r_e', '\sigma_t_r_a_n_s_v_e_r_s_e',...
                '\tau_p_l_y'}, {'Aligned ply stress in 11-direction',...
                'Aligned ply stress in 22-direction',...
                'Aligned ply stress in 12-direction'},...
                PLOT_STYLE, S_ply_aligned, z_points_norm, lineWidth,...
                nPlies, z_plies_norm, increment, fontTitle, fontTicks,...
                fontX, fontY, outputLocation, '\SP, ')

            %% CB, Optimised criterion buffer
            if isempty(CRITERION_BUFFER) == false
                % Create the figure
                f5 = figure('visible', 'off');
                figureTitle = sprintf(['Optimiser criterion for all st',...
                    'acking permutations']);

                % Plot the crierion
                histogram(CRITERION_BUFFER);
                hold on

                % Other options
                grid minor
                if (strcmpi(OUTPUT_OPTIMISED{2.0}, 'tsaih') == true ||...
                        strcmpi(OUTPUT_OPTIMISED{2.0}, 'tsaiw') == true ||...
                        strcmpi(OUTPUT_OPTIMISED{2.0}, 'azzit') == true)
                    if OUTPUT_OPTIMISED{3.0} == 1.0
                        xLabelString = sprintf('%s reserve factor',...
                            upper(OUTPUT_OPTIMISED{2.0}));
                    else
                        xLabelString = sprintf('%s value',...
                            upper(OUTPUT_OPTIMISED{2.0}));
                    end
                else
                    xLabelString = sprintf('%s value', upper(OUTPUT_OPTIMISED{2.0}));
                end
                xlabel(xLabelString, 'FontSize', fontY)
                ylabel('Number of occurrences', 'FontSize', fontX)
                title(figureTitle, 'FontSize', fontTitle)
                set(gca, 'FontSize', fontTicks)

                try
                    axis tight
                catch
                    % Don't tighten the axis
                end

                % Set the figure path
                dir = [outputLocation, '\CB, ', figureTitle];
                saveas(f5, dir, 'fig')

                % Make the figure visible
                try
                    abd.internal_makeVisible([dir, '.fig'],...
                        abd.internal_getMATLABVersion)
                catch
                    %{
                        The contents of the figure are probably too large.
                        Accept the conequences and move on
                    %}
                end
            end
        end

        %% CREATE A MATLAB FIGURE OF THE SELECTED PLOT VARIABLE
        function [] = now(figureTitle, legendStrings, plotTitle,...
                PLOT_STYLE, VARIABLE, RANGE, lineWidth, nPlies,...
                z_plies_norm, increment, fontTitle, fontTicks, fontX,...
                fontY, outputLocation, leadString)
            % Create the figure
            f = figure('visible', 'off');

            % Initialise figure handle buffer
            H = zeros(1.0, 3.0);

            % Plot the strains
            if strcmpi(PLOT_STYLE, 'compact') == true
                N = ones(1.0, 3.0);
                P = 1.0;
            else
                N = 1.0:3.0;
                P = 3.0;
            end

            for plotNumber = 1.0:3.0
                % Set the current plot space
                subplot(1.0, P, N(plotNumber))

                % Get the current domain
                DOMAIN = VARIABLE(plotNumber, :);

                hold on

                % Plot the current variable
                H(plotNumber) = plot(DOMAIN, RANGE, 'LineWidth', lineWidth);

                % Plot the ply boundaries
                abd.internal_plot.boundaries(nPlies, DOMAIN, z_plies_norm, increment)

                % Set the legend and figure title
                if (P == 1.0) && (plotNumber == 3.0)
                    % Set the legend
                    legend(H, legendStrings)

                    % Set the figure title
                    title(figureTitle, 'FontSize', fontTitle)
                    set(gca, 'FontSize', fontTicks)
                elseif P == 3.0
                    % Set the figure title
                    title(plotTitle{plotNumber}, 'FontSize', fontTitle)
                    set(gca, 'FontSize', fontTicks)
                end

                % Activate the grid
                grid minor

                % Set axis labels
                xlabel('Strain [mm/mm]', 'FontSize', fontX);
                ylabel('Thickness fraction [mm/mm]', 'FontSize', fontY);

                try
                    axis tight
                catch
                    % Don't tighten the axis
                end
            end

            % Save the MATLAB figure to a file
            abd.internal_plot.save(outputLocation, leadString,...
                figureTitle, f)
        end

        %% PLOT THE PLY BOUNDARIES
        function [] = boundaries(nPlies, DOMAIN, RANGE, increment)
            if nPlies < 51.0
                step = 1.0;

                for i = 1:nPlies + 1.0
                    line([min(min(DOMAIN)), max(max(DOMAIN))],...
                    [RANGE(i), RANGE(i)], 'Color', 'green', 'LineStyle',...
                    '--')

                    step = step + increment;
                end
            end
        end

        %% SAVE THE MATLAB FIGURE TO A FILE
        function [] = save(outputLocation, leadString, figureTitle, f)
            % Set the figure path
            dir = [outputLocation, leadString, figureTitle];
            saveas(f, dir, 'fig')

            % Make the figure visible
            abd.internal_makeVisible([dir, '.fig'],...
                abd.internal_getMATLABVersion)
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

            if (all(cellfun(@isempty, OUTPUT_FIGURE)) == true) ||...
                    (length(OUTPUT_FIGURE) ~= 2.0)
                % Incorrect number of arguments
                fprintf(['[ABD ERROR] The setting OUTPUT_FIGURE requir',...
                    'es two arguments:\n{''<mode>'', ''<layout>''}\n']);

                % Reset the error flag and RETURN
                error = true;
                return
            end

            % Process the first argument
            argument = OUTPUT_FIGURE{1.0};

            if (isempty(argument) == false) && (ischar(argument) == false)
                % Incorrect variable type
                fprintf(['[ABD ERROR] The setting OUTPUT_FIGURE(1) mus',...
                    't be a string: {''DEFAULT'' | ''SMOOTH''}\n']);

                % Reset the error flag and RETURN
                error = true;
                return
            elseif (isempty(argument) == false) &&...
                    (strcmpi(argument, 'default') == false) &&...
                    (strcmpi(argument, 'smooth') == false)
                % Incorrect mode tag
                fprintf(['[ABD ERROR] The setting OUTPUT_FIGURE(1) mus',...
                    't be one of the following: {''DEFAULT'' | ''SMOOT',...
                    'H''}\n']);

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
                fprintf(['[ABD ERROR] The setting OUTPUT_FIGURE(2) mus',...
                    't be a string: {''COMPACT'' | ''SPLIT''}\n']);

                % Reset the error flag and RETURN
                error = true;
                return
            elseif (strcmpi(argument, 'compact') == false) &&...
                    (strcmpi(argument, 'split') == false)
                % Incorrect mode tag
                fprintf(['[ABD ERROR] The setting OUTPUT_FIGURE(2) mus',...
                    't be one of the following: {''COMPACT'' | ''SPLIT',...
                    '''}\n']);

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
