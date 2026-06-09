% MC_simulation_stochastic_model.m
% =========================================================================
% Monte Carlo simulation of the three-trophic-level predator-prey model
% for transitional pipe turbulence (Eq. 1 in the paper).
%
% Species:
%   x (A): Predator = zonal flow energy
%   y (B): Prey     = turbulence energy
%   g (N): Nutrient  = baseline mean flow energy
%
% Stochastic reactions (Eq. 1):
%   A_i  --d_A-->  phi                     (death of predator)
%   B_i  --d_B-->  phi                     (death of prey)
%   A_i + B_j  --p/V-->  A_i + A_j         (nearest-neighbor predation)
%   B_i + N_j  --b/V-->  B_i + B_j         (nearest-neighbor birth)
%   B_i  --m-->  A_i                        (mutation)
%   A_i  --D_A-->  A_j                     (predator diffusion)
%   B_i  --D_B-->  B_j                     (prey diffusion)
%   phi  --g-->  N_i                        (nutrient recovery)
%   2A_i --c_A/V--> A_i                    (predator competition)
%   2B_i --c_B/V--> B_i                    (prey competition)
%   N_i  --speed=U-->  N_{i+1}             (nutrient advection)
%
% Outputs:
%   - Ensemble-averaged spatial profiles of prey (rho_B) and nutrient (rho_N)
%     for Figs. 1(a)-(c) and 2(a)-(c)
%   - Space-time contour plots for Fig. 5(a)-(c)
%   - Phase space trajectories (rho_B vs sqrt(rho_N)) for Fig. 3(a)-(c)
%   - Instantaneous density gain and loss rates for Fig. 4(a)-(i)
%   - Front positions for velocity calculation (Fig. 6)
% =========================================================================

clearvars; clc; close all;
tic

% ========== Simulation Parameters ==========
% Control parameter values (U = food advection speed)
% U = 0.0165 -> puff phase        (Re ~ 2310)
% U = 0.1125 -> weak slug phase   (Re ~ 3235)
% U = 0.525  -> strong slug phase (Re ~ 4345)
% advection_speed_range = [0.0165, 0.1125, 0.525];
advection_speed_range = [0.0165, 0.1125, 0.525];
n_ensemble = 300;        % Number of ensemble realizations
system_length = 3000;     % Lattice sites along pipe axis
system_width = 20;        % Lattice sites around pipe circumference
threshold_density = 0.5;  % Threshold to define turbulent region
time_steps = 2e4;       % Total number of time steps
plot_interval = 10;       % Interval for storing space-time data
site_capacity = 5;        % Maximum nutrient per site (N_max)
V = 5;                    % Correlation volume (controls noise amplitude)
t_nullcline = 1e4;        % Time step at which to record snapshot

n_speeds = length(advection_speed_range);
n_time_frames = time_steps / plot_interval;


% ========== Storage Arrays ==========
% Front/back tracking for Fig. 6 inset
front = zeros(n_time_frames, 2, n_ensemble, n_speeds);
back  = zeros(n_time_frames, 2, n_ensemble, n_speeds);

% Ensemble-averaged space-time data for Figs. 1, 2, 5(a)-(c)
avg_spacetime_y = zeros(n_time_frames, system_length, n_speeds);  % prey
avg_spacetime_g = zeros(n_time_frames, system_length, n_speeds);  % nutrient

% Phase space trajectories for Fig. 3 (one snapshot per U per ensemble)
nullcline = zeros(system_length, 2, n_ensemble, n_speeds);

% Energy budget for Fig. 4 (spatial profiles for each U)
y_energy = zeros(system_length, 7, n_ensemble, n_speeds);

% ========== Main Simulation Loop ==========
for i = 1:n_speeds
    food_adv_speed = advection_speed_range(i);
    fprintf('Running U = %.4f (%d/%d)...\n', food_adv_speed, i, n_speeds);

    for e = 1:n_ensemble
        if mod(e, 10) == 0
            fprintf('  Ensemble %d/%d\n', e, n_ensemble);
        end

        % Set base parameters
        D  = 0.5;
        d1 = 0.02;
        d2 = 0.02;
        p  = 0.2;
        m  = 0.0002;
        b  = 0.2;
        c  = 0.04;
        growth_rate = 20/9 * food_adv_speed^2;

        % Rescale all rates by food advection speed (co-moving frame)
        D  = D  / food_adv_speed;
        d1 = d1 / food_adv_speed;
        d2 = d2 / food_adv_speed;
        p  = p  / food_adv_speed;
        m  = m  / food_adv_speed;
        b  = b  / food_adv_speed;
        c  = c  / food_adv_speed;
        growth_rate = growth_rate / food_adv_speed;

        % Storage for space-time plots (this ensemble)
        plot_data_y = zeros(n_time_frames, system_length);
        plot_data_g = zeros(n_time_frames, system_length);

        % Initialize population grids
        x = zeros(system_width, system_length);                    % Predator A
        y = zeros(system_width, system_length);                    % Prey B
        g = site_capacity * ones(system_width, system_length);     % Nutrient N

        % Initial condition: small turbulent seed in center
        mid = system_length - 1500;
        x(:, mid-15:mid+15) = rand(system_width, 31) > 2/5;
        y(:, mid-15:mid+15) = rand(system_width, 31) > 2/5;
        g(:, mid+15:end) = 0;   % Depleted nutrient downstream of seed

        % Energy accounting arrays (accumulate over entire simulation)
        y_loss_p  = zeros(system_width, system_length);  % loss from predation
        y_loss_d  = zeros(system_width, system_length);  % loss from death
        y_loss_c  = zeros(system_width, system_length);  % loss from competition
        y_gain_b  = zeros(system_width, system_length);  % gain from birth
        y_loss_m  = zeros(system_width, system_length);  % loss from mutation
        y_diffuse = zeros(system_width, system_length);  % net from diffusion

        % Time evolution
        for t = 1:time_steps
            g(:,1) = site_capacity;  % Refill nutrient at upstream boundary

            y_before_step = y;  % Save state before this timestep

            % Each time step: attempt system_width * system_length reactions
            for k = 1:system_width * system_length
                jx = randi(system_width);
                jy = randi(system_length);
                s = rand();

                if s < 1/10
                    % --- Reaction 1: Predator diffusion ---
                    if rand() < 1 - exp(-D * x(jx, jy))
                        dir = randi(4);
                        [jx_new, jy_new] = get_neighbor(jx, jy, dir, system_width, system_length);
                        x(jx_new, jy_new) = x(jx_new, jy_new) + 1;
                        x(jx, jy) = x(jx, jy) - 1;
                    end

                elseif s < 2/10
                    % --- Reaction 2: Prey diffusion ---
                    if rand() < 1 - exp(-D * y(jx, jy))
                        dir = randi(4);
                        [jx_new, jy_new] = get_neighbor(jx, jy, dir, system_width, system_length);
                        y(jx_new, jy_new) = y(jx_new, jy_new) + 1;
                        y(jx, jy) = y(jx, jy) - 1;
                        y_diffuse(jx_new, jy_new) = y_diffuse(jx_new, jy_new) + 1;
                        y_diffuse(jx, jy) = y_diffuse(jx, jy) - 1;
                    end

                elseif s < 3/10
                    % --- Reaction 3: Prey birth (nearest-neighbor) ---
                    %     B_i + N_j -> B_i + B_j
                    dir = randi(4);
                    [jx_new, jy_new] = get_neighbor(jx, jy, dir, system_width, system_length);
                    if rand() < 1 - exp(-b * y(jx, jy) * g(jx_new, jy_new) / V)
                        y(jx_new, jy_new) = y(jx_new, jy_new) + 1;
                        g(jx_new, jy_new) = max(0, g(jx_new, jy_new) - 1);
                        y_gain_b(jx_new, jy_new) = y_gain_b(jx_new, jy_new) + 1;
                    end

                elseif s < 4/10
                    % --- Reaction 4: Predation (nearest-neighbor) ---
                    %     A_i + B_j -> A_i + A_j
                    dir = randi(4);
                    [jx_new, jy_new] = get_neighbor(jx, jy, dir, system_width, system_length);
                    if rand() < 1 - exp(-p * x(jx, jy) * y(jx_new, jy_new) / V)
                        x(jx_new, jy_new) = x(jx_new, jy_new) + 1;
                        y(jx_new, jy_new) = max(0, y(jx_new, jy_new) - 1);
                        y_loss_p(jx_new, jy_new) = y_loss_p(jx_new, jy_new) - 1;
                    end

                elseif s < 5/10
                    % --- Reaction 5: Death of predator ---
                    %     A_i -> phi
                    if rand() < 1 - exp(-d1 * x(jx, jy))
                        x(jx, jy) = max(0, x(jx, jy) - 1);
                    end

                elseif s < 6/10
                    % --- Reaction 6: Death of prey ---
                    %     B_i -> phi
                    if rand() < 1 - exp(-d2 * y(jx, jy))
                        y(jx, jy) = max(0, y(jx, jy) - 1);
                        y_loss_d(jx, jy) = y_loss_d(jx, jy) - 1;
                    end

                elseif s < 7/10
                    % --- Reaction 7: Mutation ---
                    %     B_i -> A_i
                    if rand() < 1 - exp(-m * y(jx, jy))
                        x(jx, jy) = x(jx, jy) + 1;
                        y(jx, jy) = max(0, y(jx, jy) - 1);
                        y_loss_m(jx, jy) = y_loss_m(jx, jy) - 1;
                    end

                elseif s < 8/10
                    % --- Reaction 8: Predator competition ---
                    %     2A_i -> A_i
                    if rand() < 1 - exp(-c * x(jx, jy) * (x(jx, jy) - 1) / V)
                        x(jx, jy) = max(0, x(jx, jy) - 1);
                    end

                elseif s < 9/10
                    % --- Reaction 9: Prey competition ---
                    %     2B_i -> B_i
                    if rand() < 1 - exp(-c * y(jx, jy) * (y(jx, jy) - 1) / V)
                        y(jx, jy) = max(0, y(jx, jy) - 1);
                        y_loss_c(jx, jy) = y_loss_c(jx, jy) - 1;
                    end

                else
                    % --- Reaction 10: Nutrient recovery ---
                    %     phi -> N_i
                    if rand() < 1 - exp(-growth_rate * V)
                        g(jx, jy) = min(site_capacity, g(jx, jy) + 1);
                    end
                end
            end

            % Nutrient advection: shift g one site to the right
            g = circshift(g, 1, 2);

            % Store data for space-time plots
            if mod(t, plot_interval) == 0
                idx = t / plot_interval;
                plot_data_y(idx, :) = mean(y);
                plot_data_g(idx, :) = mean(g);

                % Track front and back positions
                pos = find(mean(y) >= threshold_density);
                if length(pos) > 1
                    front(idx, :, e, i) = [t / food_adv_speed, pos(1)];
                    back(idx, :, e, i)  = [t / food_adv_speed, pos(end)];
                end
            end

            % Record snapshot for phase space and energy budget at t_nullcline
            if t == t_nullcline
                nullcline(:, 1, e, i) = mean(g)';
                nullcline(:, 2, e, i) = mean(y)';

                y_energy(:, 1, e, i) = mean(y_gain_b);   % gain from birth
                y_energy(:, 2, e, i) = mean(y_loss_p);   % loss from predation
                y_energy(:, 3, e, i) = mean(y_loss_c);   % loss from competition
                y_energy(:, 4, e, i) = mean(y_loss_d);   % loss from death
                y_energy(:, 5, e, i) = mean(y_loss_m);   % loss from mutation
                y_energy(:, 6, e, i) = mean(y_diffuse);  % diffusion
                y_energy(:, 7, e, i) = mean(y_before_step - y); % net change
            end
        end

        % Accumulate for ensemble-averaged space-time data
        avg_spacetime_y(:, :, i) = avg_spacetime_y(:, :, i) + plot_data_y;
        avg_spacetime_g(:, :, i) = avg_spacetime_g(:, :, i) + plot_data_g;
    end

    % Normalize by number of ensembles
    avg_spacetime_y(:, :, i) = avg_spacetime_y(:, :, i) / n_ensemble;
    avg_spacetime_g(:, :, i) = avg_spacetime_g(:, :, i) / n_ensemble;
end

% ========== Plotting ==========

% --- Fig. 1: Ensemble-averaged prey density snapshots ---
plot_snapshots(avg_spacetime_y, advection_speed_range, system_length, ...
    '\rho_B (prey / turbulence energy)', 'Fig. 1: Prey density snapshots');

% --- Fig. 2: Ensemble-averaged nutrient density snapshots ---
% Paper plots sqrt(rho_N) which represents mean flow speed
plot_snapshots_sqrt(avg_spacetime_g, advection_speed_range, system_length, ...
    '\surd\rho_N (mean flow speed)', 'Fig. 2: Nutrient density snapshots');

% --- Fig. 3: Phase space trajectories ---
plot_phase_space(nullcline, advection_speed_range);

% --- Fig. 4: Energy budget ---
plot_energy_budget(y_energy, advection_speed_range);

% --- Fig. 5(a)-(c): Space-time contour plots ---
plot_spacetime_contour(avg_spacetime_y, advection_speed_range, ...
    system_length, time_steps, plot_interval);

% --- Fig. 6 inset: Front velocities ---
[fv_mean, fv_std, bv_mean, bv_std] = ...
    calculate_front_velocity(front, back, advection_speed_range, plot_interval);

toc

% ========== Helper Functions ==========

function [jx_new, jy_new] = get_neighbor(jx, jy, dir, width, len)
    % Get neighbor with periodic boundary conditions on all boundaries
    switch dir
        case 1, jx_new = mod(jx - 2, width) + 1; jy_new = jy;    % Up
        case 2, jx_new = mod(jx, width) + 1;     jy_new = jy;    % Down
        case 3, jx_new = jx; jy_new = mod(jy - 2, len) + 1;      % Left
        case 4, jx_new = jx; jy_new = mod(jy, len) + 1;          % Right
    end
end

% ========== Plotting / Analysis Functions ==========

function plot_snapshots(data, speeds, sys_len, ylabel_str, title_str)
    % Plot spatial snapshots at the last time frame for each U.
    figure;
    phase_names = {'(a) Puff', '(b) Weak slug', '(c) Strong slug'};
    n = length(speeds);
    for i = 1:n
        subplot(1, n, i);
        % Use last time frame as the snapshot
        last_frame = squeeze(data(end, :, i));
        plot(1:sys_len, last_frame, 'LineWidth', 1.5);
        xlabel('Space'); ylabel(ylabel_str);
        if i <= length(phase_names)
            title(sprintf('%s (U=%.4f)', phase_names{i}, speeds(i)));
        else
            title(sprintf('U=%.4f', speeds(i)));
        end
    end
    sgtitle(title_str);
end

function plot_snapshots_sqrt(data, speeds, sys_len, ylabel_str, title_str)
    % Plot sqrt of spatial snapshots (for nutrient = mean flow speed).
    figure;
    phase_names = {'(a) Puff', '(b) Weak slug', '(c) Strong slug'};
    n = length(speeds);
    for i = 1:n
        subplot(1, n, i);
        last_frame = squeeze(data(end, :, i));
        plot(1:sys_len, sqrt(max(0, last_frame)), 'LineWidth', 1.5);
        xlabel('Space'); ylabel(ylabel_str);
        if i <= length(phase_names)
            title(sprintf('%s (U=%.4f)', phase_names{i}, speeds(i)));
        else
            title(sprintf('U=%.4f', speeds(i)));
        end
    end
    sgtitle(title_str);
end

function plot_phase_space(nullcline, speeds)
    % Fig. 3: Phase space trajectories (rho_B vs sqrt(rho_N))
    figure;
    phase_names = {'(a) Puff', '(b) Weak slug', '(c) Strong slug'};
    n = length(speeds);
    for i = 1:n
        mn = mean(nullcline(:, :, :, i), 3);  % average over ensembles
        subplot(1, n, i);
        plot(sqrt(max(0, mn(:,1))), mn(:,2), '-', 'LineWidth', 2);
        xlabel('\surd\rho_N'); ylabel('\rho_B');
        if i <= length(phase_names)
            title(sprintf('%s (U=%.4f)', phase_names{i}, speeds(i)));
        else
            title(sprintf('U=%.4f', speeds(i)));
        end
    end
    sgtitle('Fig. 3: Phase space trajectories');
end

function plot_energy_budget(y_energy, speeds)
    % Fig. 4: Energy budget for each U value
    phase_names = {'Puff', 'Weak slug', 'Strong slug'};
    row_labels = {'Rate of prey density gain', 'Rate of prey density loss', ...
                  'Rate of net density change'};
    n = length(speeds);
    figure;
    for i = 1:n
        y_avg = mean(y_energy(:, :, :, i), 3);  % (system_length x 7)
        x_axis = 1:size(y_avg, 1);

        % Row 1: Birth gain
        subplot(3, n, i);
        plot(x_axis, y_avg(:,1), 'LineWidth', 1.5);
        if i <= length(phase_names)
            title(sprintf('%s (U=%.4f)', phase_names{i}, speeds(i)));
        else
            title(sprintf('U=%.4f', speeds(i)));
        end
        if i == 1, ylabel(row_labels{1}); end

        % Row 2: Total loss
        subplot(3, n, n + i);
        plot(x_axis, sum(y_avg(:,2:5), 2), 'LineWidth', 1.5);
        if i == 1, ylabel(row_labels{2}); end

        % Row 3: Net change
        subplot(3, n, 2*n + i);
        plot(x_axis, y_avg(:,7), 'LineWidth', 1.5);
        xlabel('Space');
        if i == 1, ylabel(row_labels{3}); end
    end
    sgtitle('Fig. 4: Energy budget');
end

function plot_spacetime_contour(data, speeds, sys_len, time_steps, plot_interval)
    % Fig. 5(a)-(c): Space-time contour plots of prey density
    figure;
    phase_names = {'(a) Stochastic puff', '(b) Stochastic weak slug', ...
                   '(c) Stochastic strong slug'};
    n = length(speeds);
    for i = 1:n
        subplot(n, 1, i);
        n_frames = size(data, 1);
        t_axis = (1:n_frames) * plot_interval / speeds(i);  % physical time
        contourf(1:sys_len, t_axis, squeeze(data(:,:,i)), 'LineStyle', 'none');
        xlabel('Space'); ylabel('Time');
        if i <= length(phase_names)
            title(phase_names{i});
        else
            title(sprintf('U=%.4f', speeds(i)));
        end
        colorbar;
    end
    sgtitle('Fig. 5(a)-(c): MC space-time plots');
end

function [fv_mean, fv_std, bv_mean, bv_std] = calculate_front_velocity(front, back, speeds, dt)
    [~, ~, n_ens, n_spd] = size(front);
    fv_mean = zeros(n_spd, 1); fv_std = zeros(n_spd, 1);
    bv_mean = zeros(n_spd, 1); bv_std = zeros(n_spd, 1);
    for i = 1:n_spd
        fv = zeros(n_ens, 1); bv = zeros(n_ens, 1);
        for e = 1:n_ens
            td = diff(squeeze(front(:,1,e,i)));   
            sf = diff(squeeze(front(:,2,e,i)));
            sb = diff(squeeze(back(:,2,e,i)));
            valid = td > 0;
            if any(valid)
                fv(e) = mean(sf(valid) ./ td(valid));
                bv(e) = mean(sb(valid) ./ td(valid));
            end
        end
        fv_mean(i) = mean(fv); fv_std(i) = std(fv);
        bv_mean(i) = mean(bv); bv_std(i) = std(bv);
    end

    % ---- Save front speeds to text file for MFT comparison --------------
    % Columns: U | downstream_mean | downstream_std | upstream_mean | upstream_std
    save_dir = fileparts(mfilename('fullpath'));
    if isempty(save_dir), save_dir = pwd; end
    fname = fullfile(save_dir, 'MC_front_velocity.txt');
    fid = fopen(fname, 'w');
    fprintf(fid, 'U\t\tdownstream_mean\tdownstream_std\tupstream_mean\tupstream_std\n');
    for k = 1:n_spd
        fprintf(fid, '%.6f\t%.6f\t\t%.6f\t\t%.6f\t\t%.6f\n', ...
            speeds(k), abs(fv_mean(k)), fv_std(k), abs(bv_mean(k)), bv_std(k));
    end
    fclose(fid);
    fprintf('Saved MC front velocities to: %s\n', fname);

   
    D_phys     = 0.5;
    b_phys     = 0.2;
    d_B_phys   = 0.02;
    m_phys     = 0.0002;
    N_lam      = 5;          % site_capacity
    n_channels = 10;
    lambda_bare = b_phys * N_lam - d_B_phys - m_phys;
    c_bare = 2 * sqrt(D_phys * lambda_bare);
    c_eff  = c_bare / n_channels;

    % ---- Plot -----------------------------------------------------------
    figure;
    hold on
    
    h_mf = errorbar(speeds, abs(fv_mean), fv_std, 'ro', 'LineWidth', 1.5, 'MarkerSize', 8);
    h_mb = errorbar(speeds, abs(bv_mean), bv_std, 'bs', 'LineWidth', 1.5, 'MarkerSize', 8);

    set(gca, 'YScale', 'log');
    xlabel('U (advection speed)');
    ylabel('|Front speed|  (lattice / physical time)');
    legend([h_mf h_mb], ...
            {'Downstream front (MC)',         ...
            'Upstream front (MC)'},          ...
           'Location', 'best');
    title(sprintf('Fig. 6 : MC front velocity vs U   (c* = %.3f, c*/n_{ch} = %.3f)', ...
          c_bare, c_eff));
    grid on;
end
