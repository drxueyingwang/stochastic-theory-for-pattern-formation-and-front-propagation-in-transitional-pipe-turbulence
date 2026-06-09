% pdepe_2level_front_velocity.m
% =========================================================================
% PDE (pdepe) solver for the two-trophic-level predator-prey model of
% transitional pipe turbulence.  Computes traveling-front velocities and
% compares them to KPP-Fisher theory (Eq. 8 in the main text).
%
% Species:
%   rhoB : turbulence energy  (prey / B field)
%   rhoN : mean-flow energy   (nutrient / N field)
%
% Outputs:
%   - Diagnostic figures for each r value  (diag_r*.png)
%   - PDE front velocities saved to        PDE_front_velocity.txt
%   - KPP theory vs PDE comparison figure  front_velocity_KPP_vs_PDE.png
% =========================================================================

clear all; close all; clc;
tic

global r u01 u02;
global g_b g_p g_dB g_mmort g_cA g_cN g_DB;

% ========== Physical Parameters ==========
g_b     = 0.2;      % turbulence spontaneous decay rate              [1/s]
g_p     = 0.2;      % turbulence-nutrition interaction rate          [1/s]
g_dB    = 0.02;     % turbulence intrinsic death rate                [1/s]
g_mmort = 0.0002;   % turbulence additional mortality                [1/s]
g_cA    = 0.04;     % turbulence self-competition coefficient        [1/s]
g_cN    = 0.2;      % nutrition self-competition coefficient         [1/s]
g_DB    = 0.5;      % turbulence diffusion coefficient               [m^2/s]
% Nutrition production rate: g_func(U) = 20*U^2/9  (derived from flow speed U)
% Advection speed:           U = 0.15 * r           (set per r in the loop)

% ========== Non-Dimensionalization Scheme ==========
%   Time scale:    T    = 1/U        =>  tau    = t_phys * U
%   Density scale: rho0 = U / g_b   =>  rho_nd = rho_phys * (g_b / U)
%   Space:         xi   = x (physical; N advects at U/U = 1 in non-dim)
%
%   Non-dim coefficients (derived; see also pdefun1 below):
%     inter_nd = g_p  / g_b                    [interaction,       = g_p/g_b      ]
%     d_eff    = (g_dB + g_mmort) / U          [effective B decay, = (dB+m)/U     ]
%     c_A_nd   = g_cA / g_b                    [B self-competition, = cA/b        ]
%     alpha    = g_b  * (20/9)                 [N source, b*g/U^2 = b*20/9        ]
%     beta     = 20*U / 9                      [N decay,  g/U    = 20U/9          ]
%     c_N_nd   = g_cN / g_b                    [N self-competition, = cN/b        ]
%     D_B_nd   = g_DB / U                      [B diffusion,  = DB/U              ]
%   Rescale back: rho_phys = rho_nd * (U/g_b),  t_phys = tau / U

%% ========== Simulation Parameters ==========
ensemble = 300;           
r_all = [0.11,0.5,0.75,1,1.5,2,2.5,3,3.5];
front_v = zeros(ensemble, length(r_all), 3);

%% ========== Main Loop over r and Ensemble ==========
for ka = 1:length(r_all)
    r = r_all(ka);
    for e = 1:ensemble

        x = linspace(0, 5000, 5000);
        t = linspace(0, 300*r, 5000);

        % --- Non-dim coefficients (computed from physical parameters above) ---
        U        = 0.15 * r;                    % advection speed [m/s]
        inter_nd = g_p / g_b;                   % interaction coefficient
        d_eff    = (g_dB + g_mmort) / U;        % non-dim effective B decay = (dB+m)/U
        c_A_nd   = g_cA / g_b;                  % non-dim B self-competition = cA/b
        alpha    = g_b  * (20/9);               % non-dim N source = b*g/U^2 = b*20/9
        beta     = (20 * U) / 9;                % non-dim N decay  = g/U    = 20*U/9
        c_N_nd   = g_cN / g_b;                  % non-dim N self-competition = cN/b
        rho0     = U / g_b;                     % density scale for rescaling output

        % --- Fixed points ---

        % Laminar fixed point (rho_B = 0):
        %   Solve  alpha - beta*rhoN - c_N_nd*rhoN^2 = 0
        %   General quadratic solution:
        %     rhoN_lam = (-beta + sqrt(beta^2 + 4*c_N_nd*alpha)) / (2*c_N_nd)
        rhoN_lam = (-beta + sqrt(beta^2 + 4*c_N_nd*alpha)) / (2*c_N_nd);
        u01 = [0; rhoN_lam];

        % Invasion condition: turbulence invades laminar state if rhoN_lam > d_eff
        % (linearizing B-eq around rhoB=0: d(rhoB)/dtau = inter_nd*(rhoN_lam - d_eff)*rhoB)
        if rhoN_lam > d_eff
            % Coexistence fixed point — general quadratic.
            % From B-eq at steady state (divide by rhoB > 0):
            %   inter_nd*rhoN = d_eff + c_A_nd*rhoB  =>  rhoN = d_eff + c_A_nd*rhoB
            % Substituting w = rhoN into N-eq:
            %   (1 + c_A_nd*c_N_nd)*w^2 + (c_A_nd*beta - d_eff)*w - c_A_nd*alpha = 0
            a_q = 1 + c_A_nd * c_N_nd;
            b_q = c_A_nd * beta - d_eff;
            c_q = -c_A_nd * alpha;
            w   = (-b_q + sqrt(b_q^2 - 4*a_q*c_q)) / (2*a_q);
            rhoB_turb = max((w - d_eff) / c_A_nd, 0);
            rhoN_turb = w;
            u02 = [rhoB_turb; rhoN_turb];
        else
            % Puff regime: turbulence cannot sustain — set u02 close to laminar
            u02 = [1e-4; rhoN_lam * 0.99];
        end

        m_pde = 0;
        sol = pdepe(m_pde, @pdefun1, @pdeic1, @pdebc1, x, t);
        u1 = sol(:,:,1);
        u2 = sol(:,:,2);

        % --- Rescale non-dim solution back to physical units ---
        %   t_phys   = tau / U           = tau / (0.15*r)
        %   rho_phys = rho_nd * rho0    = rho_nd * (U / g_b)
        t  = t  / U;          % non-dim time   -> physical time
        u1 = u1 * rho0;       % non-dim rhoB   -> physical rhoB
        u2 = u2 * rho0;       % non-dim rhoN   -> physical rhoN

        % --- Track front positions at n_samples time snapshots ---
        n_samples = 10;
        dt_samp   = floor(length(t) / n_samples);
        front_x   = NaN(n_samples, 3);   % columns: [time, upstream_x, downstream_x]
        ki = 1;
        for i = dt_samp : dt_samp : length(t)
            pl = find(u1(i,1:end) < 0.5*max(u1(i,1:end)));
            if length(pl) < 2, ki = ki+1; continue; end
            gaps = pl(2:end) - pl(1:end-1);
            pl1  = find(gaps == max(gaps), 1);
            if isempty(pl1) || pl1+1 > length(pl), ki = ki+1; continue; end
            front_x(ki,:) = [t(i), x(pl(pl1)), x(pl(pl1+1))];
            ki = ki+1;
        end

        % --- Diagnostic plot for every r (first ensemble member only) ---
        if e == 1
            % (a) Turbulence profile evolution
            figure(200 + ka);
            clf;
            cmap = jet(n_samples);
            subplot(2,1,1);
            t_snap = dt_samp : dt_samp : length(t);
            for ii = 1:length(t_snap)
                plot(x, u1(t_snap(ii),:), 'Color', cmap(ii,:), 'LineWidth', 1.2);
                hold on;
            end
            xlabel('x'); ylabel('\rho_B (physical)');
            title(sprintf('Turbulence profile — r=%.2f  (U=%.4f)', r, U));
            colormap(jet); cb = colorbar;
            ylabel(cb, 'time snapshot index');
            grid on;

            % (b) Front positions vs time with polyfit
            subplot(2,1,2);
            valid = ~isnan(front_x(:,2));
            if sum(valid) >= 2
                t_fit  = front_x(valid,1);
                xup    = front_x(valid,2);   % upstream front
                xdn    = front_x(valid,3);   % downstream front
                p1_d   = polyfit(t_fit, xup, 1);
                p2_d   = polyfit(t_fit, xdn, 1);
                plot(t_fit, xup, 'bs', 'MarkerSize', 7, 'LineWidth', 1.5); hold on;
                plot(t_fit, xdn, 'r^', 'MarkerSize', 7, 'LineWidth', 1.5);
                t_line = linspace(min(t_fit), max(t_fit), 100);
                plot(t_line, polyval(p1_d, t_line), 'b--', 'LineWidth', 1.5);
                plot(t_line, polyval(p2_d, t_line), 'r--', 'LineWidth', 1.5);
                legend('Upstream front', 'Downstream front', 'Fit upstream', 'Fit downstream', ...
                       'Location', 'Best');
                xlabel('Physical time t');
                ylabel('Front position x');
                title(sprintf('Front positions — r=%.2f  v_{up}=%.4f  v_{dn}=%.4f', ...
                    r, p1_d(1), p2_d(1)));
                grid on;
            else
                text(0.5, 0.5, 'No valid front detected', 'Units', 'normalized', ...
                     'HorizontalAlignment', 'center');
            end
            drawnow;
            save_dir = fileparts(mfilename('fullpath'));
            if isempty(save_dir), save_dir = pwd; end
            diag_fname = fullfile(save_dir, sprintf('diag_r%.2f.png', r));
            saveas(figure(200 + ka), diag_fname, 'png');
            fprintf('Saved diagnostic figure: %s\n', diag_fname);
        end

        % --- Fit front position vs time to obtain front speeds ---
        valid = ~isnan(front_x(:,2));
        if sum(valid) >= 2
            p1_fit = polyfit(front_x(valid,1), front_x(valid,2), 1);  % upstream dx/dt
            p2_fit = polyfit(front_x(valid,1), front_x(valid,3), 1);  % downstream dx/dt
            v_up_lab = p1_fit(1);
            v_dn_lab = p2_fit(1);
        else
            v_up_lab = NaN;
            v_dn_lab = NaN;
        end

        front_v(e, ka, :) = [r, v_up_lab, v_dn_lab];
    end
end

%% ========== Theoretical Front Velocity (KPP theory) ==========
%
% The B equation has no advection term (effectively comoving frame), so
% the numerically measured front speeds are the INTRINSIC speeds c*.
% KPP formula at the leading edge (rhoB->0, rhoN->rhoN_lam):
%
%   lambda_phys = U * rhoN_lam_nd - (g_dB + g_mmort)     [physical growth rate]
%   c*(U) = 2 * sqrt(g_DB * lambda_phys)                  [KPP intrinsic speed]
U_per_r = 0.15;               % U = U_per_r * r

% Theory curve (uses same physical parameters as the simulation)
v_th     = (0:0.001:3.5) * U_per_r;   % advection speed array [m/s]
r_th     = v_th / U_per_r;            % corresponding r values
r_th(r_th == 0) = NaN;

alpha_val  = g_b * (20/9);
c_N_nd_val = g_cN / g_b;
beta_th    = (20 * v_th) / 9;

rhoN_lam_nd_th = (-beta_th + sqrt(beta_th.^2 + 4*c_N_nd_val*alpha_val)) ./ (2*c_N_nd_val);

lambda_th = v_th .* rhoN_lam_nd_th - (g_dB + g_mmort);
c_theory  = 2 * sqrt(g_DB .* max(lambda_th, 0));
c_theory(isnan(r_th)) = 0;

figure;
plot(v_th, c_theory, 'k-', 'LineWidth', 2);
hold on;

%% ========== Critical Advection Speed U_c ==========
%
% Condition: rhoN_lam_nd(r_c) = d_eff(r_c) = (g_dB + g_mmort) / (U_per_r * r_c)
inv_cond = @(rv) ...
    (-20*U_per_r*rv/9 + sqrt((20*U_per_r*rv/9).^2 + 4*c_N_nd_val*alpha_val)) ...
    ./ (2*c_N_nd_val) ...
    - (g_dB + g_mmort) ./ (U_per_r * rv);

r_crit = fzero(inv_cond, [0.01, 10]);
U_crit = U_per_r * r_crit;
fprintf('Turbulence onset: r_c = %.4f,  U_c = %.4f\n', r_crit, U_crit);

%% ========== Plot: KPP Theory vs PDE Numerics ==========
U_num = r_all * U_per_r;

m_front_v  = squeeze(mean(front_v, 1));    % (length(r_all), 3)
sd_front_v = squeeze(std(front_v,  0, 1)); % (length(r_all), 3)

v_upstream   = abs(m_front_v(:, 2));
v_downstream = abs(m_front_v(:, 3));
e_upstream   = sd_front_v(:, 2);
e_downstream = sd_front_v(:, 3);

% Zero out front speeds below onset (turbulent patch decays, c* = 0)
below_onset = U_num(:) < U_crit;
v_upstream(below_onset)   = 0;
v_downstream(below_onset) = 0;
e_upstream(below_onset)   = 0;
e_downstream(below_onset) = 0;

% Save PDE front speeds to text file for MFT comparison
% Columns: U | downstream_mean | downstream_std | upstream_mean | upstream_std
save_dir = fileparts(mfilename('fullpath'));
if isempty(save_dir), save_dir = pwd; end
fname = fullfile(save_dir, 'PDE_front_velocity.txt');
fid = fopen(fname, 'w');
fprintf(fid, 'U\t\tdownstream_mean\tdownstream_std\tupstream_mean\tupstream_std\n');
for k = 1:length(U_num)
    fprintf(fid, '%.6f\t%.6f\t\t%.6f\t\t%.6f\t\t%.6f\n', ...
        U_num(k), v_downstream(k), e_downstream(k), v_upstream(k), e_upstream(k));
end
fclose(fid);
fprintf('Saved PDE front velocities to: %s\n', fname);

errorbar(U_num, v_upstream,   e_upstream,   'bs', 'LineWidth', 1.5, 'MarkerSize', 8);
errorbar(U_num, v_downstream, e_downstream, 'r^', 'LineWidth', 1.5, 'MarkerSize', 8);

xline(U_crit, 'k--', sprintf('U_c = %.4f', U_crit), ...
      'LineWidth', 1.2, 'LabelVerticalAlignment', 'bottom', 'FontSize', 11);

xlabel('U (advection speed)', 'FontSize', 13);
ylabel('Intrinsic front speed c*', 'FontSize', 13);
title('Front velocity: KPP-F theory vs PDE numerics', 'FontSize', 13);
legend({'KPP-F theory', 'Upstream front (PDE)', 'Downstream front (PDE)'}, ...
       'Location', 'NorthWest', 'FontSize', 11);
set(gca, 'YScale', 'log');
grid on;

save_dir_main = fileparts(mfilename('fullpath'));
if isempty(save_dir_main), save_dir_main = pwd; end
saveas(gcf, fullfile(save_dir_main, 'front_velocity_KPP_vs_PDE.png'), 'png');
fprintf('Saved main figure: %s\n', fullfile(save_dir_main, 'front_velocity_KPP_vs_PDE.png'));

toc

%% ========== PDE Definition (pdefun1) ==========

function [c,f,s] = pdefun1(x,t,u,dudx)
global r u01 u02;
global g_b g_p g_dB g_mmort g_cA g_cN g_DB;
c = [1;1];

% Non-dim coefficients (derived from global physical parameters)
%   Changing a global at the top of the script automatically updates here.
U        = 0.15 * r;
inter_nd = g_p  / g_b;                % interaction       = p/b
d_eff_nd = (g_dB + g_mmort) / U;      % effective B decay = (dB+m)/U
c_A_nd   = g_cA / g_b;               % B self-competition = cA/b
alpha    = g_b  * (20/9);            % N source = b*g/U^2 = b*20/9  (U^2 cancels)
beta     = (20 * U) / 9;             % N decay  = g/U    = 20*U/9
c_N_nd   = g_cN / g_b;              % N self-competition = cN/b
D_B_nd   = g_DB / U;                 % B diffusion = DB/U

% PDE flux f and source s
%   rhoB: diffuses with D_B_nd; no advection (comoving frame for B)
%   rhoN: no diffusion; advects at speed 1 (= U/U) via -dudx(2) term
f = [D_B_nd * dudx(1); 0];

%   rhoB: + inter_nd*rhoN*rhoB - d_eff_nd*rhoB - c_A_nd*rhoB^2
%   rhoN: - inter_nd*rhoN*rhoB + alpha - beta*rhoN - c_N_nd*rhoN^2 - d(rhoN)/dxi
s = [ inter_nd*u(2)*u(1) - d_eff_nd*u(1) - c_A_nd*u(1)^2; ...
     -inter_nd*u(2)*u(1) + alpha          - beta*u(2)     - c_N_nd*u(2)^2 - dudx(2)];

end

%% ========== Initial Conditions (pdeic1) ==========

function u0 = pdeic1(x)
global u01 u02;

%   Turbulent seed centered at x=2500 (width ±10):
%     - inside seed: coexistence fixed point u02
%     - outside:     laminar fixed point u01
if x >= 2500-10 && x <= 2500+10
    u0 = u02;
else
    u0 = u01;
end

end

%% ========== Boundary Conditions (pdebc1) ==========

function [pl,ql,pr,qr] = pdebc1(xl,ul,xr,ur,t)
global u01 u02;

% Zero-flux (Neumann) boundary conditions on both ends
pl = [0; 0];
ql = [1; 1];
pr = [0; 0];
qr = [1; 1];
end
