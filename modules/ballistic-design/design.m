% =========================================================================
% FILE: design.m
% DESCRIPTION: SOLID ROCKET MOTOR - PRE-DESIGN
%
% Propellant : AP/HTPB composite (O/F = 6)
% Grain      : BATES (N segments, end-burning + cylindrical port)
% Nozzle     : Conical convergent-divergent with throat circular fillet
%
% Pipeline:
%   1.  Vieille Law (burn rate) fit from strand-burner data
%   2.  Preliminary sizing targets (thrust, impulse, burn time)
%   3.  Analytical BATES grain optimisation (burn time + pressure flatness)
%   4.  Conical nozzle sizing
%   5.  Internal lateral area Alat (convergent + fillet arc + divergent)
%   6.  2D proportional motor sketch
%   7.  Quasi-steady Kn/Pc evolution plots
%
% Note: SI units throughout unless stated; pressure in [bar] for Vieille Law.
% =========================================================================

clear; clc; close all;
%% =========================================================================
%%  SETUP
%% =========================================================================
project_root = fileparts(mfilename('fullpath'));
addpath(fullfile(project_root, 'engine'));
cea_work_dir = fullfile(project_root, 'cea');

% --- Common plot settings: explicit light theme, immune to user preferences ---
set(groot, ...
    'DefaultAxesFontSize',        12, ...
    'DefaultAxesFontName',        'Helvetica', ...
    'DefaultLineLineWidth',       1.8, ...
    'DefaultAxesBox',             'on', ...
    'DefaultAxesGridAlpha',       0.3, ...
    'DefaultFigureColor',         [0.97 0.97 0.97], ...  % figure background: near-white
    'DefaultAxesColor',           [1.00 1.00 1.00], ...  % axes background: white
    'DefaultAxesXColor',          [0.15 0.15 0.15], ...  % X tick/label: dark grey
    'DefaultAxesYColor',          [0.15 0.15 0.15], ...  % Y tick/label: dark grey
    'DefaultAxesZColor',          [0.15 0.15 0.15], ...  % Z tick/label
    'DefaultAxesGridColor',       [0.15 0.15 0.15], ...  % grid lines
    'DefaultAxesTitleFontWeight', 'bold', ...
    'DefaultTextColor',           [0.15 0.15 0.15]);     % all text objects

%% BALLISTIC DATA
% Strand-burner measurements (three replicates at five nominal pressures)
P = [10.1 9.7 10.3  30.2 31.0 29.8  50.2 51.0 50.3 ...
     70.2 69.0 69.8  91.2 89.1 89.0];   % Chamber pressure  [bar]
r = [4.0  3.8  4.1   5.6  6.0  5.7   7.0  7.2  7.1 ...
      8.4  8.3  8.6   8.8  9.0  9.2]; % Burning rate      [mm/s]

% Log-log linear regression: Vieille Law  r = a * P^n
%[a, Inc_a, n, Inc_n, R2] = Uncertainty_Vieille(P, r);
[a, Inc_a, n, Inc_n, R2] = Uncertainty(P, r);
Pc_target_bar = 65;                        % [bar]
rb_target_mms = a * Pc_target_bar^n;       % [mm/s]


fprintf('--- VIEILLE LAW FIT ---\n');
fprintf('Vieille Law: r = a * P^n\n');
fprintf('a      = %.6f mm/s/bar^n  (95%% CI: +/- %.6f)\n', a, Inc_a);
fprintf('n      = %.6f [-]         (95%% CI: +/- %.6f)\n', n, Inc_n);
fprintf('R^2    = %.6f [-]\n', R2);
fprintf('rb_tgt = %.6f mm/s \n\n', rb_target_mms );

% Fit plot
P_fit = linspace(min(P), max(P), 300);
r_fit = a * P_fit.^n;


figure('Name', 'Vieille Law Fit', 'NumberTitle', 'off', 'Position', [100 100 700 420]);
plot(P, r, 'ro'); hold on;
plot(P_fit, r_fit, '-', 'Color', [0.15 0.45 0.75], 'LineWidth', 2.2);
xlabel('Pressure  [bar]');
ylabel('Burning rate  [mm/s]');
title('Vieille Law  -  burn rate fit');
legend('Strand-burner data', sprintf('Fit: r = %.4f*P^{%.4f}  (R^2=%.4f)', a, n, R2), ...
    'Location', 'northwest');
grid on;

%% INPUTS & TARGETS

Pc_lim        = 70;                        % [bar]
Pc            = Pc_target_bar * 1e5;       % [Pa]
T_req         = 100e3;                     % [N]
I_tot_req     = 2.5e6;                     % [N*s]

OF              = 6;                       % [-]
ap_percentage   = OF / (OF + 1);           % [-]
htpb_percentage = 1  / (OF + 1);           % [-]
rho_htpb        = 930;                     % [kg/m^3]
rho_ap          = 1949;                    % [kg/m^3]
rho_propellant  = 1 / (ap_percentage/rho_ap + htpb_percentage/rho_htpb); % [kg/m^3]

% CEA
M_mol    = 25.955;                         % [g/mol]
gamma    = 1.1640;                         % [-]
C_star   = 1531.5;                         % [m/s]
c_ideal  = 2448.3;                         % [m/s]
epsilon  = 8.8491;                         % [-]

g0       = 9.81;                           % [m/s^2]
eta_comb = 0.95;                           % [-]

alpha    = 15;                             % [deg]
beta     = 45;                             % [deg]

lambda_conical = (1 + cosd(alpha)) / 2;    % [-]
c_eff_design   = eta_comb * c_ideal * lambda_conical; % [m/s]
Isp_design     = c_eff_design / g0;        % [s]

m_dot_target  = T_req / c_eff_design;      % [kg/s]
M_tot_target  = I_tot_req / c_eff_design;  % [kg]
t_b_target    = M_tot_target / m_dot_target; % [s]
A_t_target    = (m_dot_target * C_star) / Pc; % [m^2]
rb_target     = rb_target_mms / 1000;      % [m/s]
A_b_target    = m_dot_target / (rb_target * rho_propellant); % [m^2]

fprintf('--- PROPELLANT ---\n');
fprintf('O/F                = %.3f [-]\n',        OF);
fprintf('Pc_target          = %.2f bar\n',        Pc_target_bar);
fprintf('Pc_limit           = %.2f bar\n',        Pc_lim);
fprintf('AP mass fraction   = %.4f  (%.2f %%)\n', ap_percentage,   100*ap_percentage);
fprintf('HTPB mass fraction = %.4f  (%.2f %%)\n', htpb_percentage, 100*htpb_percentage);
fprintf('Propellant density = %.2f kg/m^3\n\n',   rho_propellant);

fprintf('--- DESIGN TARGETS ---\n');
fprintf('lambda_conical     = %.5f [-]\n',   lambda_conical);
fprintf('c_eff_design       = %.2f m/s\n',   c_eff_design);
fprintf('Isp_design         = %.2f s\n',     Isp_design);
fprintf('m_dot target       = %.4f kg/s\n',  m_dot_target);
fprintf('Propellant mass    = %.2f kg\n',    M_tot_target);
fprintf('Burn time          = %.2f s\n',     t_b_target);
fprintf('A_t target         = %.6f m^2\n',   A_t_target);
fprintf('r_b at Pc_target   = %.4f mm/s\n',  rb_target_mms);
fprintf('A_b average        = %.4f m^2\n\n', A_b_target);

%% =========================
%  OPTIMIZATION (ANALYTICAL BATES)
%  =========================
best       = struct();
best_error = inf;

w_ideal = rb_target * t_b_target;          % [m]

disp('--- RUNNING ANALYTICAL BATES OPTIMIZATION ---');

for N = 1:6
    for w_mult = linspace(0.85, 1.15, 600)
        w = w_ideal * w_mult;              % [m]

        inner_val = w^2 + 8*A_b_target / (N*pi);
        if inner_val < 0; continue; end
        D  = (3*w + sqrt(inner_val)) / 4;  % [m]
        d0 = D - 2*w;                      % [m]
        L0 = 2*d0 + 3*w;                   % [m]

        if d0 <= 0.02*D || L0 <= 2*w; continue; end

        x_vec    = linspace(0, w, 150);    % [m]
        d_x      = d0 + 2*x_vec;           % [m]
        L_x      = L0 - 2*x_vec;           % [m]
        A_b_x    = N * pi * (d_x.*L_x + 0.5*(D^2 - d_x.^2)); % [m^2]
        Kn_x     = A_b_x ./ A_t_target;    % [-]
        Pc_bar_x = (Kn_x .* a .* rho_propellant .* C_star ./ 1e8).^(1/(1-n)); % [bar]

        if any(~isreal(Pc_bar_x)) || any(~isfinite(Pc_bar_x)) || any(Pc_bar_x <= 0); continue; end

        rb_x = a .* Pc_bar_x.^n ./ 1000;   % [m/s]
        if any(~isfinite(rb_x)) || any(rb_x <= 0); continue; end

        t_est    = trapz(x_vec, 1./rb_x);  % [s]
        V_geom   = N * (pi/4) * (D^2 - d0^2) * L0; % [m^3]
        M_geom   = V_geom * rho_propellant; % [kg]
        Pc_avg   = mean(Pc_bar_x);         % [bar]
        flatness = (max(Pc_bar_x) - min(Pc_bar_x)) / Pc_avg; % [-]

        err_t = abs(t_est - t_b_target)     / t_b_target;
        err_M = abs(M_geom - M_tot_target)  / M_tot_target;
        err_P = abs(Pc_avg - Pc_target_bar) / Pc_target_bar;

        cost = 15*flatness + 10*err_P + 10*err_t + 5*err_M;
        if D > 0.95; cost = cost + (D - 0.95)*50; end

        if cost < best_error
            best_error       = cost;
            best.N           = N;
            best.D           = D;
            best.d0          = d0;
            best.L0          = L0;
            best.w           = w;
            best.Mass        = M_geom;
            best.t_est       = t_est;
            best.A_b0        = A_b_x(1);
            best.Kn_initial  = Kn_x(1);
            best.Kn_max      = max(Kn_x);
            best.Pc_mean     = Pc_avg;
            best.Pc_pp       = max(Pc_bar_x) - min(Pc_bar_x);
            best.error_Mass  = err_M;
            best.error_tb    = err_t;
            best.error_Pmean = err_P;
            best.flatness    = flatness;
            best.x_vec       = x_vec;
            best.A_b_x       = A_b_x;
            best.Kn_x        = Kn_x;
            best.Pc_bar_x    = Pc_bar_x;
            best.rb_x        = rb_x;
        end
    end
end

if isempty(fieldnames(best))
    error('No feasible BATES design found. Check input targets.');
end

%% =========================
%  RESULTS - GRAIN GEOMETRY
%  =========================
L_grain_total = best.N * best.L0;          % [m]
V_grain_total = best.N * (pi/4) * (best.D^2 - best.d0^2) * best.L0; % [m^3]

fprintf('--- BEST GRAIN CONFIGURATION ---\n');
fprintf('Number of segments N  = %d\n',       best.N);
fprintf('Outer diameter  D     = %.4f m\n',   best.D);
fprintf('Port diameter   d0    = %.4f m\n',   best.d0);
fprintf('Segment length  L0    = %.4f m\n',   best.L0);
fprintf('Total grain length    = %.4f m\n',   L_grain_total);
fprintf('Web thickness         = %.4f m\n',   best.w);
fprintf('Propellant volume     = %.6f m^3\n', V_grain_total);
fprintf('Initial burn area     = %.4f m^2\n', best.A_b0);
fprintf('Initial Kn            = %.3f\n',     best.Kn_initial);
fprintf('Maximum Kn            = %.3f\n',     best.Kn_max);
fprintf('Propellant mass       = %.2f kg  (err = %.2f %%)\n', best.Mass,    100*best.error_Mass);
fprintf('Estimated burn time   = %.2f s   (err = %.2f %%)\n', best.t_est,   100*best.error_tb);
fprintf('Average chamber Pc    = %.2f bar (err = %.2f %%)\n', best.Pc_mean, 100*best.error_Pmean);
fprintf('Pressure peak-to-peak = %.2f bar\n',   best.Pc_pp);
fprintf('Pressure flatness     = %.2f %%\n\n',  100*best.flatness);

%% =========================
%  NOZZLE SIZING
%  =========================
r_t  = sqrt(A_t_target / pi);              % [m]
d_t  = 2 * r_t;                            % [m]
A_e  = epsilon * A_t_target;               % [m^2]
r_e  = sqrt(A_e / pi);                     % [m]
d_e  = 2 * r_e;                            % [m]
R_c  = best.D / 2;                         % [m]

L_conv   = (R_c - r_t) / tand(beta);       % [m]
L_div    = (r_e - r_t) / tand(alpha);      % [m]
L_nozzle = L_conv + L_div;                 % [m]

fprintf('--- NOZZLE GEOMETRY ---\n');
fprintf('Throat radius  r_t    = %.4f m\n',   r_t);
fprintf('Throat diameter d_t   = %.4f m\n',   d_t);
fprintf('Exit radius    r_e    = %.4f m\n',   r_e);
fprintf('Exit diameter  d_e    = %.4f m\n',   d_e);
fprintf('Exit area      A_e    = %.6f m^2\n', A_e);
fprintf('Chamber radius R_c    = %.4f m\n',   R_c);
fprintf('Convergent length     = %.4f m\n',   L_conv);
fprintf('Divergent length      = %.4f m\n',   L_div);
fprintf('Total nozzle length   = %.4f m\n\n', L_nozzle);

%% =========================
%  NOZZLE LATERAL AREA (Alat)
%  =========================
rho_f     = 0.5 * r_t;                     % [m]
phi_total = beta + alpha;                  % [deg]
phi_rad   = deg2rad(phi_total);            % [rad]
arc_len   = rho_f * phi_rad;               % [m]

R_bar_fillet = r_t + rho_f - rho_f * sin(phi_rad/2) / (phi_rad/2); % [m]

r_tan_conv = r_t + rho_f * (1 - cosd(beta)); % [m]
r_tan_div  = r_t + rho_f * (1 - cosd(alpha));% [m]

s_conv = hypot(L_conv - rho_f*sind(beta),  R_c - r_tan_conv);      % [m]
s_div  = hypot(L_div  - rho_f*sind(alpha), r_e - r_tan_div);       % [m]

A_lat_conv   = pi * (R_c + r_tan_conv) * s_conv;                   % [m^2]
A_lat_fillet = 2 * pi * R_bar_fillet * arc_len;                    % [m^2]
A_lat_div    = pi * (r_tan_div + r_e) * s_div;                     % [m^2]
Alat         = A_lat_conv + A_lat_fillet + A_lat_div;              % [m^2]

fprintf('--- NOZZLE INTERNAL LATERAL AREA ---\n');
fprintf('Fillet radius  rho_f  = %.4f m\n',   rho_f);
fprintf('Fillet arc angle      = %.1f deg\n', phi_total);
fprintf('A_lat convergent      = %.6f m^2\n', A_lat_conv);
fprintf('A_lat fillet          = %.6f m^2\n', A_lat_fillet);
fprintf('A_lat divergent       = %.6f m^2\n', A_lat_div);
fprintf('A_lat TOTAL           = %.6f m^2\n\n', Alat);

%% =========================
%  2D SCALED MOTOR SKETCH  (coordinates in cm)
%  =========================
x_chamber_start = -(L_grain_total + L_conv); % [m]
x_conv_start    = -L_conv;                 % [m]
x_exit          =  L_div;                  % [m]

% Tangency points of the throat fillet
x_tan_c    = -rho_f * sind(beta);          % [m]
x_tan_d    =  rho_f * sind(alpha);         % [m]
r_tan_conv =  r_t + rho_f * (1 - cosd(beta)); % [m]
r_tan_div  =  r_t + rho_f * (1 - cosd(alpha));% [m]

% Fillet arc: circle centered at (0, r_t + rho_f)
theta    = linspace(-beta, alpha, 120);    % [deg]
x_fillet = rho_f * sind(theta);            % [m]
r_fillet = r_t + rho_f * (1 - cosd(theta));% [m]

% Nozzle upper contour without duplicated tangency points
x_noz = [x_conv_start, x_tan_c, x_fillet(2:end-1), x_tan_d, x_exit]; % [m]
r_noz = [R_c,          r_tan_conv, r_fillet(2:end-1), r_tan_div, r_e]; % [m]

% Lower contour
x_lower = fliplr(x_noz);                   % [m]
r_lower = -fliplr(r_noz);                  % [m]

t_wall = 0.05 * R_c;                       % [m]
r_cas  = R_c + t_wall;                     % [m]
cm     = 100;                              % Conversion factor to [cm]

figure('Name', '2D Motor Sketch', 'NumberTitle', 'off', 'Position', [100 100 1100 480]);
hold on;

% Outer casing envelope
fill([x_chamber_start x_exit x_exit x_chamber_start]*cm, ...
     [r_cas r_cas -r_cas -r_cas]*cm, ...
     [0.55 0.55 0.55], 'EdgeColor', 'none', 'FaceAlpha', 0.35);

% Combustion chamber free volume
fill([x_chamber_start x_conv_start x_conv_start x_chamber_start]*cm, ...
     [R_c R_c -R_c -R_c]*cm, ...
     [0.78 0.87 0.96], 'EdgeColor', [0.2 0.2 0.6], 'LineWidth', 1.2);

% Nozzle internal contour
fill([x_noz x_lower]*cm, [r_noz r_lower]*cm, ...
     [0.98 0.84 0.65], 'EdgeColor', [0.5 0.2 0.0], 'LineWidth', 1.4);

% Grain geometry
r0_grain   = best.d0 / 2;                  % [m]
seg_gap    = 0.005 * best.L0;              % [m]
L_seg_draw = best.L0 - seg_gap;            % [m]

for seg = 1:best.N
    x_offset = x_chamber_start + (seg-1) * best.L0;

    fill(([x_offset, x_offset+L_seg_draw, x_offset+L_seg_draw, x_offset])*cm, ...
         [R_c, R_c, r0_grain, r0_grain]*cm, ...
         [0.92 0.68 0.35], 'EdgeColor', [0.3 0.15 0.0], 'LineWidth', 0.8);

    fill(([x_offset, x_offset+L_seg_draw, x_offset+L_seg_draw, x_offset])*cm, ...
         [-r0_grain, -r0_grain, -R_c, -R_c]*cm, ...
         [0.92 0.68 0.35], 'EdgeColor', [0.3 0.15 0.0], 'LineWidth', 0.8);
end

% Core flow passage inside chamber
fill([x_chamber_start, x_conv_start, x_conv_start, x_chamber_start]*cm, ...
     [r0_grain, r0_grain, -r0_grain, -r0_grain]*cm, ...
     [0.95 0.95 0.95], 'EdgeColor', 'none');

% Centerline
plot([x_chamber_start x_exit]*cm, [0 0], 'k--', 'LineWidth', 0.8);

% Labels
text((x_chamber_start + x_conv_start)/2*cm, (R_c+t_wall)*cm*1.12, ...
     'Combustion Chamber', 'HorizontalAlignment', 'center', ...
     'FontSize', 9, 'Color', [0.1 0.1 0.5]);

text((x_conv_start + x_tan_c)/2*cm, r_e*cm*1.3, ...
     'Convergent', 'HorizontalAlignment', 'center', 'FontSize', 9);

text((x_tan_d + x_exit)/2*cm, r_e*cm*1.3, ...
     'Divergent', 'HorizontalAlignment', 'center', 'FontSize', 9);

text(0, -r_t*cm*2.2, 'Throat', 'HorizontalAlignment', 'center', 'FontSize', 9);

xlabel('x  [cm]');
ylabel('r  [cm]');
title(sprintf('Motor cross-section - %d BATES segments', best.N));

axis equal;
grid on;


%% =========================
%  QUASI-STEADY EVOLUTION PLOTS
%  =========================
dt_dx  = 1 ./ best.rb_x;                   % [s/m]
t_plot = cumtrapz(best.x_vec, dt_dx);      % [s]

clr_blue = [0.15 0.40 0.75];
clr_red  = [0.82 0.15 0.15];
clr_targ = [0.10 0.10 0.10];

figure('Name', 'Quasi-Steady Evolution', 'NumberTitle', 'off', 'Position', [150 150 820 540]);

subplot(2,1,1);
plot(best.x_vec*100, best.Kn_x, 'Color', clr_blue);
xlabel('Regression distance  x  [cm]');
ylabel('K_n  [-]');
title('K_n evolution vs web regression');
grid on;

subplot(2,1,2);
plot(t_plot, best.Pc_bar_x, 'Color', clr_red); hold on;
yline(Pc_target_bar, '--', 'Color', clr_targ, 'LineWidth', 1.2, ...
    'Label', sprintf('%g bar target', Pc_target_bar), 'LabelHorizontalAlignment', 'left');
xlabel('Equivalent burn time  [s]');
ylabel('Chamber pressure  [bar]');
title('Quasi-steady chamber pressure  (parabolic profile)');
ylim([0, max(best.Pc_bar_x)*1.15]);
grid on;

%% Saving

results_dir = fullfile(project_root, 'results');
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

figs = findobj('Type', 'figure');
for k = 1:numel(figs)
    fig = figs(k);

    % forza sfondo bianco figura
    set(fig, 'Color', 'w');

    % forza sfondo bianco di tutti gli assi presenti nella figura
    ax = findall(fig, 'Type', 'axes');
    set(ax, 'Color', 'w');

    fig_name = get(fig, 'Name');
    if isempty(fig_name)
        fig_name = sprintf('figure_%02d', k);
    end

    file_name = lower(regexprep(fig_name, '[^a-zA-Z0-9]+', '_'));
    file_name = regexprep(file_name, '^_|_$', '');

    exportgraphics(fig, fullfile(results_dir, [file_name '.png']), ...
        'Resolution', 200, 'BackgroundColor', 'white');
end