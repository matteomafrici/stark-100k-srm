%% EROSIVE BURNING ANALYSIS - UNIFIED CODE
%
% Erosive-burning regime analysis for the BATES motor.

% STRUCTURE:
%   Sec. 0  - Input data (single definition)
%   Sec. 1  - Port geometry and area-Mach curve
%   Sec. 2  - Rogers diagnostics (velocity-based + mass flux-based)
%   Sec. 3  - Parametric study Ap/Ath (transition thresholds)
%   Sec. 4  - Figures

% REFERENCES:
%   [1] Rogers, C.E. (2005) - "Erosive Burning Design Criteria for High
%       Power and Experimental/Amateur Rocket Motors"
%       High Power Rocketry Magazine, January 2005.
%   [2] Sutton, G.P. & Biblarz, O. (2017) - "Rocket Propulsion Elements"
%       9th Edition, Wiley. Ch. 13.
%   [3] Anderson, J.D. (2003) - "Modern Compressible Flow", 3rd Ed.,
%       McGraw-Hill. Ch. 5 (isentropic area-Mach relation).
%   [4] Lenoir & Robillard (1957), Proc. 6th IAC, pp. 663-672.
%   [5] Barrere et al. (1960), Rocket Propulsion, Elsevier.
% =========================================================================

clc; close all;
fprintf('============================================================\n');
fprintf('   EROSIVE BURNING ANALYSIS - AP/HTPB BATES MOTOR\n');
fprintf('============================================================\n\n');

% --- Figure styling: light background, dark text (matching design.m) ---
set(0, 'DefaultFigureColor',    [0.97 0.97 0.97], ...  % figure background: near-white
       'DefaultAxesColor',      [1.00 1.00 1.00], ...  % axes background: white
       'DefaultAxesXColor',     [0.15 0.15 0.15], ...  % X tick/label: dark grey
       'DefaultAxesYColor',     [0.15 0.15 0.15], ...  % Y tick/label: dark grey
       'DefaultAxesZColor',     [0.15 0.15 0.15], ...  % Z tick/label
       'DefaultAxesGridColor',  [0.15 0.15 0.15], ...  % grid lines
       'DefaultTextColor',      [0.15 0.15 0.15]);     % all text objects

%% =========================================================================
% SECTION 0 - INPUT DATA
%
% Source: BATES optimization output and Vieille fit (see design report)
% =========================================================================

% --- Vieille law: r = a * Pc^n  [mm/s, bar] ---
a_vieille  = 1.617743;      % [mm/s / bar^n]  - experimental fit
n_vieille  = 0.381551;      % [-]             - pressure exponent

% --- Non-aluminized AP/HTPB propellant ---
rho_prop   = 1685.22;       % [kg/m^3]  - propellant density
gamma_core = 1.1640;        % [-]       - specific heat ratio
                             %             (CEA, AP/HTPB, Pc=65 bar)

% --- Operating conditions ---
Pc_bar     = 65.0;          % [bar]     - chamber pressure (design target)

% --- Grain geometry (best design, N=1 BATES segment) ---
d_port     = 0.4738;        % [m]       - initial port diameter
L_grain    = 1.5441;        % [m]       - grain length
A_b0       = 2.298;        % [m^2]     - initial burning area

% --- Nozzle geometry ---
A_t        = 0.010306;      % [m^2]     - throat area

% --- Lenoir-Robillard constant ---
% Imperial: k = 0.0288 [in/s / (lb/s/in^2)^0.8]  (Barrere 1960 [5])
% SI conversion: k_SI = k_imp * 0.0254 / (703.07)^0.8
k_LR_imp   = 0.0288;
conv_G     = 703.07;                                    % [kg/m^2/s per lb/s/in^2]
k_LR       = k_LR_imp * 0.0254 / conv_G^0.8;           % [m/s / (kg/m^2/s)^0.8]

% --- Rogers thresholds [1], Table 1 (velocity-based) ---
Ma_threshold_nonerr = 0.50;   % [-]  non-erosive
Ma_threshold_maxerr = 0.70;   % [-]  max erosivity

% --- Rogers thresholds [1] at Pc = 65 bar (mass flux-based), as reported
%     in the design report (Erosive Burning Analysis section):
%       non-erosive:              G <= 703  kg/(m^2/s)
%       max recommended erosivity: G ~= 1757 kg/(m^2/s)
Pc_psia        = Pc_bar * 14.5038;    % [psia]  - for reporting only
G_thr_low_SI   = 703.0;               % [kg/m^2/s]  non-erosive limit
G_thr_hi_SI    = 1757.0;              % [kg/m^2/s]  max recommended erosivity

% --- Plot colors (unified palette) ---
col_blue   = [0.15 0.45 0.75];
col_green  = [0.10 0.65 0.10];
col_orange = [0.90 0.45 0.00];
col_red    = [0.85 0.15 0.15];
col_gray   = [0.50 0.50 0.50];

fprintf('--- SECTION 0: INPUT DATA ---\n');
fprintf('  Vieille:  a = %.6f [mm/s / bar^n],  n = %.6f [-]\n', a_vieille, n_vieille);
fprintf('  rho_prop  = %.2f kg/m^3\n',   rho_prop);
fprintf('  gamma     = %.4f [-]\n',       gamma_core);
fprintf('  Pc        = %.1f bar\n',       Pc_bar);
fprintf('  d_port    = %.4f m\n',         d_port);
fprintf('  L_grain   = %.4f m\n',         L_grain);
fprintf('  A_b0      = %.4f m^2\n',       A_b0);
fprintf('  A_t       = %.6f m^2\n',       A_t);
fprintf('  k_LR (SI) = %.6e m/s/(kg/m^2/s)^0.8\n\n', k_LR);


%% =========================================================================
% SECTION 1 - PORT GEOMETRY AND AREA-MACH CURVE
%
% Isentropic area-Mach relation (subsonic flow) [3]:
%   Ap/At = (1/Ma)*[(2/(g+1))*(1+(g-1)/2*Ma^2)]^((g+1)/(2(g-1)))
%
% The core Mach number is maximum at the aft end: worst case for
% erosive burning. Numerical inversion over a fine grid (200k points).
% =========================================================================

r_port       = d_port / 2;
A_port       = pi * r_port^2;
ratio_Ap_Ath = A_port / A_t;

% Area-Mach curve (inversion grid)
Ma_vec  = linspace(0.001, 0.999, 200000);
exp_val = (gamma_core + 1) / (2 * (gamma_core - 1));
rhs_AM  = (1 ./ Ma_vec) .* ...
          ((2/(gamma_core+1)) .* ...
          (1 + (gamma_core-1)/2 .* Ma_vec.^2)).^exp_val;

[~, idx_Ma] = min(abs(rhs_AM - ratio_Ap_Ath));
Ma_core     = Ma_vec(idx_Ma);

fprintf('--- SECTION 1: PORT GEOMETRY AND CORE MACH NUMBER ---\n');
fprintf('  Port diameter   d0    = %.4f m\n',     d_port);
fprintf('  Port area       Ap    = %.6f m^2\n',   A_port);
fprintf('  Throat area     At    = %.6f m^2\n',   A_t);
fprintf('  Ap / At               = %.4f [-]\n',   ratio_Ap_Ath);
fprintf('  Core Mach (aft end)   = %.6f [-]\n\n', Ma_core);

%% =========================================================================
% SECTION 2 - ROGERS DIAGNOSTICS
%
% 2a) Velocity-based criterion: Ma_core vs thresholds [1]
% 2b) Mass flux-based criterion: G vs thresholds [1]
%
% Burn rate and mass flow at ignition (no erosion):
%   rb0 = a * Pc^n  [mm/s -> m/s]
%   m_dot = rho_prop * rb0 * A_b0
%   G = m_dot / A_port  [kg/m^2/s]
% =========================================================================

rb0_ms    = 7.9548 / 1000;   % [m/s]
m_dot_ign = rho_prop * rb0_ms * A_b0;               % [kg/s]
G_design  = m_dot_ign / A_port;                      % [kg/m^2/s]

fprintf('--- SECTION 2a: ROGERS CRITERION - VELOCITY-BASED ---\n');
fprintf('  Non-erosive:    Ma <= %.2f\n',   Ma_threshold_nonerr);
fprintf('  Max erosivity:  Ma  = %.2f\n',   Ma_threshold_maxerr);
fprintf('  Current design: Ma  = %.6f\n',   Ma_core);
if Ma_core <= Ma_threshold_nonerr
    fprintf('  >> REGIME: NON-EROSIVE (velocity-based)\n\n');
elseif Ma_core <= Ma_threshold_maxerr
    fprintf('  >> REGIME: MODERATE EROSIVITY (velocity-based)\n\n');
else
    fprintf('  >> WARNING: BEYOND MAX RECOMMENDED EROSIVITY\n\n');
end

fprintf('--- SECTION 2b: ROGERS CRITERION - MASS FLUX-BASED ---\n');
fprintf('  rb0 (Pc=%.0f bar)     = %.4f mm/s  = %.6e m/s\n', Pc_bar, rb0_ms*1000, rb0_ms);
fprintf('  m_dot (ignition)     = %.4f kg/s\n',   m_dot_ign);
fprintf('  G design (aft end)   = %.4f kg/m^2/s\n', G_design);
fprintf('  Rogers thresholds (Pc = %.0f psia = %.0f bar):\n', Pc_psia, Pc_bar);
fprintf('    Non-erosive:     G <= %.2f kg/m^2/s\n', G_thr_low_SI);
fprintf('    Max erosivity:   G  = %.2f kg/m^2/s\n', G_thr_hi_SI);
if G_design <= G_thr_low_SI
    fprintf('  >> REGIME: NON-EROSIVE (mass flux-based)\n\n');
elseif G_design <= G_thr_hi_SI
    fprintf('  >> REGIME: MODERATE EROSIVITY (mass flux-based)\n\n');
else
    fprintf('  >> WARNING: BEYOND MAX RECOMMENDED EROSIVITY\n\n');
end

fprintf('=== DIAGNOSTIC SUMMARY ===\n');
fprintf('  Ap/At = %.2f >> 1.36   |  Ma = %.4f << %.2f\n', ...
        ratio_Ap_Ath, Ma_core, Ma_threshold_nonerr);
fprintf('  G     = %.2f << %.2f kg/m^2/s\n', G_design, G_thr_low_SI);
fprintf('  CONCLUSION: the design is largely NON-EROSIVE.\n');
fprintf('============================================================\n\n');

%% =========================================================================
% SECTION 3 - PARAMETRIC STUDY: Ap/At vs EROSIVE REGIME
%
% Vary Ap/At from 1.05 to 20 keeping constant:
%   - aft-end mass flow m_dot_ign (same rb and Ab0)
%   - throat area A_t
% For each value, Ma_core and G are computed and compared with the
% Rogers thresholds.
%
% Goal: identify the geometric threshold below which the design would
% enter the erosive regime and quantify the current margin.
% =========================================================================

fprintf('--- SECTION 3: PARAMETRIC STUDY Ap/At ---\n\n');

N_param   = 1000;
ratio_vec = linspace(1.05, 20, N_param);
Ma_param  = zeros(1, N_param);
G_param   = zeros(1, N_param);       % directly in SI [kg/m^2/s]

for ii = 1:N_param
    % Mach from the area-Mach curve
    [~, idx_i]   = min(abs(rhs_AM - ratio_vec(ii)));
    Ma_param(ii) = Ma_vec(idx_i);

    % Corresponding port area (constant A_t)
    A_port_i = ratio_vec(ii) * A_t;

    % Mass flux in SI [kg/m^2/s]
    G_param(ii) = m_dot_ign / A_port_i;
end

% Regime transition thresholds
idx_G_nonerr  = find(G_param  <= G_thr_low_SI,         1, 'first');
idx_G_maxerr  = find(G_param  <= G_thr_hi_SI,          1, 'first');
idx_Ma_nonerr = find(Ma_param <= Ma_threshold_nonerr,   1, 'first');

fprintf('  Transition thresholds:\n');
if ~isempty(idx_G_nonerr)
    fprintf('  [G]  Non-erosive for Ap/At >= %.3f  (G <= %.2f kg/m^2/s)\n', ...
        ratio_vec(idx_G_nonerr), G_thr_low_SI);
end
if ~isempty(idx_G_maxerr)
    fprintf('  [G]  Below max erosivity for Ap/At >= %.3f  (G <= %.2f kg/m^2/s)\n', ...
        ratio_vec(idx_G_maxerr), G_thr_hi_SI);
end
if ~isempty(idx_Ma_nonerr)
    fprintf('  [Ma] Non-erosive for Ap/At >= %.3f  (Ma <= %.2f)\n\n', ...
        ratio_vec(idx_Ma_nonerr), Ma_threshold_nonerr);
end

%% =========================================================================
% SECTION 4 - FIGURES
%
% FIGURE 1: Parametric diagnostics (2 subplots)
%   (a) Ma_core vs Ap/At
%   (b) G vs Ap/At  [kg/m^2/s]
% =========================================================================


% ===================== FIGURE 1 =====================
fig = figure('Name', 'Fig.1 - Parametric Study', ...
       'NumberTitle', 'off', 'Position', [80 200 1200 480]);

% --- Subplot (a): Ma_core vs Ap/At ---
subplot(1, 2, 1);
plot(ratio_vec, Ma_param, '-', 'Color', col_blue, 'LineWidth', 1.8);
hold on;
yline(Ma_threshold_nonerr, '--', 'Color', col_green, 'LineWidth', 1.4, ...
      'DisplayName', sprintf('Non-erosive (Ma = %.2f)', Ma_threshold_nonerr));
yline(Ma_threshold_maxerr, '--', 'Color', col_red, 'LineWidth', 1.4, ...
      'DisplayName', sprintf('Max erosivity (Ma = %.2f)', Ma_threshold_maxerr));
xline(ratio_Ap_Ath, ':', 'Color', col_gray, 'LineWidth', 1.3, ...
      'DisplayName', sprintf('Design A_p/A_t = %.2f', ratio_Ap_Ath));
scatter(ratio_Ap_Ath, Ma_core, 70, col_red, 'filled', ...
        'DisplayName', sprintf('Design: Ma = %.4f', Ma_core));
hold off;
xlabel('A_p / A_t  [-]', 'FontSize', 11);
ylabel('Ma_{core}  [-]', 'FontSize', 11);
title('(a) Velocity-based criterion: Ma_{core} vs A_p/A_t', 'FontSize', 11);
legend('Location', 'northeast', 'FontSize', 8);
grid on;
ylim([0, max(Ma_param)*1.15]);
xlim([ratio_vec(1), ratio_vec(end)]);

% --- Subplot (b): G vs Ap/At [kg/m^2/s] ---
subplot(1, 2, 2);
plot(ratio_vec, G_param, '-', 'Color', col_blue, 'LineWidth', 1.8);
hold on;
yline(G_thr_low_SI, '--', 'Color', col_green, 'LineWidth', 1.4, ...
      'DisplayName', sprintf('Non-erosive (G = %.1f kg/m^2/s)', G_thr_low_SI));
yline(G_thr_hi_SI, '--', 'Color', col_red, 'LineWidth', 1.4, ...
      'DisplayName', sprintf('Max erosivity (G = %.1f kg/m^2/s)', G_thr_hi_SI));
xline(ratio_Ap_Ath, ':', 'Color', col_gray, 'LineWidth', 1.3, ...
      'DisplayName', sprintf('Design A_p/A_t = %.2f', ratio_Ap_Ath));
scatter(ratio_Ap_Ath, G_design, 70, col_red, 'filled', ...
        'DisplayName', sprintf('Design: G = %.2f kg/m^2/s', G_design));
hold off;
xlabel('A_p / A_t  [-]', 'FontSize', 11);
ylabel('Mass flux  G  [kg/m^2/s]', 'FontSize', 11);
title('(b) Mass flux-based criterion: G vs A_p/A_t', 'FontSize', 11);
legend('Location', 'northeast', 'FontSize', 8);
grid on;
ylim([0, max(G_param)*1.15]);
xlim([ratio_vec(1), ratio_vec(end)]);

sgtitle(sprintf('Parametric Study - Rogers Criteria (Pc = %.0f bar)', Pc_bar), ...
        'FontSize', 13, 'FontWeight', 'bold');

% --- Save figure ----------------------------------------------------------
results_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end
fig_path = fullfile(results_dir, 'erosive_burning_parametric.png');
set(fig, 'Color', 'w');
exportgraphics(fig, fig_path, 'Resolution', 200, 'BackgroundColor', 'white');
fprintf('Figure saved: %s\n', fig_path);