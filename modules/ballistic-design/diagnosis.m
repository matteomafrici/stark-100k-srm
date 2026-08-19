% =========================================================================
% FILE: diagnosis.m
% DESCRIPTION: SRM transient diagnosis driven by quasi-steady ballistics
%              Full geometry version (chamber + convergent + nozzle)
%              CEA engine: run_cea_transient (FAC + IAC + IAC-EQ conv)
%
% Requires workspace variables from design.m (see required_vars below).
% =========================================================================

%% =========================================================================
%% PATH SETUP
%% =========================================================================
project_root = fileparts(mfilename('fullpath'));
addpath(fullfile(project_root, 'engine'));

%% =========================================================================
%% CEA EXECUTABLE - auto-detect OS
%% =========================================================================
cea_work_dir = fullfile(project_root, 'cea');

if ispc
    cea_exec = fullfile(cea_work_dir, 'FCEA2.exe');
elseif ismac
    cea_exec = fullfile(cea_work_dir, 'FCEA2_mac');
else
    cea_exec = fullfile(cea_work_dir, 'FCEA2');
end

if ~isfile(cea_exec)
    error('CEA executable not found: %s', cea_exec);
end

if ismac
    [status, ~] = system(sprintf('xattr "%s" 2>/dev/null | grep -c quarantine', cea_exec));
    if status == 0
        msg = sprintf([ ...
            'macOS has blocked the CEA executable (Gatekeeper).\n\n' ...
            'To fix this, run the following command in Terminal:\n\n' ...
            '  xattr -d com.apple.quarantine "%s"\n\n' ...
            'Then run diagnosis.m again.'], cea_exec);
        errordlg(msg, 'CEA blocked by macOS Gatekeeper', 'modal');
        error('CEA executable blocked by Gatekeeper. See the dialog for instructions.');
    end
end

%% =========================================================================
%% REQUIRED VARIABLES FROM design.m
%% =========================================================================
required_vars = { ...
    'best', 'a', 'Inc_a', 'n', 'Inc_n', 'rho_propellant', 'A_t_target', ...
    'C_star', 'Pc_target_bar', 'Pc_lim', 'T_req', 'I_tot_req', ...
    't_b_target', 'ap_percentage', 'htpb_percentage', 'lambda_conical', ...
    'g0', 'M_mol', 'gamma', 'Isp_design', 'R_c', 'x_noz', 'r_noz', ...
    'clr_blue', 'clr_red', 'clr_targ'};

for k = 1:numel(required_vars)
    assert(exist(required_vars{k}, 'var') == 1, ...
        'diagnosis.m requires "%s" from design.m.', required_vars{k});
end

%% =========================================================================
%% SETTINGS
%% =========================================================================
dt               = 0.02;       % [s]  time step (suggested: 0.01-0.05)
T_react          = 298.15;     % [K]  reactant temperature
P_amb_bar        = 1.01325;    % [bar] ambient pressure
cea_tol          = 1e-5;       % [-]  c* convergence tolerance
cea_max_iter     = 12;         % [-]  max c* iterations per time step

nozzle_station_count = 20;     % [-]  nozzle sampling stations (9-20)

freeze_cfg       = struct();
freeze_cfg.mode  = 'throat';   % 'chamber'|'throat'|'area_ratio'|'exit'
freeze_cfg.AR    = 1.02;       % used only if mode='area_ratio' (supersonic A/At)

% Combustion efficiency corrections (default = 1.0: no correction)
% Set to measured values to get realistic performance estimates
eta_cstar        = 0.95;       % [-]  c* efficiency  (typ. 0.92-0.98 AP/HTPB)
eta_Cf           = 0.98;       % [-]  Cf efficiency  (typ. 0.96-0.99)

tau_ign = 0.03;   % [s] ignition time constant
tau_bo  = 0.05;   % [s] burnout time constant

run_monte_carlo  = true;
Nmc              = 2000;       % [-]  Monte Carlo samples

cea_inp_name     = 'current_run_bates_CEA';

%% =========================================================================
%% CHAMBER REFERENCE GEOMETRY
%% =========================================================================
L_grain_total   = best.N * best.L0;        % [m]
A_c_fac         = pi * R_c^2;             % [m^2]
CR_fac          = A_c_fac / A_t_target;   % [-]
V_chamber_total = A_c_fac * L_grain_total; % [m^3]

chamber_ref = struct( ...
    'R_c',             R_c, ...
    'A_c_fac',         A_c_fac, ...
    'CR_fac',          CR_fac, ...
    'L_grain_total',   L_grain_total, ...
    'V_chamber_total', V_chamber_total, ...
    'T_react',         T_react);

%% =========================================================================
%% NOZZLE GEOMETRY (x_noz, r_noz from design.m - convergent + divergent)
%% =========================================================================
x_nozzle_raw = x_noz(:);
r_nozzle_raw = r_noz(:);

[x_nozzle_full, idx_sort]   = sort(x_nozzle_raw);
r_nozzle_full               = r_nozzle_raw(idx_sort);
[x_nozzle_full, idx_unique] = unique(x_nozzle_full, 'stable');
r_nozzle_full               = r_nozzle_full(idx_unique);

A_nozzle_full  = pi * r_nozzle_full.^2;
AR_nozzle_full = A_nozzle_full / A_t_target;

[~, throat_idx_full] = min(A_nozzle_full);
x_throat = x_nozzle_full(throat_idx_full);

AR_sup_full = AR_nozzle_full(throat_idx_full:end);
x_sup_full  = x_nozzle_full(throat_idx_full:end);

if strcmpi(freeze_cfg.mode, 'area_ratio')
    assert(freeze_cfg.AR > 1.0, 'freeze_cfg.AR must be > 1.');
    assert(freeze_cfg.AR >= min(AR_sup_full) && freeze_cfg.AR <= max(AR_sup_full), ...
        'freeze_cfg.AR is outside the supersonic branch range.');
    x_freeze = interp1(AR_sup_full, x_sup_full, freeze_cfg.AR, 'linear');
else
    x_freeze = NaN;
end

% Sample nozzle stations
x_sample = linspace(x_nozzle_full(1), x_nozzle_full(end), nozzle_station_count).';
x_extra  = x_throat;

if strcmpi(freeze_cfg.mode, 'area_ratio')
    x_extra = [x_extra; x_freeze];
elseif strcmpi(freeze_cfg.mode, 'exit')
    x_extra = [x_extra; x_nozzle_full(end)];
end

x_nozzle = unique(sort([x_sample; x_extra]));
r_nozzle = interp1(x_nozzle_full, r_nozzle_full, x_nozzle, 'linear');
A_nozzle = pi * r_nozzle.^2;
AR_nozzle = A_nozzle / A_t_target;

[~, throat_idx]  = min(A_nozzle);
AR_nozzle(throat_idx) = 1.0;

nozzle_stations = struct('x', x_nozzle, 'AR', AR_nozzle, 'throat_idx', throat_idx);

% Full axial geometry for plots
x_chamber_axis  = linspace(-(L_grain_total + abs(min(x_nozzle))), min(x_nozzle), 120).';
r_chamber_axis  = R_c       * ones(size(x_chamber_axis));
A_chamber_axis  = A_c_fac   * ones(size(x_chamber_axis));
AR_chamber_axis = CR_fac    * ones(size(x_chamber_axis));

x_axis_full  = [x_chamber_axis; x_nozzle];
r_axis_full  = [r_chamber_axis; r_nozzle];
A_axis_full  = [A_chamber_axis; A_nozzle];
AR_axis_full = [AR_chamber_axis; AR_nozzle];

axial_geom = struct( ...
    'x_chamber',      x_chamber_axis, ...
    'x_nozzle',       x_nozzle, ...
    'x_full',         x_axis_full, ...
    'r_full',         r_axis_full, ...
    'A_full',         A_axis_full, ...
    'AR_full',        AR_axis_full, ...
    'A_nozzle',       A_nozzle, ...
    'AR_nozzle',      AR_nozzle, ...
    'x_throat',       x_throat, ...
    'throat_idx',     throat_idx, ...
    'x_nozzle_full',  x_nozzle_full, ...
    'AR_nozzle_full', AR_nozzle_full);

%% =========================================================================
%% STATE TEMPLATE
%% =========================================================================
state0 = struct( ...
    't',            0.0, ...
    'd_port',       best.d0, ...
    'L_seg',        best.L0, ...
    'A_port',       pi * best.d0^2 / 4, ...
    'G_port',       0.0, ...
    'A_b',          0.0, ...
    'V_prop',       0.0, ...
    'V_free',       0.0, ...
    'm_prop',       0.0, ...
    'Kn',           0.0, ...
    'Pc_bar',       P_amb_bar, ... % P_amb_bar, Pc_target_bar
    'rb',           0, ... % a * P_amb_bar^n / 1000, a * Pc_target_bar^n / 1000, 0
    'm_dot_gen',    0.0, ...
    'm_dot_noz',    0.0, ...
    'm_dot_real',   0.0, ...
    'thrust',       0.0, ...
    'thrust_real',  0.0, ...
    'thrust_cf',    0.0, ...
    'thrust_isp',   0.0, ...
    'c_eff_exit',   Isp_design * g0, ...
    'Tc',           T_react, ...
    'R',            8314.46 / M_mol, ...
    'Mw',           M_mol, ...
    'gamma',        gamma, ...
    'cstar',        C_star, ...
    'cstar_real',   C_star, ...
    'Cf',           NaN, ...
    'Cf_real',      NaN, ...
    'Isp',          Isp_design, ...
    'Isp_real',     Isp_design, ...
    'eta_cstar',    eta_cstar, ...
    'eta_Cf',       eta_Cf, ...
    'iter_cstar',   0);

n_alloc    = max(2000, ceil(1.5 * best.t_est / dt) + 500);
state_hist = repmat(state0, n_alloc, 1);
cea_hist   = cell(n_alloc, 1);

state  = state0;
k_hist = 0;

disp('--- STARTING TRANSIENT DIAGNOSIS ---');

%% =========================================================================
%% TIME MARCHING
%% =========================================================================
t_run_start = tic;
wb = waitbar(0, 't = 0.00 s | Pc = -- bar | F = -- N', ...
    'Name', 'SRM Transient');

while (state.d_port < best.D) && (state.L_seg > 0)

    % --- Geometry update --------------------------------------------------
    state.A_port  = pi * state.d_port^2 / 4;
    state.A_b     = best.N * (pi * state.d_port * state.L_seg ...
                    + 0.5 * pi * (best.D^2 - state.d_port^2));
    state.V_prop  = best.N * (pi/4) * (best.D^2 - state.d_port^2) * state.L_seg;
    state.V_prop  = max(state.V_prop, 0.0);
    state.V_free  = max(V_chamber_total - state.V_prop, 0.0);
    state.m_prop  = rho_propellant * state.V_prop;
    state.Kn      = state.A_b / A_t_target;

    % --- c* iteration (implicit Pc) ---------------------------------------
    cstar_iter = state.cstar;
    if ~isfinite(cstar_iter) || cstar_iter <= 0
        cstar_iter = C_star;
    end

    cea_data = [];
    Pc_iter  = NaN;

    for it = 1:cea_max_iter

        % Quasi-steady ballistic pressure (uses theoretical c*)
        Pc_iter = (state.Kn * a * rho_propellant * cstar_iter / 1e8)^(1 / (1 - n));
        
        if ~isreal(Pc_iter) || ~isfinite(Pc_iter) || Pc_iter <= 0
            error('Invalid chamber pressure during c* iteration (t=%.4f s).', state.t);
        end

        cea_data = run_cea_transient( ...
            Pc_iter,        ...
            T_react,        ...
            ap_percentage,  ...
            htpb_percentage,...
            CR_fac,         ...
            nozzle_stations,...
            freeze_cfg,     ...
            cea_work_dir,   ...
            cea_inp_name,   ...
            A_t_target,     ...
            lambda_conical, ...
            eta_cstar,      ...
            eta_Cf);

        rel_err_cstar = abs(cea_data.cstar - cstar_iter) / max(abs(cea_data.cstar), 1.0);
        cstar_iter    = cea_data.cstar;   % iterate on theoretical c*

        if rel_err_cstar < cea_tol
            break;
        end
    end

    % --- State update -----------------------------------------------------
    state.iter_cstar  = it;
            alpha_ign = 1 - exp(-state.t / tau_ign);%
            state.Pc_bar = P_amb_bar + (Pc_iter - P_amb_bar) * alpha_ign;%
    state.rb          = a * state.Pc_bar^n / 1000;   % [m/s]

    % Thermodynamic (theoretical)
    state.Tc          = cea_data.Tc;
    state.R           = cea_data.R;
    state.Mw          = cea_data.Mw;
    state.gamma       = cea_data.gamma;
    state.cstar       = cea_data.cstar;
    state.Cf          = cea_data.perf.Cf_exit;
    state.Isp         = cea_data.perf.Isp_exit;
    state.c_eff_exit  = cea_data.perf.c_eff_exit;

    % Real (efficiency-corrected)
    state.eta_cstar   = cea_data.eta_cstar;
    state.eta_Cf      = cea_data.eta_Cf;
    state.cstar_real  = cea_data.cstar_real;
    state.Cf_real     = cea_data.Cf_real;
    state.Isp_real    = cea_data.Isp_real;

    % Mass flows
    state.m_dot_gen   = rho_propellant * state.A_b * state.rb;
    state.m_dot_noz   = (state.Pc_bar * 1e5 * A_t_target) / state.cstar;
    state.m_dot_real  = cea_data.perf.m_dot_real;
    state.G_port      = state.m_dot_gen / max(state.A_port, 1e-12);

    % Thrust
    state.thrust_cf   = cea_data.perf.thrust_cf;
    state.thrust_isp  = cea_data.perf.thrust_isp;
    state.thrust_real = cea_data.perf.thrust_real;
    state.thrust      = state.thrust_real;   % PRIMARY OUTPUT = real thrust

    % --- Log --------------------------------------------------------------
    k_hist             = k_hist + 1;
    state_hist(k_hist) = state;
    cea_hist{k_hist}   = cea_data;

    frac_done = min(state.t / best.t_est, 1.0);
    waitbar(frac_done, wb, sprintf('t = %.2f s | Pc = %.2f bar | F = %.1f N', ...
        state.t, state.Pc_bar, state.thrust));

    % --- Advance geometry -------------------------------------------------
    delta_web     = 2 * state.rb * dt;
    state.d_port  = min(best.D,  state.d_port + delta_web);
    state.L_seg   = max(0.0,     state.L_seg  - delta_web);
    state.t       = state.t + dt;
end

% --- BURNOUT TAIL-OFF ----------------------------------------------------
Pc_bo = state.Pc_bar;
t_bo  = state.t;

while state.Pc_bar > 1.05 * P_amb_bar
    t_tail = state.t - t_bo;
    state.Pc_bar = P_amb_bar + (Pc_bo - P_amb_bar) * exp(-t_tail / tau_bo);

    state.rb         = 0.0;
    state.m_dot_gen  = 0.0;
    state.m_dot_noz  = (state.Pc_bar * 1e5 * A_t_target) / state.cstar;
    state.thrust_cf  = state.Cf_real * state.Pc_bar * 1e5 * A_t_target * lambda_conical;
    state.thrust_isp = state.Isp_real * state.m_dot_noz * g0;
    state.thrust_real = state.thrust_cf;
    state.thrust     = state.thrust_real;

    state.t = state.t + dt;

    k_hist = k_hist + 1;
    state_hist(k_hist) = state;
    cea_hist{k_hist}   = cea_data;
end

close(wb);
t_total = toc(t_run_start);
fprintf('Done in %.2f s | %d steps | %.1f ms/step\n', ...
    t_total, k_hist, 1000*t_total/max(k_hist,1));

state_hist = state_hist(1:k_hist);
cea_hist   = cea_hist(1:k_hist);

if isempty(state_hist)
    error('Transient simulation returned no valid time steps.');
end
fprintf('\n--- TRANSIENT DIAGNOSIS COMPLETED ---\n');

%% =========================================================================
%% HISTORY VECTORS
%% =========================================================================
time_hist         = [state_hist.t].';
d_port_hist       = [state_hist.d_port].';
L_seg_hist        = [state_hist.L_seg].';
A_port_hist       = [state_hist.A_port].';
G_port_hist        = [state_hist.G_port].';
A_b_hist          = [state_hist.A_b].';
V_prop_hist       = [state_hist.V_prop].';
V_free_hist       = [state_hist.V_free].';
m_prop_hist       = [state_hist.m_prop].';
Kn_hist           = [state_hist.Kn].';
Pressure_hist     = [state_hist.Pc_bar].';
rb_hist           = [state_hist.rb].';

% Mass flows
MassFlowGen_hist  = [state_hist.m_dot_gen].';
MassFlowNoz_hist  = [state_hist.m_dot_noz].';
MassFlowReal_hist = [state_hist.m_dot_real].';

% Thrust  - theoretical vs real
ThrustCF_hist     = [state_hist.thrust_cf].';    % theoretical Cf-based
ThrustIsp_hist    = [state_hist.thrust_isp].';   % theoretical Isp-based
ThrustReal_hist   = [state_hist.thrust_real].';  % efficiency-corrected (PRIMARY)
Thrust_hist       = ThrustReal_hist;             % alias for backward compat

% Thermodynamic
Tc_hist           = [state_hist.Tc].';
R_hist            = [state_hist.R].';
Mw_hist           = [state_hist.Mw].';
gamma_hist        = [state_hist.gamma].';

% Theoretical CEA
cstar_hist        = [state_hist.cstar].';
Cf_hist           = [state_hist.Cf].';
Isp_hist          = [state_hist.Isp].';
c_eff_exit_hist   = [state_hist.c_eff_exit].';

% Real (efficiency-corrected)
cstar_real_hist   = [state_hist.cstar_real].';
Cf_real_hist      = [state_hist.Cf_real].';
Isp_real_hist     = [state_hist.Isp_real].';
eta_cstar_hist    = [state_hist.eta_cstar].';
eta_Cf_hist       = [state_hist.eta_Cf].';

iter_hist         = [state_hist.iter_cstar].';

%% =========================================================================
%% INTEGRATED PERFORMANCE (based on real thrust)
%% =========================================================================
I_tot_real   = trapz(time_hist, ThrustReal_hist);
M_prop_used  = m_prop_hist(1) - m_prop_hist(end);
M_prop_int   = trapz(time_hist, MassFlowGen_hist);
tb_real      = state.t;
MEOP         = max(Pressure_hist);
G_port_max = max(G_port_hist);
G_port_avg = mean(G_port_hist);


fprintf('--- REAL INTEGRATED PERFORMANCE (eta_cstar=%.3f, eta_Cf=%.3f) ---\n', ...
    eta_cstar, eta_Cf);
fprintf('Total impulse  = %.2f N*s  (err = %+.2f %%)\n', ...
    I_tot_real, 100*(I_tot_real - I_tot_req)/I_tot_req);
fprintf('Burn time      = %.2f s    (err = %+.2f %%)\n', ...
    tb_real, 100*(tb_real - t_b_target)/t_b_target);
fprintf('Propellant used= %.2f kg\n', M_prop_used);
fprintf('Propellant int = %.2f kg\n', M_prop_int);
fprintf('Avg thrust     = %.2f N    (err = %+.2f %%)\n', ...
    mean(ThrustReal_hist), 100*(mean(ThrustReal_hist) - T_req)/T_req);
fprintf('Avg thrust(CF) = %.2f N  [theoretical]\n', mean(ThrustCF_hist));
fprintf('MEOP           = %.2f bar\n', MEOP);
fprintf('Avg Tc         = %.2f K\n',  mean(Tc_hist));
fprintf('Max G_port = %.2f kg/(m^2*s)\n', G_port_max);
fprintf('Avg G_port = %.2f kg/(m^2*s)\n', G_port_avg);

fprintf('Avg c*         = %.2f m/s  [theoretical]\n', mean(cstar_hist));
fprintf('Avg c* real    = %.2f m/s  [corrected]\n',   mean(cstar_real_hist));
fprintf('Avg Isp        = %.2f s    [theoretical]\n', mean(Isp_hist));
fprintf('Avg Isp real   = %.2f s    [corrected]\n',   mean(Isp_real_hist));
fprintf('Avg gamma      = %.4f\n', mean(gamma_hist));
fprintf('Avg R          = %.2f J/(kg*K)\n', mean(R_hist));
fprintf('Avg Mw         = %.4f kg/kmol\n', mean(Mw_hist));
fprintf('Avg CEA iter   = %.2f\n', mean(iter_hist));
fprintf('Pc_lim exceeded= %d\n\n', any(Pressure_hist > Pc_lim));


%% =========================================================================
%% WORKSPACE EXPORTS
%% =========================================================================
diagnosis_settings = struct( ...
    'dt',                   dt, ...
    'T_react',              T_react, ...
    'P_amb_bar',            P_amb_bar, ...
    'cea_tol',              cea_tol, ...
    'cea_max_iter',         cea_max_iter, ...
    'cea_work_dir',         cea_work_dir, ...
    'cea_inp_name',         cea_inp_name, ...
    'nozzle_station_count', nozzle_station_count, ...
    'freeze_cfg',           freeze_cfg, ...
    'eta_cstar',            eta_cstar, ...
    'eta_Cf',               eta_Cf, ...
    'run_monte_carlo',      run_monte_carlo, ...
    'Nmc',                  Nmc);

diagnosis_summary = struct( ...
    'I_tot_real',           I_tot_real, ...
    'tb_real',              tb_real, ...
    'M_prop_used',          M_prop_used, ...
    'M_prop_int',           M_prop_int, ...
    'MEOP',                 MEOP, ...
    'Pc_limit_exceeded',    any(Pressure_hist > Pc_lim), ...
    'avg_Isp_real',         mean(Isp_real_hist), ...
    'avg_cstar_real',       mean(cstar_real_hist), ...
    'avg_thrust_real',      mean(ThrustReal_hist));



%% =========================================================================
%% FIGURE 1 - TEMPORAL HISTORIES
%% =========================================================================

col_eb_blue   = [0.15 0.45 0.75];
col_eb_green  = [0.10 0.65 0.10];
col_eb_orange = [0.90 0.45 0.00];
col_eb_red    = [0.85 0.15 0.15];
col_eb_gray   = [0.50 0.50 0.50];

figure('Name', 'Diagnosis - Temporal histories', ...
    'NumberTitle', 'off', 'Position', [80 40 1200 860]);

subplot(4,2,1);
plot(time_hist, ThrustReal_hist/1e3,  'Color', clr_blue,          'LineWidth', 1.8); hold on;
plot(time_hist, ThrustCF_hist/1e3,    '--', 'Color', [0.5 0.5 0.5],'LineWidth', 1.1);
yline(T_req/1e3, '--', 'Color', col_eb_orange, 'LineWidth', 1.2);
xlabel('Time [s]'); ylabel('Thrust [kN]'); title('Thrust');
legend('Real (corrected)', 'Theoretical CF', 'Target', 'Location', 'southwest'); grid on;

subplot(4,2,2);
plot(time_hist, Pressure_hist, 'Color', clr_red, 'LineWidth', 1.5); hold on;
yline(Pc_target_bar, '--', 'Color', [0.5 0.5 0.5],      'LineWidth', 1.2);
yline(Pc_lim,        ':',  'Color', col_eb_orange ,'LineWidth', 1.2);
xlabel('Time [s]'); ylabel('P_c [bar]'); title('Chamber pressure');
legend('P_c', 'Target', 'Limit', 'Location',  'southwest'); grid on;

subplot(4,2,3);
plot(time_hist, MassFlowGen_hist,  'Color', [0.15 0.60 0.30], 'LineWidth', 1.5); hold on;
plot(time_hist, MassFlowReal_hist, '--', 'Color', [0.05 0.35 0.15], 'LineWidth', 1.3);
xlabel('Time [s]'); ylabel('Mass flow [kg/s]'); title('Mass flow rates');
legend('\dot{m}_{gen}', '\dot{m}_{real}', 'Location',  'southwest','Interpreter','latex'); grid on;

subplot(4,2,4);
plot(time_hist, rb_hist*1e3, 'Color', [0.75 0.40 0.10], 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('r_b [mm/s]'); title('Burning rate'); grid on;

subplot(4,2,5);
plot(time_hist, Tc_hist, 'Color', [0.85 0.45 0.10], 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('T_c [K]'); title('Chamber temperature'); grid on;

subplot(4,2,6);
plot(time_hist, cstar_real_hist,  'Color', [0.10 0.50 0.70], 'LineWidth', 1.8); hold on;
plot(time_hist, cstar_hist, '--', 'Color', [0.50 0.70 0.85], 'LineWidth', 1.1);
yline(C_star, '--', 'Color', [0.50 0.50 0.50], 'LineWidth', 1.0);
xlabel('Time [s]'); ylabel('c* [m/s]'); title('Characteristic velocity');
legend('Real', 'Theoretical', 'Design', 'Location',  'southeast'); grid on;

subplot(4,2,7);
plot(time_hist, Isp_real_hist,  'Color', [0.45 0.20 0.80], 'LineWidth', 1.8); hold on;
plot(time_hist, Isp_hist, '--', 'Color', [0.70 0.55 0.90], 'LineWidth', 1.1);
yline(Isp_design, '--', 'Color', [0.50 0.50 0.50], 'LineWidth', 1.0);
xlabel('Time [s]'); ylabel('I_{sp} [s]'); title('Specific impulse');
legend('Real', 'Theoretical', 'Design', 'Location',  'northeast'); grid on;

subplot(4,2,8)
yyaxis left
plot(time_hist, Kn_hist, 'Color', [0.55 0.55 0.15], 'LineWidth', 1.5); hold on
ylabel('Kn [-]')

yyaxis right
plot(time_hist, G_port_hist, '--', 'Color', [0.15 0.45 0.75], 'LineWidth', 1.4);
ylabel('G_{port} [kg/(m^2 s)]')

xlabel('Time [s]');
title('Kn and port mass flux');
legend('Kn', 'G_{port}', 'Location',  'northeast');
grid on

sgtitle(sprintf('BATES SRM transient diagnosis  (\\eta_{c*}=%.2f, \\eta_{Cf}=%.2f)', ...
    eta_cstar, eta_Cf), 'FontSize', 13, 'FontWeight', 'bold', 'Color', [0 0 0]);

%% =========================================================================
%% FIGURE 2 - GEOMETRY
%% =========================================================================
figure('Name', 'Diagnosis - Chamber and nozzle geometry', ...
    'NumberTitle', 'off', 'Position', [140 80 1120 760]);

subplot(2,2,1);
plot(x_axis_full,  r_axis_full, 'Color', clr_blue, 'LineWidth', 1.8); hold on;
plot(x_axis_full, -r_axis_full, 'Color', clr_blue, 'LineWidth', 1.8);
xline(x_throat, '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.1);
xlabel('x [m]'); ylabel('r [m]'); title('Internal radius profile');
axis equal; grid on;

subplot(2,2,2);
plot(x_axis_full, AR_axis_full, 'Color', clr_red, 'LineWidth', 1.8); hold on;
xline(x_throat, '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.1);
if strcmpi(freeze_cfg.mode, 'area_ratio')
    yline(freeze_cfg.AR, ':', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.1);
end
xlabel('x [m]'); ylabel('A(x)/A_t [-]'); title('Axial area ratio profile'); grid on;

subplot(2,2,3);
plot(time_hist, d_port_hist, 'Color', [0.15 0.45 0.75], 'LineWidth', 1.5); hold on;
plot(time_hist, L_seg_hist,  'Color', [0.75 0.30 0.20], 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('Geometry [m]');
title('Port diameter and segment length');
legend('d_{port}', 'L_{seg}', 'Location', 'best'); grid on;

subplot(2,2,4);
plot(time_hist, V_free_hist*1e3, 'Color', [0.10 0.55 0.70], 'LineWidth', 1.5); hold on;
plot(time_hist, A_port_hist,     '--', 'Color', [0.40 0.40 0.40], 'LineWidth', 1.3);
xlabel('Time [s]'); ylabel('V_f [L] / A_{port} [m^2]');
title('Chamber free volume and port area');
legend('V_f', 'A_{port}', 'Location', 'best'); grid on;

sgtitle('Chamber lumped data and nozzle axial geometry', ...
    'FontSize', 13, 'FontWeight', 'bold', 'Color', [0 0 0]);

%% =========================================================================
%% FIGURE 3 - CEA NOZZLE SNAPSHOTS
%% =========================================================================
snapshot_ids = unique([1, round(numel(cea_hist)/2), numel(cea_hist)]);
has_nozzle   = all(cellfun(@(c) isstruct(c) && isfield(c,'nozzle') ...
                   && ~isempty(c.nozzle), cea_hist(snapshot_ids)));

if has_nozzle
    figure('Name', 'Diagnosis - CEA nozzle snapshots', ...
        'NumberTitle', 'off', 'Position', [200 90 1200 780]);

    snap_colors = {clr_blue, [0.85 0.45 0.10], clr_red};
    snap_labels = {'t_{ini}', 't_{mid}', 't_{end}'};

    for s = 1:numel(snapshot_ids)
        i   = snapshot_ids(s);
        noz = cea_hist{i}.nozzle;
        lbl = snap_labels{s};
        col = snap_colors{s};

        subplot(2,2,1); hold on;
        if isfield(noz,'Mach')
            plot(noz.x, noz.Mach, 'Color', col, 'LineWidth', 1.4, 'DisplayName', lbl);
        end
        subplot(2,2,2); hold on;
        if isfield(noz,'P')
            plot(noz.x, noz.P, 'Color', col, 'LineWidth', 1.4, 'DisplayName', lbl);
        end
        subplot(2,2,3); hold on;
        if isfield(noz,'T')
            plot(noz.x, noz.T, 'Color', col, 'LineWidth', 1.4, 'DisplayName', lbl);
        end
        subplot(2,2,4); hold on;
        plot(noz.x, noz.AR, 'Color', col, 'LineWidth', 1.4, 'DisplayName', lbl);
    end

    subplot(2,2,1); xlabel('x [m]'); ylabel('Mach [-]'); title('Mach number');
    grid on; legend('Location','best');
    subplot(2,2,2); xlabel('x [m]'); ylabel('P [bar]'); title('Static pressure');
    grid on; legend('Location','best');
    subplot(2,2,3); xlabel('x [m]'); ylabel('T [K]'); title('Static temperature');
    grid on; legend('Location','best');
    subplot(2,2,4); xlabel('x [m]'); ylabel('A/A_t [-]'); title('Area ratio');
    grid on; legend('Location','best');

    sgtitle('CEA nozzle snapshots', 'FontSize', 13, 'FontWeight', 'bold','Color', [0 0 0]);
end

%% =========================================================================
%% FIGURE 4 - COMPACT PERFORMANCE (real vs theoretical)
%% =========================================================================
figure('Name', 'Real Motor Performance', 'NumberTitle', 'off', ...
    'Position', [200 100 900 620]);

subplot(2,2,1);
plot(time_hist, ThrustReal_hist/1e3, 'Color', clr_blue, 'LineWidth', 1.8); hold on;
plot(time_hist, ThrustCF_hist/1e3,   '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.1);
yline(T_req/1e3, '--', 'Color', col_eb_orange, 'LineWidth', 1.2);
xlabel('Time [s]'); ylabel('Thrust [kN]'); title('Thrust profile');
legend('Real', 'Theoretical', 'Target', 'Location', 'best');
xlim([0 max(time_hist)]); ylim([0, max(ThrustReal_hist/1000)*1.15]); grid on;

subplot(2,2,2);
plot(time_hist, Pressure_hist, 'Color', clr_red, 'LineWidth', 1.5); hold on;
yline(Pc_target_bar, '--', 'Color', clr_targ, 'LineWidth', 1.2);
xlabel('Time [s]'); ylabel('P_c [bar]'); title('Chamber pressure');
xlim([0 max(time_hist)]); ylim([0, max(Pressure_hist)*1.15]); grid on;

subplot(2,2,3);
plot(time_hist, Isp_real_hist, 'Color', [0.45 0.20 0.80], 'LineWidth', 1.8); hold on;
plot(time_hist, Isp_hist, '--', 'Color', [0.70 0.55 0.90], 'LineWidth', 1.1);
yline(Isp_design, '--', 'Color', [0.30 0.30 0.30], 'LineWidth', 1.0);
xlabel('Time [s]'); ylabel('I_{sp} [s]'); title('Specific impulse');
legend('Real', 'Theoretical', 'Design', 'Location', 'best'); grid on;

subplot(2,2,4);
plot(time_hist, cstar_real_hist, 'Color', [0.10 0.50 0.70], 'LineWidth', 1.8); hold on;
plot(time_hist, cstar_hist, '--', 'Color', [0.50 0.70 0.85], 'LineWidth', 1.1);
yline(C_star, '--', 'Color', [0.30 0.30 0.30], 'LineWidth', 1.0);
xlabel('Time [s]'); ylabel('c* [m/s]'); title('Characteristic velocity');
legend('Real', 'Theoretical', 'Design', 'Location', 'best'); grid on;

sgtitle(sprintf('BATES SRM - Real vs Theoretical  (\\eta_{c*}=%.2f, \\eta_{Cf}=%.2f)', ...
    eta_cstar, eta_Cf), 'FontSize', 13, 'FontWeight', 'bold', 'Color', [0 0 0]);

%% =========================================================================
%% MONTE CARLO UNCERTAINTY
%% =========================================================================
if run_monte_carlo
    rng(0);  % fixed seed for reproducible Monte Carlo
    Pc_meop_mc = zeros(Nmc, 1);
    tb_mc      = zeros(Nmc, 1);
    Itot_mc    = zeros(Nmc, 1);

    sigma_a     = Inc_a / 1.96;
    sigma_n     = Inc_n / 1.96;
    sigma_Cstar = 0.01 * C_star;
    sigma_rho   = 0.01 * rho_propellant;

    c_eff_mc  = mean(cstar_real_hist);  % use real c* for MC
    Cf_mc     = mean(Cf_real_hist);     % use real Cf for MC

    for k = 1:Nmc
        a_k     = max(1e-6, a     + sigma_a     * randn);
        n_k     =           n     + sigma_n     * randn;
        Cstar_k = max(1,    C_star+ sigma_Cstar * randn) * eta_cstar;
        rho_k   = max(1, rho_propellant + sigma_rho * randn);

        d_c = best.d0;
        L_c = best.L0;
        t   = 0;

        t_hist_k = [];
        T_hist_k = [];
        P_hist_k = [];

        while (d_c < best.D) && (L_c > 0)
            A_b_c = best.N * (pi*d_c*L_c + (pi/2)*(best.D^2 - d_c^2));
            Kn_c  = A_b_c / A_t_target;
            Pc_c  = (Kn_c * a_k * rho_k * Cstar_k / 1e8)^(1/(1-n_k));

            if ~isfinite(Pc_c) || Pc_c <= 0, break; end

            rb_c = (a_k * Pc_c^n_k) / 1000;
            if ~isfinite(rb_c) || rb_c <= 0, break; end

            T_c = Cf_mc * eta_Cf * Pc_c * 1e5 * A_t_target * lambda_conical;

            t_hist_k(end+1) = t;       
            T_hist_k(end+1) = T_c;     
            P_hist_k(end+1) = Pc_c;    

            d_c = d_c + 2*rb_c*dt;
            L_c = L_c - 2*rb_c*dt;
            t   = t   + dt;
        end

        if numel(t_hist_k) >= 2
            Pc_meop_mc(k) = max(P_hist_k);
            tb_mc(k)      = t;
            Itot_mc(k)    = trapz(t_hist_k, T_hist_k);
        else
            Pc_meop_mc(k) = NaN;
            tb_mc(k)      = NaN;
            Itot_mc(k)    = NaN;
        end
    end

    valid_mc   = isfinite(Pc_meop_mc) & isfinite(tb_mc) & isfinite(Itot_mc);
    Pc_meop_mc = Pc_meop_mc(valid_mc);
    tb_mc      = tb_mc(valid_mc);
    Itot_mc    = Itot_mc(valid_mc);

    fprintf('--- MONTE CARLO UNCERTAINTY (%d samples) ---\n', Nmc);
    fprintf('Valid samples  = %d / %d\n',        numel(Pc_meop_mc), Nmc);
    fprintf('MEOP           = %.2f +/- %.2f bar\n', mean(Pc_meop_mc), std(Pc_meop_mc));
    fprintf('Burn time      = %.2f +/- %.2f s\n',   mean(tb_mc),      std(tb_mc));
    fprintf('Total impulse  = %.2f +/- %.2f N*s\n\n', mean(Itot_mc),  std(Itot_mc));

    figure('Name', 'Monte Carlo Uncertainty', 'NumberTitle', 'off', ...
        'Position', [250 100 900 620]);

    clr_h1   = [0.20 0.55 0.85];
    clr_h2   = [0.85 0.30 0.20];
    clr_conv = [0.15 0.60 0.35];

    subplot(2,2,1);
    histogram(Pc_meop_mc, 30, 'FaceColor', clr_h1, 'EdgeColor', 'w', 'FaceAlpha', 0.8);
    xline(mean(Pc_meop_mc), '--k', 'LineWidth', 1.3);
    xlabel('MEOP [bar]'); ylabel('Count'); title('MEOP distribution'); grid on;

    subplot(2,2,2);
    histogram(tb_mc, 30, 'FaceColor', clr_h2, 'EdgeColor', 'w', 'FaceAlpha', 0.8);
    xline(mean(tb_mc), '--k', 'LineWidth', 1.3);
    xlabel('Burn time [s]'); ylabel('Count'); title('Burn time distribution'); grid on;

    subplot(2,2,3);
    running_mean_meop = arrayfun(@(i) mean(Pc_meop_mc(1:i)), 1:numel(Pc_meop_mc));
    plot(1:numel(Pc_meop_mc), running_mean_meop, 'Color', clr_conv);
    yline(mean(Pc_meop_mc), '--k', 'LineWidth', 1.0);
    xlabel('Samples'); ylabel('Running mean [bar]'); title('MEOP convergence'); grid on;

    subplot(2,2,4);
    if numel(Pc_meop_mc) >= 2
        running_std_meop = arrayfun(@(i) std(Pc_meop_mc(1:i)), 2:numel(Pc_meop_mc));
        plot(2:numel(Pc_meop_mc), running_std_meop, 'Color', [0.65 0.35 0.75]);
    end
    xlabel('Samples'); ylabel('Running std [bar]'); title('MEOP std convergence'); grid on;

    sgtitle('Monte Carlo uncertainty', 'FontSize', 13, 'FontWeight', 'bold', 'Color', [0 0 0]);
end

%% =========================================================================
%% SAVE .mat + FIGURES
%% =========================================================================
results_dir = fullfile(project_root, 'results');
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

% --- .mat -----------------------------------------------------------------
save_path = fullfile(results_dir, 'diagnosis_results.mat');
save(save_path, ...
    'time_hist', 'Pressure_hist', 'ThrustReal_hist', 'ThrustCF_hist', 'ThrustIsp_hist', ...
    'cstar_hist', 'cstar_real_hist', 'Isp_hist', 'Isp_real_hist', ...
    'Cf_hist', 'Cf_real_hist', 'eta_cstar_hist', 'eta_Cf_hist', ...
    'Tc_hist', 'gamma_hist', 'Mw_hist', 'R_hist', ...
    'MassFlowGen_hist', 'MassFlowNoz_hist', 'MassFlowReal_hist', ...
    'rb_hist', 'Kn_hist', 'A_b_hist', 'd_port_hist', 'L_seg_hist', ...
    'V_free_hist', 'A_port_hist', 'm_prop_hist', ...
    'state_hist', 'cea_hist', ...
    'diagnosis_settings', 'diagnosis_summary', 'chamber_ref', 'axial_geom');
fprintf('Results saved to: %s\n', save_path);

% --- Figures --------------------------------------------------------------
fig_specs = { ...
    'Diagnosis - Temporal histories',          'temporal_histories'; ...
    'Diagnosis - Chamber and nozzle geometry', 'chamber_nozzle_geometry'; ...
    'Diagnosis - CEA nozzle snapshots',        'cea_nozzle_snapshots'; ...
    'Real Motor Performance',                  'real_motor_performance'; ...
    'Monte Carlo Uncertainty',                 'monte_carlo_uncertainty'};

saved_figs = 0;
for fi = 1:size(fig_specs, 1)
    fig_h = findobj('Type', 'figure', 'Name', fig_specs{fi, 1});
    if isempty(fig_h)
        continue;
    end
    fig_h = fig_h(1);   % to avoid duplications
    set(fig_h, 'Color', 'w');
    set(findall(fig_h, 'Type', 'axes'), 'Color', 'w');
    fig_path = fullfile(results_dir, [fig_specs{fi, 2}, '.png']);
    exportgraphics(fig_h, fig_path, 'Resolution', 200);
    fprintf('Figure saved : %s\n', fig_path);
    saved_figs = saved_figs + 1;
end
fprintf('Total figures saved: %d\n\n', saved_figs);