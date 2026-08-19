function cea_out = run_cea_transient( ...
        Pc_bar, T_react_K, ap_frac, htpb_frac, CR_fac, nozzle_stations, ...
        freeze_cfg, work_dir, inp_name, A_t_m2, lambda_nozzle, eta_cstar, eta_Cf)
% RUN_CEA_TRANSIENT
%   FAC (Finite Area Combustor) equilibrium call for chamber thermodynamics
%   + IAC (Infinite Area Combustor) rocket call for nozzle performance.
%
%   INPUTS
%     Pc_bar           Chamber pressure [bar]
%     T_react_K        Reactant temperature [K]
%     ap_frac          AP mass fraction [-]
%     htpb_frac        HTPB mass fraction [-]
%     CR_fac           Chamber-to-throat area ratio Ac/At [-]
%     nozzle_stations  Struct with fields: x, AR, throat_idx
%                        x          : axial coordinates [m], full nozzle
%                        AR         : A(x)/At [-], convergent + divergent
%                        throat_idx : index of throat in x/AR arrays
%     freeze_cfg       Struct with fields: mode, AR
%                        mode : 'chamber' | 'throat' | 'area_ratio' | 'exit'
%                        AR   : freeze area ratio (only for mode='area_ratio')
%     work_dir         Path to FCEA2 executable directory
%     inp_name         Base name for .inp/.out files (unique per call)
%     A_t_m2           Throat area [m^2] - optional, enables thrust/mdot
%     lambda_nozzle    Nozzle efficiency factor [-] - optional, default 1.0
%     eta_cstar        Combustion efficiency on c* [-] - optional, default 1.0
%     eta_Cf           Thrust coefficient efficiency [-] - optional, default 1.0
%
%   OUTPUTS  (cea_out struct)
%     .Tc, .Pc, .Mw, .gamma, .cp, .R, .cstar   - chamber scalars (theoretical)
%     .Pt, .Tt, .Macht                          - throat scalars
%     .Pe, .Te, .Mache, .Cf, .Isp, .c_eff       - exit scalars (theoretical)
%     .cstar_real, .Cf_real, .Isp_real           - corrected scalars
%     .perf.{cstar_fac, Cf_exit, Isp_exit, c_eff_exit,
%            m_dot, thrust_cf, thrust_isp,
%            cstar_real, Cf_real, Isp_real,
%            m_dot_real, thrust_real, Isp_si_real}
%     .nozzle    - spatial profile struct (x, AR, T, P, Mach, Cf, Isp, ...)
%     .fac, .iac - raw parsed structs from CEA output

if nargin < 10 || isempty(A_t_m2)
    A_t_m2 = [];
end
if nargin < 11 || isempty(lambda_nozzle)
    lambda_nozzle = 1.0;
end
if nargin < 12 || isempty(eta_cstar)
    eta_cstar = 1.0;
end
if nargin < 13 || isempty(eta_Cf)
    eta_Cf = 1.0;
end

if ispc
    FCEA2 = fullfile(work_dir, 'FCEA2.exe');
elseif ismac
    FCEA2 = fullfile(work_dir, 'FCEA2_mac');
else
    FCEA2 = fullfile(work_dir, 'FCEA2');
end
if ~isfile(FCEA2)
    error('CEA executable not found: %s', FCEA2);
end

nozzle = normalize_nozzle_stations(nozzle_stations);
freeze = normalize_freeze_cfg(freeze_cfg, nozzle);

fac_case = [inp_name, '_fac'];
iac_case = [inp_name, '_iac'];

% ------------------------------------------------------------------
% FAC call - equilibrium, Finite Area Combustor
% Provides: Tc, Pc, Mw, gamma, cp, R, cstar, Pt, Tt
% ------------------------------------------------------------------
fac_map = build_fac_station_map(nozzle);
fac_inp = build_fac_input(Pc_bar, T_react_K, ap_frac, htpb_frac, CR_fac, fac_map, fac_case);
fac_files = run_fcea_case(work_dir, FCEA2, fac_case, fac_inp);
fac = parse_cea_result(fac_files.out, 'fac', fac_map);

% ------------------------------------------------------------------
% IAC call - frozen (or equilibrium if mode='exit')
% Provides: nozzle supersonic profile, Cf, Isp at exit
% ------------------------------------------------------------------
[iac_map, iac_mode, iac_nfz, freeze] = build_iac_station_map(nozzle, freeze);
iac_inp = build_iac_input(Pc_bar, T_react_K, ap_frac, htpb_frac, iac_map, iac_case, iac_mode, iac_nfz);
iac_files = run_fcea_case(work_dir, FCEA2, iac_case, iac_inp);
iac = parse_cea_result(iac_files.out, 'iac', iac_map);

% ------------------------------------------------------------------
% IAC-EQ call - equilibrium, subsonic convergent stations only
% Provides: T, P, Mach, gamma along the convergent (sub branch)
% ------------------------------------------------------------------
conv_map = build_conv_station_map(nozzle);
if ~isempty(conv_map.sub_AR)
    conv_case  = [inp_name, '_conv'];
    conv_inp   = build_iac_eq_input( ...
        Pc_bar, T_react_K, ap_frac, htpb_frac, conv_map, conv_case);
    conv_files = run_fcea_case(work_dir, FCEA2, conv_case, conv_inp);
    conv_eq    = parse_cea_result(conv_files.out, 'iac_eq', conv_map);
else
    conv_eq    = [];
    conv_files = struct('inp','','out','','plt','','csv','');
end

% ------------------------------------------------------------------
% Merge profiles into unified nozzle spatial struct
% ------------------------------------------------------------------
nozzle_profile = merge_nozzle_profiles(nozzle, fac, iac, conv_eq, freeze);

% ------------------------------------------------------------------
% Build output struct - scalars
% ------------------------------------------------------------------
cea_out = struct();
cea_out.fac    = fac;
cea_out.iac    = iac;
cea_out.freeze = freeze;
cea_out.nozzle = nozzle_profile;

cea_out.Tc    = fac.chamber.T;
cea_out.Pc    = fac.chamber.P;
cea_out.Mw    = fac.chamber.Mw;
cea_out.gamma = fac.chamber.gamma;
cea_out.cp    = fac.chamber.cp;
cea_out.R     = fac.chamber.R;
cea_out.cstar = fac.chamber.cstar;

cea_out.Pt    = nozzle_profile.throat.P;
cea_out.Tt    = nozzle_profile.throat.T;
cea_out.Macht = nozzle_profile.throat.Mach;

cea_out.Pe    = iac.exit.P;
cea_out.Te    = iac.exit.T;
cea_out.Mache = iac.exit.Mach;
cea_out.Cf    = iac.exit.Cf;
cea_out.Isp   = iac.exit.Isp;
cea_out.c_eff = iac.exit.c_eff;

% ------------------------------------------------------------------
% Performance struct - theoretical CEA
% ------------------------------------------------------------------
cea_out.perf = struct();
cea_out.perf.cstar_fac   = fac.chamber.cstar;
cea_out.perf.Cf_exit     = iac.exit.Cf;
cea_out.perf.Isp_exit    = iac.exit.Isp;
cea_out.perf.c_eff_exit  = iac.exit.c_eff;

if ~isempty(A_t_m2)
    m_dot = (Pc_bar * 1e5 * A_t_m2) / max(fac.chamber.cstar, 1.0);
    cea_out.perf.m_dot      = m_dot;
    cea_out.perf.thrust_cf  = iac.exit.Cf * Pc_bar * 1e5 * A_t_m2 * lambda_nozzle;
    cea_out.perf.thrust_isp = m_dot * iac.exit.Isp * 9.80665 * lambda_nozzle;
else
    cea_out.perf.m_dot      = NaN;
    cea_out.perf.thrust_cf  = NaN;
    cea_out.perf.thrust_isp = NaN;
end

% ------------------------------------------------------------------
% Combustion efficiency corrections (eta_cstar, eta_Cf)
% ------------------------------------------------------------------
cea_out.eta_cstar  = eta_cstar;
cea_out.eta_Cf     = eta_Cf;
cea_out.cstar_real = cea_out.cstar * eta_cstar;
cea_out.Cf_real    = cea_out.Cf    * eta_Cf;
cea_out.Isp_real   = cea_out.Isp   * eta_cstar * eta_Cf;

cea_out.perf.cstar_real  = cea_out.cstar_real;
cea_out.perf.Cf_real     = cea_out.Cf_real;
cea_out.perf.Isp_real    = cea_out.Isp_real;

if ~isempty(A_t_m2)
    m_dot_real               = (Pc_bar * 1e5 * A_t_m2) / max(cea_out.cstar_real, 1.0);
    cea_out.perf.m_dot_real  = m_dot_real;
    cea_out.perf.thrust_real = cea_out.Cf_real * Pc_bar * 1e5 * A_t_m2 * lambda_nozzle;
    cea_out.perf.Isp_si_real = cea_out.Isp_real * 9.80665;
else
    cea_out.perf.m_dot_real  = NaN;
    cea_out.perf.thrust_real = NaN;
    cea_out.perf.Isp_si_real = NaN;
end

% ------------------------------------------------------------------
cea_out.files      = struct();
cea_out.files.fac  = fac_files;
cea_out.files.iac  = iac_files;
cea_out.files.conv = conv_files;

end

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function nozzle = normalize_nozzle_stations(nozzle_stations)

assert(isstruct(nozzle_stations),          'nozzle_stations must be a struct.');
assert(isfield(nozzle_stations, 'x'),      'nozzle_stations.x is required.');
assert(isfield(nozzle_stations, 'AR'),     'nozzle_stations.AR is required.');

x  = nozzle_stations.x(:);
AR = nozzle_stations.AR(:);

assert(numel(x) == numel(AR), 'x and AR must have the same length.');

[x, idx] = sort(x);
AR = AR(idx);

[x, iu] = unique(x, 'stable');
AR = AR(iu);

if isfield(nozzle_stations, 'throat_idx') && ~isempty(nozzle_stations.throat_idx)
    [~, throat_idx] = min(abs(x - x(nozzle_stations.throat_idx)));
else
    [~, throat_idx] = min(AR);
end

AR(throat_idx) = 1.0;

nozzle            = struct();
nozzle.x          = x;
nozzle.AR         = AR;
nozzle.N          = numel(x);
nozzle.throat_idx = throat_idx;
nozzle.x_throat   = x(throat_idx);
nozzle.idx_sub    = (1:throat_idx-1).';
nozzle.idx_throat = throat_idx;
nozzle.idx_sup    = (throat_idx+1:numel(x)).';

end

% -------------------------------------------------------------------------

function freeze = normalize_freeze_cfg(freeze_cfg, nozzle)

freeze      = struct();
freeze.mode = lower(string(freeze_cfg.mode));

if freeze.mode == "area_ratio"
    assert(isfield(freeze_cfg, 'AR') && ~isempty(freeze_cfg.AR), ...
        'freeze_cfg.AR is required when mode = area_ratio.');
    freeze.AR = double(freeze_cfg.AR);
else
    freeze.AR = NaN;
end

if freeze.mode == "exit"
    if isempty(nozzle.idx_sup)
        error('No supersonic stations available for mode = exit.');
    end
    freeze.AR = nozzle.AR(nozzle.idx_sup(end));
end

if freeze.mode == "area_ratio"
    if freeze.AR <= 1.0
        error('freeze_cfg.AR must be > 1 for a nozzle freeze station.');
    end
end

if freeze.mode == "area_ratio" || freeze.mode == "exit"
    nozzle = ensure_station_at_AR(nozzle, freeze.AR);
    freeze.inserted_station = true;
else
    freeze.inserted_station = false;
end

freeze.nozzle = nozzle;

switch freeze.mode
    case "chamber",    freeze.label = 'Frozen after chamber';
    case "throat",     freeze.label = 'Frozen after throat';
    case "area_ratio", freeze.label = sprintf('Frozen after A/At = %.6f', freeze.AR);
    case "exit",       freeze.label = 'Equilibrium to exit';
    otherwise,         error('Unsupported freeze mode: %s', freeze.mode);
end

end

% -------------------------------------------------------------------------

function nozzle = ensure_station_at_AR(nozzle, AR_target)

if isempty(nozzle.idx_sup)
    error('Cannot insert a supersonic freeze station without supersonic nodes.');
end

AR_sup = nozzle.AR(nozzle.idx_sup);
x_sup  = nozzle.x(nozzle.idx_sup);

if any(abs(AR_sup - AR_target) < 1e-10)
    return;
end

if AR_target < min(AR_sup) || AR_target > max(AR_sup)
    error('freeze_cfg.AR is outside the supersonic branch range.');
end

x_new  = interp1(AR_sup, x_sup, AR_target, 'linear');
x_all  = [nozzle.x;  x_new];
AR_all = [nozzle.AR; AR_target];

[x_all,  idx] = sort(x_all);
AR_all = AR_all(idx);

[x_all, iu] = unique(x_all, 'stable');
AR_all = AR_all(iu);

[~, throat_idx] = min(AR_all);

nozzle.x          = x_all;
nozzle.AR         = AR_all;
nozzle.N          = numel(x_all);
nozzle.throat_idx = throat_idx;
nozzle.x_throat   = x_all(throat_idx);
nozzle.idx_sub    = (1:throat_idx-1).';
nozzle.idx_throat = throat_idx;
nozzle.idx_sup    = (throat_idx+1:numel(x_all)).';

end

% -------------------------------------------------------------------------

function fac_map = build_fac_station_map(nozzle)

fac_map            = struct();
fac_map.x          = nozzle.x;
fac_map.AR         = nozzle.AR;
fac_map.throat_idx = nozzle.throat_idx;
fac_map.sub_idx    = nozzle.idx_sub;
fac_map.sup_idx    = nozzle.idx_sup;
fac_map.sub_AR     = nozzle.AR(fac_map.sub_idx).';
fac_map.sup_AR     = nozzle.AR(fac_map.sup_idx).';
fac_map.sub_x      = nozzle.x(fac_map.sub_idx).';
fac_map.sup_x      = nozzle.x(fac_map.sup_idx).';

end

% -------------------------------------------------------------------------

function [iac_map, iac_mode, nfz, freeze] = build_iac_station_map(nozzle_in, freeze)

nozzle = freeze.nozzle;

iac_map            = struct();
iac_map.x          = nozzle.x;
iac_map.AR         = nozzle.AR;
iac_map.throat_idx = nozzle.throat_idx;
iac_map.freeze_mode = freeze.mode;
iac_map.freeze_AR  = freeze.AR;

switch freeze.mode
    case "chamber"
        iac_mode        = 'frozen';
        nfz             = 1;
        iac_map.sub_idx = nozzle.idx_sub;
        iac_map.sup_idx = nozzle.idx_sup;

    case "throat"
        iac_mode        = 'frozen';
        nfz             = 2;
        iac_map.sub_idx = nozzle.idx_sub;
        iac_map.sup_idx = nozzle.idx_sup;

    case "area_ratio"
        iac_mode    = 'frozen';
        nfz         = 3;
        sup_idx     = nozzle.idx_sup;
        AR_sup      = nozzle.AR(sup_idx);
        kf          = find(abs(AR_sup - freeze.AR) < 1e-10, 1, 'first');
        if isempty(kf)
            error('Freeze area ratio not found in supersonic branch.');
        end
        freeze.sup_local_idx  = kf;
        freeze.sup_global_idx = sup_idx(kf);
        iac_map.sub_idx       = zeros(0,1);
        iac_map.sup_idx       = sup_idx(kf:end);

    case "exit"
        iac_mode        = 'equilibrium';
        nfz             = [];
        iac_map.sub_idx = nozzle.idx_sub;
        iac_map.sup_idx = nozzle.idx_sup;

    otherwise
        error('Unsupported freeze mode: %s', freeze.mode);
end

iac_map.sub_AR = nozzle.AR(iac_map.sub_idx).';
iac_map.sup_AR = nozzle.AR(iac_map.sup_idx).';
iac_map.sub_x  = nozzle.x(iac_map.sub_idx).';
iac_map.sup_x  = nozzle.x(iac_map.sup_idx).';

freeze.nozzle           = nozzle_in;
freeze.nozzle_augmented = nozzle;

if freeze.mode == "area_ratio"
    freeze.x = nozzle.x(freeze.sup_global_idx);
else
    freeze.x = NaN;
end

end

% -------------------------------------------------------------------------

function txt = build_fac_input(Pc_bar, T_react_K, ap_frac, htpb_frac, CR_fac, fac_map, case_name)

lines = {};
lines{end+1} = sprintf('problem case=%s rocket fac equilibrium p,bar=%.8f acat=%.10f', ...
    case_name, Pc_bar, CR_fac);

if ~isempty(fac_map.sub_AR)
    lines{end+1} = sprintf(' sub,ae/at=%s', num_list(fac_map.sub_AR));
end
if ~isempty(fac_map.sup_AR)
    lines{end+1} = sprintf(' sup,ae/at=%s', num_list(fac_map.sup_AR));
end

lines{end+1} = 'reac';
lines{end+1} = sprintf(' oxid=NH4CLO4(I) wt%%=%.8f t(k)=%.4f', ap_frac * 100,  T_react_K);
lines{end+1} = sprintf([' fuel=HTPB h,kj/mol=-58.0 C 7.075 H 10.65 O 0.223 N 0.063 ' ...
    'wt%%=%.8f t(k)=%.4f'], htpb_frac * 100, T_react_K);
lines{end+1} = 'output siunits short';
lines{end+1} = 'end';

txt = sprintf('%s\n', lines{:});

end

% -------------------------------------------------------------------------

function txt = build_iac_input(Pc_bar, T_react_K, ap_frac, htpb_frac, iac_map, case_name, iac_mode, nfz)

lines = {};

if strcmpi(iac_mode, 'equilibrium')
    lines{end+1} = sprintf('problem case=%s rocket equilibrium p,bar=%.8f', ...
        case_name, Pc_bar);
else
    lines{end+1} = sprintf('problem case=%s rocket frozen nfz=%d p,bar=%.8f', ...
        case_name, nfz, Pc_bar);
end

if ~isempty(iac_map.sub_AR)
    lines{end+1} = sprintf(' sub,ae/at=%s', num_list(iac_map.sub_AR));
end
if ~isempty(iac_map.sup_AR)
    lines{end+1} = sprintf(' sup,ae/at=%s', num_list(iac_map.sup_AR));
end

lines{end+1} = 'reac';
lines{end+1} = sprintf(' oxid=NH4CLO4(I) wt%%=%.8f t(k)=%.4f', ap_frac * 100,  T_react_K);
lines{end+1} = sprintf([' fuel=HTPB h,kj/mol=-58.0 C 7.075 H 10.65 O 0.223 N 0.063 ' ...
    'wt%%=%.8f t(k)=%.4f'], htpb_frac * 100, T_react_K);
lines{end+1} = 'output siunits short';
lines{end+1} = 'end';

txt = sprintf('%s\n', lines{:});

end

% -------------------------------------------------------------------------

function files = run_fcea_case(work_dir, FCEA2, case_name, input_text)

inp_file = fullfile(work_dir, [case_name, '.inp']);
out_file = fullfile(work_dir, [case_name, '.out']);
plt_file = fullfile(work_dir, [case_name, '.plt']);
csv_file = fullfile(work_dir, [case_name, '.csv']);

if isfile(inp_file), delete(inp_file); end
write_text_file(inp_file, input_text);

if isfile(out_file), delete(out_file); end
if isfile(plt_file), delete(plt_file); end
if isfile(csv_file), delete(csv_file); end

% Write the case name to a temporary input file for FCEA2
% (avoids echo/pipe issues on Windows)
stdin_file = fullfile(work_dir, [case_name, '_stdin.txt']);
fid = fopen(stdin_file, 'w');
fprintf(fid, '%s\n', case_name);
fclose(fid);

if ispc
    cmd = sprintf('cd /d "%s" && "%s" < "%s"', work_dir, FCEA2, stdin_file);
else
    cmd = sprintf('cd "%s" && "%s" < "%s"', work_dir, FCEA2, stdin_file);
end

[status, cmdout] = system(cmd);

if isfile(stdin_file), delete(stdin_file); end

if status ~= 0
    error('FCEA2 failed for case "%s":\n%s', case_name, cmdout);
end

if ~isfile(out_file)
    error('CEA output file not found for case "%s".', case_name);
end

files     = struct();
files.inp = inp_file;
files.out = out_file;
files.plt = plt_file;
files.csv = csv_file;

end

% -------------------------------------------------------------------------

function parsed = parse_cea_result(out_file, model, station_map)

raw = fileread(out_file);

T     = extract_cea_row(raw, {'T, K'});
P     = extract_cea_row(raw, {'P, BAR'});
Mw    = extract_cea_row(raw, {'M, (1/n)', 'M, (1/N)'});
gamma = extract_cea_row(raw, {'GAMMAs', 'GAMMA'});
cp    = extract_cea_row(raw, {'Cp, KJ/(KG)(K)', 'CP, KJ/(KG)(K)'});
Mach  = extract_cea_row(raw, {'MACH NUMBER'});
cstar = extract_cea_row(raw, {'CSTAR, M/SEC'});
cf    = extract_cea_row(raw, {'CF'});
Isp   = extract_cea_row(raw, {'Isp, M/SEC', 'ISP, M/SEC'});

if isempty(T) || isempty(P) || isempty(Mw) || isempty(gamma)
    error('CEA parsing failed: key rows not found in %s', out_file);
end

switch lower(model)
    case 'fac'
        idx.injector = 1;
        idx.chamber  = 2;
        idx.throat   = 3;
        idx.exit0    = 4;
    case 'iac'
        idx.injector = [];
        idx.chamber  = 1;
        idx.throat   = 2;
        idx.exit0    = 3;
    case 'iac_eq'
        % IAC equilibrium convergent-only: chamber | throat | sub(1..Ns)
        idx.injector = [];
        idx.chamber  = 1;
        idx.throat   = 2;
        idx.exit0    = 3;
    otherwise
        error('Unknown model type: %s', model);
end

parsed           = struct();
parsed.model     = lower(model);
parsed.raw_file  = out_file;

if ~isempty(idx.injector) && numel(T) >= idx.injector
    parsed.injector = build_station_state(idx.injector, T, P, Mw, gamma, cp, Mach, cstar, cf, Isp);
else
    parsed.injector = [];
end

parsed.chamber = build_station_state(idx.chamber, T, P, Mw, gamma, cp, Mach, cstar, cf, Isp);
parsed.throat  = build_station_state(idx.throat,  T, P, Mw, gamma, cp, Mach, cstar, cf, Isp);

n_sub_req = numel(station_map.sub_idx);
n_sup_req = numel(station_map.sup_idx);

parsed.sub = init_branch(0);
parsed.sup = init_branch(0);

% FAC: sub/sup branches empty - chamber scalars only
if strcmpi(model, 'fac')
    parsed.exit = parsed.throat;
    return;
end

exit_idx0       = idx.exit0;
n_cols_avail    = max(numel(T) - exit_idx0 + 1, 0);

if strcmpi(model, 'iac_eq')
    n_sub = min(n_sub_req, n_cols_avail);
    n_sup = 0;
else
    n_sub = 0;
    n_sup = min(n_sup_req, n_cols_avail);
end

parsed.sub = init_branch(n_sub);
parsed.sup = init_branch(n_sup);

for k = 1:n_sub
    s = build_station_state(exit_idx0 + k - 1, T, P, Mw, gamma, cp, Mach, cstar, cf, Isp);
    parsed.sub = assign_branch_state(parsed.sub, k, station_map.sub_x(k), station_map.sub_AR(k), s);
end

for k = 1:n_sup
    s = build_station_state(exit_idx0 + n_sub + k - 1, T, P, Mw, gamma, cp, Mach, cstar, cf, Isp);
    parsed.sup = assign_branch_state(parsed.sup, k, station_map.sup_x(k), station_map.sup_AR(k), s);
end

if n_sup > 0
    parsed.exit = get_branch_state(parsed.sup, n_sup);
elseif n_sub > 0
    parsed.exit = get_branch_state(parsed.sub, n_sub);
else
    parsed.exit = parsed.throat;
end

end

% -------------------------------------------------------------------------

function nozzle_out = merge_nozzle_profiles(nozzle_ref, fac, iac, conv_eq, freeze)

noz = freeze.nozzle_augmented;
N   = noz.N;

nozzle_out         = struct();
nozzle_out.x       = noz.x(:);
nozzle_out.AR      = noz.AR(:);
nozzle_out.T       = NaN(N,1);
nozzle_out.P       = NaN(N,1);
nozzle_out.Mw      = NaN(N,1);
nozzle_out.gamma   = NaN(N,1);
nozzle_out.cp      = NaN(N,1);
nozzle_out.R       = NaN(N,1);
nozzle_out.Mach    = NaN(N,1);
nozzle_out.velocity = NaN(N,1);
nozzle_out.Cf      = NaN(N,1);
nozzle_out.Isp     = NaN(N,1);
nozzle_out.c_eff   = NaN(N,1);
nozzle_out.region  = strings(N,1);
nozzle_out.mode    = strings(N,1);
nozzle_out.is_freeze_point = false(N,1);
nozzle_out.is_exit         = false(N,1);

% POINT 4 - Use IAC-EQ convergent data if available.
if ~isempty(conv_eq) && ~isempty(conv_eq.sub.x)
    for k = 1:numel(conv_eq.sub.x)
        idx_k = find(abs(nozzle_out.x - conv_eq.sub.x(k)) < 1e-12, 1);
        nozzle_out = assign_nozzle_row(nozzle_out, idx_k, conv_eq.sub, k, "sub", "equilibrium");
    end
end

if ~isempty(fac.sub.x)
    for k = 1:numel(fac.sub.x)
        idx_k = find(abs(nozzle_out.x - fac.sub.x(k)) < 1e-12, 1);
        nozzle_out = assign_nozzle_row(nozzle_out, idx_k, fac.sub, k, "sub", "equilibrium");
    end
end

idx_t = noz.throat_idx;
nozzle_out = assign_nozzle_state(nozzle_out, idx_t, fac.throat, "throat", "equilibrium");

switch freeze.mode

    case "throat"
        nozzle_out.is_freeze_point(idx_t) = true;
        if ~isempty(iac.sup.x)
            for k = 1:numel(iac.sup.x)
                idx_k = find(abs(nozzle_out.x - iac.sup.x(k)) < 1e-12, 1);
                nozzle_out = assign_nozzle_row(nozzle_out, idx_k, iac.sup, k, "sup", "frozen");
            end
        end

    case "chamber"
        nozzle_out = assign_nozzle_state(nozzle_out, idx_t, iac.throat, "throat", "frozen");
        if ~isempty(iac.sub.x)
            for k = 1:numel(iac.sub.x)
                idx_k = find(abs(nozzle_out.x - iac.sub.x(k)) < 1e-12, 1);
                nozzle_out = assign_nozzle_row(nozzle_out, idx_k, iac.sub, k, "sub", "frozen");
            end
        end
        if ~isempty(iac.sup.x)
            for k = 1:numel(iac.sup.x)
                idx_k = find(abs(nozzle_out.x - iac.sup.x(k)) < 1e-12, 1);
                nozzle_out = assign_nozzle_row(nozzle_out, idx_k, iac.sup, k, "sup", "frozen");
            end
        end

    case "area_ratio"
        if ~isempty(fac.sup.x)
            for k = 1:numel(fac.sup.x)
                idx_k = find(abs(nozzle_out.x - fac.sup.x(k)) < 1e-12, 1);
                if fac.sup.AR(k) <= freeze.AR + 1e-10
                    nozzle_out = assign_nozzle_row(nozzle_out, idx_k, fac.sup, k, "sup", "equilibrium");
                end
            end
        end
        if ~isempty(iac.sup.x)
            for k = 1:numel(iac.sup.x)
                idx_k = find(abs(nozzle_out.x - iac.sup.x(k)) < 1e-12, 1);
                if iac.sup.AR(k) > freeze.AR + 1e-10
                    nozzle_out = assign_nozzle_row(nozzle_out, idx_k, iac.sup, k, "sup", "frozen");
                end
            end
        end
        idx_f = find(abs(nozzle_out.AR - freeze.AR) < 1e-10 & nozzle_out.x > noz.x_throat, 1);
        if ~isempty(idx_f)
            nozzle_out.is_freeze_point(idx_f) = true;
        end

    case "exit"
        if ~isempty(iac.sup.x)
            for k = 1:numel(iac.sup.x)
                idx_k = find(abs(nozzle_out.x - iac.sup.x(k)) < 1e-12, 1);
                nozzle_out = assign_nozzle_row(nozzle_out, idx_k, iac.sup, k, "sup", "equilibrium");
            end
        end
        if ~isempty(noz.idx_sup)
            nozzle_out.is_freeze_point(noz.idx_sup(end)) = true;
        end

    otherwise
        error('Unsupported freeze mode in merge_nozzle_profiles: %s', freeze.mode);
end

if ~isempty(noz.idx_sup)
    nozzle_out.is_exit(noz.idx_sup(end)) = true;
else
    nozzle_out.is_exit(noz.idx_throat) = true;
end

nozzle_out.throat = row_to_state(nozzle_out, idx_t);

if any(nozzle_out.is_exit)
    idx_e = find(nozzle_out.is_exit, 1);
else
    idx_e = idx_t;
end
nozzle_out.exit = row_to_state(nozzle_out, idx_e);

nozzle_out.reference = nozzle_ref;
nozzle_out.freeze    = freeze;

end

% =========================================================================
% BRANCH HELPERS
% =========================================================================

function branch = init_branch(n)
branch          = struct();
branch.x        = NaN(n,1);
branch.AR       = NaN(n,1);
branch.T        = NaN(n,1);
branch.P        = NaN(n,1);
branch.Mw       = NaN(n,1);
branch.gamma    = NaN(n,1);
branch.cp       = NaN(n,1);
branch.R        = NaN(n,1);
branch.Mach     = NaN(n,1);
branch.velocity = NaN(n,1);
branch.cstar    = NaN(n,1);
branch.Cf       = NaN(n,1);
branch.Isp      = NaN(n,1);
branch.c_eff    = NaN(n,1);
end

function branch = assign_branch_state(branch, k, x, AR, s)
branch.x(k)        = x;
branch.AR(k)       = AR;
branch.T(k)        = s.T;
branch.P(k)        = s.P;
branch.Mw(k)       = s.Mw;
branch.gamma(k)    = s.gamma;
branch.cp(k)       = s.cp;
branch.R(k)        = s.R;
branch.Mach(k)     = s.Mach;
branch.velocity(k) = s.velocity;
branch.cstar(k)    = s.cstar;
branch.Cf(k)       = s.Cf;
branch.Isp(k)      = s.Isp;
branch.c_eff(k)    = s.c_eff;
end

function s = get_branch_state(branch, k)
s          = struct();
s.x        = branch.x(k);
s.AR       = branch.AR(k);
s.T        = branch.T(k);
s.P        = branch.P(k);
s.Mw       = branch.Mw(k);
s.gamma    = branch.gamma(k);
s.cp       = branch.cp(k);
s.R        = branch.R(k);
s.Mach     = branch.Mach(k);
s.velocity = branch.velocity(k);
s.cstar    = branch.cstar(k);
s.Cf       = branch.Cf(k);
s.Isp      = branch.Isp(k);
s.c_eff    = branch.c_eff(k);
end

% =========================================================================
% NOZZLE PROFILE HELPERS
% =========================================================================

function nozzle_out = assign_nozzle_row(nozzle_out, idx, branch, k, region, mode)
s = get_branch_state(branch, k);
nozzle_out = assign_nozzle_state(nozzle_out, idx, s, region, mode);
end

function nozzle_out = assign_nozzle_state(nozzle_out, idx, s, region, mode)
if isempty(idx) || any(isnan(idx))
    return;
end
nozzle_out.T(idx)        = s.T;
nozzle_out.P(idx)        = s.P;
nozzle_out.Mw(idx)       = s.Mw;
nozzle_out.gamma(idx)    = s.gamma;
nozzle_out.cp(idx)       = s.cp;
nozzle_out.R(idx)        = s.R;
nozzle_out.Mach(idx)     = s.Mach;
nozzle_out.velocity(idx) = s.velocity;
nozzle_out.Cf(idx)       = s.Cf;
nozzle_out.Isp(idx)      = s.Isp;
nozzle_out.c_eff(idx)    = s.c_eff;
nozzle_out.region(idx)   = region;
nozzle_out.mode(idx)     = mode;
end

function s = row_to_state(nozzle_out, idx)
s          = struct();
s.x        = nozzle_out.x(idx);
s.AR       = nozzle_out.AR(idx);
s.T        = nozzle_out.T(idx);
s.P        = nozzle_out.P(idx);
s.Mw       = nozzle_out.Mw(idx);
s.gamma    = nozzle_out.gamma(idx);
s.cp       = nozzle_out.cp(idx);
s.R        = nozzle_out.R(idx);
s.Mach     = nozzle_out.Mach(idx);
s.velocity = nozzle_out.velocity(idx);
s.Cf       = nozzle_out.Cf(idx);
s.Isp      = nozzle_out.Isp(idx);
s.c_eff    = nozzle_out.c_eff(idx);
s.region   = nozzle_out.region(idx);
s.mode     = nozzle_out.mode(idx);
end

% =========================================================================
% CEA OUTPUT PARSER
% =========================================================================

function s = build_station_state(i, T, P, Mw, gamma, cp, Mach, cstar, cf, Isp)
s         = struct();
s.T       = pick_value(T,     i);
s.P       = pick_value(P,     i);
s.Mw      = pick_value(Mw,    i);
s.gamma   = pick_value(gamma, i);
s.cp      = 1000 * pick_value(cp, i);
s.R       = 8314.462618 / s.Mw;
s.Mach    = pick_value(Mach,  i);
s.velocity = s.Mach * sqrt(max(s.gamma * s.R * s.T, 0.0));

if isempty(cstar)
    s.cstar = NaN;
else
    s.cstar = pick_value(cstar, min(i, numel(cstar)));
end

if isempty(cf)
    s.Cf = NaN;
else
    s.Cf = pick_value(cf, min(i, numel(cf)));
end

if isempty(Isp)
    s.Isp   = NaN;
    s.c_eff = NaN;
else
    s.Isp   = pick_value(Isp, min(i, numel(Isp))) / 9.80665;
    s.c_eff = s.Isp * 9.80665;
end
end

function v = pick_value(vec, i)
if isempty(vec)
    v = NaN;
elseif i <= numel(vec)
    v = vec(i);
else
    v = vec(end);
end
end

function vals = extract_cea_row(raw, labels)
if ischar(labels) || isstring(labels)
    labels = cellstr(labels);
end

vals  = [];
lines = regexp(raw, '\r\n|\n|\r', 'split');

for j = 1:numel(labels)
    label = char(labels{j});
    patt  = ['^\s*', regexptranslate('escape', label), '\s+(.+?)\s*$'];

    tmp    = [];
    nmatch = 0;

    for i = 1:numel(lines)
        tok = regexp(lines{i}, patt, 'tokens', 'once', 'ignorecase');
        if isempty(tok), continue; end

        row_vals = parse_cea_numeric_tokens(tok{1});
        if isempty(row_vals), continue; end

        nmatch = nmatch + 1;

        if nmatch == 1
            tmp = row_vals;
        else
            if strcmpi(strtrim(label), 'CF') || contains(upper(label), 'ISP')
                if numel(row_vals) >= 2
                    row_vals = row_vals(2:end);
                else
                    row_vals = [];
                end
            else
                if numel(row_vals) >= 3
                    row_vals = row_vals(3:end);
                else
                    row_vals = [];
                end
            end
            tmp = [tmp, row_vals]; %#ok
        end
    end

    if ~isempty(tmp)
        vals = tmp;
        return;
    end
end
end

function vals = parse_cea_numeric_tokens(str_line)
tokens = regexp(str_line, '\S+', 'match');
vals   = [];
for k = 1:numel(tokens)
    val = convert_cea_number(tokens{k});
    if ~isnan(val)
        vals(end+1) = val; %#ok
    end
end
end

function val = convert_cea_number(tok)
tok = strtrim(tok);
tok = strrep(tok, 'D', 'E');
tok = strrep(tok, 'd', 'E');
tok = strrep(tok, '*', '');

if isempty(tok)
    val = NaN; return;
end

% Standard scientific notation: 1.23E+04
if ~isempty(regexp(tok, '^[+-]?\d*\.?\d+[Ee][+-]?\d+$', 'once'))
    val = str2double(tok); return;
end

% Fortran-style: 1.23+04 or 1.23-04
if ~isempty(regexp(tok, '^[+-]?\d*\.?\d+[+-]\d+$', 'once'))
    sp = regexp(tok(2:end), '[+-]', 'once', 'end');
    if ~isempty(sp)
        sp  = sp + 1;
        tok = [tok(1:sp-1), 'E', tok(sp:end)];
        val = str2double(tok); return;
    end
end

% Plain decimal
if ~isempty(regexp(tok, '^[+-]?\d*\.?\d+$', 'once'))
    val = str2double(tok); return;
end

val = NaN;
end

function write_text_file(file_path, txt)
fid = fopen(file_path, 'w');
if fid < 0
    error('Cannot open file for writing: %s', file_path);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s', txt);
end

function conv_map = build_conv_station_map(nozzle)
conv_map            = struct();
conv_map.x          = nozzle.x;
conv_map.AR         = nozzle.AR;
conv_map.throat_idx = nozzle.throat_idx;
conv_map.sub_idx    = nozzle.idx_sub;
conv_map.sup_idx    = zeros(0,1);
conv_map.sub_AR     = nozzle.AR(nozzle.idx_sub).';
conv_map.sup_AR     = zeros(1,0);
conv_map.sub_x      = nozzle.x(nozzle.idx_sub).';
conv_map.sup_x      = zeros(1,0);
end

function txt = build_iac_eq_input(Pc_bar, T_react_K, ap_frac, htpb_frac, conv_map, case_name)
lines = {};
lines{end+1} = sprintf('problem case=%s rocket equilibrium p,bar=%.8f', ...
    case_name, Pc_bar);
if ~isempty(conv_map.sub_AR)
    lines{end+1} = sprintf(' sub,ae/at=%s', num_list(conv_map.sub_AR));
end
lines{end+1} = 'reac';
lines{end+1} = sprintf(' oxid=NH4CLO4(I) wt%%=%.8f t(k)=%.4f', ap_frac*100, T_react_K);
lines{end+1} = sprintf([' fuel=HTPB h,kj/mol=-58.0 C 7.075 H 10.65 O 0.223 N 0.063 ' ...
    'wt%%=%.8f t(k)=%.4f'], htpb_frac*100, T_react_K);
lines{end+1} = 'output siunits short';
lines{end+1} = 'end';
txt = sprintf('%s\n', lines{:});
end

function s = num_list(vec)
if isempty(vec)
    s = ''; return;
end
parts = arrayfun(@(v) sprintf('%.8f', v), vec, 'UniformOutput', false);
s = strjoin(parts, ',');
end
