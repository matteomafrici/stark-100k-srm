% =========================================================
% COOLING JACKET THERMAL ANALYSIS - Stark-100k Motor
% Nozzle throat region: area ratio 2 (sub) to area ratio 2 (sup)
% Wall material: CuCrZr (4 mm)
% Coolant: tap water, co-flow, pressurised
% TBC: YSZ (if needed)
% =========================================================
clc; clear all; close all;

%% ============================================================
%  SECTION 1 - GAS PROPERTIES AT THROAT (from RPA / CEA output)
% ============================================================
% All values at throat conditions (M=1), non-aluminized AP/HTPB
% at p_c = 70 bar. Source: RPA thermochemical code output
% (consistent with assignment data).

p_c     = 65e5;             % Chamber pressure [Pa]
Tc      = 2801.53;          % Adiabatic flame temperature [K]
gamma   = 1.1640;           % Ratio of specific heats [-] (frozen, at throat)
R_gas   = 8314 / 25.0;      % Specific gas constant [J/kgK]
                            % Molar mass ~25 g/mol for AP/HTPB non-aluminized
                            % (Sutton & Biblarz, Table 5-6, ~24-26 g/mol)
cp_g    = 2745.9;           % Specific heat gas [J/kgK]
mu_g    = 9.3556e-5;        % Dynamic viscosity [Pa*s]
lambda_g= 0.42194;          % Thermal conductivity gas [W/mK]
Pr_g    = mu_g * cp_g / lambda_g;  % Prandtl number gas [-]

% Throat temperature (isentropic relation)
T_t = Tc / (1 + (gamma-1)/2);      % [K]  (M=1 -> T_t = T0/(1+(g-1)/2))

% Speed of sound and velocity at throat
a_t  = sqrt(gamma * R_gas * T_t);  % [m/s]
v_t  = a_t;                        % M=1 at throat

% Throat density (ideal gas)
p_t  = p_c * (2/(gamma+1))^(gamma/(gamma-1));  % [Pa]
rho_t= p_t / (R_gas * T_t);                    % [kg/m^3]

fprintf('=== GAS PROPERTIES AT THROAT ===\n');
fprintf('T_throat     = %.1f K\n', T_t);
fprintf('p_throat     = %.4f bar\n', p_t/1e5);
fprintf('rho_throat   = %.4f kg/m^3\n', rho_t);
fprintf('v_throat     = %.2f m/s\n', v_t);
fprintf('Pr_gas       = %.4f\n', Pr_g);

%% ============================================================
%  SECTION 2 - NOZZLE GEOMETRY
% ============================================================
% Jacket covers from area ratio epsilon=2 (subsonic) to epsilon=2 (supersonic)
% Throat diameter from assignment: D_t = 0.1067 m

D_t  = 0.114551328;              % Throat diameter [m]
r_t  = D_t/2;               % Throat radius [m]
A_t  = pi/4 * D_t^2;        % Throat area [m^2]

% Area ratio epsilon = A/A_t = 2 -> radius at epsilon=2
A_eps2 = 2 * A_t;
r_eps2 = sqrt(A_eps2 / pi); % [m]

% Half-angles (from assignment)
alpha_div = 15;             % Divergence half-angle [deg]
beta_conv = 45;             % Convergence half-angle [deg]

% Axial lengths of convergent and divergent sections
L_conv = (r_eps2 - r_t) / tand(beta_conv);   % [m]
L_div  = (r_eps2 - r_t) / tand(alpha_div);   % [m]

% Slant lengths (frustum lateral surface)
s_conv = sqrt(L_conv^2 + (r_eps2 - r_t)^2);
s_div  = sqrt(L_div^2  + (r_eps2 - r_t)^2);

% Lateral (wetted) area of the jacket - frustum formula
A_conv = pi * (r_eps2 + r_t) * s_conv;       % [m^2]
A_div  = pi * (r_eps2 + r_t) * s_div;        % [m^2]
A_jacket = A_conv + A_div;                   % [m^2] total wetted area

% Total axial length of jacket
L_jacket = L_conv + L_div;                   % [m]
x = linspace(-L_conv, L_div, 200);
r = zeros(size(x));
for i = 1:length(x)
    if x(i) < 0
        r(i) = r_t - x(i)*tand(beta_conv);
    else
        r(i) = r_t + x(i)*tand(alpha_div);
    end
end
A = pi*r.^2;
AR = A / A_t;

fprintf('\n=== NOZZLE GEOMETRY ===\n');
fprintf('Throat diameter D_t   = %.4f m\n', D_t);
fprintf('Radius at eps=2       = %.4f m\n', r_eps2);
fprintf('L_convergent          = %.4f m\n', L_conv);
fprintf('L_divergent           = %.4f m\n', L_div);
fprintf('Total jacket length   = %.4f m\n', L_jacket);
fprintf('Wetted area A_jacket  = %.5f m^2\n', A_jacket);

%% ============================================================
%  SECTION 3.1 - GAS-SIDE HEAT TRANSFER COEFFICIENT (Dittus equation)
% ============================================================

m_dot=43.7396;
vc=m_dot./(rho_t.*A);
Re=(rho_t*vc*2.*r)./mu_g;
Nu=0.0265*(Re.^0.8)*(Pr_g^0.3);
hg=Nu.*lambda_g./(2.*r);


%% ============================================================
%  SECTION 3.2 - GAS-SIDE HEAT TRANSFER COEFFICIENT (Bartz equation)
% ============================================================
% The Bartz equation is the standard method for rocket nozzle
% convective heat transfer (Bartz, 1957; Sutton & Biblarz Section 8).
% It accounts for the nozzle geometry via the throat curvature
% radius R_c and the local area ratio.
%
%   h_g = (0.026/D_t^0.2) * (mu_g^0.2 * cp_g / Pr_g^0.6)
%         * (p_c * g0 / c_star)^0.8 * (D_t/R_c)^0.1 * sigma
%
% where sigma is the correction factor for property variation
% across the boundary layer (wall-to-gas temperature ratio).
% At the throat (worst case), sigma ~ 0.5-0.9.
% We use the simplified Bartz form evaluated at the throat.

% Characteristic velocity c* [m/s]
% c* = p_c * A_t / m_dot
% m_dot from assignment: F=100kN, Isp~235s -> m_dot = F/(Isp*g0)
g0      = 9.81;
Isp_est = 235;              % [s] estimated (AP/HTPB non-aluminized, lit.)
m_dot   = 100000 / (Isp_est * g0);  % [kg/s]
c_star  = p_c * A_t / m_dot;        % [m/s]

fprintf('\n=== PROPULSION PARAMETERS ===\n');
fprintf('m_dot        = %.3f kg/s\n', m_dot);
fprintf('c*           = %.1f m/s\n', c_star);

% Throat radius of curvature R_c
% Standard assumption: R_c = 0.5 * r_t (Sutton Section 8, typical for solid motors)
% This is conservative (smaller R_c -> higher h_g)
R_c = 0.5 * r_t;

% Bartz simplified equation at throat (epsilon=1, local correction = 1)
% sigma: boundary layer correction. For T_wall/T_gas ~ 0.3-0.5, sigma ~ 0.6
% Conservative assumption: sigma = 0.6 (Huzel & Huang, Modern Engineering
% for Design of Liquid-Propellant Rocket Engines, AIAA 1992, p.82)
sigma = 0.6;

%convective heat transfer using Bartz model
hg_b = (0.026 ./(2*r.^0.2)) .* ...
     (mu_g^0.2 * cp_g / Pr_g^0.6) * ...
     (p_c / c_star)^0.8 .* ...
     (2.*r./ R_c).^0.1 * sigma;        % [W/m^2K]

fprintf('\n=== GAS-SIDE HEAT TRANSFER ===\n');
fprintf('c*           = %.1f m/s\n', c_star);
fprintf('R_c (throat) = %.4f m\n', R_c);
fprintf('sigma        = %.2f\n', sigma);


figure;
plot(x,hg, 'LineWidth', 2); hold on
plot(x,hg_b, 'LineWidth', 2);
xlabel('x [m]');
ylabel('hg [W/m^2K]')
legend('hg Dittus [W/m^2K]', 'hg Bartz [W/m^2K]');
title('Profil between A/A* = 2 boundaries');
grid on

%% ============================================================
%  SECTION 4 - ADIABATIC WALL TEMPERATURE
% ============================================================
% Recovery factor for turbulent flow: r = Pr^(1/3)
% (White, Viscous Fluid Flow, 3rd ed., Section 7.4)
% T_aw = T_t * (1 + r*(gamma-1)/2 * M^2) evaluated at M=1
% Since T_t is already the stagnation temperature at throat
% and M=1: T_aw = T_t * (1 + r*(gamma-1)/2)

%area Mach function
area_mach = @(M)(1./M).*((2/(gamma+1).*(1 + (gamma-1)/2.*M.^2)).^((gamma+1)/(2*(gamma-1))));

%solve Mach
M = zeros(size(x));

for i = 1:length(x)
    if abs(AR(i)-1) < 1e-6
        M(i) = 1;   
    elseif x(i) < 0
        % subsonic branch
        M(i) = fsolve(@(M) area_mach(M)-AR(i), 0.3);    
    else
        % supersonic branch
        M(i) = fsolve(@(M) area_mach(M)-AR(i), 2);
    end    
end

Recov_fact_t=(1+(Pr_g^(1/3)).*((gamma-1)/2).*M.^2)./(1+((gamma-1)/2).*M.^2);
T_aw_t=Recov_fact_t.*Tc; 
figure;
plot(x,T_aw_t, 'LineWidth', 2);
xlabel('x [m]');
ylabel('T_aw [K]');
title('T_aw Profil between A/A* = 2 boundaries');
grid on


fprintf('\n=== ADIABATIC WALL TEMPERATURE ===\n');
fprintf('Recovery factor r     = %.4f\n', mean(Recov_fact_t));
fprintf('T_aw (mean)         = %.1f K\n', mean(T_aw_t));

%% ============================================================
%  SECTION 5 - WALL MATERIAL: CuCrZr
% ============================================================
% CuCrZr (C18150) is a precipitation-hardened copper alloy used
% in rocket combustion chambers (Vulcain, RS-68, etc.).
% Key properties (Aurubis datasheet, ITER MPH database):
%   - Thermal conductivity: ~320 W/mK (vs 390 for pure Cu)
%   - Max service temperature: ~300 degC (573 K) to retain
%     precipitation hardening (Cr2Zr precipitates dissolve above
%     ~400 degC, causing rapid strength loss)
%   - Melting point: ~1075 degC (1348 K)
% Conservative structural limit: T_wall_max = 573 K (300 degC)


t_wall  = 0.00219;          % Wall thickness [m] (2.1 mm, from literature/geometry)
k_wall  = 320;              % CuCrZr thermal conductivity [W/mK]
T_wall_max = 573;           % Max allowable wall temp (hot side) [K]
                            % (precipitation hardening limit, Aurubis DS)

%% ============================================================
%  SECTION 6 - COOLANT (TAP WATER) PROPERTIES
% ============================================================
% Tap water at bulk temperature T_bulk = (T_in + T_out)/2
% Assignment: T_in = 18 degC, T_out must not boil during 29 s burn.
% Constraint: water must NOT change phase -> T_out < T_sat(p_water)
%
% Cooling jacket pressure:
% Literature recommendation: p_coolant = 5-25% of p_chamber
% (Huzel & Huang, AIAA 1992; Sutton Section 8)
% Conservative choice: p_coolant = 0.15 * p_c 
% -> T_sat(10.5 bar) ~ 182 degC (455 K)  [steam tables]
% We set T_out_max = 50 degC (353 K) as design limit, well below
% saturation, providing a safety margin of >100 degC against boiling.

p_water  = 0.15 * p_c;     % Coolant pressure 
T_in_w   = 18 + 273.15;    % Inlet water temperature [K]
T_out_w  = 50 + 273.15;    % Max outlet temperature [K] (design limit)
T_bulk_w = (T_in_w + T_out_w) / 2;  % Bulk temperature [K]

% Water properties at T_bulk ~ 322 K (49 degC)
% Source: NIST WebBook / Incropera Table A.6
rho_w   = 988;              % Density [kg/m^3]
cp_w    = 4181;             % Specific heat [J/kgK]
mu_w    = 5.5e-4;           % Dynamic viscosity [Pa*s]
lambda_w= 0.644;            % Thermal conductivity [W/mK]
Pr_w    = mu_w * cp_w / lambda_w;  % Prandtl number

T_sat   = 182+273.15;              % Saturation temperature at 9.75 bar [K]
fprintf('\n=== COOLANT PROPERTIES ===\n');
fprintf('p_coolant    = %.2f bar\n', p_water/1e5);
fprintf('T_sat at p   = %.1f K (%.1f C)\n', T_sat, T_sat-273.15);
fprintf('T_out design = %.1f K (%.1f C) -- margin: %.1f K\n', ...
        T_out_w, T_out_w-273.15, T_sat-T_out_w);
fprintf('Pr_water     = %.4f\n', Pr_w);

%% ============================================================
%  SECTION 7 - THERMAL RESISTANCE NETWORK (1D steady-state)
% ============================================================
% The heat path from hot gas to coolant is modelled as series
% resistances per unit area [m^2K/W]:
%
%   Gas (convection) -> TBC (conduction) -> Wall (conduction) -> Water
%
%   R_gas  = 1/h_g
%   R_tbc  = t_tbc / k_tbc
%   R_wall = t_wall / k_wall
%   R_w    = 1/h_w
%   R_tot  = R_gas + R_tbc + R_wall + R_w
%
%   q = (T_aw - T_bulk_w) / R_tot
%
% The temperature at each interface is then:
%   T_tbc_hot  = T_aw  - q * R_gas          (TBC hot face)
%   T_tbc_cold = T_tbc_hot - q * R_tbc      (TBC cold face = wall hot face)
%   T_wall_cold= T_tbc_cold - q * R_wall    (wall cold face)
%   T_water    = T_bulk_w                   (bulk coolant)

R_gas_res  = 1 ./ hg;%using Dittus-Boethler correlation
R_gas_res_b  = 1 ./ hg_b; %using Bartz correlation
R_wall_res = t_wall / k_wall; %wall resistance plane model
R_wall_res_cyl=r.*log((r+t_wall)./r)./(k_wall); %wall resistance cylindical model

%% ============================================================
%  SECTION 9 - CASE A: WITH TBC plane (YSZ)
% ============================================================
% YSZ (8 wt% yttria-stabilized zirconia) properties:
%   k_tbc = 0.9-2.5 W/mK depending on porosity and deposition.
%   For APS-deposited YSZ (standard): k_tbc ~ 1.0-1.5 W/mK
%   Conservative (lower conductivity = better insulation): 0.9 W/mK
%   Source: Clarke et al., MRS Bulletin 2012; Marple et al. 2007
%   (documents in literature review)

k_tbc = 0.90;
t_tbc=linspace(0,450,451); %Thickness range tbc [micrometer]
R_tbc_1=(t_tbc*(10^-6))./k_tbc;%thermic resistance of the tbc

%plane model
R_s_tbc_CuCrZr=R_wall_res+min(R_gas_res)+R_tbc_1;
q_TBC_CuCrZr=(max(T_aw_t)-T_bulk_w)./R_s_tbc_CuCrZr;
T_wh_TBC_CuCrZr=max(T_aw_t)-q_TBC_CuCrZr.*(min(R_gas_res)+R_tbc_1);

R_s_notbc_CuCrZr=R_wall_res+(R_gas_res);
q_noTBC=(T_aw_t-T_bulk_w)./R_s_notbc_CuCrZr;
T_wh_noTBC_CuCrZr=T_aw_t-q_noTBC.*(R_gas_res);

valid = (T_wh_TBC_CuCrZr < T_wall_max/1.2);
if ~any(valid)
    error('No TBC thickness satisfies both thermal constraints!');
end
% Take minimum thickness satisfying both constraints
idx = find(valid, 1, 'first');
t_tbc_opt = t_tbc(idx);
fprintf('Optimal TBC thickness = %.4f micrometers\n', t_tbc_opt);

%Tube inner diameter
D=0.002;

%V required for minimal tbc thickness find in previous part
R_tbc_opt = (t_tbc_opt*10^-6) / k_tbc;
R_tot_tbc_CuCrZr_opt=R_wall_res+min(R_gas_res)+R_tbc_opt;
q_TBC_CuCrZr_opt=(max((T_aw_t))-T_bulk_w)./R_tot_tbc_CuCrZr_opt;
Q_TBC_CuCrZr_opt=q_TBC_CuCrZr_opt*A_jacket;
%Using dittus correlation to recover the velocity
h_req_TBC_CuCrZr_opt = q_TBC_CuCrZr_opt / (T_sat/1.1 - T_bulk_w);
Nu_TBC_CuCrZr_opt = h_req_TBC_CuCrZr_opt * D / lambda_w;
Re_TBC_CuCrZr_opt = (Nu_TBC_CuCrZr_opt/(0.0243*Pr_w^0.4))^(1/0.8);
v_req_TBC_CuCrZr_opt = Re_TBC_CuCrZr_opt * mu_w / (rho_w * D);

t = t_tbc_opt*1e-6;
t_max = 400e-6;
dt = 5e-6;
v_req = 100;
t_hist = [];
v_hist = [];

%give couple verifying that 2 constraints
while (v_req > 15) && (t < t_max)

    R_tbc = t / k_tbc;
    R_tot = min(R_gas_res) + R_tbc + R_wall_res;

    q = (max(T_aw_t) - T_bulk_w) / R_tot;
    Q=q*A_jacket;

    h_req = q / (T_sat/1.2 - T_bulk_w);

    Nu = h_req * D / lambda_w;
    Re = (Nu/(0.0243*Pr_w^0.4))^(1/0.8);

    v_req = Re * mu_w / (rho_w * D);

    t = t + dt;
    t_hist(end+1) = t*1e6; % en um
    v_hist(end+1) = v_req;
end

t_opt = t;
v_final_req = v_req;

m_dot = Q/(cp_w*(T_out_w - T_in_w));
A = pi*(D/2)^2;

N_hist = [];
v_real_hist = [];

tol = 1e-2; % relative criterion
N_range = 1:0.1:60;
for i = 1:length(N_range)
    N = N_range(i);
    v_real = m_dot / (rho_w * N * A);
    % stockage
    N_hist(end+1) = N;
    v_real_hist(end+1) = v_real;
   
    if abs(v_real - v_req)/v_req < tol
        N_sol = ceil(N);
        v_real_sol = v_real;
        break
    end
end


Re_real = rho_w*v_real*D/mu_w;
Nu_real=0.0243*(Re_real^(0.8))*Pr_w^0.4;
h_w_real=Nu_real*lambda_w/D;
f = 0.316*Re_real^(-0.25);
DeltaP = f*(L_jacket/D)*(rho_w*v_real^2/2);

T_wc=T_bulk_w+q/h_w_real;  %Temperature at interface wall/coolant
R_tbc_f=t_opt/k_tbc;
q_TBC=(T_aw_t-T_bulk_w)./(R_tbc_f+R_wall_res+R_gas_res+1/h_w_real);
Q_TBC=0; %initialization 
for i = 1:length(x)-1
    dx = x(i+1) - x(i);
    % slope
    if x(i) < 0
        drdx = -tan(beta_conv);
    else
        drdx = tan(alpha_div);
    end
    ds = sqrt(1 + drdx^2) * dx;
    dA = 2*pi*r(i)*ds;
    Q_TBC = Q_TBC + q_TBC(i)*dA;
end
T_wh_tbc=T_aw_t-q_TBC.*(R_gas_res+R_tbc_f); %Temperature at interface wall/TBC
T_tbc_gas=T_aw_t-q_TBC.*(R_gas_res); %Temperature at interface gas/TBC
fprintf('T_wc = %.4f K\n', max(T_wc));
fprintf('T_wh_tbc = %.4f K\n', max(T_wh_tbc));
fprintf('T_tbc_gas = %.2f K\n', max(T_tbc_gas));
fprintf('v_final_req = %.4f m/s\n', v_final_req);
fprintf('v_real = %.4f m/s\n', v_real);
fprintf('m_dot = %.4f kg/s\n', m_dot);


%% ============================================================
%  SECTION 10 - CASE B: WITH TBC (YSZ) cyl model
%  Follow the same workflow as CASE A but using a cylindrical model for the resistance of the wall, comparison is made in the plots
%  section
% ============================================================
R_s_tbc_CuCrZr_cyl=min(R_wall_res_cyl)+min(R_gas_res)+R_tbc_1;
q_TBC_CuCrZr_cyl=(max(T_aw_t)-T_bulk_w)./R_s_tbc_CuCrZr_cyl;
T_wh_TBC_CuCrZr_cyl=max(T_aw_t)-q_TBC_CuCrZr_cyl.*(min(R_gas_res)+R_tbc_1);

R_s_notbc_CuCrZr_cyl=R_wall_res_cyl+(R_gas_res);
q_noTBC_cyl=(T_aw_t-T_bulk_w)./R_s_notbc_CuCrZr_cyl;
T_wh_noTBC_CuCrZr_cyl=T_aw_t-q_noTBC_cyl.*(R_gas_res);

valid_cyl = (T_wh_TBC_CuCrZr_cyl < T_wall_max/1.2);
if ~any(valid_cyl)
    error('No TBC thickness satisfies both thermal constraints!');
end
idx_cyl = find(valid_cyl, 1, 'first');
t_tbc_opt_cyl = t_tbc(idx_cyl);
fprintf('Optimal TBC thickness cyl = %.4f micrometers\n', t_tbc_opt_cyl);

R_tbc_opt_cyl = (t_tbc_opt_cyl*10^-6) / k_tbc;
R_tot_tbc_CuCrZr_opt_cyl=min(R_wall_res_cyl)+min(R_gas_res)+R_tbc_opt_cyl;
q_TBC_CuCrZr_opt_cyl=(max((T_aw_t))-T_bulk_w)./R_tot_tbc_CuCrZr_opt_cyl;
Q_cyl=q_TBC_CuCrZr_opt_cyl*A_jacket;

t_cyl = t_tbc_opt_cyl*1e-6;
v_req_cyl=100;
t_hist_cyl = [];
v_hist_cyl = [];
while (v_req_cyl > 15) && (t_cyl< t_max)

    R_tbc_cyl = t_cyl / k_tbc;
    R_tot_cyl = min(R_gas_res) + R_tbc_cyl + mean(R_wall_res_cyl);

    q_cyl = (max(T_aw_t) - T_bulk_w) / R_tot_cyl;
    Q_cyl=q_cyl*A_jacket;

    h_req_cyl = q_cyl / (T_sat/1.2 - T_bulk_w);

    Nu_cyl = h_req_cyl * D / lambda_w;
    Re_cyl = (Nu_cyl/(0.0243*Pr_w^0.4))^(1/0.8);

    v_req_cyl = Re_cyl * mu_w / (rho_w * D);

    t_cyl = t_cyl + dt;
    
    t_hist_cyl(end+1) = t_cyl*1e6; % en um
    v_hist_cyl(end+1) = v_req_cyl;

end

t_opt_cyl = t_cyl;
v_final_req_cyl = v_req_cyl;

m_dot_cyl = Q_cyl/(cp_w*(T_out_w - T_in_w));
A = pi*(D/2)^2;
N_hist_cyl = [];
v_real_hist_cyl = [];

for i = 1:length(N_range)
    N_cyl = N_range(i);
    v_real_cyl = m_dot_cyl / (rho_w * N_cyl * A);
    N_hist_cyl(end+1) = N;
    v_real_hist_cyl(end+1) = v_real;
    if abs(v_real_cyl - v_final_req_cyl)/v_final_req_cyl < tol
        N_sol_cyl = N_cyl;
        v_real_sol_cyl = v_real_cyl;
        break
    end
end

Re_real_cyl = rho_w*v_real_cyl*D/mu_w;
Nu_real_cyl=0.0243*(Re_real_cyl^(0.8))*Pr_w^0.4;
h_w_real_cyl=Nu_real_cyl*lambda_w/D;
f_cyl = 0.316*Re_real_cyl^(-0.25);
DeltaP_cyl = f_cyl*(L_jacket/D)*(rho_w*v_real_cyl^2/2);

T_wc_cyl=T_bulk_w+q_cyl/h_w_real_cyl;
R_tbc_f_cyl=t_opt_cyl/k_tbc;
q_TBC_cyl=(T_aw_t-T_bulk_w)./(R_tbc_f_cyl+R_wall_res_cyl+R_gas_res+1/h_w_real_cyl);
T_wh_tbc_cyl=T_aw_t-q_TBC_cyl.*(R_gas_res+R_tbc_f_cyl);
T_tbc_gas_cyl=T_aw_t-q_TBC_cyl.*(R_gas_res);
fprintf('T_wc_cyl = %.4f K\n', max(T_wc_cyl));
fprintf('T_wh_tbc_cyl = %.4f K\n', max(T_wh_tbc_cyl));
fprintf('T_tbc_gas_cyl = %.2f K\n', max(T_tbc_gas_cyl));
fprintf('v_final_req_cyl = %.4f m/s\n', v_final_req_cyl);
fprintf('v_real_cyl = %.4f m/s\n', v_real_cyl);
fprintf('m_dot_cyl = %.4f kg/s\n', m_dot_cyl);

%% ============================================================
%  SECTION 11 - CASE C: WITH TBC (YSZ) model, Bartz
%  Follow the same workflow as CASE A/B but using Bartz model to determine the convective heat transfer on gas side, 
%  Comparison is made in the plots section
% ============================================================
R_s_tbc_CuCrZr_b=min(R_wall_res)+max(R_gas_res_b)+R_tbc_1;
q_TBC_CuCrZr_b=(max(T_aw_t)-T_bulk_w)./R_s_tbc_CuCrZr_b;
T_wh_TBC_CuCrZr_b=max(T_aw_t)-q_TBC_CuCrZr_b.*(min(R_gas_res_b)+R_tbc_1);

R_s_notbc_CuCrZr_b=min(R_wall_res)+max(R_gas_res_b);
q_noTBC_b=(T_aw_t-T_bulk_w)./R_s_notbc_CuCrZr_b;
T_wh_noTBC_CuCrZr_b=T_aw_t-q_noTBC_b.*(R_gas_res_b);

valid_b = (T_wh_TBC_CuCrZr_b < T_wall_max/1.2);
if ~any(valid_b)
    error('No TBC thickness satisfies both thermal constraints!');
end
idx_b = find(valid_b, 1, 'first');
t_tbc_opt_b = t_tbc(idx_b);
fprintf('Optimal TBC thickness b = %.4f micrometers\n', t_tbc_opt_b);

R_tbc_opt_b = (t_tbc_opt_b*10^-6) / k_tbc;
R_tot_tbc_CuCrZr_opt_b=min(R_wall_res)+min(R_gas_res_b)+R_tbc_opt_b;
q_TBC_CuCrZr_opt_b=(max((T_aw_t))-T_bulk_w)./R_tot_tbc_CuCrZr_opt_b;
Q_b=q_TBC_CuCrZr_opt_b*A_jacket;

t_b = t_tbc_opt_b*1e-6;
v_req_b=100;
t_hist_b = [];
v_hist_b = [];
while (v_req_b > 15) && (t_b< t_max)

    R_tbc_b = t_b / k_tbc;
    R_tot_b = min(R_gas_res_b) + R_tbc_b + mean(R_wall_res);

    q_b = (max(T_aw_t) - T_bulk_w) / R_tot_b;
    Q_b=q_b*A_jacket;

    h_req_b = q_b / (T_sat/1.2 - T_bulk_w);

    Nu_b = h_req_b * D / lambda_w;
    Re_b = (Nu_b/(0.0243*Pr_w^0.4))^(1/0.8);

    v_req_b = Re_b * mu_w / (rho_w * D);

    t_b = t_b + dt;
    
    t_hist_b(end+1) = t_b*1e6; % en um
    v_hist_b(end+1) = v_req_b;

end

t_opt_b = t_b;
v_final_req_b = v_req_b;

m_dot_b = Q_b/(cp_w*(T_out_w - T_in_w));
A = pi*(D/2)^2;
N_hist_b = [];
v_real_hist_b = [];

for i = 1:length(N_range)
    N_b = N_range(i);
    v_real_b = m_dot_b / (rho_w * N_b * A);
    N_hist_b(end+1) = N;
    v_real_hist_b(end+1) = v_real;
    if abs(v_real_b - v_final_req_b)/v_final_req_b < tol
        N_sol_b = N_b;
        v_real_sol_b = v_real_b;
        break
    end
end

Re_real_b = rho_w*v_real_b*D/mu_w;
Nu_real_b=0.0243*(Re_real_b^(0.8))*Pr_w^0.4;
h_w_real_b=Nu_real_b*lambda_w/D;
f_b = 0.316*Re_real_b^(-0.25);
DeltaP_b = f_b*(L_jacket/D)*(rho_w*v_real_b^2/2);

T_wc_b=T_bulk_w+q_b/h_w_real_b;
R_tbc_f_b=t_opt_b/k_tbc;
q_TBC_b=(T_aw_t-T_bulk_w)./(R_tbc_f_b+R_wall_res+R_gas_res_b+1/h_w_real_b);
T_wh_tbc_b=T_aw_t-q_TBC_b.*(R_gas_res_b+R_tbc_f_b);
T_tbc_gas_b=T_aw_t-q_TBC_b.*(R_gas_res_b);
fprintf('T_wc_b = %.4f K\n', max(T_wc_b));
fprintf('T_wh_tbc_b = %.4f K\n', max(T_wh_tbc_b));
fprintf('T_tbc_gas_b = %.2f K\n', max(T_tbc_gas_b));
fprintf('v_final_req_b = %.4f m/s\n', v_final_req_b);
fprintf('v_real_b = %.4f m/s\n', v_real_b);
fprintf('m_dot_b = %.4f kg/s\n', m_dot_b);



%% ============================================================
%  SECTION 12 - PLOTS
% ============================================================

figure;
plot(t_tbc, T_wh_TBC_CuCrZr, 'LineWidth', 2); hold on
plot(t_tbc, T_wh_TBC_CuCrZr_b, 'LineWidth', 2); hold on
plot(t_tbc, T_wh_TBC_CuCrZr_cyl, 'LineWidth', 2); 
legend('Tinterface CuCrZr/TBC', 'Tinterface CuCrZr/TBC Bartz', 'Tinterface CuCrZr/TBC cylindrical model'); 
yline(T_wall_max/1.2, '--', 'T hot limit');
xlabel('TBC thickness [um]')
ylabel('Temperature [K]');
title('Hot Wall temperatures vs TBC thickness');
grid on

%Temperature tbc/gas cyl/plane/bartz
figure;
plot(x, T_tbc_gas, 'LineWidth', 2); hold on
plot(x, T_tbc_gas_b, 'LineWidth', 2); hold on
plot(x, T_tbc_gas_cyl, 'LineWidth', 2); 
legend('Tinterface gas/TBC', 'Tinterface gas/TBC Bartz','Tinterface gas/TBC cylindrical model'); 
yline(1600, '--', 'T hot limit');
xlabel('x [m]');
ylabel('Temperature [K]');
title('tbc Wall temperatures over the domain');
grid on

%Profil 
figure;
plot(x, r, 'LineWidth',2);
xlabel('x [m]');
ylabel('r [m]');
title('Profil between A/A* = 2 boundaries');
grid on

%Temperature wall_hot profil no TBC
figure;
plot(x, T_wh_noTBC_CuCrZr, 'LineWidth',2); hold on
plot(x, T_wh_noTBC_CuCrZr_b, 'LineWidth',2); hold on
plot(x, T_wh_noTBC_CuCrZr_cyl, 'LineWidth',2);
yline(T_wall_max/1.2, '--', 'T hot limit');
legend('Tinterface gas/metal','Tinterface gas/metal Bartz', 'Tinterface gas/metal cylindrical model');
xlabel('x [m]');
ylabel('T [k]');
title('Temperature wall hot no TBC profil between A/A* = 2 boundaries');
grid on

%Temperature wall_hot profil TBC plane/cylindrical/bartz
figure;
plot(x, T_wh_tbc, 'LineWidth',2); hold on
plot(x, T_wh_tbc_b, 'LineWidth',2); hold on
plot(x, T_wh_tbc_cyl, 'LineWidth',2);
legend('Tinterface TBC/metal','Tinterface TBC/metal Bartz', 'Tinterface TBC/metal cylindrical model');
xlabel('x [m]');
ylabel('T [k]');
grid on

%Mach number
figure;
plot(x, M, 'LineWidth',2);
xlabel('x [m]');
ylabel('Mach number');
title('Mach distribution between A/A* = 2 boundaries');
grid on

%Heat flux with TBC CuCrZr wall plane/cyl model/bartz
figure;
%plot(x, q_noTBC./10^6, 'LineWidth',2); hold on
%plot(x, q_noTBC_b./10^6, 'LineWidth',2); hold on
%plot(x, q_noTBC_cyl./10^6, 'LineWidth',2); hold on
plot(x, q_TBC./10^6, 'LineWidth',2); hold on 
plot(x, q_TBC_b./10^6, 'LineWidth',2); hold on 
plot(x, q_TBC_cyl./10^6, 'LineWidth',2);
legend('heat flux plane TBC','heat flux plane TBC Bartz','heat flux cylindrical TBC')
xlabel('x [m]');
ylabel('Convective Heat flux [MW/m^2]');
grid on

%Convergence of required coolant velocity v_req as a function of TBC thickness
figure;
plot(t_hist, v_hist, 'LineWidth', 2); hold on
plot(t_hist_b, v_hist_b, 'LineWidth', 2); hold on
plot(t_hist_cyl, v_hist_cyl, 'LineWidth', 2); hold on
yline(15, '--r')
title('Convergence of required coolant velocity $v_{req}$ as a function of TBC thickness', ...
      'Interpreter','latex')
xlabel('TBC thickness [\mum]')
ylabel('v_{req} [m/s]')
legend({'$v_{req, plane}$','$v_{req, Bartz}$', '$v_{req, cyl}$','$v_{lim}=15$'},'Interpreter','latex', 'Location','best')
grid on

%Convergence of v_real toward v_{req} with tube count
figure;
plot(N_hist, v_real_hist, 'LineWidth', 2) ;hold on
yline(v_final_req_cyl, '--r')
xlabel('Number of tubes $N$', 'Interpreter','latex')
ylabel('Velocity [m/s]', 'Interpreter','latex')
title('Convergence of $v_{real}$ toward $v_{req}$ with tube count', ...
      'Interpreter','latex')
legend({'$v_{real, plane}$',}, ...
       'Interpreter','latex', 'Location','best')
grid on
