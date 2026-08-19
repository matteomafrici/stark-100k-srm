%% DEF
% BILAYER INCONEL718 USING THE RESISTANCE ANALOGY 

clc; clear all; close all;

%% SECTION 1 - GAS PROPERTIES (Input from RPA/CEA)
p_c     = 6.5e6;             % Chamber pressure [Pa]
Tc      = 2836.99;           % Chamber temperature [K]
gamma   = 1.1871;           
R_gas   = 8314 / 25.0;   

% Cp, thermal conductivity and Prandtl number refer to
% equilibrium-reaction properties
cp_g    = 2388;                        
mu_g    = 9.2063e-5;        
lambda_g=        0.38519;         
Pr_g    = mu_g * cp_g / lambda_g;  
M       = 0.313;            % Local Mach number in the convergent part
r_t  = 0.1145 /2 ;               % Throat radius [m]
D_t  = 2 * r_t;              % Throat diameter [m]
A_t  = pi/4 * D_t^2;        % Throat area [m^2]

% The system is oversized considering the area where A/A_t = 2, 
% which experiences the highest heat flux

% Local gas properties
T_gas = Tc / (1 + (gamma-1)/2 * M^2);     
% Density (ideal gas)
p  = p_c / (1 + (gamma-1)/2 * M^2)^(gamma/(gamma-1));  % [Pa]
rho= p / (R_gas * T_gas);                    % [kg/m^3]
mu_0 = 9.2624e-5;    % IN COMBUSTION CHAMBER, AT T0=Tc, NOT OF THE GAS IN THE POINT WITH EPS = 2
c_star = 1509.9;
At_A = 0.5;
rho_e = 6.4843;
S = 115;
Te = 2812.44;
r = (1+Pr_g^(1/3) * (gamma-1)/2 * M^2)/(1+ (gamma-1)/2 * M^2);
Taw = Te*(1+r * (gamma-1)/2 * M^2);

Tw_Cu = 573;    % Assumption
Tw_In = 1300;   % Assumption
T_ref_Cu = 0.5 * Tw_Cu + 0.5 * Te + 0.22 * (Taw - Te);
rho_ref_Cu = rho_e / (Te / T_ref_Cu);

mu_ref_Cu = mu_0 * (T_ref_Cu/Tc)^1.5 * (Tc+S)/(T_ref_Cu+S);

h_g_Cu = 0.026 / (2*r_t)^0.2 * (2).^0.1 * mu_0^0.2 * cp_g / Pr_g^0.6 * (p_c / c_star)^0.8 * (At_A)^0.9 * (rho_ref_Cu/rho_e)^0.8 * (mu_ref_Cu/mu_0)^0.2;

R_gas_res_Cu = 1 / h_g_Cu;

T_ref_In = 0.5 * Tw_In + 0.5 * Te + 0.22 * (Taw - Te);
mu_ref_In = mu_0 * (T_ref_In/Tc)^1.5 * (Tc+S)/(T_ref_In+S);
rho_ref_In = rho_e / (Te / T_ref_In);

h_g_In = 0.026 / (2*r_t)^0.2 * (2).^0.1 * mu_0^0.2 * cp_g / Pr_g^0.6 * (p_c / c_star)^0.8 * (At_A)^0.9 * (rho_ref_In/rho_e)^0.8 * (mu_ref_In/mu_0)^0.2 ;clc

R_gas_res_In = 1 / h_g_In;

%% SECTION 2 - GEOMETRY (Entrance Section -> AR=2)

% Area ratio epsilon = A/A_t = 2 -> radius at epsilon=2
A_eps2 = 2 * A_t;
r_eps2 = sqrt(A_eps2 / pi); % [m]

% Half-angles (from assignment)
alpha_div = 15;             % Divergence half-angle [deg]
beta_conv = 45;             % Convergence half-angle [deg]

% Only the convergent-section length from the inlet to A/A_t = 2 is needed:

L_tot_conv = 0.3784;
arearatio = 2;
L_in2 = L_tot_conv - r_t * (sqrt(arearatio)-1)/tan(beta_conv);
r_ing = r_t + L_tot_conv * tan(beta_conv);

% Axial lengths of convergent and divergent sections
L_conv = (r_eps2 - r_t) / tand(beta_conv);   % [m]
L_div  = (r_eps2 - r_t) / tand(alpha_div);   % [m]

% Slant lengths (frustum lateral surface)
s_conv = sqrt(L_conv^2 + (r_eps2 - r_t)^2);
s_div  = sqrt(L_div^2  + (r_eps2 - r_t)^2);

% Lateral (wetted) area of the cooling jacket (frustum geometry)
A_conv = pi * (r_eps2 + r_t) * s_conv;       % [m^2]
A_div  = pi * (r_eps2 + r_t) * s_div;        % [m^2]
A_jacket = A_conv + A_div;                   % [m^2] total wetted area

L_sez = (r_eps2 - r_t) / tand(beta_conv);
s_conv = sqrt(L_sez^2 + (r_eps2 - r_t)^2);
A_scambio = pi * (r_eps2 + r_t) * s_conv; 

%% SECTION 3 - WALL MATERIALS (Updated for Bilayer)
t_wall = 0.0068;           % 6.8 mm
k_rsz  = 0.6;               % RSZ (Top coat)
k_ysz  = 0.9;               % YSZ (Bond coat/Intermediate)

% CuCrZr 
k_wall_Cu = 320;
T_target_Cu = 600; 
R_wall_Cu = t_wall / k_wall_Cu;

% Inconel 718 
k_wall_In = 27;   
T_target_In = 1300;   
R_wall_In = t_wall / k_wall_In;

%% SECTION 4 - CASE A: BILAYER TBC (AIR/PASSIVE)
h_aria_ext = 500;           
T_air = 288;
R_air_res = 1 / h_aria_ext;

% --- LOGICA BILAYER ---
% Since we have two unknowns (t_rsz and t_ysz) but only one equation (T_target), 
% we must impose an additional condition. For example: we can fix the YSZ thickness 
% (e.g., 0.10 mm) and solve for the required RSZ thickness.
t_ysz_fisso = 0.00010; % 0.10 mm of YSZ (Bond coat standard)
R_ysz = t_ysz_fisso / k_ysz;

% New function: Gas-Metal Interface Temperature for a Bilayer (Inconel example)  
% T_int = T_aw - q * (R_gas + R_rsz + R_ysz)
T_int_bilayer_In = @(R_rsz) Taw - (((Taw - T_air) / ...
    (R_gas_res_In + R_rsz + R_ysz + R_wall_In + R_air_res)) * (R_gas_res_In + R_rsz + R_ysz));

%% SOLUTION FOR INCONEL 718
if T_int_bilayer_In(0) < T_target_In
    t_rsz_necessario = 0;
else
    % Find the required RSZ resistance
    R_rsz_opt = fsolve(@(R_rsz) T_int_bilayer_In(R_rsz) - T_target_In, 0.0002);
    t_rsz_necessario = R_rsz_opt * k_rsz; % [m]
end

% Calculation of Heat Flux and Layer Temperatures
R_tbc_tot = R_rsz_opt + R_ysz;
q_bilayer = (Taw - T_air) / (R_gas_res_In + R_tbc_tot + R_wall_In + R_air_res);

T_superficie_RSZ = Taw - q_bilayer * R_gas_res_In;
T_interfaccia_RSZ_YSZ = T_superficie_RSZ - q_bilayer * R_rsz_opt;
T_interfaccia_YSZ_metallo = T_interfaccia_RSZ_YSZ - q_bilayer * R_ysz;

%% RESULTS
fprintf('=== RESULT: BILAYER TBC ON INCONEL 718 ===\n');
fprintf('YSZ thickness (fixed):    %.3f mm\n', t_ysz_fisso * 1000);
fprintf('Required RSZ thickness:   %.3f mm\n', t_rsz_necessario * 1000);
fprintf('Total TBC thickness:       %.3f mm\n', (t_ysz_fisso + t_rsz_necessario) * 1000);
fprintf('Heat flux:            %.3f MW/m^2\n', q_bilayer/1e6);
fprintf('-------------------------------------------\n');
fprintf('Gas-side surface temperature:    %.1f K\n', T_superficie_RSZ);
fprintf('RSZ/YSZ interface temperature: %.1f K\n', T_interfaccia_RSZ_YSZ);
fprintf('Metal-side interface temperature: %.1f K\n', T_interfaccia_YSZ_metallo);

% --- INTEGRITY CHECK ---
T_max_TBC = 1600; 
if T_superficie_RSZ > T_max_TBC
    fprintf('STATUS: DEGRADATION. Surface temperature > %.0f K\n', T_max_TBC);
else
    fprintf('STATUS: VERIFIED.\n');
end