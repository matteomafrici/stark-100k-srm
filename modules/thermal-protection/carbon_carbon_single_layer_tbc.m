% DEFINITIVO
% VARIOUS CARBON-COMPOSITE CONFIGURATIONS WITH A SINGLE TBC LAYER

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

k_wall_CC = [2, 5, 10, 11.5, 21, 54,80, 121.5, 180, 233];
T_max_CC = 2200;

t_wall_CC = 0.029;
R_CC = t_wall_CC ./ k_wall_CC;

k_tbc = 0.6;
T_max_TBC = 1600;

%% SECTION 4 - CASE A: ONLY AIR (PASSIVE)
% To obtain a TBC thickness of a few mm without water:
h_aria_ext = 500;           % Forced convection + radiation [W/m^2K]
T_air = 288;
R_air_res = 1 / h_aria_ext;


for i = 1:length(k_wall_CC)
    
    % Extraction of current values from vectors
    k_w = k_wall_CC(i);
    R_w = R_CC(i);
    
    % Definition of the TBC-metal interface temperature function
    % Taw - q * (R_gas + R_tbc)
    T_int_func = @(R_tbc) Taw - (((Taw - T_air) / ...
        (R_gas_res_Cu + R_tbc + R_w + R_air_res)) * (R_gas_res_Cu + R_tbc));


    % Check whether the TBC is required
    if T_int_func(0) < T_max_CC
        t_tbc_necessario = 0;
        R_tbc_sol = 0;
        fprintf('Status: TBC not required. The metal withstands passively.\n');
    else
        % Solving for the required thickness using fsolve
        % I start from a trial resistance (es. 100 micron / k_tbc)
        options = optimset('Display','off');
        R_tbc_sol = fsolve(@(R_tbc) T_int_func(R_tbc) - T_max_CC, 0.0002, options);
        t_tbc_necessario = R_tbc_sol * k_tbc * 1000; % Convertion in mm
    end

    % Post-solution parameter calculation
    q_final = (Taw - T_air) / (R_gas_res_Cu + R_tbc_sol + R_w + R_air_res);
    T_s_tbc = Taw - (q_final * R_gas_res_Cu);
    T_int_metallo = Taw - q_final * (R_gas_res_Cu + R_tbc_sol);

    % Results
    fprintf('Required TBC thickness: %.3f mm\n', t_tbc_necessario);
    fprintf('TBC surface temperature: %.1f K\n', T_s_tbc);
    fprintf('TBC-metal interface temperature: %.1f K\n', T_int_metallo);

    % TBC integrity check
    if T_s_tbc > T_max_TBC
        fprintf('WARNING: TBC is degrading! T_sup (%.1f K) > Limit (%.0f K)\n', T_s_tbc, T_max_TBC);
    else
        fprintf('TBC status: integrity verified.\n');
    end
    fprintf('-------------------------------------------\n');
end


