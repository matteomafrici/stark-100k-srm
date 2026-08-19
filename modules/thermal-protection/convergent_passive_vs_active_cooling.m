%% DEF
%  ============================================================
%  STARK-100K MOTOR THERMAL ANALYSIS - ENTRANCE SECTION -> AR 2
% ============================================================

%  === RESULT: AIR-ONLY CASE (PASSIVE) ===
% Required TBC thickness (on CuCrZr): 8.668 mm
% === RESULT: AIR-ONLY CASE (PASSIVE) ===
% Required TBC thickness (on Inconel 718): 1.972 mm
% === RESULT: CASE COOLING JACKET (ACTIVE) ===
% Required TBC thickness (on CuCrZr): 0.293 mm
% === RESULT: CASE COOLING JACKET (ACTIVE) ===
% Required TBC thickness (on Inconel 718): 0.197 mm
% IT IS CLEARLY EVIDENT THAT, WITH THE COOLING JACKET, THE REQUIRED TBC THICKNESS 
% REMAINS WELL WITHIN THE MAXIMUM ALLOWABLE LIMIT OF 1 MM

clc; clear all; close all;

%% SECTION 1 - GAS PROPERTIES (Input from CEA)
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

h_g_In = 0.026 / (2*r_t)^0.2 * (2).^0.1 * mu_0^0.2 * cp_g / Pr_g^0.6 * (p_c / c_star)^0.8 * (At_A)^0.9 * (rho_ref_In/rho_e)^0.8 * (mu_ref_In/mu_0)^0.2 ;

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

% Axial lengths of the convergent and divergent sections
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

%% SECTION 3 - WALL MATERIALS
t_wall = 0.0068;             % 6.8 mm
t_wall_1 = 0.0117;          % 11.7 mm (Copper thickness required to withstand chamber pressure)
k_tbc  = 0.6;               % RSZ coating (value from course reference material)

% CuCrZr 
k_wall_Cu = 320;
T_target_Cu = 600; 
R_wall_Cu = t_wall_1 / k_wall_Cu;

% Inconel 718 
k_wall_In = 27;    
T_target_In = 1300;  
R_wall_In = t_wall / k_wall_In;

%% SECTION 4 - CASE A: ONLY AIR (PASSIVE)
% To obtain a TBC thickness of a few mm without water:
h_aria_ext = 500;           % Forced convection + radiation [W/m^2K]
T_air = 288;
R_air_res = 1 / h_aria_ext;

% Gas-metal interface temperature function for passive (air-only) cooling
T_int_aria_Cu = @(R_tbc) Taw - (((Taw - T_air) / ...
    (R_gas_res_Cu + R_tbc + R_wall_Cu + R_air_res)) * (R_gas_res_Cu + R_tbc));

% Solution TBC in air
if T_int_aria_Cu(0) < T_target_Cu     % Check whether the material withstands 
    % the temperature without TBC
    t_tbc_aria_Cu = 0;
else
    % Determine the thickness required to prevent copper from melting.
    % 0.0002 = t_TBC_min / k_TBC = 100 micrometres / 0.5 = minimum thermal resistance,
    % which is used as the starting point for the iterations.
    R_tbc_aria_Cu = fsolve(@(R_tbc) T_int_aria_Cu(R_tbc) - T_target_Cu, 0.0002);
    t_tbc_aria_Cu = R_tbc_aria_Cu * k_tbc * 1000;
end

q_tbc_Cu = (Taw - T_air) / (R_gas_res_Cu + R_tbc_aria_Cu + R_wall_Cu + R_air_res);

T_interfaccia_metallo = Taw - q_tbc_Cu * (R_gas_res_Cu + R_tbc_aria_Cu);
T_superficie_TBC = Taw - q_tbc_Cu * R_gas_res_Cu;
T_parete_esterna = T_interfaccia_metallo - q_tbc_Cu * R_wall_Cu;

% fprintf('Gas-side external temperature: %.1f K\n', T_superficie_TBC);
% fprintf('TBC-metal interface temperature: %.1f K\n', T_interfaccia_metallo);

fprintf('=== RESULT: AIR-ONLY CASE WITH CuCrZr (PASSIVE) ===\n');
fprintf('Required TBC thickness (on CuCrZr): %.3f mm\n', t_tbc_aria_Cu);
% fprintf('Heat flux (on CuCrZr): %.3f MW/m^2\n', q_tbc_Cu/1e6);

% --- INTEGRITY VERIFICATION OF THE TBC WITH CuCrZr (PASSIVE CASE) ---
T_max_TBC = 1600; % file Durable RSZ/YSZ coatings with 1600  C thermal 
% shock resistance for rocket engine thrust chamber applications

%  Gas-side surface temperature
T_s_tbc_aria = Taw - (q_tbc_Cu * R_gas_res_Cu);

fprintf('--- MATERIAL VERIFICATION WITH CuCrZr (AIR) ---\n');
fprintf('TBC surface temperature: %.1f K\n', T_s_tbc_aria);

if T_s_tbc_aria > T_max_TBC
    fprintf('STATUS: TBC DEGRADES. The surface temperature (%.1f K) ', T_s_tbc_aria);
    fprintf('exceeds the material limit (%.0f K). The %.2f-mm thickness cannot be used without TBC degradation.\n', T_max_TBC, t_tbc_aria_Cu);
else
    fprintf('STATUS: VERIFIED.\n');
end

% Temperature function at the gas-metal interface in air
T_int_aria_In = @(R_tbc) Taw - (((Taw - T_air) / ...
    (R_gas_res_In + R_tbc + R_wall_In + R_air_res)) * (R_gas_res_In + R_tbc));

% Solution TBC in air
if T_int_aria_In(0) < T_target_In % I verify whether, without any 
    % thermal protection (i.e., with a TBC thickness of zero), the material 
    % is already capable of withstanding the gas temperature
    t_tbc_aria_In = 0;
else
    % Determine the thickness required to prevent Inconel from melting
    R_tbc_aria_In = fsolve(@(R_tbc) T_int_aria_In(R_tbc) - T_target_In, 0.0002);
    t_tbc_aria_In = R_tbc_aria_In * k_tbc * 1000;
end

q_tbc_In = (Taw - T_air) / (R_gas_res_In + R_tbc_aria_In + R_wall_In + R_air_res);

fprintf('=== RESULT: AIR-ONLY CASE WITH INCONEL (PASSIVE) ===\n');
fprintf('Required TBC thickness (on Inconel 718): %.3f mm\n', t_tbc_aria_In);
% fprintf('Heat flux (on Inconel 718): %.3f MW/m^2\n', q_tbc_In/1e6);

T_interfaccia_metallo_In = Taw - q_tbc_In * (R_gas_res_In + R_tbc_aria_In);
T_superficie_TBC_In = Taw - q_tbc_In * R_gas_res_In;
T_parete_esterna_In = T_interfaccia_metallo_In - q_tbc_In * R_wall_In;


% --- INTEGRITY VERIFICATION OF THE TBC WITH INCONEL (PASSIVE CASE)  ---

% Gas-side surface temperature
T_s_tbc_aria_In = Taw - (q_tbc_In * R_gas_res_In);

fprintf('--- MATERIAL VERIFICATION (AIR) WITH INCONEL 718 ---\n');
fprintf('TBC surface temperature: %.1f K\n', T_s_tbc_aria_In);

if T_s_tbc_aria_In > T_max_TBC
    fprintf('STATUS: TBC DEGRADES. The surface temperature (%.1f K) ', T_s_tbc_aria_In);
    fprintf('exceeds the TBC limit (%.0f K). The %.2f-mm thickness cannot be used without TBC degradation.\n', T_max_TBC, t_tbc_aria_In);
else
    fprintf('STATUS: VERIFIED.\n');
end

fprintf('TBC-metal interface temperature: %.2f K\n', T_interfaccia_metallo_In);

%% TBC THICKNESS OPTIMIZATION
t_vettore = linspace(0, 0.008, 500); % Scan from 0 to 8 mm
R_tbc_vettore = t_vettore / k_tbc;

% Temperature storage vectors
T_metallo = zeros(size(t_vettore));
T_superficie_TBC = zeros(size(t_vettore));

for i = 1:length(R_tbc_vettore)
    R_tot_Cu = R_gas_res_Cu + R_tbc_vettore(i) + R_wall_Cu + R_air_res;
    q = (Taw - T_air) / R_tot_Cu;
    
    T_metallo(i) = Taw - q * (R_gas_res_Cu + R_tbc_vettore(i));
    T_superficie_TBC(i) = Taw - q * R_gas_res_Cu;
end

% Find the indices that satisfy both constraints
indici_validi = find(T_metallo < T_target_Cu & T_superficie_TBC < T_max_TBC);

if isempty(indici_validi)
    fprintf('WARNING: No thickness satisfies both constraints!\n');
    % If this occurs, h_w (water velocity) must be increased or the material must be changed
else
    t_ottimale_mm = t_vettore(indici_validi(1)) * 1000; % I take the minimum required
    fprintf('Suggested TBC thickness: %.3f mm\n', t_ottimale_mm);
end 

%% SECTION 5 - COOLANT (TAP WATER) PROPERTIES

fprintf ('\n\nFROM HERE ON WE ASSUME THAT WE ALSO HAVE COOLING JACKET IN THE CONVERGENT PART')

p_water  = 0.15 * p_c;     
T_in_w   = 18 + 273.15;    
T_out_w  = 80 + 273.15;    
T_bulk_w = (T_in_w + T_out_w) / 2;  

rho_w   = 988;              
cp_w    = 4181;             
mu_w    = 5.5e-4;           
lambda_w= 0.644;            
Pr_w    = mu_w * cp_w / lambda_w;  

%% SECTION 6 - WATER-SIDE HEAT TRANSFER (h_w)
gap_w   = 0.005;            
D_h_w   = 2 * gap_w;       
v_w_guess = 5.0;            % Water speed [m/s]

Re_w = rho_w * v_w_guess * D_h_w / mu_w;
Nu_w = 0.023 * Re_w^0.8 * Pr_w^0.4;
h_w  = Nu_w * lambda_w / D_h_w;    
R_w_res = 1 / h_w;

%% SECTION 7 - CASE B: COOLING JACKET (ACTIVE)

p_water  = 0.15 * p_c;     
T_in_w   = 18 + 273.15;    
T_out_w  = 80 + 273.15;    
T_bulk_w = (T_in_w + T_out_w) / 2;  

% Calculation of the active TBC for CuCrZr
T_int_Cu_w = @(R_tbc) Taw - ((Taw - T_bulk_w) / (R_gas_res_Cu + R_tbc + R_wall_Cu + R_w_res)) * (R_gas_res_Cu + R_tbc);
if T_int_Cu_w(0) < T_target_Cu
    R_sol_Cu_w = 0; 
    t_tbc_Cu_mm = 0;
else 
    R_sol_Cu_w = fsolve(@(R) T_int_Cu_w(R) - T_target_Cu, 0.00002); 
    t_tbc_Cu_mm = R_sol_Cu_w * k_tbc * 1000;
end

% Gas-metal interface temperature function with water cooling
T_int_acqua_Cu = @(R_tbc) Taw - (((Taw - T_bulk_w) / ...
    (R_gas_res_Cu + R_tbc + R_wall_Cu + R_w_res)) * (R_gas_res_Cu + R_tbc));

% TBC solution with cooling jacket
if T_int_acqua_Cu(0) < T_target_Cu % Verify whether, without any thermal protection 
% (i.e., with a TBC thickness equal to zero), the material is already able 
% to withstand the gas temperature
    t_tbc_acqua_Cu = 0;
else
    % Find the thickness required to prevent the Inconel from melting
    R_tbc_acqua_Cu = fsolve(@(R_tbc) T_int_acqua_Cu(R_tbc) - T_target_Cu, 0.0002);
    t_tbc_acqua_Cu = R_tbc_acqua_Cu * k_tbc * 1000;
end

q_tbc_Cu = (Taw - T_bulk_w) / (R_gas_res_Cu + R_tbc_acqua_Cu + R_wall_Cu + R_w_res);
% --- TBC INTEGRITY VERFICATION WITH INCONEL AND COOLING JACKET ---

% Gas-side surface temperature calculation
T_s_tbc_acqua_Cu = Taw - (q_tbc_Cu * R_gas_res_Cu);

fprintf('--- MATERIAL VERIFICATION (WATER) WITH INCONEL 718 ---\n');
fprintf('TBC surface temperature: %.1f K\n', T_s_tbc_acqua_Cu);

if T_s_tbc_acqua_Cu > T_max_TBC
    fprintf('STATUS: TBC DEGRADES. The surface temperature (%.1f K) ', T_s_tbc_acqua_Cu);
     fprintf('exceeds the TBC limit (%.0f K). The %.2f mm thickness cannot be used without TBC degradation.\n', T_max_TBC, t_tbc_acqua_Cu);
else
     fprintf('STATUS: VERIFIED.\n');
end


 fprintf('=== RESULT: COOLING JACKET CASE (ACTIVE) ===\n');
 fprintf('Required TBC thickness (on CuCrZr): %.3f mm\n', t_tbc_Cu_mm);

% Calculation of the active TBC for Inconel
T_int_In_w = @(R_tbc) Taw - ((Taw - T_bulk_w) / (R_gas_res_In + R_tbc + R_wall_In + R_w_res)) * (R_gas_res_In + R_tbc);
if T_int_In_w(0) < T_target_In
    R_sol_In_w = 0; 
    t_tbc_In_mm = 0;

else 
    R_sol_In_w = fsolve(@(R) T_int_In_w(R) - T_target_In, 0.00002); 
    t_tbc_In_mm = R_sol_In_w * k_tbc * 1000; 

end

% Gas-metal interface temperature function with water cooling
T_int_acqua_In = @(R_tbc) Taw - (((Taw - T_bulk_w) / ...
    (R_gas_res_In + R_tbc + R_wall_In + R_w_res)) * (R_gas_res_In + R_tbc));

% TBC solution with cooling jacket
if T_int_acqua_In(0) < T_target_In % Verify whether, without any thermal protection 
% (i.e., with a TBC thickness equal to zero), the material is already able 
% to withstand the gas temperature
    t_tbc_acqua_In = 0;
else
    % Find the thickness required to prevent the Inconel from melting
    R_tbc_acqua_In = fsolve(@(R_tbc) T_int_acqua_In(R_tbc) - T_target_In, 0.0002);
    t_tbc_acqua_In = R_tbc_acqua_In * k_tbc * 1000;
end

q_tbc_In = (Taw - T_bulk_w) / (R_gas_res_In + R_tbc_acqua_In + R_wall_In + R_w_res);
% --- TBC INTEGRITY VERFICATION WITH INCONEL AND COOLING JACKET ---

% Gas-side surface temperature calculation
T_s_tbc_acqua_In = Taw - (q_tbc_In * R_gas_res_In);

fprintf('--- MATERIAL VERIFICATION (WATER) WITH INCONEL 718 ---\n');
fprintf('TBC surface temperature: %.1f K\n', T_s_tbc_acqua_In);

if T_s_tbc_acqua_In > T_max_TBC
    fprintf('STATUS: TBC DEGRADES. The surface temperature (%.1f K) ', T_s_tbc_acqua_In);
     fprintf('exceeds the TBC limit (%.0f K). The %.2f mm thickness cannot be used without TBC degradation.\n', T_max_TBC, t_tbc_acqua_In);
else
     fprintf('STATUS: VERIFIED.\n');
end

fprintf('=== RESULT: COOLING JACKET CASE (ACTIVE) ===\n');
fprintf('Required TBC thickness (on Inconel 718): %.3f mm\n', t_tbc_In_mm);