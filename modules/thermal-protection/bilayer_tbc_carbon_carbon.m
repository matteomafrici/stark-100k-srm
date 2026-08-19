% BILAYER TBC WITH C/C - Resistance model (DEFINITIVO)
% ===============================================================
clear; clc; close all;

% --- Figure styling: light background, dark text (matching design.m) ---
set(0, 'DefaultFigureColor',    [0.97 0.97 0.97], ...  % figure background: near-white
       'DefaultAxesColor',      [1.00 1.00 1.00], ...  % axes background: white
       'DefaultAxesXColor',     [0.15 0.15 0.15], ...  % X tick/label: dark grey
       'DefaultAxesYColor',     [0.15 0.15 0.15], ...  % Y tick/label: dark grey
       'DefaultAxesZColor',     [0.15 0.15 0.15], ...  % Z tick/label
       'DefaultAxesGridColor',  [0.15 0.15 0.15], ...  % grid lines
       'DefaultTextColor',      [0.15 0.15 0.15]);     % all text objects

%% FIXED INPUT 
p_c     = 6.5e6;             
Tc      = 2836.99;           
gamma   = 1.1871;           
R_gas   = 8314 / 25.0;   
cp_g    = 2388;                        
mu_g    = 9.2063e-5;        
lambda_g= 0.38519;        
Pr_g    = mu_g * cp_g / lambda_g;  
M       = 0.313;            
r_t     = 0.1145 / 2;               
c_star  = 1509.9;
At_A    = 0.5;
rho_e   = 6.4843;
S       = 115;
Taw     = 2.8382e+03;     
Tw_CC   = 2200;   
Te      = 2812.44;

% Calculation h_g (Bartz)
T_ref_CC = 0.5 * Tw_CC + 0.5 * Te + 0.22 * (Taw - Te);
rho_ref_CC = rho_e / (Te / T_ref_CC);
mu_0 = 9.2624e-5;
mu_ref_CC = mu_0 * (T_ref_CC/Tc)^1.5 * (Tc+S)/(T_ref_CC+S);
h_g_CC = 0.026 / (2*r_t)^0.2 * (2).^0.1 * mu_0^0.2 * cp_g / Pr_g^0.6 * (p_c / c_star)^0.8 * (At_A)^0.9 * (rho_ref_CC/rho_e)^0.8 * (mu_ref_CC/mu_0)^0.2;

R_gas_res = 1 / h_g_CC;
h_air    = 500;   
T_air    = 288;    
R_air_res = 1 / h_air;

% Material properties [k, rho, cp]
mat.rsz = [0.6,   6500,  620];
mat.ysz = [0.9,   6000,  640]; 

L_cc   = 0.029;    % Calculated with sigma_yield,CC = 150 MPa, reduced to 100 
% MPa with FS = 1.5
k_wall_CC_vec = [2, 5, 10,11.5, 21, 54, 80, 121.5, 180, 233];

% Thickness search range[m]
L_rsz_search = 0.0001:0.00005:0.0030; 
L_ysz_search = 0.0001:0.00005:0.0030; 
T_limit  = 2200;       

%% Parametric loop over k
results = []; 
fprintf('%-10s | %-10s | %-10s | %-10s\n', 'k_CC', 'RSZ [mm]', 'YSZ [mm]', 'TOT [mm]');
fprintf('------------------------------------------------------------\n');

for k_idx = 1:length(k_wall_CC_vec)
    current_k_cc = k_wall_CC_vec(k_idx);
    R_cc_res = L_cc / current_k_cc;
    
    min_total = inf;
    best_combo = [0, 0];
    found = false;
    
    for L1 = L_rsz_search
        R_rsz = L1 / mat.rsz(1);
        for L2 = L_ysz_search
            R_ysz = L2 / mat.ysz(1);
            
            % Total resistance calculation
            R_tot = R_gas_res + R_rsz + R_ysz + R_cc_res + R_air_res;
            
            % Heat flux q
            q = (Taw - T_air) / R_tot;
            
            % TBC/metal interface temperature
            Tint = Taw - q * (R_gas_res + R_rsz + R_ysz);
            
            if Tint <= T_limit
                if (L1 + L2) < min_total
                    min_total = L1 + L2;
                    best_combo = [L1, L2];
                    found = true;
                end
                break; 
            end
        end
    end
    
    if found
        fprintf('%-10.1f | %-10.2f | %-10.2f | %-10.2f\n', ...
            current_k_cc, best_combo(1)*1e3, best_combo(2)*1e3, min_total*1e3);
        results = [results; current_k_cc, best_combo(1)*1e3, best_combo(2)*1e3, min_total*1e3];
    else
        fprintf('%-10.1f | %-33s\n', current_k_cc, 'NO SOLUTION');
    end
end

%% PLOT
fig1 = figure('Color', 'w');
plot(results(:,1), results(:,4), '-sr', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on;
xlabel('Thermal conductivity C/C [W/mK]');
ylabel('Optimized thickness TBC (RSZ+YSZ) [mm]');
title('Stationary optimization (Resistance model)');

%% Save figure
results_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end
set(fig1, 'Color', 'w');
exportgraphics(fig1, fullfile(results_dir, 'bilayer_cc_thickness_vs_conductivity.png'), ...
    'Resolution', 200, 'BackgroundColor', 'white');
fprintf('Figure saved to: %s\n', results_dir);