clc; clear all;

%% 1. Physical data
T_c = 2892.5;                      
T_target_interfaccia = 1300; % Maximum temperature at the interface between RSZ/Metal
t_burn = 29;            

% RSZ properties
k_rsz = 0.6; 
rho_rsz = 6601; 
cp_rsz = 569.8;
alpha_rsz = k_rsz / (rho_rsz * cp_rsz);

% Inconel 718 properties
t_alloy = 0.0068;
k_alloy = 27;   
rho_alloy = 8193.25; 
cp_alloy = 435;
alpha_alloy = k_alloy / (rho_alloy * cp_alloy);

%% Thermal convective heat transfer
cp_g    = 2745.9;           % [J/kgK]
mu_g    = 9.3556e-5;        % [Pa*s]
lambda_g= 0.42194;          % [W/mK]
Pr_g    = mu_g * cp_g / lambda_g;  

D_p = 0.4738;               % Port Diameter [m]
m_dot = 43.7396;            % Mass flow rate [kg/s]

Re_g = (4 * m_dot)/(pi * D_p * mu_g); 
f = (0.790 * log(Re_g) - 1.64)^(-2); 
Nu_g = ((f/8) * (Re_g - 1000) * Pr_g)/(1 + 1.27 * (f/8)^(1/2) * (Pr_g^(2/3) - 1)); 
h_g = Nu_g * lambda_g/D_p; 

%% 2. Iterative cycle for the thickness of the thermal protection
s_min = 0.0001; s_max = 0.002;
s_rsz = (s_min + s_max) / 2;
tolleranza = 1e-3; 
diff = 100;

fprintf('Optimal thickness for T_target_interfac = %d K...\n', T_target_interfaccia);

while abs(diff) > tolleranza
    nx_rsz = 25; nx_alloy = 25;
    dx_rsz = s_rsz / (nx_rsz - 1);
    dx_alloy = t_alloy / (nx_alloy - 1);
    
    % Von Neumann stability criterion
    dt = 0.8 * min(dx_rsz^2/(2*alpha_rsz), dx_alloy^2/(2*alpha_alloy));
    nt = ceil(t_burn / dt);
    
    T = ones(1, nx_rsz + nx_alloy) * 293.15;
    
    for n = 1:nt
        T_old = T;
        % 1. Hot face RSZ
        T(1) = T_old(1) + 2*alpha_rsz*dt/dx_rsz^2 * (T_old(2) - T_old(1) + (dx_rsz*h_g/k_rsz)*(T_c - T_old(1)));
        
        % 2. Internal RSZ
        for i = 2:nx_rsz-1
            T(i) = T_old(i) + alpha_rsz*dt/dx_rsz^2 * (T_old(i+1) - 2*T_old(i) + T_old(i-1));
        end
        
       % 3. Interface
T(nx_rsz) = (k_rsz/dx_rsz * T_old(nx_rsz-1) + k_alloy/dx_alloy * T_old(nx_rsz+1)) / (k_rsz/dx_rsz + k_alloy/dx_alloy);

        % 4. Internal Inconel 718
        for i = nx_rsz+1 : nx_rsz+nx_alloy-1
            T(i) = T_old(i) + alpha_alloy*dt/dx_alloy^2 * (T_old(i+1) - 2*T_old(i) + T_old(i-1));
        end
        
        % 5. External face Inconel 718
        T(end) = T_old(end) + 2*alpha_alloy*dt/dx_alloy^2 * (T_old(end-1) - T_old(end));
    end
    
    % TARGET: temperature at the interface
    T_int_effettiva = T(nx_rsz);
    diff = T_int_effettiva - T_target_interfaccia;
    
    if diff > 0
        s_min = s_rsz;
    else
        s_max = s_rsz;
    end
    s_rsz = (s_min + s_max) / 2;
    
    if (s_max - s_min) < 1e-8, break; end
end

%% 3. Results
fprintf('\n================================================\n');
fprintf('         Results TBC\n');
fprintf('================================================\n');
fprintf('RSZ thickness:      %.4f mm\n', s_rsz * 1000);
fprintf('T top surface RSZ:      %.2f K\n', T(1));
fprintf('T interface (Target):      %.2f K\n', T(nx_rsz));
fprintf('T external surface Alloy:      %.2f K\n', T(end));
fprintf('Convective coefficient h_g: %.2f W/m^2K\n', h_g);
fprintf('================================================\n');


nx_rsz = 25; nx_alloy = 25;
dx_rsz = s_rsz / (nx_rsz - 1);
dx_alloy = t_alloy / (nx_alloy - 1);

dt = 0.8 * min(dx_rsz^2/(2*alpha_rsz), dx_alloy^2/(2*alpha_alloy));
nt = ceil(t_burn / dt);

T = ones(1, nx_rsz + nx_alloy) * 293.15;

T_int_time = zeros(1, nt);
T_surf_time = zeros(1, nt);
T_outer_time = zeros(1, nt);
time = (0:nt-1)*dt;

for n = 1:nt
    T_old = T;

    T(1) = T_old(1) + 2*alpha_rsz*dt/dx_rsz^2 * ...
        (T_old(2) - T_old(1) + (dx_rsz*h_g/k_rsz)*(T_c - T_old(1)));

    for i = 2:nx_rsz-1
        T(i) = T_old(i) + alpha_rsz*dt/dx_rsz^2 * ...
            (T_old(i+1) - 2*T_old(i) + T_old(i-1));
    end

    T(nx_rsz) = (k_rsz/dx_rsz * T(nx_rsz-1) + ...
                 k_alloy/dx_alloy * T(nx_rsz+1)) / ...
                (k_rsz/dx_rsz + k_alloy/dx_alloy);

    for i = nx_rsz+1 : nx_rsz+nx_alloy-1
        T(i) = T_old(i) + alpha_alloy*dt/dx_alloy^2 * ...
            (T_old(i+1) - 2*T_old(i) + T_old(i-1));
    end

    T(end) = T_old(end) + 2*alpha_alloy*dt/dx_alloy^2 * ...
             (T_old(end-1) - T_old(end));

    T_int_time(n) = T(nx_rsz);
    T_surf_time(n) = T(1);
    T_outer_time(n) = T(end); 
end

%% PLOT 1 - Final spatial profile
x_rsz   = linspace(0, s_rsz*1000, nx_rsz);
x_alloy = linspace(s_rsz*1000, (s_rsz+t_alloy)*1000, nx_alloy);
x_tot   = [x_rsz, x_alloy];

fig1 = figure('Color','white','Position',[100 100 900 400]);
hold on;
patch([x_rsz, fliplr(x_rsz)], [293.15*ones(1,nx_rsz), fliplr(T(1:nx_rsz))], ...
    [0.85 0.92 0.98], 'EdgeColor','none','FaceAlpha',0.4);
patch([x_alloy, fliplr(x_alloy)], [293.15*ones(1,nx_alloy), fliplr(T(nx_rsz+1:end))], ...
    [0.98 0.88 0.82], 'EdgeColor','none','FaceAlpha',0.4);
plot(x_tot, T, 'Color',[0.12 0.37 0.65], 'LineWidth', 2.2);
xline(s_rsz*1000, '--', 'Color',[0.5 0.5 0.5], 'LineWidth', 1.2, ...
    'Label','RSZ/Alloy Interface', 'LabelOrientation','horizontal','FontSize',9);
yline(T_target_interfaccia, ':', 'Color',[0.85 0.35 0.18], 'LineWidth', 1.5, ...
    'Label','T_{target} = 1300 K', 'LabelOrientation','horizontal','FontSize',9);
text(s_rsz*1000/2, T(1)*0.97, ['RSZ'], ...
    'HorizontalAlignment','center','FontSize',10,'Color',[0.12 0.37 0.65]);
text((s_rsz + t_alloy/2)*1000, T(end)+120, 'Alloy 718', ...
    'HorizontalAlignment','center','FontSize',10,'Color',[0.85 0.35 0.18]);
xlabel('Thickness of TBC [mm]','FontSize',11);
ylabel('Temperature [K]','FontSize',11);
title('Temperature profile at t = 25 s','FontSize',13,'FontWeight','normal');
grid on; box on;
set(gca,'FontSize',10,'GridAlpha',0.25,'LineWidth',0.8);
xlim([0, (s_rsz+t_alloy)*1000]);

%% PLOT 2 - Temperature transient profile
fig2 = figure('Color','white','Position',[100 100 900 420]);
hold on;
plot(time, T_surf_time,  'Color',[0.85 0.35 0.18], 'LineWidth', 2);
plot(time, T_int_time,   'Color',[0.12 0.37 0.65], 'LineWidth', 2);
plot(time, T_outer_time, 'Color',[0.11 0.62 0.46], 'LineWidth', 2);
yline(T_target_interfaccia, '--', 'Color',[0.85 0.35 0.18], ...
    'LineWidth',1.2,'Alpha',0.5, ...
    'Label','T_{target} = 1300 K','LabelOrientation','horizontal','FontSize',9);
xlabel('Time [s]','FontSize',11);
ylabel('Temperature [K]','FontSize',11);
title('Transient Temperature','FontSize',13,'FontWeight','normal');
legend({'Hot wall RSZ','Interface RSZ/Alloy','External face Alloy'}, ...
    'Location','northwest','FontSize',10,'Box','off');
grid on; box on;
set(gca,'FontSize',10,'GridAlpha',0.25,'LineWidth',0.8);
xlim([0, time(end)]);


%% PLOT 3 - Bisection method convergence
s_lo2  = 0.0001; s_hi2 = 0.002;
s_try2 = (s_lo2 + s_hi2) / 2;
s_history = [];
diff2 = 100;

while abs(diff2) > 1e-3
    s_history(end+1) = s_try2; %#ok<AGROW>
    nx2   = 25;
    dx_r2 = s_try2 / (nx2-1);
    dx_a2 = t_alloy / (nx2-1);
    dt2   = 0.8 * min(dx_r2^2/(2*alpha_rsz), dx_a2^2/(2*alpha_alloy));
    nt2   = ceil(t_burn / dt2);
    T2    = ones(1, 2*nx2) * 293.15;
    for n = 1:nt2
        To2 = T2;
        T2(1) = To2(1) + 2*alpha_rsz*dt2/dx_r2^2 * ...
            (To2(2)-To2(1) + (dx_r2*h_g/k_rsz)*(T_c-To2(1)));
        for i = 2:nx2-1
            T2(i) = To2(i) + alpha_rsz*dt2/dx_r2^2*(To2(i+1)-2*To2(i)+To2(i-1));
        end
        T2(nx2) = (k_rsz/dx_r2 * To2(nx2-1) + k_alloy/dx_a2 * To2(nx2+1)) / ...
                  (k_rsz/dx_r2 + k_alloy/dx_a2);
        for i = nx2+1:2*nx2-1
            T2(i) = To2(i) + alpha_alloy*dt2/dx_a2^2*(To2(i+1)-2*To2(i)+To2(i-1));
        end
        T2(end) = To2(end) + 2*alpha_alloy*dt2/dx_a2^2*(To2(end-1)-To2(end));
    end
    diff2 = T2(nx2) - T_target_interfaccia;
    if diff2 > 0, s_lo2 = s_try2; else, s_hi2 = s_try2; end
    s_try2 = (s_lo2 + s_hi2) / 2;
    if (s_hi2 - s_lo2) < 1e-8, break; end
end

fig3 = figure('Color','white','Position',[100 100 700 380]);
hold on;
plot(1:length(s_history), s_history*1000, 'o-', ...
    'Color',[0.12 0.37 0.65], 'LineWidth', 1.8, ...
    'MarkerFaceColor','white', 'MarkerSize', 6);
yline(s_rsz*1000, '--', 'Color',[0.85 0.35 0.18], 'LineWidth', 1.2, ...
    'Label',sprintf('Convergence: %.3f mm', s_rsz*1000), ...
    'LabelOrientation','horizontal','FontSize',9);
xlabel('Iteration','FontSize',11);
ylabel('Thickness RSZ [mm]','FontSize',11);
title('Bisection Method Convergence','FontSize',13,'FontWeight','normal');
grid on; box on;
set(gca,'FontSize',10,'GridAlpha',0.25,'LineWidth',0.8);
xlim([1, length(s_history)]);

%% Thermal shock analysis and stress resistance TBC

Bi = (s_rsz) * h_g/k_rsz
A4 = 5.6;
Kic = 1.5e6;    % 1.5 MPa*sqrt(m)
E = 25e9;       % 25 GPa
alpha = 10.1e-6; % CTE
DeltaT = (A4 * Kic) / (E * alpha * sqrt(pi * s_rsz))
t = 0.115+0.45/(1+2.25 * Bi)

%% Thermal rate analysis
dt_sample = 1; 
t_punti = 0:dt_sample:t_burn; 

fprintf('\n==========================================================\n');
fprintf('   Dynamic analysis\n');
fprintf('   Time saved every %.1f s\n', dt_sample);
fprintf('==========================================================\n');
fprintf('%-10s | %-15s | %-15s\n', 'Time [s]', 'T_surf [K]', 'Temperature increase [K]');
fprintf('----------------------------------------------------------\n');

T_precedente = T_surf_time(1); % Temperature at the starting time

for tp = t_punti(2:end)
    [~, idx] = min(abs(time - tp));
    
    T_attuale = T_surf_time(idx);
    Aumento_termico = T_attuale - T_precedente;
    
    fprintf('%-10.1f | %-15.2f | +%-15.2f\n', time(idx), T_attuale, Aumento_termico);
    
    T_precedente = T_attuale;
end
fprintf('==========================================================\n');

%% Save figures
results_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end
figs = [fig1, fig2, fig3];
names = {'chamber_wall_profile', 'chamber_wall_transient', 'chamber_wall_bisection'};
for k = 1:numel(figs)
    set(figs(k), 'Color', 'w');
    exportgraphics(figs(k), fullfile(results_dir, [names{k} '.png']), ...
        'Resolution', 200, 'BackgroundColor', 'white');
end
fprintf('Figures saved to: %s\n', results_dir);

