%% Greenhouse Simulation Parameters

clear
clc

% ---------------- Simulation Settings ----------------
T_sample= 1;                % Time step [s]

% ---------------- Initial Conditions ----------------
T0 = 20;                   % Initial temperature [°C]
L0 = 1;                  % Initial water level [L]
N0 = 1;                  % Initial nutrient concentration [mg/L]

% ---------------- Temperature Dynamics ----------------
C_T = 12000;               % Thermal capacity [J/°C]
k_T = 20;                  % Heat transfer coefficient [W/°C]
T_opt = 22;                % Optimal temperature for nutrient uptake [°C]

% ---------------- Water Level Dynamics ----------------
q_p = 0.002;               % Plant water consumption rate [L/s]
epsilon = 5e-5;            % Evaporation coefficient [L/(s·°C)]
L_opt = 10;                % Optimal water level for nutrient uptake[L]

% ---------------- Nutrient Dynamics ----------------
R = 10;                    % Nutrient injection rate [mg/s]
k_N = 2e-4;                % Base nutrient uptake rate [1/s]
alpha = 0.05;              % Temperature sensitivity of uptake [1/°C]
N_opt = 150;

% ---------------- External Profiles ----------------
tempo=[0:60*60:24*60*60]';

T_ext = [18; 18; 20; 22; 22; 22; 24; 24; 25; 25; 26; 28; 30; 30; 30; 32; 34; 35; 36; 33; 30; 26; 24; 22; 17];
T_ext2sim=[tempo T_ext];

%% Graphs

% Plot T(t)
figure();
plot(T2sim.Time, T2sim.Data);
xlabel('Time (s)');
ylabel('Temperature (°C)');
title('Temperature T(t)');
grid on;

% Plot N(t)
figure();
plot(N2sim.Time, N2sim.Data);
xlabel('Time (s)');
ylabel('Nutrient concentration (mg/L)');
title('Nutrient concentration N(t)');
grid on;

% Plot L(t)
figure();
plot(L2sim.Time, L2sim.Data);
xlabel('Time (s)');
ylabel('Water level (L)');
title('Water level L(t)');
grid on;

%% Calculations

% Water used
water_used = trapz(uL2sim.Time, uL2sim.Data);
fprintf('Water Used (in Litre): %.4f\n', water_used);

% Heating/Coolig Energy
area_under_Q_curve = trapz(Q2sim.Time, abs(Q2sim.Data));
%fprintf('Heating/Coolig Energy (in Joules): %.4f\n', area_under_Q_curve);
fprintf('Heating/Coolig Energy (in kJ): %.4f\n', area_under_Q_curve/1000);

% Nutrient Added
nutrient_added = trapz(Ninject2sim.Time, Ninject2sim.Data);
fprintf('Nutrient added (in mg): %.4f\n', nutrient_added);

% Nutrient absorption
absorption_rate = NAbsorbRate2sim.Data .* L2sim.Data;
plant_N_absorption = trapz(NAbsorbRate2sim.Time, absorption_rate);
fprintf('Plant Nutrient Absorption (in mg): %.4f\n', plant_N_absorption);

% crop yield
nu_N = 0.8;
crop_yield = nu_N * plant_N_absorption;
%fprintf('Crop yield (in mg): %.4f\n', crop_yield);
fprintf('Crop yield (in kg): %.4f\n', crop_yield/1000000);

% water efficiency
fprintf('Water Efficiency (kg/L): %.7f\n', crop_yield/1000000/water_used);

% Energy efficiency
fprintf('Energy Efficiency (kg/kWh): %.4f\n', crop_yield/1000000/(area_under_Q_curve/3600000));

% water efficiency
fprintf('Nutrient Efficiency (mg/mg): %.4f\n', plant_N_absorption/nutrient_added);