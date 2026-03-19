%% =======================================================================
% MAIN SCRIPT: Heavy Vehicle Crosswind Dynamic Simulation
% DESCRIPTION: Simulates the dynamic response, trajectory deviation, and 
%              rollover risk of a heavy goods vehicle subjected to 
%              turbulent crosswinds. Uses a 4-contact tire model and a 
%              PID-based driver model for trajectory tracking.
%
% AUTHORS: Francisco Javier Martin Lopez
%          Alejandro Rivera Miguez
%          Mikel Segovia Diaz
% COURSE: Wind Engineering - Politecnico di Milano (A.A. 2024-2025)
%% =======================================================================

close all; clear all; clc;
addpath("Auxiliary functions\")

%% 1. GLOBAL VARIABLES & INITIAL CONFIGURATION
global g dt
global m Jz Jx Fza Fzp
global a b p c zR zG xP zP Rr
global c_roA c_roP k_roA k_roP dNA dNP
global BF BR C E q s N0
global kp kd ki L x_ref y_ref ErrInt
global Af Al den fid hr
global FyW MzW MxW Wind_Coefficients vento
global indice xp_prec iii xpp
global delta d
global vx

fprintf('--- Starting Heavy Vehicle Crosswind Simulation ---\n');

% Load wind data (Ensure the 'Vento' folder and files exist in the path)
vento = struct();
vento.v = load('Vento/ventoV25_Iu07Lu150.dat');
% Alternative wind profiles:
% vento.v = load('Vento/ventoV25_Iu14Lu150.dat');
% vento.v = load('Vento/ventoV30_Iu25Lu150.dat');

% Wind data structuring
vento.t = vento.v(:,1);       % Wind time vector
vento.v = vento.v(:,2:end);   % Wind velocity matrix
vento.x = (0:5:5*(size(vento.v,2)-1))'; % Spatial x-coordinates of the wind

% Simulation Parameters
V = 112;            % Vehicle speed [km/h]
dt = 0.03;          % Integration time step [s]
tfin = 100;         % Final simulation time [s]
lambda = 0.9;       % Adherence factor (Friction)
tau_roll = 0.644;   % Roll stiffness distribution front-to-total [0.25 to 0.75]
road_select = 3;    % Road type: 1 = snow, 2 = wet, 3 = dry

% Variables and arrays initialization
v_roll_vec = [];
indice = 0;
iii = 2;
xpp(1,:) = zeros(1,7); % State derivatives
delta(1,1) = 0;        % Steering angle
dNA(1,1) = 0;          % Front load transfer
dNP(1,1) = 0;          % Rear load transfer

%% 2. VEHICLE DATA & DRIVER CONTROLLER
% Call external function to load dynamic truck data
Truck_rollover_vehicle_data_evo_fun(lambda, road_select, tau_roll);

% Driver Model (PID Controller Gains)
kp = 5;                             % Proportional gain [rad/m]
kd = 0.05;                          % Derivative gain
ki = 0.1;                           % Integral gain
L_lookahead = 3;                    % Look-ahead prediction distance [m] (Mapped as L globally)
L = L_lookahead;

% Reference Path Definition (Straight line simulation)
x_path = [0 20 25 70 75 1000];      % Key points in X [m]
y_path = [0 0  1  1  0  0] * 0;     % Key points in Y [m] (Currently set to straight line Y=0)
x_ref = 0:1:1000;                   % Dense X vector for interpolation
y_ref = interp1(x_path, y_path, x_ref); % Path interpolation
y_ref = smooth(y_ref, 20);          % Path smoothing to avoid sharp steering inputs
ErrInt = 0;                         % Initial integral error

%% 3. AERODYNAMIC PARAMETERS
Af = 6.6;               % Frontal area [m^2]
Al = 18.9;              % Lateral area [m^2]
hr = 2.62;              % Reference height [m]
den = 1.204;            % Air density [kg/m^3]
xP = p/2 - b;           % Aerodynamic center of pressure X-coordinate [m]
zP = 0;                 % Aerodynamic center of pressure Z-coordinate [m]

% Load aerodynamic coefficients matrix
load Wind_Coefficients;  

%% 4. NUMERICAL INTEGRATION (DYNAMIC SIMULATION)
vx = V / 3.6;           % Convert vehicle speed to m/s

% State vector (x) definition:
% 1: vy (lateral vel), 2: psip (yaw rate), 3: rhop (roll rate)
% 4: X (global pos X), 5: Y (global pos Y), 6: psi (yaw angle), 7: rho (roll angle)
xp_prec = zeros(7,1); 
x0 = zeros(7,1);        % Null initial conditions
t0 = 0;

fprintf('Solving differential equations using custom ode45c solver...\n');
% Solve differential equations using custom solver ode45c
[t, x] = ode45c('Truck_rollover_Equations_fun', t0, tfin, dt, x0);

% Extract and calculate state variables after simulation
vx   = vx * ones(size(t)); % Longitudinal velocity is constant
vy   = x(:,1);
psip = x(:,2);
rhop = x(:,3);
X    = x(:,4);
Y    = x(:,5);
psi  = x(:,6);
rho  = x(:,7);

% Extract accelerations (state derivatives calculated inside the solver)
vyp   = xpp(:,1);
psipp = xpp(:,2);
rhopp = xpp(:,3);

ay = vyp + vx .* psip;     % True vehicle lateral acceleration
Xp = xpp(:,4);
Yp = xpp(:,5);
beta = atan2(vy, vx);      % Sideslip angle

%% 5. WIND FORCES & VERTICAL LOADS CALCULATION
% Recalculate aerodynamic forces at each time step
for ii = 1:length(t)
    [FyW(ii,1), MxW(ii,1), MzW(ii,1), alfa(ii,1), Uy(ii,1), Vr(ii,1)] = ...
        aerodynamics_force(Xp(ii), Yp(ii), psi(ii), X(ii), Y(ii), t(ii));
end

% Roll Load Transfer
MA = k_roA * rho + c_roA * rhop;    % Roll moment at front axle
MP = k_roP * rho + c_roP * rhop;    % Roll moment at rear axle
dNa = MA / c;                       % Load variation at front axle
dNp = MP / c;                       % Load variation at rear axle

% Vertical (Normal) load on each tire
NFL = Fza/2 - dNa; % Front Left
NFR = Fza/2 + dNa; % Front Right
NRL = Fzp/2 - dNp; % Rear Left
NRR = Fzp/2 + dNp; % Rear Right

% Rollover Warning System
for i = 1:length(t)
    % If any wheel's load drops below 10% of its static load, trigger warning
    if NFL(i) < 0.1*Fza/2 || NFR(i) < 0.1*Fza/2 || NRL(i) < 0.1*Fzp/2 || NRR(i) < 0.1*Fzp/2
        warning('⚠️ Rollover danger detected at time t = %.2f s', t(i));
        break % Exit loop upon first detection
    end
end

% Tire slip angles
alphaF = delta - atan((vy + psip*a) ./ vx);
alphaR = -atan((vy - psip*b) ./ vx);

% Pacejka Tire Model (Magic Formula) - Coefficient D calculation
DFL = (q + s .* (NFL - N0) ./ N0) .* NFL; 
DFR = (q + s .* (NFR - N0) ./ N0) .* NFR; 
DRL = (q + s .* (NRL - N0) ./ N0) .* NRL; 
DRR = (q + s .* (NRR - N0) ./ N0) .* NRR; 

% Lateral forces on each tire
FyFL = DFL .* sin(C .* atan(BF .* alphaF - E .* (BF .* alphaF - atan(BF .* alphaF))));
FyFR = DFR .* sin(C .* atan(BF .* alphaF - E .* (BF .* alphaF - atan(BF .* alphaF))));
FyF  = (FyFL + FyFR) .* cos(delta); % Total lateral force on front axle

FyRL = DRL .* sin(C .* atan(BR .* alphaR - E .* (BR .* alphaR - atan(BR .* alphaR))));
FyRR = DRR .* sin(C .* atan(BR .* alphaR - E .* (BR .* alphaR - atan(BR .* alphaR))));
FyR  = (FyRL + FyRR);               % Total lateral force on rear axle

fprintf('Simulation complete. Generating figures...\n');

%% 6. RESULTS PLOTTING
% Extract data for Figure 1 from Wind_Coefficients
AOA = Wind_Coefficients(1, :);
Cy  = Wind_Coefficients(3, :);
Cmz = Wind_Coefficients(4, :);
Cmx = Wind_Coefficients(5, :);

% --- FIGURE 1: Basic Aerodynamic Coefficients ---
figure('Name', 'Aerodynamic Coefficients', 'Color', 'w');
subplot(1,3,1);
plot(AOA, Cy, '*', 'LineWidth', 1.5);
xlabel('$\alpha$ [deg]', 'Interpreter', 'latex', 'FontSize', 12); ylabel('$C_y$', 'Interpreter', 'latex', 'FontSize', 12);
title('Lateral Force', 'Interpreter', 'latex', 'FontSize', 14); grid on; set(gca, 'FontSize', 12);
subplot(1,3,2);
plot(AOA, Cmz, '*', 'LineWidth', 1.5);
xlabel('$\alpha$ [deg]', 'Interpreter', 'latex', 'FontSize', 12); ylabel('$C_{mz}$', 'Interpreter', 'latex', 'FontSize', 12);
title('Yaw Moment', 'Interpreter', 'latex', 'FontSize', 14); grid on; set(gca, 'FontSize', 12);
subplot(1,3,3);
plot(AOA, Cmx, '*', 'LineWidth', 1.5);
xlabel('$\alpha$ [deg]', 'Interpreter', 'latex', 'FontSize', 12); ylabel('$C_{mx}$', 'Interpreter', 'latex', 'FontSize', 12);
title('Roll Moment', 'Interpreter', 'latex', 'FontSize', 14); grid on; set(gca, 'FontSize', 12);

% --- FIGURE 2: Vehicle Dynamics ---
figure('Name', 'Vehicle Dynamics Overview', 'Color', 'w');
subplot(231); plot(t, delta*180/pi, 'LineWidth', 1.5); xlabel('Time [s]'); ylabel('[deg]'); grid on; title('Steer Angle');
subplot(232); plot(t, ay, 'LineWidth', 1.5); xlabel('Time [s]'); ylabel('[m/s^2]'); grid on; title('Lateral Acceleration');
subplot(233); plot(t, vx*3.6, t, vy*3.6, 'LineWidth', 1.5); xlabel('Time [s]'); ylabel('[km/h]'); grid on; legend('v_x','v_y'); title('Speed');
subplot(234); plot(t, beta*180/pi, 'LineWidth', 1.5); xlabel('Time [s]'); ylabel('[deg]'); grid on; title('Sideslip Angle');
subplot(235); plot(t, psip, 'LineWidth', 1.5); xlabel('Time [s]'); ylabel('[rad/s]'); grid on; title('Yaw Rate');
subplot(236); plot(t, rhop, 'LineWidth', 1.5); xlabel('Time [s]'); ylabel('[rad/s]'); grid on; title('Roll Rate');

% --- FIGURE 3: Wind Analysis & Trajectory ---
figure('Name', 'Wind Forces and Trajectory', 'Color', 'w');
subplot(231); plot(t, alfa*180/pi, t, Uy, t, Vr, 'LineWidth', 1.5); grid on; legend('Angle','Wind Speed','Rel Speed'); xlabel('Time [s]'); title('Wind Yaw Angle & Speed');
subplot(232); plot(t, FyW, t, MxW, t, MzW, 'LineWidth', 1.5); grid on; xlabel('Time [s]'); ylabel('[N], [Nm]'); legend('Fy','Mx','Mz'); title('Wind Forces & Moments');
subplot(233); plot(t, NFL/Fza*2, t, NFR/Fza*2, t, NRL/Fzp*2, t, NRR/Fzp*2, 'LineWidth', 1.5); grid on; xlabel('Time [s]'); title('Load Ratio N/N_{static}'); legend('FL','FR','RL','RR');
subplot(234); plot(X, Y, 'LineWidth', 1.5); grid on; hold on; plot(x_ref, y_ref, 'r--', 'LineWidth', 1.5); xlabel('X [m]'); ylabel('Y [m]'); legend('CoG','Ref'); xlim([0 max(X)]); title('Trajectory');
subplot(235); plot(t, psi*(180/pi), 'LineWidth', 1.5); grid on; xlabel('Time [s]'); ylabel('[deg]'); title('Yaw Angle');
subplot(236); plot(t, rho*(180/pi), 'LineWidth', 1.5); grid on; xlabel('Time [s]'); ylabel('[deg]'); title('Roll Angle');

% --- FIGURE 4: Trajectory Animation (Quiver) ---
figure('Name', 'Vehicle Trajectory and Orientation', 'Color', 'w');
plot(X, Y, 'LineWidth', 1.5); hold on; title('Vehicle Trajectory and Orientation', 'Interpreter', 'latex', 'FontSize', 14);
for ii = 1:10:length(t) % Downsampled loop for faster quiver plotting
    Xg = interp1(t, X, t(ii));
    Yg = interp1(t, Y, t(ii));
    Psi = interp1(t, psi, t(ii));
    Xpg = 5 * cos(Psi);
    Ypg = 5 * sin(Psi);
    quiver(Xg, Yg, Xpg, Ypg, 'r', 'MaxHeadSize', 0.5); 
end
axis equal; grid on; xlabel('X [m]'); ylabel('Y [m]');

% --- FIGURE 6: Velocity and Turbulence Relationship ---
iu = [7 14 25]';
v25 = [140 130 115]';
v30 = [140 125 80]';
figure('Name', 'Turbulence Effects', 'Color', 'w');
plot(v25, iu, '*-', v30, iu, '*-', 'LineWidth', 1.5, 'MarkerSize', 8); grid on; grid minor;
legend('U = 25 m/s','U = 30 m/s', 'Location', 'best');
xlabel('Vehicle Critical Speed [km/h]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Turbulence Intensity $I_u$ [\%]', 'Interpreter', 'latex', 'FontSize', 12);
title('Critical Speed vs. Turbulence', 'Interpreter', 'latex', 'FontSize', 14);

% --- FIGURE 7: Absolute Load on Wheels ---
figure('Name', 'Absolute Wheel Loads', 'Color', 'w');
plot(t, NFL, t, NFR, t, NRL, t, NRR, 'LineWidth', 1.5);
grid on; grid minor; xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 12); ylabel('Load [N]', 'Interpreter', 'latex', 'FontSize', 12);
title('Absolute Load on Wheels', 'Interpreter', 'latex', 'FontSize', 14);
legend('Front Left (FL)', 'Front Right (FR)', 'Rear Left (RL)', 'Rear Right (RR)', 'Location', 'best');

% --- FIGURE 8: Normalized Load & Rollover Limit ---
figure('Name', 'Rollover Threshold Check', 'Color', 'w');
plot(t, NFL/Fza*2, 'LineWidth', 1.5); hold on; grid on; grid minor;
plot(t, NFR/Fza*2, 'LineWidth', 1.5);
plot(t, NRL/Fzp*2, 'LineWidth', 1.5);
plot(t, NRR/Fzp*2, 'LineWidth', 1.5);
plot(t, 0.1 * ones(size(t)), '--r', 'LineWidth', 2); % 10% Rollover Limit
xlabel('Time $t$ $[s]$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Normalized Vertical Load $N/N_{static}$ $[-]$', 'Interpreter', 'latex', 'FontSize', 14);
title('Wheel Load Ratio and Rollover Limit', 'Interpreter', 'latex', 'FontSize', 16);
xlim([0 50]);
legend('Front Left', 'Front Right', 'Rear Left', 'Rear Right', 'Rollover Limit', 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'best');