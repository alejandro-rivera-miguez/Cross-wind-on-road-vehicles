%% =======================================================================
% SCRIPT: Plot_Crosswind_Parametric_Results
% DESCRIPTION: Visualizes the results of the parametric study on heavy 
%              vehicle stability under crosswinds. Generates 3D wind 
%              profiles, and rollover speed boundaries depending on 
%              payload, road conditions, and roll stiffness.
%
% AUTHORS: Francisco Javier Martin Lopez
%          Alejandro Rivera Miguez
%          Mikel Segovia Diaz
% COURSE: Wind Engineering - Politecnico di Milano
%% =======================================================================

clear all; close all; clc;

%% 1. 3D WIND FIELD SURFACE PLOT
fprintf('--- Generating Wind Field Surface Plot ---\n');
% Load wind data matrix
data = load('Vento/ventoV30_Iu14Lu150.dat');

% Extract time and velocities
t = data(:, 1);             % Time column vector (2001x1)
x = 0:5:995;                % Position row vector (1x200)
v = data(:, 2:end);         % Wind velocities matrix (2001x200)

% Create meshgrid for surface plotting
[T, X] = meshgrid(t, x);    % T: (200x2001), X: (200x2001)

figure('Name', 'Wind Field Surface', 'Color', 'w');
surf(T', X', v, 'EdgeColor', 'none'); 
xlabel('Time $t$ $[s]$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Distance $x$ $[m]$', 'Interpreter', 'latex', 'FontSize', 14);
zlabel('Wind velocity $U$ $[m/s]$', 'Interpreter', 'latex', 'FontSize', 14);
title('Turbulent Wind Field Profile', 'Interpreter', 'latex', 'FontSize', 16);

colormap(jet);
c = colorbar;
c.Label.String = 'Wind speed $[m/s]$';
c.Label.Interpreter = 'latex';
c.Label.FontSize = 14;
view(2); % Top-down view (2D projection)
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');

%% 2. ROLLOVER SPEED VS PAYLOAD FACTOR
fprintf('--- Generating Payload Factor Plot ---\n');
% Data
lambda = 0:0.1:1;
V = [0 0 0 0 32 61 92 109 109 110 108];

figure('Name', 'Payload Factor Influence', 'Color', 'w');
plot(lambda, V, '-o', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'b');
grid on; grid minor; box on;
xlabel('Payload factor $\lambda$ $[-]$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Rollover speed $V_{Rollover}$ $[km/h]$', 'Interpreter', 'latex', 'FontSize', 14);
title('Influence of Payload Factor on Rollover Speed', 'Interpreter', 'latex', 'FontSize', 16);
xlim([0 1]);
ylim([0 120]);
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');

%% 3. CHARACTERISTIC WIND CURVES (ROAD CONDITIONS & TURBULENCE)
fprintf('--- Generating Characteristic Wind Curves ---\n');
W_speed = [25, 30]; % Wind speeds [m/s]

% Rollover speeds [km/h] for each condition and turbulence index
V_dry  = [234, 172, 114, 132, 108, 0];      % Dry road
V_wet  = [253, 178, 164, 170, 119, 62];     % Wet road
V_snow = [392, 291, 280, 294, 129, 89];     % Snow road

% Group values based on turbulence index positions:
% Iu = 7%  -> indices 1 & 4
% Iu = 14% -> indices 2 & 5
% Iu = 25% -> indices 3 & 6
V_dry_7   = [V_dry(1), V_dry(4)];
V_dry_14  = [V_dry(2), V_dry(5)];
V_dry_25  = [V_dry(3), V_dry(6)];

V_wet_7   = [V_wet(1), V_wet(4)];
V_wet_14  = [V_wet(2), V_wet(5)];
V_wet_25  = [V_wet(3), V_wet(6)];

V_snow_7  = [V_snow(1), V_snow(4)];
V_snow_14 = [V_snow(2), V_snow(5)];
V_snow_25 = [V_snow(3), V_snow(6)];

% --- Dry Conditions ---
figure('Name', 'Wind Curve - Dry', 'Color', 'w');
plot(V_dry_7, W_speed, '-or', 'LineWidth', 1.5, 'MarkerFaceColor', 'r'); hold on;
plot(V_dry_14, W_speed, '-og', 'LineWidth', 1.5, 'MarkerFaceColor', 'g');
plot(V_dry_25, W_speed, '-ob', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
legend('$I_u = 7\%$', '$I_u = 14\%$', '$I_u = 25\%$', 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'southwest');
xlabel('Rollover speed $[km/h]$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Wind speed $[m/s]$', 'Interpreter', 'latex', 'FontSize', 14);
title('Characteristic Wind Curve (Dry Conditions)', 'Interpreter', 'latex', 'FontSize', 16);
ylim([24 31]); grid on; grid minor; set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');

% --- Wet Conditions ---
figure('Name', 'Wind Curve - Wet', 'Color', 'w');
plot(V_wet_7, W_speed, '-or', 'LineWidth', 1.5, 'MarkerFaceColor', 'r'); hold on;
plot(V_wet_14, W_speed, '-og', 'LineWidth', 1.5, 'MarkerFaceColor', 'g');
plot(V_wet_25, W_speed, '-ob', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
legend('$I_u = 7\%$', '$I_u = 14\%$', '$I_u = 25\%$', 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'southwest');
xlabel('Rollover speed $[km/h]$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Wind speed $[m/s]$', 'Interpreter', 'latex', 'FontSize', 14);
title('Characteristic Wind Curve (Wet Conditions)', 'Interpreter', 'latex', 'FontSize', 16);
ylim([24 31]); grid on; grid minor; set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');

% --- Snow Conditions ---
figure('Name', 'Wind Curve - Snow', 'Color', 'w');
plot(V_snow_7, W_speed, '-or', 'LineWidth', 1.5, 'MarkerFaceColor', 'r'); hold on;
plot(V_snow_14, W_speed, '-og', 'LineWidth', 1.5, 'MarkerFaceColor', 'g');
plot(V_snow_25, W_speed, '-ob', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
legend('$I_u = 7\%$', '$I_u = 14\%$', '$I_u = 25\%$', 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'southwest');
xlabel('Rollover speed $[km/h]$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Wind speed $[m/s]$', 'Interpreter', 'latex', 'FontSize', 14);
title('Characteristic Wind Curve (Snow Conditions)', 'Interpreter', 'latex', 'FontSize', 16);
ylim([24 31]); grid on; grid minor; set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');

%% 4. ROLLOVER SPEED VS ROLL STIFFNESS DISTRIBUTION
fprintf('--- Generating Roll Stiffness Plot ---\n');
% Data
tau_roll = 0.25:0.05:0.75;
V_tau = [0 0 0 0 49 71 98 118 107 93 71];

figure('Name', 'Roll Stiffness Influence', 'Color', 'w');
plot(tau_roll, V_tau, '-o', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'b');
grid on; grid minor; box on;
xlabel('Ratio of roll stiffness $\tau_{Roll}$ $[-]$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Rollover speed $V_{Rollover}$ $[km/h]$', 'Interpreter', 'latex', 'FontSize', 14);
title('Influence of Roll Stiffness Ratio on Rollover Speed', 'Interpreter', 'latex', 'FontSize', 16);
xlim([0.25 0.75]);
ylim([0 120]);
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');

fprintf('--- All plots generated successfully! ---\n');