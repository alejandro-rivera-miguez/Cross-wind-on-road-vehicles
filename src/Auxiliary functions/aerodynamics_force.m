%% =======================================================================
% FUNCTION: forze_aerodinamiche
% DESCRIPTION: Computes the instantaneous aerodynamic forces and moments 
%              acting on the heavy vehicle based on its kinematics and 
%              the turbulent wind field profile.
%
% AUTHORS: Francisco Javier Martin Lopez
%          Alejandro Rivera Miguez
%          Mikel Segovia Diaz
% COURSE: Wind Engineering - Politecnico di Milano
%
% INPUTS:
%   Xp    - Vehicle velocity along global X-axis [m/s]
%   Yp    - Vehicle velocity along global Y-axis [m/s]
%   psi   - Vehicle yaw/heading angle w.r.t global X-axis [rad]
%   X     - Vehicle global X-coordinate [m]
%   Y     - Vehicle global Y-coordinate [m]
%   time  - Current simulation time [s]
%
% OUTPUTS:
%   Fy    - Aerodynamic lateral force along vehicle local y-axis [N]
%   Mx    - Aerodynamic roll moment along vehicle local x-axis [Nm]
%   Mz    - Aerodynamic yaw moment along vehicle local z-axis [Nm]
%   alfa  - Aerodynamic angle of attack / wind yaw angle [rad]
%   Uy    - Absolute wind speed seen by the vehicle [m/s]
%   Vr    - Relative wind speed magnitude [m/s]
%% =======================================================================
function [Fy, Mx, Mz, alfa, Uy, Vr] = aerodynamics_force(Xp, Yp, psi, X, Y, time)

    % Global parameters defined in the main script
    global Wind_Coefficients vento Al hr den

    %% 1. WIND SPEED INTERPOLATION
    % Interpolate absolute wind speed from 3D wind data field at current time and position
    Uy = interp2(vento.x, vento.t, vento.v, X, time);
    
    % Smooth ramp-up to avoid numerical instability (shock loads) at simulation start
    if time < 2
        Uy = Uy * time / 2;
    end

    %% 2. RELATIVE VELOCITY KINEMATICS
    % Relative Wind speed in the absolute reference frame
    Vrx = -Xp;       % Vehicle velocity component in X (opposite to relative wind)
    Vry = Uy - Yp;   % Wind speed minus vehicle lateral velocity
    
    % Relative wind speed modulus
    Vr  = sqrt(Vrx^2 + Vry^2); 

    % Relative Wind speed transformed into the local (vehicle) reference frame
    vrx = Vrx * cos(psi) + Vry * sin(psi);
    vry = -Vrx * sin(psi) + Vry * cos(psi);

    %% 3. AERODYNAMIC ANGLES & COEFFICIENTS
    % Wind yaw angle (angle between relative wind and vehicle's longitudinal axis)
    % Using atan2 ensures proper quadrant handling
    alfa = atan2(vry, -vrx);  
    
    % Interpolate aerodynamic coefficients based on the current wind angle of attack
    % Wind_Coefficients matrix format: [Alpha; Cx; Cy; Mz; Mx]
    C = interp1(Wind_Coefficients(1,:), Wind_Coefficients(2:5,:).', alfa * 180/pi, 'linear');
    
    % Extract relevant lateral and moment coefficients
    Cy  = C(2);
    Cmz = C(3);
    Cmx = C(4);

    %% 4. AERODYNAMIC FORCES & MOMENTS CALCULATION
    % Calculate lateral force (Fy), roll moment (Mx), and yaw moment (Mz)
    Fy = 0.5 * den * Al * Cy * Vr^2;
    Mx = 0.5 * den * Al * hr * Cmx * Vr^2;
    Mz = 0.5 * den * Al * hr * Cmz * Vr^2;

end