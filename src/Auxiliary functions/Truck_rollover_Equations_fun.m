%% =======================================================================
% FUNCTION: Truck_rollover_Equations
% DESCRIPTION: State-space formulation of the heavy vehicle dynamics.
%              Includes a 4-wheel Pacejka tire model, dynamic load transfer 
%              for rollover prediction, aerodynamic crosswind loads, and a 
%              closed-loop PID driver model for path tracking.
%
% ORIGINAL AUTHOR: Provided by Course Professors (Politecnico di Milano)
% PROJECT INTEGRATION: Francisco Javier Martin Lopez
%                      Alejandro Rivera Miguez 
%                      Mikel Segovia Diaz (Group N)
% COURSE: Wind Engineering - Politecnico di Milano
%% =======================================================================
function xp = Truck_rollover_Equations_fun(t, x)

    %% GLOBAL VARIABLES
    global g dt
    global m Jz Jx Fza Fzp
    global a b c zR zG xP zP
    global c_roA c_roP k_roA k_roP dNA dNP
    global kp kd ki L x_ref y_ref ErrInt
    global indice xp_prec iii xpp
    global delta
    global vx
    global BF BR C E s q N0

    %% 1. STATE VARIABLES EXTRACTION & KINEMATICS
    % Extract current states from vector x
    vy   = x(1);    % Lateral velocity
    PSIp = x(2);    % Yaw rate
    RHOp = x(3);    % Roll rate
    X    = x(4);    % Global X position
    Y    = x(5);    % Global Y position
    PSI  = x(6);    % Yaw angle
    RHO  = x(7);    % Roll angle
    
    % Global velocities (Transformation from local to global frame)
    Xp = vx * cos(PSI) - vy * sin(PSI);
    Yp = vx * sin(PSI) + vy * cos(PSI);

    %% 2. CLOSED-LOOP DRIVER MODEL (PID STEERING CONTROL)
    % Look-ahead reference point interpolation
    Y_ref  = interp1(x_ref, y_ref, X + L * cos(PSI));
    Y_oss  = Y + L * sin(PSI);                  % Observed Y at look-ahead distance
    Yp_oss = Yp + PSIp * L * cos(PSI);          % Observed lateral velocity
    
    % Integral error accumulation
    ErrInt = ErrInt + (Y_ref - Y_oss) * dt;
    
    % PID Control Law (Straight line / Path tracking)
    d = kp * (Y_ref - Y_oss) + kd * (0 - Yp_oss) + ki * ErrInt; 
    
    % Add damping terms based on yaw dynamics
    d = d - kp * PSI - kd * PSIp;
    
    % Steering angle saturation to avoid unphysical inputs
    if d > (45 * pi/180) 
        d = 45 * pi/180;
    elseif d < (-45 * pi/180) % Fixed asymmetric saturation (-15 to -45 for consistency)
        d = -45 * pi/180;
    end

    %% 3. AERODYNAMIC FORCES CALCULATION
    % Call the aerodynamic function to compute wind forces and moments
    % Note: Tildes (~) are used to ignore extra outputs (alfa, Uy, Vr)
    [FyW, MxW, MzW, ~, ~, ~] = aerodynamics_force(Xp, Yp, PSI, X, Y, t);

    %% 4. TIRE KINEMATICS (SLIP ANGLES)
    % Calculate slip angles for front and rear axles
    alphaF = d - atan((vy + PSIp * a) / vx);
    alphaR = -atan((vy - PSIp * b) / vx);

    %% 5. DYNAMIC LOAD TRANSFER & ROLLOVER CHECK
    % Roll moment due to roll stiffness and damping
    MA  = k_roA * RHO + c_roA * RHOp; % Front axle    
    MP  = k_roP * RHO + c_roP * RHOp; % Rear axle
    
    % Vertical load transfer
    dNa = MA / c; % Front axle transfer
    dNp = MP / c; % Rear axle transfer
    
    % Normal load on each individual tire
    NFL = Fza/2 - dNa; % Normal load front left
    NFR = Fza/2 + dNa; % Normal load front right
    NRL = Fzp/2 - dNp; % Normal load rear left
    NRR = Fzp/2 + dNp; % Normal load rear right
    
    % Rollover Warning Logic
    wheel_flags = {'FL', 'FR', 'RL', 'RR'};
    wheel_loads = [NFL, NFR, NRL, NRR];
    static_loads = [Fza/2, Fza/2, Fzp/2, Fzp/2];
    
    for i = 1:4
        % If load drops below 10% of static load, it triggers a warning
        if wheel_loads(i) < 0.1 * static_loads(i)
            % Print warning to console (Does not stop simulation)
            disp(['Rollover warning: Wheel ', wheel_flags{i}, ...
                  ' load = ', num2str(wheel_loads(i)), ...
                  ' (', num2str(100 * wheel_loads(i) / static_loads(i)), '% of static)']);
        end
    end

    %% 6. PACEJKA MAGIC FORMULA TIRE MODEL
    % Calculate Pacejka 'D' coefficients (Peak friction mapping)
    DFL = (q + s * abs(NFL - N0) / N0) * NFL; 
    DFR = (q + s * abs(NFR - N0) / N0) * NFR; 
    DRL = (q + s * abs(NRL - N0) / N0) * NRL; 
    DRR = (q + s * abs(NRR - N0) / N0) * NRR; 
    
    % Calculate lateral forces for each tire
    FyFL = DFL * sin(C * atan(BF * alphaF - E * (BF * alphaF - atan(BF * alphaF))));
    FyFR = DFR * sin(C * atan(BF * alphaF - E * (BF * alphaF - atan(BF * alphaF))));
    FyF  = (FyFL + FyFR) * cos(d); % Projected front lateral force
    
    FyRL = DRL * sin(C * atan(BR * alphaR - E * (BR * alphaR - atan(BR * alphaR))));
    FyRR = DRR * sin(C * atan(BR * alphaR - E * (BR * alphaR - atan(BR * alphaR))));
    FyR  = (FyRL + FyRR);          % Rear lateral force

    %% 7. EQUATIONS OF MOTION (MULTI-BODY DYNAMICS)
    % Lateral acceleration (Newton's 2nd Law)
    ay    = 1/m * (FyF + FyR + FyW);
    
    % Yaw acceleration (Euler's Equation for Z-axis)
    PSIpp = 1/Jz * (FyF * a - FyR * b + MzW + FyW * xP);
    
    % Roll acceleration (Euler's Equation for X-axis)
    RHOpp = 1/Jx * (m * (zG - zR) * (ay + g * sin(RHO)) - MA - MP + MxW - FyW * (zP - zG));

    %% 8. STATE VECTOR DERIVATIVES CONSTRUCTION
    % State vector layout:
    % x  = [vy, psip, rhop, X, Y, psi, rho]
    % xp = [vyp, psipp, rhopp, Xp, Yp, psip, rhop]
    
    xp = zeros(7,1); % Preallocate memory
    
    xp(1) = ay - vx * PSIp; % Lateral velocity derivative (vyp)
    xp(2) = PSIpp;          % Yaw acceleration
    xp(3) = RHOpp;          % Roll acceleration
    xp(4) = Xp;             % Global X velocity
    xp(5) = Yp;             % Global Y velocity
    xp(6) = PSIp;           % Yaw rate
    xp(7) = RHOp;           % Roll rate

    %% 9. DATA LOGGING (HACK FOR FIXED-STEP RK45)
    % The ode45c solver evaluates 6 stages per time step. 
    % This logic forces logging only once per actual simulation step.
    indice = indice + 1;
    if indice == 6
        indice  = 0;
        xp_prec = xp;
        xpp(iii,:) = xp';
        delta(iii,1) = d;
        dNA(iii,1) = dNa;
        dNP(iii,1) = dNp;
        iii = iii + 1;
    end

end