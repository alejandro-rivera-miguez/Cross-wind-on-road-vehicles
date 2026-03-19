%% =======================================================================
% FUNCTION: Truck_rollover_vehicle_data_evo_fun
% DESCRIPTION: Initializes the heavy vehicle parameters (mass, inertia,
%              geometry) and tire model (Pacejka Magic Formula) 
%              coefficients based on the selected payload, road 
%              condition, and roll stiffness distribution.
%
% AUTHORS: Francisco Javier Martin Lopez
%          Alejandro Rivera Miguez
%          Mikel Segovia Diaz
% COURSE: Wind Engineering - Politecnico di Milano
%% =======================================================================
function Truck_rollover_vehicle_data_evo_fun(lambda, road_select, tau_roll)

    %% GLOBAL VARIABLES DECLARATION
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
    global rollover_flag

    %% 1. VEHICLE DATA: GEOMETRY, MASS, STIFFNESS AND DAMPING
    g       = 9.80665;      % [m/s^2]   Gravity
    p       = 4.42;         % [m]       Wheel-base
    c       = 1.74;         % [m]       Track width
    Rr      = 0.4;          % [m]       Wheel rolling radius
    
    mb      = 4842;         % [kg]      Vehicle body mass
    Jxb     = 1016;         % [kg*m^2]  Vehicle body moment of inertia (x-axis)
    Jzb     = 4795;         % [kg*m^2]  Vehicle body moment of inertia (z-axis)
    ab      = 0.534;        % [m]       Distance of the body CoG from front axle
    zb      = 1.058;        % [m]       Height of the body CoG from ground
    
    al      = 3.35;         % [m]       Distance of the lumped mass CoG from front axle
    ll      = 6.16;         % [m]       Length of cargo bed
    wl      = 2.5;          % [m]       Width of cargo bed
    hl      = 2.62;         % [m]       Height of cargo bed
    
    zl_min  = 0.88;         % [m]       Minimum height of the lumped mass CoG
    zl_max  = zl_min + hl/2;% [m]       Maximum height of the lumped mass CoG
    zl      = zl_min + lambda*(zl_max - zl_min); % [m] Lumped mass CoG height
    
    ml_max  = 8600 - mb;    % [kg]      Maximum lumped mass
    ml      = ml_max * lambda; % [kg]   Lumped mass simulating the load (Payload)
    
    Jxl     = ml/12 * (wl^2 + (lambda*hl)^2); % [kg*m^2] Lumped mass inertia (x-axis)
    Jzl     = ml/12 * (ll^2 + wl^2);          % [kg*m^2] Lumped mass inertia (z-axis)
    
    m       = mb + ml;                        % [kg] Total sprung mass
    a       = (al*ml + ab*mb) / m;            % [m]  Distance of sprung mass CoG from front axle
    b       = p - a;                          % [m]  Distance of sprung mass CoG from rear axle
    zG      = (zb*mb + zl*ml) / (mb + ml);    % [m]  Height of sprung mass CoG from ground
    
    Jz      = Jzb + mb*(a - ab)^2 + Jzl + ml*(a - al)^2; % Total moment of inertia (z-axis)
    Jx      = Jxb + mb*(zG - zb)^2 + Jxl + ml*(zG - zl)^2; % Total moment of inertia (x-axis)
    
    Fza     = m * b / p * g;    % [N]       Static load on front axle
    Fzp     = m * a / p * g;    % [N]       Static load on rear axle
    zR      = Rr * 1.25;        % [m]       Roll centre height from ground
    
    c_roA   = 4 * 0.5 * 10.36 * 1e3; % [Nms/rad] Roll damping front axle
    c_roP   = 2 * 0.5 * 6.28 * 1e3;  % [Nms/rad] Roll damping rear axle
    
    k_roTOT = 1.5 * (252 + 139) * 1e3; % [Nm/rad]  Total vehicle roll stiffness
    k_roA   = k_roTOT * tau_roll;      % [Nm/rad]  Roll stiffness front axle
    k_roP   = k_roTOT * (1 - tau_roll);% [Nm/rad]  Roll stiffness rear axle

    %% 2. PACEJKA TIRE COEFFICIENTS AND ROAD CONDITIONS
    C = 1.5;
    E = -0.5;
    
    % By selecting the road surface condition, the tire-road contact forces 
    % parameters get updated accordingly.
    % 1: Snow | 2: Wet | 3: Dry
    if road_select == 1
        q = 0.3;
    elseif road_select == 2
        q = 0.6;
    elseif road_select == 3
        q = 0.9;
    else
        warning('Error in road condition selection. Ensure road_select is 1, 2, or 3.');
    end
    
    s = 0; % -0.01;
    N0 = 20000; % [N] Nominal vertical load
    
    BF = 2.7502; % Extracted from lateral stiffness calculation models
    BR = 2.5 * BF;
    
end