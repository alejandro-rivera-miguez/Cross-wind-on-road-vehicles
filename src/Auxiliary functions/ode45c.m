%% =======================================================================
% FUNCTION: ode45c
% DESCRIPTION: Custom ordinary differential equation solver using 4th/5th 
%              order Runge-Kutta-Fehlberg formulas. 
%              NOTE: The adaptive step-size mechanism has been intentionally 
%              disabled in this version to enforce a fixed-step integration (h).
%              Fixed-step solvers are strictly required for vehicle dynamic 
%              simulations involving discrete driver/PID controllers.
%
% ORIGINAL AUTHOR: C.B. Moler (The MathWorks, Inc., 1984-1994)
% PROJECT INTEGRATION: Francisco Javier Martin Lopez
%                      Alejandro Rivera Miguez
%                      Mikel Segovia Diaz (Group N)
% COURSE: Wind Engineering - Politecnico di Milano
%% =======================================================================
function [tout, yout] = ode45c(ypfun, t0, tfinal, h, y0, tol, trace)
%	INPUT:
%	ypfun - String/Handle containing name of user-supplied problem description.
%	        Call: yprime = fun(t,y) where ypfun = 'fun'.
%	t0    - Initial value of time t.
%	tfinal- Final value of time t.
%   h     - Fixed integration time step.
%	y0    - Initial value column-vector.
%	tol   - The desired accuracy (Default: tol = 1.e-6). Not used in fixed-step.
%	trace - If nonzero, each step is printed. (Default: trace = 0).
%
%	OUTPUT:
%	tout  - Returned integration time points (column-vector).
%	yout  - Returned solution, one solution column-vector per tout-value.

    % The Fehlberg coefficients (Runge-Kutta Butcher Tableau):
    alpha = [1/4  3/8  12/13  1  1/2]';
    beta  = [ [    1      0      0     0      0    0]/4
              [    3      9      0     0      0    0]/32
              [ 1932  -7200   7296     0      0    0]/2197
              [ 8341 -32832  29440  -845      0    0]/4104
              [-6080  41040 -28352  9295  -5643    0]/20520 ]';
    gamma = [ [902880  0  3953664  3855735  -1371249  277020]/7618050
              [ -2090  0    22528    21970    -15048  -27360]/752400 ]';
    
    pow = 1/5;
    
    % Default argument handling
    if nargin < 6, tol = 1.e-6; end
    if nargin < 7, trace = 0; end
    
    % --- Initialization ---
    t = t0;
    hmax = (tfinal - t)/16; % Max step limit (legacy from adaptive version)
    
    y = y0(:);
    f = zeros(length(y), 6);
    
    % Preallocate output arrays in chunks to drastically improve speed
    chunk = 128;
    tout = zeros(chunk, 1);
    yout = zeros(chunk, length(y));
    k = 1;
    tout(k) = t;
    yout(k,:) = y.';
    
    if trace
       clc, t, h, y
    end
    
    % Initialize GUI progress bar
    att = waitbar(0, 'Wait, numeric integration in progress...', 'Color', [0, 1, 0.8]);
    
    %% --- Main Integration Loop ---
    while (t < tfinal) && (t + h > t)
       
        % Ensure we don't overshoot the final time
        if t + h > tfinal
            h = tfinal - t; 
        end
        
        % 1. Compute the RK slopes (k1 to k6)
        temp = feval(ypfun, t, y);
        f(:,1) = temp(:);
        for j = 1:5
            temp = feval(ypfun, t + alpha(j)*h, y + h*f*beta(:,j));
            f(:,j+1) = temp(:);
        end
       
        % 2. Estimate the error and the acceptable error 
        % (Legacy Fehlberg logic - computed but ignored for fixed-step)
        delta = norm(h * f * gamma(:,2), 'inf');
        tau = tol * max(norm(y, 'inf'), 1.0);
       
        % 3. Update the solution (Fixed-step logic enforced)
        % Original adaptive condition 'if delta <= tau' is bypassed.
        t = t + h;
        y = y + h * f * gamma(:,1);
        k = k + 1;
        
        % Expand chunk memory if limits are exceeded
        if k > length(tout)
            tout = [tout; zeros(chunk, 1)];
            yout = [yout; zeros(chunk, length(y))];
        end
        
        % Store results
        tout(k) = t;
        yout(k,:) = y.';
       
        if trace
            home, t, h, y
        end
       
        % 4. Update the step size (Adaptive logic disabled)
        % if delta ~= 0.0
        %    h = min(hmax, 0.8 * h * (tau/delta)^pow);
        % end
        
        % Update progress bar
        waitbar(t/tfinal, att);
    end
    
    % --- Cleanup and Output Formatting ---
    close(att); % Close progress bar
    
    if (t < tfinal)
       disp('Singularity likely.')
       t
    end
    
    % Truncate preallocated arrays to actual size
    tout = tout(1:k);
    yout = yout(1:k, :);

end