close all; clear; clc; 

%% Global variables
%% Physical parameters
N = 50; % Number of vertices 

rho_alum = 2700; % Density of Aluminum
BeamLength = 1; % Beam length (meter)
deltaL = BeamLength / (N - 1); % Discrete length

totalTime = 1; % seconds (Total Time)
dt = 0.01; % % Time step size 

R_outer = 0.013; % Beam outer radius (meter)
R_inner = 0.011; % Beam inner radius (meter)

g = 9.8; % Gravity (m/s^2)
Y = 70e9; % Young's modulus 
% Using Y instead of E to avoid ambiguity
% Utility quantities
ne = N - 1; % Number of edges
EI = Y * pi * (R_outer^4 - R_inner^4) / 4;
EA = Y * pi * (R_outer^2 - R_inner^2);

% Mass matrix
mass = pi*(R_outer^2 - R_inner^2)*deltaL*rho_alum;
M = zeros(2 * N, 2 * N); 
for i = 1:N % Assign value to the diagonal matrix
    mass = pi*(R_outer^2 - R_inner^2)*deltaL*rho_alum; 
    M(2 * i - 1, 2 * i - 1) = mass; 
    M(2 * i, 2 * i) = mass; 
end 

% Gravity (only on y direction)
W = zeros(2 * N, 1);
for i = 1:N 
    W(2 * i) = mass * g;  
end 

% Initial DOF vector
q0 = zeros(2 * N, 1);
for i = 1:N
    q0 (2 * i - 1) = deltaL * (i - 1);% x coordinate
    q0 (2 * i) = 0; % y coordinate
end

% New position and velocity
q = q0; % DOF vector
u = (q - q0) / dt; % Velocity vector

% Number of time steps
Nsteps = round( totalTime / dt );
all_y = zeros(Nsteps, N); % y-position (mid sphere)
all_v = zeros(Nsteps, N); % y-velocity (mid sphere)
all_y(1,:) = q(2:2:end); 
all_v(1,:) = u(2:2:end); 

f_applyPoint = round(0.75 / (BeamLength / (N - 1)));


% Tolerance
tol = EI / BeamLength^2 * 1e-3;

free_index = [3:(N*2 - 1)]; 

force_vs_ymax = zeros(10, 3);

for j = 1:10
    force = 2000 + 1000 * (j - 1); 
    fprintf('Force = %f\n', force);
    % Time marching scheme
    for c = 2:Nsteps 
        %fprintf('Steps = %f\n', (c - 1) ); 

        q = q0; % Guess
        q_free = q(free_index); 
        % Newton Raphson
        err = 10 * tol;
        while err > tol
            % Inertia
            f = M / dt * ( (q-q0) / dt - u );
            J = M / dt^2;
            
            f(f_applyPoint * 2) = f(f_applyPoint * 2) + force; 

            for k = 2:N 
                xk = q(2 * (k - 1) - 1); 
                yk = q(2 * (k - 1)); 
                xkp1 = q(2 * k - 1); 
                ykp1 = q(2 * k); 
                dF = gradEs(xk, yk, xkp1, ykp1, deltaL, EA);
                dJ = hessEs(xk, yk, xkp1, ykp1, deltaL, EA);

                bg = 2 * (k - 1) - 1; % begin
                fl = 2 * k; % end 
                f(bg:fl) = f(bg:fl) + dF;
                J(bg:fl,bg:fl) = J(bg:fl,bg:fl) + dJ;
            end


           
            for k = 3:N
                xkm1 = q((k - 2) * 2 - 1); 
                ykm1 = q((k - 2) * 2); 
                xk = q((k - 1) * 2 - 1); 
                yk = q((k - 1) * 2); 
                xkp1 = q(2 * k - 1); 
                ykp1 = q(2 * k); 
                curvature0 = 0;
                dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, ...
                    curvature0, deltaL, EI);
                dJ = hessEb(xkm1, ykm1, xk, yk, xkp1, ykp1, ...
                    curvature0, deltaL, EI);

                bg = 2 * (k - 2) - 1; % start 
                fl = 2 * k; % end
                f(bg:fl) = f(bg:fl) + dF;
                J(bg:fl,bg:fl) = J(bg:fl,bg:fl) + dJ;
            end
            
            f = f - W; 
            % Update (make sure both end-points = 0)
            f_free = f(free_index); 
            J_free = J(free_index, free_index); 
            q_free = q_free - J_free \ f_free; 

            q(free_index) = q_free;
            
            err = sum( abs(f(free_index)) );
        end

        % Update
        u = (q - q0) / dt; % Velocity
        q0 = q; % Old position
        
        
        figure(1);
        plot(q(1:2:end), q(2:2:end), 'ro-');
        axis equal
        drawnow
        
        
        % Store
        all_mid_y(c,:) = q(2:2:end);
        all_mid_v(c,:) = u(2:2:end);
    end
    tempt = -1 * all_mid_y(Nsteps,:); 
    force_vs_ymax(j, 1) = force;
    force_vs_ymax(j, 2) = max(tempt);
    force_vs_ymax(j, 3) = beams_theory(force, BeamLength, min(BeamLength, BeamLength - 0.75), EI, BeamLength);
end


figure(2);
plot(force_vs_ymax(:,1) / 1000, force_vs_ymax(:,2), 'r-');
xlabel('Force, [1000 N]');
ylabel('Y Max, [Meters]');
hold on; 
plot(force_vs_ymax(:,1) / 1000, force_vs_ymax(:,3), 'b-');



function ymax = beams_theory(P, L, c, EI, l)
    ymax = (P*c*(L^2 - c^2)^1.5) / (9*sqrt(3)*EI*l);
end


















