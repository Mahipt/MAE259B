close all; clear; clc; 


%% Global Variables 
% Physical Parameters
g = 9.8; % Gravity 
visc = 1000; % Viscosity (Pa-S) 

N = 3; % Number of vertice 

rho_metal = 7000; % Density of metal spehere 
rho_fluid = 1000; % Density of the fluid  
rho = rho_metal - rho_fluid; % Relative Density of the sphere

R1 = 0.025; % Radius of sphere 1 
R2 = 0.025; % Radius of sphere 2
R3 = 0.025; % Radius of sphere 3 
RodLength = 0.1;  % Total rod length 
deltaL = RodLength / (N - 1); % Rod length in each section 
r0 = 0.001; % Rod radius

% pre-cal parameters
% Mass matrix 
M = zeros(N * 2, N * 2); % 2 DOF for x and y
M(1,1) = 4/3*pi*R1^3*rho_metal;
M(2,2) = 4/3*pi*R1^3*rho_metal;
M(3,3) = 4/3*pi*R2^3*rho_metal;
M(4,4) = 4/3*pi*R2^3*rho_metal;
M(5,5) = 4/3*pi*R3^3*rho_metal;
M(6,6) = 4/3*pi*R3^3*rho_metal;

% Viscous damping matrix 
C = zeros(N * 2, N * 2); % viscosity
C1 = 6*pi*visc*R1;
C2 = 6*pi*visc*R2;
C3 = 6*pi*visc*R3;
C(1,1) = C1;
C(2,2) = C1;
C(3,3) = C2;
C(4,4) = C2;
C(5,5) = C3;
C(6,6) = C3;

% Gravity (only for y direction) 
W = zeros(N * 2, 1); 
W(2) = -4/3*pi*R1^3*rho*g;
W(4) = -4/3*pi*R2^3*rho*g;
W(6) = -4/3*pi*R3^3*rho*g;

% Time step size 
totalTime = 10; 
dt = 0.001; 

% Utility quantities 
Y = 1e9; % Young's modulus
ne = N - 1; % Number of edges 
EI = Y * pi * r0^4 / 4; 
EA = Y * pi * r0^2; 

% Initial Condition and update variable
q0 = zeros(2 * N, 1);% Initial DOF vector 
for i = 1:N % Giving q0 the initial value
	q0(2 * i - 1) = (i - 1) * deltaL;% init x position 
	q0(2 * i) = 0;% initial y position 
end 
q = q0; % DOF vector (initial with init-cond)
u = (q - q0) / dt; % updated velocity (initial with init-cond) 
Nsteps = round(totalTime / dt); % Calculate the total steps

N_mid = N - 2; % Number of the midpoint
all_mid_y = zeros(Nsteps, N_mid); 
all_mid_v = zeros(Nsteps, N_mid); 
for i = i:N_mid
	all_mid_y = q((i + 1) * 2); 
	all_mid_v = u((i + 1) * 2); 
end 

% Calculation start from here 
%{ 
	Calculate effect of each segment in x and y direction, 
	and update them to the node that it affect directly. 
%}
% Hyper-parameter for calculation 
tol = EI / RodLength^2 * 1e-3; 
err = 10 * tol; 
for i = 2:Nsteps
	% initialize hyper-parameters for each cal in each steps 
	q = q0; % initial position (from guessing) 
	err = 10 * tol; % Start with large error 
	while err > tol % (Newton Raphson) 
		% Intertia
		f = M / dt * ((q - q0) / dt - u); % F = m * a
		J = M / (dt^2); 
		%
        % Elastic forces
        %
        % Linear spring 1 between nodes 1 and 2
        xk = q(1);
        yk = q(2);
        xkp1 = q(3);
        ykp1 = q(4);
        dF = gradEs(xk, yk, xkp1, ykp1, deltaL, EA);
        dJ = hessEs(xk, yk, xkp1, ykp1, deltaL, EA);
        f(1:4) = f(1:4) + dF;
        J(1:4,1:4) = J(1:4,1:4) + dJ;
        
        % Linear spring 2 between nodes 2 and 3
        xk = q(3);
        yk = q(4);
        xkp1 = q(5);
        ykp1 = q(6);
        dF = gradEs(xk, yk, xkp1, ykp1, deltaL, EA);
        dJ = hessEs(xk, yk, xkp1, ykp1, deltaL, EA);
        f(3:6) = f(3:6) + dF;
        J(3:6,3:6) = J(3:6,3:6) + dJ;
        
        % Bending spring between nodes 1, 2, and 3
        xkm1 = q(1);
        ykm1 = q(2);
        xk = q(3);
        yk = q(4);
        xkp1 = q(5);
        ykp1 = q(6);
        curvature0 = 0;
        dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, ...
            curvature0, deltaL, EI);
        dJ = hessEb(xkm1, ykm1, xk, yk, xkp1, ykp1, ...
            curvature0, deltaL, EI);
        f(1:6) = f(1:6) + dF;
        J(1:6,1:6) = J(1:6,1:6) + dJ;
        
        % Viscous force
        f = f + C * ( q - q0 ) / dt;
        J = J + C / dt;
        
        % Weight
        f = f - W;
        
        % Update
        q = q - J \ f;
        
        err = sum( abs(f) );
	end 
	% update 
	u = (q - q0) / dt; 
	q0 = q; 

	% Plot position 
	figure(1); 
	plot(q(1:2:end), q(2:2:end), 'ro-'); 
	axis equal; 
	drawnow; 

	% Store
	for j = 1:N_mid
	    all_mid_y(i, j) = q((j + 1) * 2);
	    all_mid_v(i, j) = u((j + 1) * 2);
	end
end 

figure(2);
timeArray = (1:Nsteps) * dt;
plot(timeArray, all_mid_v, 'k-');
xlabel('Time, t [sec]');
ylabel('Velocity of mid-node, v [meter/sec]');


