close all; clear; clc; 


%% Global Variables 
% Physical Parameters
g = 9.8; % Gravity 
visc = 1000; % Viscosity (Pa-S) 

N = 1; % Number of vertice 
R = 0.025; % Radius of spheres
rho_metal = 7000; % Density of metal spehere 
rho_fluid = 1000; % Density of the fluid  
rho = rho_metal - rho_fluid; % Relative Density of the sphere

% pre-cal parameters 
M = (4/3) * pi() * R^3 * rho_metal; 
C = visc; % viscosity
W = (-4/3) * pi() * R^3 * rho * g; % Gravity on y direction 

% Time step size 
totalTime = 1000; 
dt = 0.01; 

% Initial Condition and update variable
q0 = zeros(2, 1); 
q0(1) = 0.1; %% (X direction) Middle of the pool  
q0(2) = 0; 

q = q0; % updated position of the sphere (inital with init-cond) 
u = (q - q0) / dt; % updated velocity (initial with init-cond) 
steps = round(totalTime / dt); % Calculate the total steps


% every positions 
% x direction, y direction, and time steps
position = zeros(steps, 3);
position(:, 3) = (1:steps) * dt; % writing all time steps

% Calculation start from here 
%{ 
	Calculate effect of each segment in x and y direction, 
	and update them to the node that it affect directly. 
%}
% Hyper-parameter for calculation 
tol = 1e-3; 
err = 10 * tol; 
for i = 1:steps
	% initialize hyper-parameters for each cal in each steps 
	q = q0; % initial position  
	err = 10 * tol; % Start with large error 
	while err > tol
			
	end 

end 




