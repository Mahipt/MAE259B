%% Time setup
total_Time = 5;
dt = 0.01;
N_step = total_Time / dt;

%% Beam parameters 
N = 50; % Number of nodes 
l = 0.5; % beam length
R = 0.002; % beam outer radius 
r = 0.001; % beam inner radius
rho_beam = 980; % PDMS density
E = 1.47e6; % Youngs Modulus

dl = l / (N - 1); % length of each segment 
A = pi*(R^2-r^2); % Beam cross-section 

m = A*l*rho_beam/(N-1); % mass of each node  
I = pi/4*(R^4 - r^4); % Moment of inerita 

EA = E*A; % Stretching stiffness 
EI = E*I; % Bending stiffness 

% Beam Mass Matrix 
M = zeros(2*N,2*N);
for i = 1 : 2*N
    M(i,i) = m;
end 
%% Beam initial configuration 
% Initial position 
nodes_I = zeros(N,2);
for i = 1:N 
    nodes_I(i,1) = (1-i) * dl + init_pos;
    nodes_I(i,2) = init_pos;
end 

% Initial DOF vector
q0 = zeros(2*N,1);
for i = 1:N
    q0 (2*i-1) = nodes_I(i,1); % x coordinate 
    q0 (2*i) = nodes_I(i,2); % y coordinate
end

% Initial position and velo 
q = q0;
u = (q - q0) / dt;

% tolerance 
tol = EI / l^2 * 1e-3; 
