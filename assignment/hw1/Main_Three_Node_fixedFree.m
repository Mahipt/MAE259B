%% Global variables

%% Physical parameters
% Number of vertices
N = 3;

% Time step size
dt = 0.01; % second

% Rod length
RodLength = 0.1; % meter

% Discrete length
deltaL = RodLength / (N-1);

% Radius of spheres
R1 = 0.005;
R2 = 0.025;
R3 = 0.005;

% Density
rho_metal = 7000;
rho_f = 1000;
rho = rho_metal - rho_f;

% Rod radius
r0 = 0.001;

% Young's modulus
Y = 1e9; % Using Y instead of E to avoid ambiguity

% Gravity
g = 9.8; % m/s^2

% Viscosity
visc = 1000; % Pa-s

% Total time
totalTime = 10; % seconds

% Utility quantities
ne = N - 1; % Number of edges
EI = Y * pi * r0^4 / 4;
EA = Y * pi * r0^2;

% Geometry
nodes = zeros(N, 2);
for c = 1:N
    nodes(c,1) = (c-1) * deltaL;
%     nodes(c,2) = 0;
end

% Mass matrix
M = zeros(2*N,2*N);
M(1,1) = 4/3*pi*R1^3*rho_metal;
M(2,2) = 4/3*pi*R1^3*rho_metal;
M(3,3) = 4/3*pi*R2^3*rho_metal;
M(4,4) = 4/3*pi*R2^3*rho_metal;
M(5,5) = 4/3*pi*R3^3*rho_metal;
M(6,6) = 4/3*pi*R3^3*rho_metal;

% Viscous damping matrix
C = zeros(6,6);
C1 = 6*pi*visc*R1;
C2 = 6*pi*visc*R2;
C3 = 6*pi*visc*R3;
C(1,1) = C1;
C(2,2) = C1;
C(3,3) = C2;
C(4,4) = C2;
C(5,5) = C3;
C(6,6) = C3;

% Gravity
W = zeros(2*N,1);
W(2) = -4/3*pi*R1^3*rho*g;
W(4) = -4/3*pi*R2^3*rho*g;
W(6) = -4/3*pi*R3^3*rho*g;

% Initial DOF vector
q0 = zeros(2*N,1);
for c=1:N
    q0 ( 2*c - 1 ) = nodes(c,1); % x coordinate
    q0 ( 2*c ) = nodes(c,2); % y coordinate
end

% New position and velocity
q = q0; % DOF vector
u = (q - q0) / dt; % Velocity vector

% Number of time steps
Nsteps = round( totalTime / dt );
all_mid_y = zeros( Nsteps, 1); % y-position of R2
all_mid_v = zeros( Nsteps, 1); % y-velocity of R2

all_mid_y(1) = q(4);
all_mid_v(1) = u(4);

% Tolerance
tol = EI / RodLength^2 * 1e-3;

% Free and fixed index
free_index = [1,2,5,6];
fixed_index = [3,4];
imposedV = -0.01 * 10;

% Time marching scheme
for c=2:Nsteps
    
    fprintf('Time = %f\n', (c-1) * dt );
    
    q_fixed = [ nodes(2,1);
        imposedV * (c-1) * dt];
    
    q = q0; % Guess
    
    % Updating the fixed index
    q(fixed_index) = q_fixed;
    q_free = q(free_index);
    
    % Newton Raphson
    err = 10 * tol;
    while err > tol        
        % Inertia
        f = M / dt * ( (q-q0) / dt - u );
        J = M / dt^2;
        
        %
        % Elastic forces
        %
        % Linear spring 1 between nodes 1 and 2
        for k=1:N-1
            xk = q(2*k-1);
            yk = q(2*k);
            xkp1 = q(2*k+1);
            ykp1 = q(2*k+2);
            dF = gradEs(xk, yk, xkp1, ykp1, deltaL, EA);
            dJ = hessEs(xk, yk, xkp1, ykp1, deltaL, EA);
            f(2*k-1:2*k+2) = f(2*k-1:2*k+2) + dF;
            J(2*k-1:2*k+2,2*k-1:2*k+2) = ...
                J(2*k-1:2*k+2,2*k-1:2*k+2) + dJ;
        end
        
        % Bending spring between nodes 1, 2, and 3
        for k=2:N-1
            xkm1 = q(2*k-3);
            ykm1 = q(2*k-2);
            xk = q(2*k-1);
            yk = q(2*k);
            xkp1 = q(2*k+1);
            ykp1 = q(2*k+2);
            curvature0 = 0;
            dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, ...
                curvature0, deltaL, EI);
            dJ = hessEb(xkm1, ykm1, xk, yk, xkp1, ykp1, ...
                curvature0, deltaL, EI);
            f(2*k-3:2*k+2) = f(2*k-3:2*k+2) + dF;
            J(2*k-3:2*k+2,2*k-3:2*k+2) = ...
                J(2*k-3:2*k+2,2*k-3:2*k+2) + dJ;
        end
        
        % Viscous force
        f = f + C * ( q - q0 ) / dt;
        J = J + C / dt;
        
        % Weight
        f = f - W;
        
        % Update
        f_free = f(free_index);
        J_free = J(free_index, free_index);
        q_free = q_free - J_free \ f_free;
        
        err = sum( abs(f_free) );
        
        % Plug free DOFs back into the full DOF vector
        q(free_index) = q_free;

    end
    
    % Update
    u = (q - q0) / dt; % Velocity
    q0 = q; % Old position
    
    figure(1);
    plot( q(1:2:end), q(2:2:end), 'ro-');
    axis equal
    drawnow
    
    % Store
    all_mid_y(c) = q(4);
    all_mid_v(c) = u(4);
end

figure(2);
timeArray = (1:Nsteps) * dt;
plot(timeArray, all_mid_v, 'k-');
xlabel('Time, t [sec]');
ylabel('Velocity of mid-node, v [meter/sec]');
