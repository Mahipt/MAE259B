clc; clear; close all;

%% Time setup
total_Time = 5;
dt = 0.001;
N_step = total_Time / dt;
counter = 0;

%% Load Maze  
rowDim = 4;
colDim = 4;
[lines,board,hist,result] = ...
    maze(rowDim,colDim,[1,1],[rowDim,colDim],true,true);
init_pos = 1/(2*rowDim);

%% Controller Gain value
K_p = 100;
K_d = 10;

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
% q = q0;
% u = (q - q0) / dt;

% tolerance 
tol = EI / l^2 * 1e-3; 

%% Blood parameters 
rho_blood = 1025; % Blood density 
Cd = 1.05; % Coefficient of drag
Area = pi*R^2; % Front space area

%% Trajectory config 
total_Length=0;
for i = 2:length(result)
    differ = abs(result(i,:) - result(i-1,:));
    total_Length = total_Length + sum(differ); % traj length 
end

%% Controller Initialization 
% Position 
position_hist = result(1,:);
step_time_index = round( total_Time/(length(result)-1) * (1/dt)+0.5 );

% Filling the Trace
for i = 2:length(result)
    dx = linspace( result(i-1,1) , result(i,1) , step_time_index);
    dy = linspace( result(i-1,2) , result(i,2) , step_time_index);
    
    if i == 2
        start_idx = 1;
        end_idx   = (i-1)*step_time_index;
        position_hist(start_idx:end_idx, :) = [dx;dy]';
    else
        start_idx = (i-2)*step_time_index;
        end_idx   = (i-1)*step_time_index;
        position_hist(start_idx-(i-3):end_idx-(i-2), :) = [dx;dy]';
    end
end

% Reference Force & Velo
refer_Velo_mag = total_Length / total_Time;
refer_Force_mag = 1/2 * rho_blood * refer_Velo_mag^2 * Area * Cd;

% Initialize
unit_vector = zeros(length(position_hist),2);
refer_Velo  = zeros(length(position_hist),2);
refer_Force = zeros(length(position_hist),2);

% Initialize First
unit_vector(1,:)= (position_hist(2,:) - position_hist(1,:)) / abs(sum(position_hist(2,:) - position_hist(1,:)));
refer_Velo(1,:) = refer_Velo_mag * unit_vector(1,:);
refer_Force(1,:) = refer_Force_mag * unit_vector(1,:);

% Force and Velocity at Ideal Situation
for i = 2:length(position_hist)
    unit_vector(i,:)= (position_hist(i,:) - position_hist(i-1,:)) / abs(sum(position_hist(i,:) - position_hist(i-1,:)));
    refer_Velo(i,:) = refer_Velo_mag * unit_vector(i,:);
    refer_Force(i,:) = refer_Force_mag * unit_vector(i,:);
end

%% Simulink
position_hist = timetable(position_hist,'TimeStep',seconds(dt));
refer_Velo = timetable(refer_Velo,'TimeStep',seconds(dt));
refer_Force = timetable(refer_Force,'TimeStep',seconds(dt));
unit_vector = timetable(unit_vector,'TimeStep',seconds(dt));

end_time = length(position_hist.position_hist)/(1/dt);

real_robot_sim = sim('Control_Block_3.slx');

%% Get data from Simulation
real_Posi = squeeze(real_robot_sim.real_Posi.Data)';
real_Velo = squeeze(real_robot_sim.real_Velo.Data)';
Total_Force = squeeze(real_robot_sim.Total_Force.Data)';
q = squeeze(real_robot_sim.q.Data)';
time_hist = real_robot_sim.tout;

%% DER Simulation
% Ex_force = zeros(2*N,1);
% 
% % time marching scheme
% for i = 2 : N_step
%     
%     fprintf('Time = %f\n',(i-1) * dt);
% 
%     q = q0; % Guess
% 
%     err = 10 * tol;
%     % Newton Raphson
%     while err > tol
% 
%     % Inertia 
%     f = M/dt * ((q - q0) / dt - u);
%     J = M/ dt^2;
% 
%     % Elastic forces 
%     % Linear Spring 1 between node 1 and 2
%     for k = 1: N-1
%         xk = q(2*k-1);
%         yk = q(2*k);
%         xkp1 = q(2*k+1);
%         ykp1 = q(2*k+2); 
%         dF = gradEs(xk, yk, xkp1, ykp1, dl, EA);
%         dJ = hessEs(xk, yk, xkp1, ykp1, dl, EA);
%         f(2*k-1:2*k+2) = f(2*k-1:2*k+2) + dF;
%         J(2*k-1:2*k+2, 2*k-1:2*k+2) = ...
%         J(2*k-1:2*k+2, 2*k-1:2*k+2) + dJ;
%     end
% 
%     % Bending spring 
%     for k = 2:N-1
%         xkm1 = q(2*k-3);
%         ykm1 = q(2*k-2);
%         xk = q(2*k-1);
%         yk = q(2*k);
%         xkp1 = q(2*k+1);
%         ykp1 = q(2*k+2);
%         curvature0 = 0;
%         dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, dl, EI);
%         dJ = hessEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, dl, EI);
%         f(2*k-3:2*k+2) = f(2*k-3:2*k+2) + dF;
%         J(2*k-3:2*k+2, 2*k-3:2*k+2) = ...
%             J(2*k-3:2*k+2, 2*k-3:2*k+2) + dJ;
%     end
%     
%     % External force 
%     Ex_force(1) = Total_Force(i,1);
%     Ex_force(2) = Total_Force(i,2);
% 
%     % Weight 
%     f = f + Ex_force;
% 
%     % Update 
%     q(1:(2*N)) = q(1:(2*N)) - J(1:(2*N),1:(2*N)) \ f(1:(2*N));
%     err = sum(abs(f(1:(2*N))));
%     end
% 
%     % Update
%     u = (q - q0) / dt;
%     q0 = q;
%     
%     % plot of sphere position
%     figure(1);
%     plot(result(:,1),result(:,2), 'y');
%     hold on
%     plot(lines(:,1),lines(:,2), 'b');
%     hold on
%     plot(q(1:2:end),q(2:2:end),'ro-');
%     plot(real_Posi(i,1), real_Posi(i,2),'ro');
% %   axis equal
%     axis([-0.5,1.1,-0.1,1.1])
%     hold off 
%     drawnow
% end 

%% Plot 
for i = 1:length(real_Posi)-1
    figure(5)
    plot(result(:,1),result(:,2), 'y');
    hold on
    plot(lines(:,1),lines(:,2), 'b');
    hold on
    plot(q(i,1:2:end),q(i,2:2:end),'ro-');
    hold on
    plot(real_Posi(i,1), real_Posi(i,2),'ro');
    axis([-0.5,1.1,-0.1,1.1])

    % quiver(real_Posi(i,1),real_Posi(i,2),Total_Force(i,1)*1000,Total_Force(i,2)*1000, 5);
    hold off
end

