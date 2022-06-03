clc;clear;close all;

rowDim = 10;
colDim = 10;
total_Time = 10;
resistance = 10;
dt = 0.01;
K_p = 100;
K_d = 10;

r0 = 0.001;             % rod radius cross section
Cd = 0.42;              % Half-Sphere Head

% Parameters
rho = 1000;
Area = pi() * r0^2;

[lines,board,hist,result] = maze(rowDim,colDim,[1,1],[rowDim,colDim],true,true);

total_Length=0;
for i = 2:length(result)
    differ = abs(result(i,:) - result(i-1,:));
    total_Length = total_Length + sum(differ);
end

%% Position
time_hist_idx = 0:1:total_Time*(1/dt);
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

% figure()
% plot(position_hist(:,1), position_hist(:,2))
% for i = 2:length(position_hist)
%     velo_error(i,:) = position_hist(i,:) - position_hist(i-1,:);
% end

%% Reference Force & Velo
refer_Velo_mag = total_Length / total_Time;
refer_Force_mag = 1/2 * rho* refer_Velo_mag^2 * Area * Cd;

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

real_robot_sim = sim('Control_Block_2.slx');

%% Plot Simulation

real_Posi = squeeze(real_robot_sim.real_Posi.Data)';
real_Velo = squeeze(real_robot_sim.real_Velo.Data)';
Total_Force = squeeze(real_robot_sim.Total_Force.Data)';
time_hist = real_robot_sim.tout;

v = VideoWriter('Test','MPEG-4');
open(v)

for i = 1:length(real_Posi)-1
    figure(5)
    set(gcf, 'Position', [1000 0 1600 1600]);

    plot(result(:,1),result(:,2), 'r'); hold on
    plot(lines(:,1),lines(:,2), 'b');
    plot(real_Posi(i,1), real_Posi(i,2),'ro','MarkerSize',20);
    title('One Node Robot Trajectory Following');
    xlabel('X');ylabel('Y');
    axis([-0.1,1.1,-0.1,1.1])
    
%     F_arraw_x = [position_hist(i,1), position_hist(i,1) + Total_Force(i,1)*1000];
%     F_arraw_y = [position_hist(i,2), position_hist(i,2) + Total_Force(i,2)*1000];
    quiver(real_Posi(i,1),real_Posi(i,2),Total_Force(i,1)*500,Total_Force(i,2)*500, 5,'m');
    hold off

    M = getframe(gcf);
    writeVideo(v,M);
end
close(v)

