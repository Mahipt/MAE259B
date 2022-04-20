close all; clear; clc; 


vel_timeStep = zeros(19, 2); 

for j = 1:5
    fprintf('Node = %f\n', 0.1^j);
    %% Global variables

    %% Physical parameters
    % Number of vertices (must be an odd number) 
    N = 21;

    % Time step size1
    % Total time
    totalTime = 50; % seconds
    dt = 0.1^(j); % second

    % Rod length
    RodLength = 0.1; % meter

    % Discrete length
    deltaL = RodLength / (N-1);

    % Radius of spheres
    R = zeros(N, 1); 
    for i = 1:N % Initialize the sphere 
    	if i ~= ((N + 1) / 2) % if not the middle sphere
    		R(i) = deltaL / 10; 
    	else
    		R(i) = 0.025; % middle sphere 
    	end
    end 

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


    % Utility quantities
    ne = N - 1; % Number of edges
    EI = Y * pi * r0^4 / 4;
    EA = Y * pi * r0^2;

    % Mass matrix
    M = zeros(2 * N, 2 * N); 
    for i = 1:N % Assign value to the diagonal matrix
    	M(2 * i - 1, 2 * i - 1) = (4/3) * pi * R(i)^3 * rho_metal; 
    	M(2 * i, 2 * i) = (4/3) * pi * R(i)^3 * rho_metal; 
    end 

    % Viscous damping matrix
    C = zeros(2 * N, 2 * N); 
    for i = 1:N % Assign value to the diagonal matrix
    	C(2 * i - 1, 2 * i - 1) = 6 * pi * visc * R(i); 
    	C(2 * i, 2 * i) = 6 * pi * visc * R(i); 
    end 

    % Gravity
    W = zeros(2 * N, 1);
    for i = 1:N 
    	W(2 * i) = (-4/3) * pi * R(i)^3 * rho * g;  
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
    all_mid_y = zeros(Nsteps, 1); % y-position (mid sphere)
    all_mid_v = zeros(Nsteps, 1); % y-velocity (mid sphere)

    all_mid_y(1) = q((N + 1));
    all_mid_v(1) = u((N + 1));

    % Tolerance
    tol = EI / RodLength^2 * 1e-3;

    % Time marching scheme
    for c = 2:Nsteps 
        %fprintf('Time = %f\n', (c - 1) * dt );
        
        q = q0; % Guess
        % Newton Raphson
        err = 10 * tol;
        while err > tol
            % Inertia
            f = M / dt * ( (q-q0) / dt - u );
            J = M / dt^2;
            
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
        
    	        
            % Viscous force
            f = f + C * ( q - q0 ) / dt;
            J = J + C / dt;
            
            % Weight
            f = f - W;   

            %tempt = J \ f  
            % Update
            q = q - J \ f;
            
            err = sum( abs(f) );
        end

        % Update
        u = (q - q0) / dt; % Velocity
        q0 = q; % Old position
        
        %{
        figure(1);
        plot( q(1:2:end), q(2:2:end), 'ro-');
        axis equal
        drawnow
        %}
        
        % Store
        all_mid_y(c) = q((N + 1));
        all_mid_v(c) = u((N + 1));
    end

    %{
    figure(2);
    timeArray = (1:Nsteps) * dt;
    plot(timeArray, all_mid_v, 'k-');
    xlabel('Time, t [sec]');
    ylabel('Velocity of mid-node, v [meter/sec]');


    figure(3); 
    timeArray = (1:Nsteps) * dt;
    plot(timeArray, all_mid_y, 'k-');
    xlabel('Time, t [sec]');
    ylabel('Positio of mid-node, y[meter]');
    %}

    vel_nodes(j, 2) = u(end); 
    vel_nodes(j, 1) = dt; 

end

figure(1); 
plot(vel_nodes(:, 1), vel_nodes(:, 2), 'k-'); 
xlabel('Time steps, [Steps size]');
ylabel('Terminal Velocity, v[meter/s]');




