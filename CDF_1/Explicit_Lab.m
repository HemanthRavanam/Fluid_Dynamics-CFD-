clear all
close all
clc

%% Defining the mesh
Lx = 1;       
Ly = 1;
Nx = 41;
Ny = 41;
h = 1 / 40;
dt=0.0001;
lambda = dt/(h*h);
T = 0.16;
Nt = T / dt;

x = linspace(0, Lx, Nx);        
y = linspace(0, Ly, Ny);

% Initialize solution matrix with initial condition
u = zeros(Nx, Ny);

%% Initialising the problem
u(:, Ny) = 1 - sin((pi * flipud(y') / 2));   
u(:, 1) = 1 - flipud(y').^3;
u(1, :) = 0;
u(Nx, :) = 1;

error_mag = 1;
error_req = 1e-6;
iterations = 0;

%Tracking the error magnitude
error_track = 0;

%Visualization
x_dom = ((1:Nx)-1).*h;
y_dom = 1-((1:Ny)-1).*h;
[X,Y] = meshgrid(x_dom,y_dom);

%% Calculations
while error_mag > error_req
    for t = 1:Nt
    % Initialize temporary variable for updated solution
    u_new = u;
        for i = 2:(Nx-1)
            for j = 2:(Ny-1)
                u_new(i,j) = u(i,j) + lambda * (u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1) - 4 * u(i,j));
            end
        end
        % Apply boundary conditions
        u_new(:, Ny) = 1 - sin((pi * flipud(y') / 2));
        u_new(:, 1) = 1 - flipud(y').^3;
        u_new(1, :) = 0;                         
        u_new(Nx, :) =1;
    end


    iterations = iterations + 1;
    % Calculation of error magnitude
    error_mag = 0;
    for i = 2:(Nx-1)
        for j = 2:(Ny-1)
            error_mag = error_mag + abs(u(i,j) - u_new(i,j));
            error_track(iterations) = error_mag;
        end
    end

    if rem(iterations, 1000) == 0
        iterations;
        error_mag;
    end

    % Assigning new to be old
    u = u_new;

    %Update visualization every few iterations to avoid flicker
    if mod(iterations, 10) == 0
        contourf(X, Y, u, 30);
        colorbar
        colormap(jet)
        set(gca, 'YDir', 'normal');
        xlabel('X Axis');
        ylabel('Y Axis');
        title({'2D  Heat Conduction in Transient State for Explicit Method.'});
        pause(0.01);
    end
end


