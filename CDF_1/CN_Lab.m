clear all
close all
clc

%% Defining Parameters and Mesh
Lx = 1;       
Ly = 1;       
Nx = 41;
Ny = 41;
dt = 0.0001;
T = 0.16;     
Nt = T / dt;

x = linspace(0, Lx, Nx);        
y = linspace(0, Ly, Ny);   
dx= abs(x(2)-x(1));

% Discretization
h = dx;
lambda = dt / (2*h*h);
iteration=0;

% Initialize solution matrix with initial condition
u = zeros(Nx, Ny);

% Boundary Conditions
u(:, Ny) = 1 - sin((pi * flipud(y') / 2));   
u(:, 1) = 1 - flipud(y').^3;
u(1, :) = 0;
u(Nx, :) = 1; 

% Set up figure for dynamic plotting
figure;
colorbar;
colormap(jet);
title('2D Heat Conduction in Transient State (Crank-Nicolson Method)');
xlabel('x');
ylabel('y');

%Visualization
x_dom = ((1:Nx)-1).*h;
y_dom = 1-((1:Ny)-1).*h;
[X,Y] = meshgrid(x_dom,y_dom);

% Time-Stepping using Crank-Nicolson Method
for t = 1:Nt
    % Initialize temporary variable for updated solution
    u_new = u;
    
    for i = 2:Nx-1
        for j = 2:Ny-1
            % Calculate the right-hand side for Crank-Nicolson
            u_new(i, j) = (lambda) * (u(i+1, j) + u(i-1, j) + u(i, j+1) + u(i, j-1)) ...
                          + (1 - 4* lambda) * u(i, j);
        end
    end
    
    % Apply boundary conditions
    u_new(:, Ny) = 1 - sin((pi * flipud(y') / 2));
    u_new(:, 1) = 1 - flipud(y').^3;
    u_new(1, :) = 0;                         
    u_new(Nx, :) =1;

    % Update solution for the next time step
    u = u_new;

    iteration = iteration + 1;
    
    % Clear the current plot and redraw with updated data
    cla;
    contourf(X, Y, u, 30);
    colorbar;
    colormap(jet);
    title({'2D  Heat Conduction in Transient State for CN Method.'})
    xlabel('x');
    ylabel('y')
    drawnow;
    pause(0.01);
end
