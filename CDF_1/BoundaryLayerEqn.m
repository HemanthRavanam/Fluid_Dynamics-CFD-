clc;
clear;
close all;

%% Parameters
Lx = 1; 
Ly = 0.2; %Ly = 2*delta_99, delta_99= 4.91*x/sqrt(Re)
Re = 1e4; 
kinematic_viscosity = 1.516e-05;
U_inf = 20;  
nu = U_inf * Lx / Re;

% Grid parameters
Nx = 1000; 
Ny = 100;  
dx = Lx / Nx;  
dy = Ly / Ny;  

% Initialize u and v velocity fields (zero at the start)
u = zeros(Ny, Nx); 
v = zeros(Ny, Nx);  

% Boundary Conditions
u(:,1) = 0;
u(:, 1) = U_inf;  
u(end, :) = U_inf; 
u(1,1) = 0;

v(1, :) = 0;
v(:, 1) = 0;
v(1, 1) = 0;

% Loop through the grid and update u and v velocities
for i = 1:Nx-1
    for j = 2:Ny-1
        % Discretize x-momentum equation (solve for u at the next point)
        u_next = u(j+1, i);  % u(i+1,j-1)
        u_current = u(j, i);  % u(i,j-1)
        u_prev = u(j-1, i);  % u(i-1,j-1)

        % Update u
        u(j, i+1) = u_current + nu * (u_next - 2 * u_current + u_prev) / u_current * dx / (dy^2) ...
                   - (u_next - u_prev) / 2 * v(j, i) / u_current * dx / dy;
    end
end

% Loop to update v velocity components
for i = 1:Nx-1 
    for j = 1:Ny-1
        % Update v 
        v(j+1, i) = v(j, i) - (dy/dx) * (u(j, i+1) - u(j,i));
    end
end

% x and y grid arrays for plotting
x_arr = linspace(0, Lx, Nx);
y_arr = linspace(0, Ly, Ny);
[X,Y] = meshgrid(x_arr,y_arr);

% Compute theoretical boundary layer thickness (delta_x)
delta_x = (4.91 ./ sqrt(Re)) .* sqrt(x_arr);

%% Plot 1: u-Velocity Contour
figure;

% Contour plot for u-velocity
contourf(X, Y, u, 100, 'LineStyle', 'none');  % Contour of u-velocity
colorbarHandle = colorbar;  % Create colorbar and get handle
colormap(jet);
ylabel(colorbarHandle, 'Velocity (m/s)', 'FontSize', 12);  % Add units to colorbar
hold on;

% Overlay theoretical boundary layer thickness
plot(x_arr, delta_x, 'r-', 'LineWidth', 2, 'DisplayName', 'Theoretical \delta(x)');

% Labels and title
xlabel('x (m)');
ylabel('y (m)');
title('u-Velocity Contour');
legend('Location', 'northeast');
grid on;

%% Plot 2: v-Velocity Profiles
figure;
hold on;

% Select x positions where you want to extract the v-velocity profile
x_positions = 0.1:0.1:Lx;  % x = 0.1, 0.2, ..., 1

% Plot v-velocity profiles at selected x positions
for x_pos = x_positions
    % Find the closest index for the current x position
    [~, idx_x] = min(abs(x_arr - x_pos));
    
    % Extract the v-velocity profile at the selected x position
    v_profile = v(:, idx_x);  % All y-values for the selected x position
    
    % Plot the v-velocity profile
    plot(v_profile, y_arr, 'LineWidth', 2, 'DisplayName', ['x = ', num2str(x_pos)]);
end

% Labels and title
xlabel('v-Velocity (m/s)');
ylabel('y (m)');
title('v-Velocity Profiles at Various x Positions');
legend('Location', 'best');
grid on;

%% Plot 3: u-Velocity Profiles
figure;
hold on;

% Plot u-velocity profiles at selected x positions
for x_pos = x_positions
    % Find the closest index for the current x position
    [~, idx_x] = min(abs(x_arr - x_pos));
    
    % Extract the u-velocity profile at the selected x position
    u_profile = u(:, idx_x);  % All y-values for the selected x position
    
    % Plot the u-velocity profile
    plot(u_profile, y_arr, 'LineWidth', 2, 'DisplayName', ['x = ', num2str(x_pos)]);
end

% Labels and title
xlabel('u-Velocity (m/s)');
ylabel('y (m)');
title('u-Velocity Profiles at Various x Positions');
legend('Location', 'best');
grid on;