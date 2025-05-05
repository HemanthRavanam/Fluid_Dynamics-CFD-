clear all
close all
clc

%% Defining the Grid parameters
i_max = 51;  
j_max = 51;
dom_length = 1;
h = dom_length/(i_max-1);
x = 0:h:dom_length; 
y = 0:h:dom_length; 
Re = 1; 
nu = 1/Re;
% Under-relaxation factors
alpha = 0.01;
alpha_p = 0.8;

%% Initializing the variables
u_final(i_max,j_max)=0;
v_final(i_max,j_max)=0;
p_final(i_max,j_max)=0;
vorticity_final(i_max,j_max)=0;
u_final(1,:) = 1;

%Staggered variables
u_star(i_max+1,j_max)=0;
d_e(i_max+1,j_max)=0;

v_star(i_max,j_max+1)=0;
d_n(i_max,j_max+1)=0;

p_star(i_max+1,j_max+1)=0;
pc(i_max+1,j_max+1)=0;

q(i_max+1,j_max+1)=0;

u_star(1,:)=2;

u_new(i_max+1,j_max)=0;
v_new(i_max,j_max+1)=0;
p_new(i_max+1,j_max+1)=0;
vorticity_new(i_max+1,j_max+1)=0;
u_new(1,:)=1;

%% Solving the governing equations with initial guess u* v* p*
error = 1;
iterations = 0;
error_req = 1e-7;

while error > error_req
    % x-momentum eq. - Interior Nodes
    for i = 2:i_max
        for j = 2:j_max - 1
            u_E = 0.5*(u_star(i,j) + u_star(i,j+1));
            u_W = 0.5*(u_star(i,j) + u_star(i,j-1));
            v_N = 0.5*(v_star(i-1,j) + v_star(i-1,j+1));
            v_S = 0.5*(v_star(i,j) + v_star(i,j+1));

            a_E = -0.5*u_E*h + nu;
            a_W = 0.5*u_W*h + nu;
            a_N = -0.5*v_N*h + nu;
            a_S = 0.5*v_S*h + nu;

            a_e = 0.5*u_E*h - 0.5*u_W*h + 0.5*v_N*h - 0.5*v_S*h + 4*nu;

            A_e = -h;
            d_e(i,j) = A_e/a_e;

            u_star_m = (a_E*u_star(i,j+1) + a_W*u_star(i,j-1) + a_N*u_star(i-1,j) + a_S*u_star(i+1,j))/a_e + d_e(i,j)*(p_star(i,j+1) - p_star(i,j));
            u_star(i,j)=(1-alpha)*u_star(i,j)+alpha*u_star_m;
        end
    end

    % x-momentum eq. - Boundaries
    u_star(1,:) = 2 - u_star(2,:);
    u_star(i_max + 1,:) = -u_star(i_max,:);
    u_star(2:i_max,1) = 0;
    u_star(2:i_max,j_max) = 0;

    % y-momentum eq. - Interior Nodes
    for i = 2:i_max - 1
        for j = 2:j_max
            u_E = 0.5*(u_star(i,j) + u_star(i+1,j));
            u_W = 0.5*(u_star(i,j-1) + u_star(i+1,j-1));
            v_N = 0.5*(v_star(i-1,j) + v_star(i,j));
            v_S = 0.5*(v_star(i,j) + v_star(i+1,j));

            a_E = -0.5*u_E*h + nu;
            a_W = 0.5*u_W*h + nu;
            a_N = -0.5*v_N*h + nu;
            a_S = 0.5*v_S*h + nu;

            a_n = 0.5*u_E*h - 0.5*u_W*h + 0.5*v_N*h - 0.5*v_S*h + 4*nu;

            A_n = -h;
            d_n(i,j) = A_n/a_n;

            v_star_m = (a_E*v_star(i,j+1) + a_W*v_star(i,j-1) + a_N*v_star(i-1,j) + a_S*v_star(i+1,j))/a_n + d_n(i,j)*(p_star(i,j) - p_star(i+1,j));
            v_star(i,j)=(1-alpha)*v_star(i,j)+alpha*v_star_m;
        end
    end

    % y-momentum eq. - Boundaries
    v_star(:,1) = -v_star(:,2);
    v_star(:,j_max + 1) = -v_star(:,j_max);
    v_star(1,2:j_max) = 0;
    v_star(i_max,2:j_max) = 0;

    % Initialising the pressure correction
    pc(1:i_max+1,1:j_max+1)=0;

    % Continuity equation for pressure correction
    for i = 2:i_max
        for j = 2:j_max
            a_E = -d_e(i,j)*h;
            a_W = -d_e(i,j-1)*h;
            a_N = -d_n(i-1,j)*h;
            a_S = -d_n(i,j)*h;
            a_P = a_E + a_W + a_N + a_S;
            q(i,j) = -(u_star(i,j) - u_star(i,j-1))*h + (v_star(i,j) - v_star(i-1,j))*h;

            pc(i,j) = (a_E*pc(i,j+1) + a_W*pc(i,j-1) + a_N*pc(i-1,j) + a_S*pc(i+1,j) + q(i,j))/a_P;
        end
    end

    % Correcting the pressure field
    for i = 2:i_max
        for j = 2:j_max
            p_new(i,j) = p_star(i,j) + alpha_p*pc(i,j);
        end
    end

    % Continuity eq. - Boundary Nodes
    p_new(1,:) = p_new(2,:);
    p_new(i_max + 1,:) = p_new(i_max,:);
    p_new(:,1) = p_new(:,2);
    p_new(:,j_max + 1) = p_new(:,j_max);

    % Correcting the velocities
    for i = 2:i_max
        for j = 2:j_max - 1
            u_new(i,j) = u_star(i,j) + d_e(i,j)*(pc(i,j+1) - pc(i,j));
        end
    end

    % x-momentum eq. - Boundary Nodes
    u_new(1,:) = 2 - u_new(2,:);
    u_new(i_max + 1,:) = -u_new(i_max,:);
    u_new(2:i_max,1) = 0;
    u_new(2:i_max,j_max) = 0;

    for i = 2:i_max - 1
        for j = 2:j_max
            v_new(i,j) = v_star(i,j) + d_n(i,j)*(pc(i,j) - pc(i+1,j));
        end
    end

    % y-momentum eq. - Boundary
    v_new(:,1) = -v_new(:,2);
    v_new(:,j_max + 1) = -v_new(:,j_max);
    v_new(1,2:j_max) = 0;
    v_new(i_max,2:j_max) = 0;

    % Continuity residual as error measure
    error = 0;
    for i = 2:i_max
        for j = 2:j_max
            error = error + abs(q(i,j));
        end
    end

    % Finishing the iteration
    u_star = u_new;
    v_star = v_new;
    p_star = p_new;
    iterations = iterations + 1;

end

% After the results converged, Mapping the staggered variables to
for i = 1:i_max
    for j = 1:j_max
        u_final(i,j) = 0.5*(u_star(i,j) + u_star(i+1,j));
        v_final(i,j) = 0.5*(v_star(i,j) + v_star(i,j+1));
        p_final(i,j) = 0.25*(p_star(i,j) + p_star(i,j+1) + p_star(i+1,j) + p_star(i+1,j+1));
    end
end

% Compute vorticity in the interior nodes using central difference
for i = 2:i_max-1
    for j = 2:j_max-1
        vorticity(i,j) = (v_final(i,j+1) - v_final(i,j-1))/(2*h) - (u_final(i+1,j) - u_final(i-1,j))/(2*h);
    end
end

% Compute vorticity at boundaries using one-sided differences
for i = 2:i_max-1
    % Bottom boundary (j = 1)
    vorticity(i,1) = (v_final(i,2) - v_final(i,1))/h - (u_final(i+1,1) - u_final(i-1,1))/(2*h);

    % Top boundary (j = j_max)
    vorticity(i,j_max) = (v_final(i,j_max) - v_final(i,j_max-1))/h - (u_final(i+1,j_max) - u_final(i-1,j_max))/(2*h);
end

for j = 2:j_max-1
    % Left boundary (i = 1)
    vorticity(1,j) = (v_final(1,j+1) - v_final(1,j-1))/(2*h) - (u_final(2,j) - u_final(1,j))/h;

    % Right boundary (i = i_max)
    vorticity(i_max,j) = (v_final(i_max,j+1) - v_final(i_max,j-1))/(2*h) - (u_final(i_max,j) - u_final(i_max-1,j))/h;
end

% Corner points (can be set as average of neighboring values)
vorticity(1,1) = 0.5 * (vorticity(2,1) + vorticity(1,2));
vorticity(1,j_max) = 0.5 * (vorticity(2,j_max) + vorticity(1,j_max-1));
vorticity(i_max,1) = 0.5 * (vorticity(i_max-1,1) + vorticity(i_max,2));
vorticity(i_max,j_max) = 0.5 * (vorticity(i_max-1,j_max) + vorticity(i_max,j_max-1));
% Interior points for x- and y- momentum equations
du_dx = zeros(i_max, j_max);
du_dy = zeros(i_max, j_max);
dv_dx = zeros(i_max, j_max);
dv_dy = zeros(i_max, j_max);

for i = 2:i_max-1
    for j = 2:j_max-1
        du_dx(i,j) = (u_final(i,j+1) - u_final(i,j-1)) / (2 * h);
        du_dy(i,j) = (u_final(i+1,j) - u_final(i-1,j)) / (2 * h);

        dv_dx(i,j) = (v_final(i,j+1) - v_final(i,j-1)) / (2 * h);
        dv_dy(i,j) = (v_final(i+1,j) - v_final(i-1,j)) / (2 * h);
    end
end

% Compute LHS (convective terms)
lhs_x_momentum = u_final .* du_dx + v_final .* du_dy;
lhs_y_momentum = u_final .* dv_dx + v_final .* dv_dy;

% Compute RHS (pressure gradient and viscous terms)
dp_dx = zeros(i_max, j_max);
dp_dy = zeros(i_max, j_max);
d2u_dx2 = zeros(i_max, j_max);
d2u_dy2 = zeros(i_max, j_max);
d2v_dx2 = zeros(i_max, j_max);
d2v_dy2 = zeros(i_max, j_max);

for i = 2:i_max-1
    for j = 2:j_max-1
        dp_dx(i,j) = (p_final(i,j+1) - p_final(i,j-1)) / (2 * h);
        dp_dy(i,j) = (p_final(i+1,j) - p_final(i-1,j)) / (2 * h);

        d2u_dx2(i,j) = (u_final(i,j+1) - 2 * u_final(i,j) + u_final(i,j-1)) / h^2;
        d2u_dy2(i,j) = (u_final(i+1,j) - 2 * u_final(i,j) + u_final(i-1,j)) / h^2;

        d2v_dx2(i,j) = (v_final(i,j+1) - 2 * v_final(i,j) + v_final(i,j-1)) / h^2;
        d2v_dy2(i,j) = (v_final(i+1,j) - 2 * v_final(i,j) + v_final(i-1,j)) / h^2;
    end
end

rhs_x_momentum = -dp_dx + nu * (d2u_dx2 + d2u_dy2);
rhs_y_momentum = -dp_dy + nu * (d2v_dx2 + d2v_dy2);

% Compute the difference between LHS and RHS
momentum_x_diff = lhs_x_momentum - rhs_x_momentum;
momentum_y_diff = lhs_y_momentum - rhs_y_momentum;

% Interior points for x- and y- momentum equations
du_dx = zeros(i_max, j_max);
du_dy = zeros(i_max, j_max);
dv_dx = zeros(i_max, j_max);
dv_dy = zeros(i_max, j_max);

for i = 2:i_max-1
    for j = 2:j_max-1
        du_dx(i,j) = (u_final(i,j+1) - u_final(i,j-1)) / (2 * h);
        du_dy(i,j) = (u_final(i+1,j) - u_final(i-1,j)) / (2 * h);

        dv_dx(i,j) = (v_final(i,j+1) - v_final(i,j-1)) / (2 * h);
        dv_dy(i,j) = (v_final(i+1,j) - v_final(i-1,j)) / (2 * h);
    end
end

% Compute LHS (convective terms)
lhs_x_momentum = u_final .* du_dx + v_final .* du_dy;
lhs_y_momentum = u_final .* dv_dx + v_final .* dv_dy;

% Compute RHS (pressure gradient and viscous terms)
dp_dx = zeros(i_max, j_max);
dp_dy = zeros(i_max, j_max);
d2u_dx2 = zeros(i_max, j_max);
d2u_dy2 = zeros(i_max, j_max);
d2v_dx2 = zeros(i_max, j_max);
d2v_dy2 = zeros(i_max, j_max);

for i = 2:i_max-1
    for j = 2:j_max-1
        dp_dx(i,j) = (p_final(i,j+1) - p_final(i,j-1)) / (2 * h);
        dp_dy(i,j) = (p_final(i+1,j) - p_final(i-1,j)) / (2 * h);

        d2u_dx2(i,j) = (u_final(i,j+1) - 2 * u_final(i,j) + u_final(i,j-1)) / h^2;
        d2u_dy2(i,j) = (u_final(i+1,j) - 2 * u_final(i,j) + u_final(i-1,j)) / h^2;

        d2v_dx2(i,j) = (v_final(i,j+1) - 2 * v_final(i,j) + v_final(i,j-1)) / h^2;
        d2v_dy2(i,j) = (v_final(i+1,j) - 2 * v_final(i,j) + v_final(i-1,j)) / h^2;
    end
end

rhs_x_momentum = -dp_dx + nu * (d2u_dx2 + d2u_dy2);
rhs_y_momentum = -dp_dy + nu * (d2v_dx2 + d2v_dy2);

% Compute the difference between LHS and RHS
momentum_x_diff = lhs_x_momentum - rhs_x_momentum;
momentum_y_diff = lhs_y_momentum - rhs_y_momentum;

% Compute continuity equation (divergence of velocity field)
continuity_eqn = du_dx + dv_dy;  

%% Contour's
x_dom = ((1:i_max)-1).*h;
y_dom = 1-((1:i_max)-1).*h;
[X,Y] = meshgrid(x_dom,y_dom);
figure(1);
contourf(X,Y,u_final, 30, 'LineStyle', 'none')
colorbar
colormap('jet')
xlabel('x')
ylabel('y')
title({['U Velocity - 2D Cavity Flow with Re = ', num2str(Re)]});
% print(gcf,'u-velocity Re1000','-dpng','-r300');

figure(2);
contourf(X,Y,v_final, 30, 'LineStyle', 'none')
colorbar
colormap('jet')
xlabel('x')
ylabel('y')
title({['V Velocity - 2D Cavity Flow with Re = ', num2str(Re)]});
% print(gcf,'v-velocity Re1','-dpng','-r300');

figure(3);
contourf(X,Y,p_final, 30, 'LineStyle', 'none')
colorbar
colormap('jet')
xlabel('x')
ylabel('y')
title({['Pressure - 2D Cavity Flow with Re = ', num2str(Re)]});
% print(gcf,'pressure Re1','-dpng','-r300');

figure(4);
hold on
quiver(X, Y, u_final, v_final, 3, 'k');
xlabel('x');
ylabel('y');
title({['Stream Lines - 2D Cavity Flow with Re = ', num2str(Re)]});
% print(gcf,'Stream line Re1','-dpng','-r300');

% Plot vorticity contour
figure(5);
contourf(X, Y, vorticity, 30, 'LineStyle', 'none');
colorbar;
colormap('jet');
xlabel('x');
ylabel('y');
title({['Vorticity Contour - 2D Cavity Flow with Re = ', num2str(Re)]});
% print(gcf,'Vorticity Re 1','-dpng','-r300');

figure(6);
subplot(1,2,1);
pcolor(y, x, flipud(fliplr(momentum_x_diff)));
shading interp;
colorbar;
title({['X-Momentum Equation Difference Re = ', num2str(Re)]});
xlabel('x');
ylabel('y');
hold on
figure(6);
subplot(1,2,2);
pcolor(x, y, flipud(fliplr(momentum_y_diff)));
shading interp;
colorbar;
title({['Y-Momentum Equation Difference Re = ', num2str(Re)]});
xlabel('x');
ylabel('y');
% print(gcf,'Momentum eqn Re1','-dpng','-r300');

sgtitle('Momentum Equation Difference');

figure(7);
pcolor(y, x, flipud(fliplr(continuity_eqn)));
shading interp;
colorbar;
title({['Continuity Equation Residual  Re = ', num2str(Re)]});
xlabel('x');
ylabel('y');
% print(gcf,'Continuity Eqn Residual Re1','-dpng','-r300');
