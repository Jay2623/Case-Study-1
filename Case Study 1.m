Lx = 50;
Ly = 50;
Nx = 50;
Ny = 10;
dx = Lx / (Nx - 1);
dy = Ly / (Ny - 1);
alpha = 1e-4;
T_initial = 100;
T_left = 200;
T_right = 50;
T_top = 75;
T_bottom = 125;
t_end = 100;
dt = 0.1;
Nt = ceil(t_end / dt);
x = linspace(0, Lx, Nx);
y = linspace(0, Ly, Ny);
[X, Y] = meshgrid(x, y);
T = T_initial * ones(Ny, Nx);
T(:, 1) = T_left;
T(:, end) = T_right;
T(1, :) = T_top;
T(end, :) = T_bottom;
q_left = 0;
q_right = 0;
q_top = 0;
q_bottom = 0;
T_history = zeros(Ny, Nx, Nt);
for t = 1:Nt
T_old = T;
for i = 2:Nx-1
for j = 2:Ny-1
% Neumann boundary conditions
if i == 2 % Left boundary
dT_dx = (T_old(j, i+1) - T_old(j, i)) / dx;
elseif i == Nx - 1 % Right boundary
dT_dx = (T_old(j, i) - T_old(j, i-1)) / dx;
else
dT_dx = (T_old(j, i+1) - T_old(j, i-1)) / (2*dx);
end
if j == 2 % Top boundary
dT_dy = (T_old(j+1, i) - T_old(j, i)) / dy;
elseif j == Ny - 1 % Bottom boundary
dT_dy = (T_old(j, i) - T_old(j-1, i)) / dy;
else
dT_dy = (T_old(j+1, i) - T_old(j-1, i)) / (2*dy);
end
% Update temperature using Neumann boundary conditions
T(j, i) = T_old(j, i) + alpha * dt * (...
(T_old(j, i+1) - 2*T_old(j, i) + T_old(j, i-1))/(dx^2) +...
(T_old(j+1, i) - 2*T_old(j, i) + T_old(j-1, i))/(dy^2) +...
dT_dx^2 + dT_dy^2);
end
end
T_history(:, :, t) = T;
end
figure;
for t = 1:10:Nt
contourf(X, Y, T_history(:, :, t), 'LineStyle', 'none');
colorbar;
xlabel('Distance (m)');
ylabel('Distance (m)');
title(['Temperature Distribution at t = ' num2str(t*dt) ' s']);
colormap(jet);
caxis([min(T(:)), max(T(:))]);
colorbar('Ticks', linspace(min(T(:)), max(T(:)), 5), 'TickLabels', {'50°C',
'75°C', '100°C', '125°C', '200°C'});
pause(0.1);
end