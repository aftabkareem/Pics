% Parameters
L = 1;                  % Length of the domain
Nx = 50;                % Number of spatial points
dx = L/Nx;              % Spatial step size
CFL = 1;              % Courant number
dt = CFL*dx/abs(c);     % Time step size, c is the velocity
Nt = 100;               % Number of time steps
x = linspace(0, L, Nx); % Spatial grid

% Initial condition (assuming a cosine wave)
u0 = cos(2*pi*x);
u_explicit = u0;

% Time integration using explicit upwind scheme
for n = 1:Nt
    % Vectorized calculation for interior points
    u_explicit(2:Nx) = u_explicit(2:Nx) - (c*dt/dx)*(u_explicit(2:Nx) - u_explicit(1:Nx-1));
end


% Plot the results
plot(x, u_explicit);
xlabel('Spatial coordinate');
ylabel('Scalar value');
title('Explicit Upwind');