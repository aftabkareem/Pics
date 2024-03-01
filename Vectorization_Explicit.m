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

-----------new code----------------------------------------
----------------------------------------------------------------------

% Parameters
L = 1;                  % Length of the domain
Nx = 50;                % Number of spatial points
dx = L / Nx;            % Spatial step size
c = 1;                  % Velocity
CFL = 1;                % Courant number
dt = CFL * dx / abs(c); % Time step size
Nt = 100;               % Number of time steps
x = linspace(0, L, Nx); % Spatial grid

% Initial condition (assuming a cosine wave)
u0 = cos(2 * pi * x);
u_explicit = u0;

% Initialize a matrix to store wave profiles at different time instances
wave_profiles_explicit = zeros(Nx, Nt);
wave_profiles_explicit(:, 1) = u0; % Store the initial condition

% Time integration using explicit upwind scheme
for n = 1:Nt
    u_next = zeros(size(u_explicit)); % Initialize the next time step
    % Update the solution using the explicit upwind scheme
    for i = 2:(Nx - 1) % Avoiding the boundaries
        u_next(i) = u_explicit(i) - (c * dt / dx) * (u_explicit(i) - u_explicit(i - 1));
    end
    
    % Apply boundary conditions
    u_next(1) = u_explicit(1); % Assuming a fixed value at the left boundary
    % For a more realistic model, consider implementing open or periodic boundary conditions
    
    % Update the current wave profile
    u_explicit = u_next;
    
    % Store the current wave profile for visualization
    wave_profiles_explicit(:, n) = u_explicit;
end

% Plotting the final wave profile
figure;
surf(x, dt * (0:Nt-1), wave_profiles_explicit, 'EdgeColor', 'none');
xlabel('Spatial coordinate');
ylabel('Time');
zlabel('Wave amplitude');
title('Wave Propagation Over Time');
view(2); % Adjust view to top-down for better visualization
colorbar;
