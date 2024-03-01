% Parameters
L = 1;                  % Length of the domain
Nx = 50;                % Number of spatial points
dx = L/Nx;              % Spatial step size
c = 1;
CFL = 1;                % Courant number
dt = CFL*dx/abs(c);     % Time step size, c is the velocity
Nt = 100;               % Number of time steps
x = linspace(0, L, Nx); % Spatial grid

% Initial condition (assuming a cosine wave)
u0 = cos(2*pi*x);
u_explicit = u0;

% Initialize time variable to store the time when wave leaves the domain
time_wave_leaves = 0;

% Time integration using explicit upwind scheme
for n = 1:Nt
    % Update the solution using the explicit upwind scheme
    u_explicit(2:Nx) = u_explicit(2:Nx) - (c*dt/dx)*(u_explicit(2:Nx) - u_explicit(1:Nx-1));
    
    % Check if the wave has reached the right boundary
    if u_explicit(end) <= 0.5  % Assuming wavefront is at 0.5
        time_wave_leaves = n*dt;
        break;  % Exit the loop once wave leaves the domain
    end
end

% Display the time taken for the wave to leave the domain
fprintf('Time taken for the wave to leave the domain: %.2f\n', time_wave_leaves);

