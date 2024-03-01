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



-----------new code--------------------
-------------------------------------

% Parameters
L = 1;                  % Length of the domain
Nx = 50;                % Number of spatial points
dx = L/Nx;              % Spatial step size
c = 1;                  % Velocity
CFL = 1;                % Courant number
dt = CFL*dx/abs(c);     % Time step size
Nt = 100;               % Number of time steps
x = linspace(0, L, Nx); % Spatial grid

% Initial condition (assuming a cosine wave)
u0 = cos(2*pi*x);
u_explicit = u0;

% Initialize a matrix to store wave profiles at different time instances
wave_profiles_explicit = zeros(Nx, Nt);
wave_profiles_explicit(:, 1) = u0;

% Initialize time variable to store the time when wave leaves the domain
time_wave_leaves = 0;

% Time integration using explicit upwind scheme
for n = 1:Nt
    u_next = u_explicit; % Temporary array for the new time step
    u_next(2:Nx) = u_explicit(2:Nx) - (c*dt/dx)*(u_explicit(2:Nx) - u_explicit(1:Nx-1));
    
    % Update the wave profiles storage
    wave_profiles_explicit(:, n) = u_next;
    
    % Update the current profile for the next iteration
    u_explicit = u_next;
    
    % Check if the wave has reached the right boundary
    % Here, using a more sophisticated criterion based on the wave's maximum amplitude
    if max(abs(u_explicit(Nx-5:Nx))) < 0.05 % Checking the last few points for a significant decrease
        time_wave_leaves = n*dt;
        break; % Exit the loop once the amplitude significantly decreases
    end
end

% Display the time taken for the wave to leave the domain
if time_wave_leaves > 0
    fprintf('Time taken for the wave to leave the domain: %.2f seconds\n', time_wave_leaves);
else
    fprintf('The wave has not left the domain within the given time frames.\n');
end

% Plotting the wave profiles at various time instances
surf(x, dt*(1:Nt), wave_profiles_explicit', 'EdgeColor', 'none');
xlabel('Spatial coordinate');
ylabel('Time');
zlabel('Wave Amplitude');
title('Wave Propagation Over Time');
view([0 90]); % View from top
colorbar;


% Display the time taken for the wave to leave the domain
fprintf('Time taken for the wave to leave the domain: %.2f\n', time_wave_leaves);

