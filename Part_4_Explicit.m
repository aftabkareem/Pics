% Parameters
L = 1;                  % Length of the domain
Nx = 50;                % Number of spatial points
dx = L/Nx;              % Spatial step size
CFL = 1;                % Courant number
c = 1;                  % Velocity
dt = CFL*dx/abs(c);     % Time step size
Nt = 100;               % Number of time steps
x = linspace(0, L, Nx); % Spatial grid

% Initial condition (assuming a cosine wave)
u0 = cos(2*pi*x);
u_explicit = u0;

% Initialize a matrix to store wave profiles at different time instances
wave_profiles_explicit = zeros(Nx, Nt);

% Time integration using explicit upwind scheme
for n = 1:Nt
    for i = 2:Nx
        if c > 0
            u_explicit(i) = u_explicit(i) - (c*dt/dx)*(u_explicit(i) - u_explicit(i-1));
        else
            u_explicit(i) = u_explicit(i) - (c*dt/dx)*(u_explicit(i+1) - u_explicit(i));
        end
    end

    % Store the current wave profile
    wave_profiles_explicit(:, n) = u_explicit;
end

% Plot the results
plot(x, wave_profiles_explicit);


--------------new code--------------------------------
-----------------------------------------------------------------------

% Parameters
L = 1;
Nx = 50;
dx = L/Nx;
CFL = 1;
c = 1;
dt = CFL*dx/abs(c);
Nt = 100;
x = linspace(0, L, Nx);

% Initial condition
u0 = cos(2*pi*x);
u_explicit = u0;

% Initialize a matrix to store wave profiles at different time instances
wave_profiles_explicit = zeros(Nx, Nt);
wave_profiles_explicit(:, 1) = u0; % Store initial condition

% Time integration using explicit upwind scheme
for n = 1:Nt-1 % Adjust loop to fill wave_profiles_explicit correctly
    u_next = u_explicit;
    for i = 2:Nx
        u_next(i) = u_explicit(i) - (c*dt/dx)*(u_explicit(i) - u_explicit(i-1));
    end
    u_explicit = u_next;

    % Store the current wave profile
    wave_profiles_explicit(:, n+1) = u_explicit;
end

% Plot the results using a surface plot for clarity
surf(linspace(0, Nt*dt, Nt), x, wave_profiles_explicit, 'EdgeColor', 'none');
xlabel('Time');
ylabel('Spatial coordinate');
zlabel('Scalar value');
title('Wave Profiles at Various Time Instances');
colorbar;
view(2); % Adjust the view to look directly down the 'Time' axis

xlabel('Spatial coordinate');
ylabel('Scalar value');
title('Wave Profiles at Various Time Instances');
legend(arrayfun(@(n) sprintf('t = %.2f', n*dt), 1:Nt, 'UniformOutput', false));
