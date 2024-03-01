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
u_implicit = u0;

% Initialize a matrix to store wave profiles at different time instances
wave_profiles_implicit = zeros(Nx, Nt);

% Time integration using implicit upwind scheme
for n = 1:Nt
    % Construct coefficient matrix and right-hand side vector
    A = zeros(Nx, Nx);
    B = u_implicit;
    for i = 2:Nx-1
        if c > 0
            A(i, i) = 1 + (c * dt / dx);
            A(i, i-1) = -(c * dt / dx);
            B(i) = u_implicit(i);
        else
            A(i, i) = 1 - (c * dt / dx);
            A(i, i+1) = (c * dt / dx);
            B(i) = u_implicit(i);
        end
    end

    % Boundary conditions
    A(1, 1) = 1;        % Dirichlet boundary condition at left boundary
    A(Nx, Nx) = 1;      % Dirichlet boundary condition at right boundary
    B(1) = u0(1);       % Dirichlet boundary condition at left boundary
    B(Nx) = u0(Nx);     % Dirichlet boundary condition at right boundary
    
    % Solve the linear system
    u_implicit = A \ B;
    
    % Store the current wave profile
    wave_profiles_implicit(:, n) = u_implicit;
end

% Plot the results
plot(x, wave_profiles_implicit);
xlabel('Spatial coordinate');
ylabel('Scalar value');
title('Wave Profiles at Various Time Instances for Implicit Upwind');


-----------------------new code-------------------------
----------------------------------------------------------------

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
u_implicit = u0;

% Initialize a matrix to store wave profiles at different time instances
wave_profiles_implicit = zeros(Nx, Nt);

% Time integration using implicit upwind scheme
for n = 1:Nt
    % Construct coefficient matrix A for implicit upwind scheme
    A = eye(Nx); % Start with an identity matrix for diagonal
    for i = 2:Nx
        A(i, i) = 1 + CFL; % Main diagonal
        A(i, i-1) = -CFL;  % Sub-diagonal for backward difference
    end
    
    % Boundary conditions
    % Dirichlet boundary condition at left boundary
    A(1, 1) = 1;
    % Dirichlet boundary condition at the right boundary is naturally handled by the loop
    
    % Right-hand side vector
    B = u_implicit;
    B(1) = u0(1); % Applying Dirichlet condition at the left boundary for all timesteps
    B(Nx) = u0(Nx); % The right boundary condition is implicitly handled

    % Solve the linear system
    u_implicit = A \ B;
    
    % Store the current wave profile
    wave_profiles_implicit(:, n) = u_implicit;
end

% Plot the final wave profile
surf(linspace(0, Nt*dt, Nt), x, wave_profiles_implicit, 'EdgeColor', 'none');
xlabel('Time');
ylabel('Spatial coordinate');
zlabel('Scalar value');
title('Wave Profiles at Various Time Instances for Implicit Upwind');
view(2);
colorbar;

