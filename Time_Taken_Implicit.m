% Parameters
L = 1;                  % Length of the domain
Nx = 50;                % Number of spatial points
dx = L / Nx;            % Spatial step size
CFL = 0.8;              % Adjusted CFL value for stability
c = 1;                  % Velocity
dt = CFL * dx / abs(c); % Time step size
Nt = 100;               % Number of time steps
x = linspace(0, L, Nx); % Spatial grid

% Initial condition (assuming a cosine wave)
u0 = cos(2 * pi * x);
u_implicit = u0;

% Initialize time variable to store the time when wave leaves the domain
time_wave_leaves_implicit = 0;

% Time integration using implicit upwind scheme with iterative solver
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
    
    % Solve the linear system using BiCGStab iterative solver
    [u_implicit, ~, ~, ~] = bicgstab(A, B);
    
    % Check if the wave has reached the right boundary
    if u_implicit(end) <= 0.5  % Assuming wavefront is at 0.5
        time_wave_leaves_implicit = n * dt;
        break;  % Exit the loop once wave leaves the domain
    end
end

% Display the time taken for the wave to leave the domain
fprintf('Time taken for the wave to leave the domain (implicit): %.2f\n', time_wave_leaves_implicit);
