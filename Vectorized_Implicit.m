% Parameters
L = 1;                  % Length of the domain
Nx = 50;                % Number of spatial points
dx = L/Nx;              % Spatial step size
CFL = 0.5;              % Courant number
c = 1;
dt = CFL*dx/abs(c);     % Time step size, c is the velocity
Nt = 100;               % Number of time steps
x = linspace(-1, L, Nx); % Spatial grid

% Initial condition (assuming a cosine wave)
u0 = cos(2*pi*x);
u_implicit = u0;

% Time integration using implicit upwind scheme
for n = 1:Nt
    % Construct coefficient matrix
    A = zeros(Nx, Nx);
    A(2:Nx, 2:Nx) = eye(Nx-1) - (c*dt/dx)*diag(ones(Nx-2, 1), 1) + (c*dt/dx)*diag(ones(Nx-2, 1), -1);
    
    % Boundary conditions
    A(1,1) = 1;
    A(1,2) = 0;  % Neumann boundary condition at left boundary
    A(Nx,Nx) = 1;
    A(Nx,Nx-1) = 0;  % Neumann boundary condition at right boundary
    
    % Right-hand side adjustments for boundary conditions
    B = u_implicit;
    B(1) = u0(1);  % Dirichlet boundary condition at left boundary
    B(Nx) = u0(Nx);  % Dirichlet boundary condition at right boundary
    
    % Solve the linear system
    u_implicit = A\B;
end

% Plot the results
plot(x, u_implicit);
xlabel('Spatial coordinate');
ylabel('Scalar value');
title('Implicit Upwind');

---new modified code-----------------------
----------------------------------------


% Parameters
L = 1;                  % Length of the domain
Nx = 50;                % Number of spatial points
dx = L/Nx;              % Spatial step size
CFL = 0.5;              % Courant number
c = 1;                  % Velocity
dt = CFL*dx/abs(c);     % Time step size
Nt = 100;               % Number of time steps
x = linspace(0, L, Nx); % Spatial grid, corrected to start from 0

% Initial condition (assuming a cosine wave)
u0 = cos(2*pi*x);
u_implicit = u0;

% Time integration using implicit upwind scheme
for n = 1:Nt
    % Construct coefficient matrix A for implicit upwind scheme
    A = eye(Nx); % Start with an identity matrix
    for i = 2:Nx
        A(i, i) = 1 + CFL;    % Main diagonal
        A(i, i-1) = -CFL;     % Sub-diagonal for backward difference
    end
    
    % Boundary conditions
    % Assuming Dirichlet condition at both ends for simplicity
    % Left boundary (first row of A is already [1, 0, 0, ..., 0])
    % Right boundary (last row of A ensures u[Nx] = u0[Nx])
    A(Nx, Nx-1) = 0;  % Ensuring no flux from outside the domain
    
    % Right-hand side vector
    B = u_implicit; % In implicit schemes, the current timestep is used as the RHS
    B(1) = u0(1);    % Applying Dirichlet condition at the left boundary for all timesteps
    B(Nx) = u0(Nx);  % Applying Dirichlet condition at the right boundary for all timesteps
    
    % Solve the linear system
    u_implicit = A\B;
end

% Plot the results
plot(x, u_implicit, 'LineWidth', 2);
xlabel('Spatial coordinate');
ylabel('Scalar value');
title('Solution using Implicit Upwind Scheme');
grid on;


