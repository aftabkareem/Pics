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
u_implicit = u0;

% Time integration using implicit upwind scheme
for n = 1:Nt
    A = zeros(Nx,Nx); % Coefficient matrix
    B = zeros(Nx,1);  % Right-hand side
    for i = 2:Nx-1
        if c > 0
            A(i,i) = 1 + (c*dt/dx);
            A(i,i-1) = -(c*dt/dx);
            B(i) = u_implicit(i);
        else
            A(i,i) = 1 - (c*dt/dx);
            A(i,i+1) = (c*dt/dx);
            B(i) = u_implicit(i);
        end
    end
    % Boundary conditions
    A(1,1) = 1; B(1) = u_implicit(1);
    A(Nx,Nx) = 1; B(Nx) = u_implicit(Nx);
    
    % Solve the linear system
    u_implicit = A\B;
end

% Plot the results
plot(x, u_implicit);
xlabel('Spatial coordinate');
ylabel('Scalar value');
title('Implicit Upwind');
