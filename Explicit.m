% Parameters
L = 1;                  % Length of the domain
Nx = 50;                % Number of spatial points
dx = L/Nx;              % Spatial step size
c=1;
CFL = 1;              % Courant number
dt = CFL*dx/abs(c);     % Time step size, c is the velocity
Nt = 100;               % Number of time steps
x = linspace(0, L, Nx); % Spatial grid

% Initial condition (assuming a cosine wave)
u0 = cos(2*pi*x);
u_explicit = u0;

% Time integration using explicit upwind scheme
for n = 1:Nt
    for i = 2:Nx
        if c > 0
            u_explicit(i) = u_explicit(i) - (c*dt/dx)*(u_explicit(i) - u_explicit(i-1));
        else
            u_explicit(i) = u_explicit(i) - (c*dt/dx)*(u_explicit(i+1) - u_explicit(i));
        end
    end
end

% Plot the results
plot(x, u_explicit);
xlabel('Spatial coordinate');
ylabel('Scalar value');
title('Explicit Upwind');


-----------------new code------------------------------------
--------------------------------------------------------------

% Parameters
L = 1;
Nx = 50;
dx = L/Nx;
c = 1;
CFL = 1;
dt = CFL*dx/abs(c);
Nt = 100;
x = linspace(0, L, Nx);

% Initial condition
u0 = cos(2*pi*x);
u_explicit = u0;
u_next = u0; % Temporary vector for the next state

% Time integration using explicit upwind scheme
for n = 1:Nt
    for i = 2:Nx
        % Apply upwind differencing for positive velocity
        u_next(i) = u_explicit(i) - (c*dt/dx)*(u_explicit(i) - u_explicit(i-1));
    end
    u_explicit = u_next; % Update the wave state
end

% Plot the results
plot(x, u_explicit);
xlabel('Spatial coordinate');
ylabel('Scalar value');
title('Explicit Upwind');
