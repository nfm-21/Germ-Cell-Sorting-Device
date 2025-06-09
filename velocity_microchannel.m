%% Velocity profile in y-z plane for rectangular microchannel

% Parameters
h = 0.62e-3;           % Channel height (m)
w = 1e-3;              % Channel width (m)
mu = 0.001;            % Dynamic viscosity (Pa.s)
Q_in = 3e-9;           % Volumetric flow rate (m^3/s)
L = 5e-2;              % Channel length (m)

% Compute pressure drop using Fourier approximation
sum_term = 0;
for n = 1:2:99 % Odd n terms
    sum_term = sum_term + (1/n^5) * (192/pi^5) * (h/w) * tanh(n*pi*w/(2*h));
end
DeltaP = (12 * mu * L * Q_in) / (h^3 * w * (1 - sum_term)); % Pressure drop

% Grid points (y across width, z along height)
y = linspace(-w/2, w/2, 100);
z = linspace(0, h, 100);

% Initialize velocity matrix
u = zeros(length(y), length(z));

% Compute velocity at each (y,z) position using Fourier series
for i = 1:length(y)
    for j = 1:length(z)
        sum_u = 0;
        for n = 1:2:99
            sum_u = sum_u + (1/n^3) * ...
                (1 - cosh(n*pi*y(i)/h)/cosh(n*pi*w/(2*h))) * ...
                sin(n*pi*z(j)/h);
        end
        u(i,j) = (4*h^2)/(pi^3 * mu) * (DeltaP/L) * sum_u;
    end
end

u = u * 1e3;  % Convert from m/s to mm/s

% Convert coordinates to micrometers for plotting
y_um = y * 1e6;
z_um = z * 1e6;

% Plot velocity profile as a contour plot
figure;
contourf(y_um, z_um, u', 20, 'LineColor', 'none');
xlabel('Width (µm)');
ylabel('Height (µm)');
colorbar;
ylabel(colorbar, 'Velocity (mm/s)');
colormap(pmkmp(128));
caxis([0, 10]);
axis equal

% Set font to Times New Roman and size 15 for axes ticks
set(gca, 'FontName', 'Aptos', 'FontSize', 15);

% Set font for colorbar ticks and label
cb = colorbar;
set(cb, 'FontName', 'Aptos', 'FontSize', 15);
ylabel(cb, 'Velocity (mm/s)', 'FontName', 'Aptos', 'FontSize', 15);

%% Calculating average velocity to validate formula

% -------------------------
% Velocity profile in y-z plane for rectangular microchannel
% -------------------------

% Parameters
h = 0.62e-3;           % Channel height (m)
w = 1e-3;              % Channel width (m)
mu = 0.001;            % Dynamic viscosity (Pa.s)
Q_in = 3e-9;           % Volumetric flow rate (m^3/s)
L = 5e-2;              % Channel length (m)

% Compute pressure drop using Fourier approximation
sum_term = 0;
for n = 1:2:99 % Odd n terms
    sum_term = sum_term + (1/n^5) * (192/pi^5) * (h/w) * tanh(n*pi*w/(2*h));
end
DeltaP = (12 * mu * L * Q_in) / (h^3 * w * (1 - sum_term)); % Pressure drop

% Grid points (y across width, z along height)
y = linspace(-w/2, w/2, 100);
z = linspace(0, h, 100);

% Initialize velocity matrix
u = zeros(length(y), length(z));

% Compute velocity at each (y,z) position using Fourier series
for i = 1:length(y)
    for j = 1:length(z)
        sum_u = 0;
        for n = 1:2:99
            sum_u = sum_u + (1/n^3) * ...
                (1 - cosh(n*pi*y(i)/h)/cosh(n*pi*w/(2*h))) * ...
                sin(n*pi*z(j)/h);
        end
        u(i,j) = (4*h^2)/(pi^3 * mu) * (DeltaP/L) * sum_u;
    end
end

% u = u * 1e3;  % Convert from m/s to mm/s

% Calculate average velocity
% Calculate the spacing between the grid points (in mm)
dy = (y(end) - y(1)) / (length(y) - 1);  % Grid spacing in y-direction (mm)
dz = (z(end) - z(1)) / (length(z) - 1);  % Grid spacing in z-direction (mm)

% Calculate the area of each grid cell (in mm²)
cell_area = dy * dz;  % Area of each grid cell (in mm²)

% Calculate the total velocity, summed over the grid
total_velocity = sum(u(:)) * cell_area;

% Area of the channel cross-section in mm²
% Convert width and height from meters to mm
area = w*h;  % Area in mm²

% Average velocity in mm/s
u_avg = total_velocity / area;  % in mm/s

% Display average velocity
fprintf('Average velocity: %.3f mm/s\n', u_avg);
