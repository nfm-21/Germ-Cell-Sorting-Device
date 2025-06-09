% Parameters
h = 0.62e-3; % Channel height (m)
w = 1e-3; % Channel width (m)
mu = 0.001; % Dynamic viscosity (Pa.s)
Q_in = 3e-9; % Volumetric flow rate (m^3/s)
L = 5e-2; % Length of the channel (m)

% Calculate DeltaP from Q
sum_term = 0;
for n = 1:2:99 % Sum over odd n
    sum_term = sum_term + (1/n^5) * (192/pi^5) * (h/w) * tanh(n*pi*w/(2*h));
end
DeltaP = (12 * mu * L * Q_in) / (h^3 * w * (1 - sum_term));

% Grid points
y = linspace(-w/2, w/2, 100);
z = linspace(0, h, 100);

% Initialize shear stress matrices
tau_y = zeros(length(y), length(z)); % miu*du/dy
tau_z = zeros(length(y), length(z)); % miu*du/dz
tau_total = zeros(length(y), length(z)); % Total shear stress

% Calculate shear stress in y and z directions
for i = 1:length(y)
    for j = 1:length(z)
        sum_tau_y = 0;
        sum_tau_z = 0;
        for n = 1:2:99 % Sum over odd n
            sum_tau_y = sum_tau_y + (1/n^2) * (sinh(n*pi*y(i)/h) / cosh(n*pi*w/(2*h))) * sin(n*pi*z(j)/h);
            sum_tau_z = sum_tau_z + (1/n^2) * (1 - cosh(n*pi*y(i)/h) / cosh(n*pi*w/(2*h))) * cos(n*pi*z(j)/h);
        end
        tau_y(i,j) = (4*h/pi^2) * (DeltaP/L) * abs(sum_tau_y);
        tau_z(i,j) = (4*h/pi^2) * (DeltaP/L) * abs(sum_tau_z);
        tau_total(i,j) = sqrt(tau_y(i,j)^2 + tau_z(i,j)^2); % Total shear stress
    end
end
% Convert coordinates to micrometers for plotting
y_um = y * 1e6; % Convert from meters to µm
z_um = z * 1e6; % Convert from meters to µm

% Plot total shear stress as a contour plot
figure;
contourf(y_um, z_um, tau_total', 20, 'LineColor', 'none');
xlabel('Width (µm)');
ylabel('Height (µm)');
% Set y-ticks at 200 µm intervals up to channel height
ytick_positions = 0:200:round(h*1e6);
yticks(ytick_positions);
yticklabels(string(ytick_positions));

colorbar;
ylabel(colorbar, 'Shear Stress (Pa)');

% Adjust x-ticks (in micrometers)
% Set x-ticks at 100 µm intervals
tick_positions = -500:100:500;
xticks(tick_positions);

% Create labels: empty for most, values for ends
xtick_labels = repmat({''}, size(tick_positions));
xtick_labels(1) = {'-500'};
xtick_labels(end) = {'500'};

% Apply the custom labels
xticklabels(xtick_labels);
colormap(pmkmp(128));
% caxis([0, 0.07]);
axis equal

% Set font to Times New Roman and size 15 for axes ticks
set(gca, 'FontName', 'Aptos', 'FontSize', 15);

% Set font for colorbar ticks and label
cb = colorbar;
set(cb, 'FontName', 'Aptos', 'FontSize', 15);
ylabel(cb, 'Shear Stress (Pa)', 'FontName', 'Aptos', 'FontSize', 15);
