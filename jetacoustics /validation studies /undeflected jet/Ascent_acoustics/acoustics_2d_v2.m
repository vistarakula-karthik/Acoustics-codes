% Parameters for the polar grid
r_vals = linspace(0, 250, 6);             % 6 concentric circles (0 to 250 mm)
theta_vals = linspace(0, 2*pi, 36);       % 36 angular positions (every 10 degrees)

% Create meshgrid in polar coordinates
[R, TH] = meshgrid(r_vals, theta_vals);

% Convert to Cartesian coordinates
X = R .* cos(TH);
Y = R .* sin(TH);
Z = zeros(size(X));                       % All mics lie in the XY plane (z = 0)

% Convert meshgrid to mic_locations list (optional: reshape into N x 3 array)
mic_locations = [X(:), Y(:), Z(:)];

% Plot mic positions
figure;
plot(X, Y, 'bo');                         % Plot all mic positions
hold on;
plot(0, 0, 'r*', 'MarkerSize', 10);       % Nozzle at origin
axis equal;
grid on;
title('Microphone Locations in Concentric Circular Grid');
xlabel('X-axis (mm)');
ylabel('Y-axis (mm)');
