num_rings = 5;                      % number of concentric circles
radius_step = 50;                   % distance between rings
mics_per_ring = [8 12 16 24 32];    % number of mics per ring (can vary per ring)

mic_locations = [];                 % initialize
nozzle_location = [0 0 0];          % origin

for r = 1:num_rings
    radius = r * radius_step;
    num_mics = mics_per_ring(r);
    theta = linspace(0, 2*pi, num_mics + 1);  % +1 to close the circle
    theta(end) = [];  % remove duplicate at 2*pi

    x = radius * cos(theta);
    y = radius * sin(theta);
    z = zeros(size(x));  % flat in z=0 plane

    ring_locations = [x' y' z'];
    mic_locations = [mic_locations; ring_locations];
end

% Optional: visualize
figure;
plot(mic_locations(:,1), mic_locations(:,2), 'bo');
hold on;
plot(0, 0, 'r*', 'MarkerSize', 10);  % nozzle
axis equal;
grid on;
title('Microphones in Concentric Circles');
xlabel('X-axis');
ylabel('Y-axis');
