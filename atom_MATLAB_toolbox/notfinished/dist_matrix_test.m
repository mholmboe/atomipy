% Example atom struct with N atoms
N = 5; % Number of atoms
atom(1:N).x = rand(1, N) * 10; % Random x coordinates
atom(1:N).y = rand(1, N) * 10; % Random y coordinates
atom(1:N).z = rand(1, N) * 10; % Random z coordinates

% Triclinic cell parameters
a = 5; % Length of vector a
b = 6; % Length of vector b
c = 7; % Length of vector c
alpha = 90; % Angle between b and c in degrees
beta = 90; % Angle between a and c in degrees
gamma = 90; % Angle between a and b in degrees

% Convert angles from degrees to radians
alpha = deg2rad(alpha);
beta = deg2rad(beta);
gamma = deg2rad(gamma);

% Compute the transformation matrix from Cartesian to fractional coordinates
V = a * b * c * sqrt(1 - cos(alpha)^2 - cos(beta)^2 - cos(gamma)^2 + 2 * cos(alpha) * cos(beta) * cos(gamma));
T = [a, b * cos(gamma), c * cos(beta);
     0, b * sin(gamma), c * (cos(alpha) - cos(beta) * cos(gamma)) / sin(gamma);
     0, 0, V / (a * b * sin(gamma))];

% Initialize Cartesian coordinates matrix
cartesian_coords = zeros(N, 3);

% Compute Cartesian coordinates
for i = 1:N
    cartesian_coords(i, :) = [atom(i).x, atom(i).y, atom(i).z];
end

% Initialize distance matrix
distance_matrix = zeros(N, N);

% Compute pairwise distances
for i = 1:N
    for j = 1:N
        % Compute the distance vector
        delta = cartesian_coords(i, :) - cartesian_coords(j, :);

        % Apply periodic boundary conditions
        delta(1) = delta(1) - round2dec(delta(1) / a) * a;
        delta(2) = delta(2) - round2dec(delta(2) / b) * b;
        delta(3) = delta(3) - round2dec(delta(3) / c) * c;

        % Compute the distance
        distance_matrix(i, j) = norm(delta);
    end
end

% Display the result
disp('Distance Matrix:');
disp(distance_matrix);
