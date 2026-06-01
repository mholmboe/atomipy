function atom = frac_atom_test(atom,Box_dim)

if numel(Box_dim)==9
    Cell=Box_dim2Cell(Box_dim);
elseif numel(Box_dim)==6
    Cell=Box_dim;
elseif numel(Box_dim)==3
    Cell=Box_dim;
end

a=Cell(1);
b=Cell(2);
c=Cell(3);

if numel(Cell)==3
    % Convert angles from degrees to radians
    alpha_rad = deg2rad(90);
    beta_rad = deg2rad(90);
    gamma_rad = deg2rad(90);
else
    % Convert angles from degrees to radians
    alpha_rad = deg2rad(Cell(4));
    beta_rad = deg2rad(Cell(5));
    gamma_rad = deg2rad(Cell(6));
end

% Compute the transformation matrix from Cartesian to fractional coordinates
V = a * b * c * sqrt(1 - cos(alpha_rad)^2 - cos(beta_rad)^2 - cos(gamma_rad)^2 + 2 * cos(alpha_rad) * cos(beta_rad) * cos(gamma_rad));
T = [a, b * cos(gamma_rad), c * cos(beta_rad);
    0, b * sin(gamma_rad), c * (cos(alpha_rad) - cos(beta_rad) * cos(gamma_rad)) / sin(gamma_rad);
    0, 0, V / (a * b * sin(gamma_rad))];

% Initialize fractional coordinates
fractional_coords = zeros(N, 3);

% Compute fractional coordinates
for i = 1:N
    cartesian_coords = [atom(i).x; atom(i).y; atom(i).z];
    fractional_coords(i, :) = (T \ cartesian_coords)';
end

atom

% Display the result
disp('Fractional Coordinates:');
disp(fractional_coords);
