clear; close all; clc;

%% Creation of sample synthetic data on a unit cube.
% Number of nodes in each direction.
num_x = 100;
num_y = 100;
num_z = 100;
% Step between nodes in each direction.
hx = 1 / (num_x - 1);
hy = 1 / (num_y - 1);
hz = 1 / (num_z - 1);
% Vectors representing axis ticks along each direction.
x = linspace(0, 1, num_x);
y = linspace(0, 1, num_y);
z = linspace(0, 1, num_z);


%% Set up data.
% u_x, u_y, u_z represent the x,y,z components of a vector field. For
% example, at (u_x(i,j,k), u_y(i,j,k), u_z(i,j,k)) = u(i,j,k).
R = 0.25; % radius of cylinder
u_max = 5; % maximum speed
parabola = @(r) -1 * u_max / R^2 * (r - R) * (r + R);

% % Sanity check for parabola.
% testx = linspace(-1, 1);
% testy = zeros(100);
% for i = 1 : 100
%     r = testx(i);
%     testy(i) = -1 * 0.1 / R^2 * (r - R) * (r + R);
% end
% plot(testx, testy);

% Speed in each direction.
u1 = zeros(num_x, num_y, num_z);
u2 = zeros(num_x, num_y, num_z);
u3 = zeros(num_x, num_y, num_z);
for i = linspace(1, num_x, num_x)
    for j = linspace(1, num_y, num_y)
        for k = linspace(1, num_z, num_z)
            % For a cylinder centered around x = 0.5, y = 0.5,
            % distance from center = (x - 0.5)^2 + (y - 0.5)^2
            d_sq = (x(i) - 0.5)^2 + (y(j) - 0.5)^2;
            if (d_sq == 0)
                u3(i, j, k) = u_max;
            elseif (d_sq <= R^2)
                d = sqrt(d_sq);
                u3(i, j, k) = parabola(d);
%                 disp(parabola(d))
            else
                u3(i, j, k) = 0;
            end
        end
    end
end

figure(1);
imshow3D(u3);


%% Calculate gradients
% Calculate partial derivatives
du1dx = zeros(num_x - 2, num_y - 2, num_z - 2);
du1dy = zeros(num_x - 2, num_y - 2, num_z - 2);
du1dz = zeros(num_x - 2, num_y - 2, num_z - 2);
du2dx = zeros(num_x - 2, num_y - 2, num_z - 2);
du2dy = zeros(num_x - 2, num_y - 2, num_z - 2);
du2dz = zeros(num_x - 2, num_y - 2, num_z - 2);
du3dx = zeros(num_x - 2, num_y - 2, num_z - 2);
du3dy = zeros(num_x - 2, num_y - 2, num_z - 2);
du3dz = zeros(num_x - 2, num_y - 2, num_z - 2);
for i = linspace(1, num_x - 2, num_x - 2) % NOTE: different bounds from building data.
    for j = linspace(1, num_y - 2, num_y - 2)
        for k = linspace(1, num_z - 2, num_z - 2)
            du1dx(i, j, k) = (u1(i+2, j+1, k+1) - u1(i, j+1, k+1)) / (2 * hx);
            du1dy(i, j, k) = (u1(i+1, j+2, k+1) - u1(i+1, j, k+1)) / (2 * hy);
            du1dz(i, j, k) = (u1(i+1, j+1, k+2) - u1(i+1, j+1, k)) / (2 * hz);
            du2dx(i, j, k) = (u2(i+2, j+1, k+1) - u2(i, j+1, k+1)) / (2 * hx);
            du2dy(i, j, k) = (u2(i+1, j+2, k+1) - u2(i+1, j, k+1)) / (2 * hy);
            du2dz(i, j, k) = (u2(i+1, j+1, k+2) - u2(i+1, j+1, k)) / (2 * hz);
            du3dx(i, j, k) = (u3(i+2, j+1, k+1) - u3(i, j+1, k+1)) / (2 * hx);
            du3dy(i, j, k) = (u3(i+1, j+2, k+1) - u3(i+1, j, k+1)) / (2 * hy);
            du3dz(i, j, k) = (u3(i+1, j+1, k+2) - u3(i+1, j+1, k)) / (2 * hz);
        end
    end
end

% Get magnitude of gradients.
% Each element of grad stores [du1dx, du1dy, du1dz
%                              du2dx, du2dy, du2dz
%                              du3dx, du3dy, du3dz].
gradients = cell(num_x - 2, num_y - 2, num_z - 2);
% grad_test = zeros(num_x - 2, num_y - 2, num_z - 2);
% Iterate through i, j, k (each point along the field).
for i = linspace(1, num_x - 2, num_x - 2) % NOTE: different bounds from building data.
    for j = linspace(1, num_y - 2, num_y - 2)
        for k = linspace(1, num_z - 2, num_z - 2)
            gradients{i, j, k} = zeros(3, 3);
            % Fill in indices of the 9 partial derivatives.
            gradients{i, j, k}(1, 1) = du1dx(i, j, k);
            gradients{i, j, k}(1, 2) = du1dy(i, j, k);
            gradients{i, j, k}(1, 3) = du1dz(i, j, k);
            gradients{i, j, k}(2, 1) = du2dx(i, j, k);
            gradients{i, j, k}(2, 2) = du2dy(i, j, k);
            gradients{i, j, k}(2, 3) = du2dz(i, j, k);
            gradients{i, j, k}(3, 1) = du3dx(i, j, k);
            gradients{i, j, k}(3, 2) = du3dy(i, j, k);
            gradients{i, j, k}(3, 3) = du3dz(i, j, k);
            % Sanity check
%             grad_test(i, j, k) = sqrt(du3dx(i, j, k)^2 + du3dy(i, j, k)^2 + du3dz(i, j, k)^2);
        end
    end
end
% celldisp(gradients);
% figure(6);
% imshow3D(grad_test);

%% Compute shear stress
% Compute E
mu = 1; % Viscosity.
E = cell(num_x - 2, num_y - 2, num_z - 2); % E is nxnxn, each elem storing a 3x3 -- 5-D.
for i = linspace(1, num_x - 2, num_x - 2) % NOTE: same bounds as gradient.
    for j = linspace(1, num_y - 2, num_y - 2)
        for k = linspace(1, num_z - 2, num_z - 2)
            E{i, j, k} = mu * gradients{i, j, k} + (gradients{i, j, k})'; % 3x3.
        end
    end
end

% Compute normal vector at each point.
% NOTE: the cylinder is along the z-axis
%   --> axis of rotation for cylinder is x = 0.5, y = 0.5.
% At each element, store the unit normal vector in a column vector.
normals = cell(num_x - 2, num_y - 2, num_z - 2);
for i = linspace(1, num_x - 2, num_x - 2)
    for j = linspace(1, num_y - 2, num_y - 2)
        for k = linspace(1, num_z - 2, num_z - 2)
            vector = [x(i) - 0.5; y(j) - 0.5; 0];
            normals{i, j, k} = vector/norm(vector);
        end
    end
end
% celldisp(normals);

% Compute w.
w = cell(num_x - 2, num_y - 2, num_z - 2);
for i = linspace(1, num_x - 2, num_x - 2)
    for j = linspace(1, num_y - 2, num_y - 2)
        for k = linspace(1, num_z - 2, num_z - 2)
            w{i, j, k} = E{i, j, k} * normals{i, j, k};
        end
    end
end

% Compute magnitude of w
mag_w = zeros(num_x - 2, num_y - 2, num_z - 2);
for i = linspace(1, num_x - 2, num_x - 2)
    for j = linspace(1, num_y - 2, num_y - 2)
        for k = linspace(1, num_z - 2, num_z - 2)
            mag_w(i, j, k) = norm(w{i, j, k});
%             rdotw(i, j, k) = dot(normals{i, j, k}, w{i, j, k});
        end
    end
end






%% Plot the values of velocity field and magnitude of w at one line across the center

plot_x = x(2 : size(x, 2) - 1);
plot_velocity = zeros(size(x, 2));
plot_w = zeros(size(x, 2));
for i = linspace(2, num_x - 1, num_x - 2)
    plot_velocity(i-1) = u3(i, 50, 50);
    plot_w(i-1) = norm(w{i-1, 49, 49});
end
figure(9);
plot(x, plot_velocity);
hold on
plot(x, plot_w);
hold off

% Check if the other components of shear stress are zero in this example...

% This particular test is for shear stress in the direction of parallel to
% the cylinder axis.

% % Compute unit vectors for this.
% intos = cell(num_x - 2, num_y - 2, num_z - 2); % "into" as in "into the page"
% for i = linspace(1, num_x - 2, num_x - 2)
%     for j = linspace(1, num_y - 2, num_y - 2)
%         for k = linspace(1, num_z - 2, num_z - 2)
%             vector = [x(i); y(j) - 0.5; z(k+1) - z(k)];
%             intos{i, j, k} = vector/norm(vector);
%         end
%     end
% end
%
% % Compute r dot w
% qdotw = zeros(num_x - 2, num_y - 2, num_z - 2);
% for i = linspace(1, num_x - 2, num_x - 2)
%     for j = linspace(1, num_y - 2, num_y - 2)
%         for k = linspace(1, num_z - 2, num_z - 2)
%             qdotw(i, j, k) = dot(intos{i, j, k}, w{i, j, k});
%         end
%     end
% end
%
% figure(8);
% imshow3D(qdotw);


% % This test is for tangential vectors.
% % Compute unit vectors for this.
% tangents = cell(num_x - 2, num_y - 2, num_z - 2); % "into" as in "into the page"
% for i = linspace(1, num_x - 2, num_x - 2)
%     for j = linspace(1, num_y - 2, num_y - 2)
%         for k = linspace(1, num_z - 2, num_z - 2)
%             vector = [x(i); y(j) - 0.5; z(k+1) - z(k)];
%             tangents{i, j, k} = vector/norm(vector);
%         end
%     end
% end
%
% % Compute w.
% w = cell(num_x - 2, num_y - 2, num_z - 2);
% for i = linspace(1, num_x - 2, num_x - 2)
%     for j = linspace(1, num_y - 2, num_y - 2)
%         for k = linspace(1, num_z - 2, num_z - 2)
%             w{i, j, k} = E{i, j, k} * tangents{i, j, k};
%         end
%     end
% end
%
% % Compute r dot w
% pdotw = zeros(num_x - 2, num_y - 2, num_z - 2);
% for i = linspace(1, num_x - 2, num_x - 2)
%     for j = linspace(1, num_y - 2, num_y - 2)
%         for k = linspace(1, num_z - 2, num_z - 2)
%             pdotw(i, j, k) = dot(tangents{i, j, k}, w{i, j, k});
%         end
%     end
% end
%
% figure(8);
% imshow3D(pdotw);
% % Compute r dot w
% qdotw = zeros(num_x - 2, num_y - 2, num_z - 2);
% for i = linspace(1, num_x - 2, num_x - 2)
%     for j = linspace(1, num_y - 2, num_y - 2)
%         for k = linspace(1, num_z - 2, num_z - 2)
%             qdotw(i, j, k) = dot(intos{i, j, k}, w{i, j, k});
%         end
%     end
% end
%
% figure(8);
% imshow3D(qdotw);
% disp(rdotw);

figure(7);
imshow3D(mag_w);


%% Read in a .stl file and operate on it.
% Read in the file.
filename = "test_cylinder3.stl";
radius = 0.23;
cylinder = stlread(filename);
points = cylinder.Points;
connectivity = cylinder.ConnectivityList;

% Transform the file to be what we want the cylinder to look like.
x_min = min(points(:, 1));
x_max = max(points(:, 1));
y_min = min(points(:, 2));
y_max = max(points(:, 2));
z_min = min(points(:, 3));
z_max = max(points(:, 3));
z_lower_lim = z(2);
z_upper_lim = z(num_z - 1) - 0.0001;
for idx = 1 : size(points, 1)
    points(idx, 1) = radius * 2 * ((points(idx, 1) - x_min) / (x_max - x_min) - 0.5) + 0.5;
    points(idx, 2) = radius * 2 * ((points(idx, 2) - y_min) / (y_max - y_min) - 0.5) + 0.5;
    points(idx, 3) = (z_upper_lim - z_lower_lim) * ((points(idx, 3) - z_min) / (z_max - z_min) - 0.5) + 0.5;
end
TR = triangulation(connectivity, points);
stlwrite(TR, "test_cylinder_stl.stl", 'text')

% Determine which voxel a point on the cylinder lies and give it the value
% of w at that voxel.
% At each point (x,y,z), find maximum (i,j,k) less than (x,y,z).
% NOTE: We don't have access to the largest and smallest values of x,y,z.
% We only have access to values between x(2) and x(num_nodes-1) for a total
% of num_nodes-2 voxels.
w_cylinder_u3_component = zeros(size(points, 1), 1); % Stores the scalar w values associated with each point on the cylinder.
for point = 1 : size(points, 1)
    x0 = points(point, 1);
    y0 = points(point, 2);
    z0 = points(point, 3);
    x_idx = binary_search_voxel(x0, x) - 1; % We subtract 1 because these are indices in w, which is smaller than the original.
    y_idx = binary_search_voxel(y0, y) - 1;
    z_idx = binary_search_voxel(z0, z) - 1;
    % Beginning of special cases.
    if (x0 == x(x_idx))
        if (y0 == y(y_idx))
            if (z0 == z(z_idx))
                % Most special case: the corner case.
                w_cylinder_u3_component(point) = mag_w(x_idx, y_idx, z_idx);
            else
                % On a x-y edge
                w_cylinder_u3_component(point) = 1/2 * (mag_w(x_idx, y_idx, z_idx) + ...
                    mag_w(x_idx, y_idx, z_idx+1));
            end
        else
            if (z0 == z(z_idx))
                % On a x-z edge
                w_cylinder_u3_component(point) = 1/2 * (mag_w(x_idx, y_idx, z_idx) + ...
                    mag_w(x_idx, y_idx+1, z_idx));
            else
                % On a x face
                w_cylinder_u3_component(point) = 1/4 * (mag_w(x_idx, y_idx, z_idx) + ...
                    mag_w(x_idx, y_idx+1, z_idx) + ...
                    mag_w(x_idx, y_idx, z_idx+1) + ...
                    mag_w(x_idx, y_idx+1, z_idx+1));
            end
        end
    else
        if (y0 == y(y_idx))
            if (z0 == z(z_idx))
                % On a y-z edge
                w_cylinder_u3_component(point) = 1/2 * (mag_w(x_idx, y_idx, z_idx) + ...
                    mag_w(x_idx+1, y_idx, z_idx));
            else
                % On a y face
                w_cylinder_u3_component(point) = 1/4 * (mag_w(x_idx, y_idx, z_idx) + ...
                    mag_w(x_idx+1, y_idx, z_idx) + ...
                    mag_w(x_idx, y_idx, z_idx+1) + ...
                    mag_w(x_idx+1, y_idx, z_idx+1));
            end
        else
            if (z0 == z(z_idx))
                % On a z face
                w_cylinder_u3_component(point) = 1/4 * (mag_w(x_idx, y_idx, z_idx) + ...
                    mag_w(x_idx+1, y_idx, z_idx) + ...
                    mag_w(x_idx, y_idx+1, z_idx) + ...
                    mag_w(x_idx+1, y_idx+1, z_idx));
            else
                % End of special cases, only generic case left: average all 8 points around this point.
                w_cylinder_u3_component(point) = 1/8 * (mag_w(x_idx, y_idx, z_idx) + ...
                    mag_w(x_idx+1, y_idx, z_idx) + ...
                    mag_w(x_idx, y_idx+1, z_idx) + ...
                    mag_w(x_idx, y_idx, z_idx+1) + ...
                    mag_w(x_idx+1, y_idx+1, z_idx) + ...
                    mag_w(x_idx+1, y_idx, z_idx+1) + ...
                    mag_w(x_idx, y_idx+1, z_idx+1) + ...
                    mag_w(x_idx+1, y_idx+1, z_idx+1));
            end
        end
    end
end


%% In general, compute an orthonormal frame at each point.
% Compute a normal vector.
radial = zeros(size(points, 1), 3); % The matrix of radial vectors at each point.
% Iterate through the triangles and store the normal vector at each of its
% triangles to find the area-weighted normal vector.
r_triangle_aggregate = zeros(size(points, 1), 4); % (sum_x, sum_y, sum_z, num)
F = faceNormal(TR);
for triangle = 1 : size(connectivity, 1)
    % Compute area of the triangle for weighting purposes.
    point1 = points(connectivity(triangle, 1), 1:3);
    point2 = points(connectivity(triangle, 2), 1:3);
    point3 = points(connectivity(triangle, 3), 1:3);
    area = area_of_triangle(point1, point2, point3);
    % Weighting process.
    for point = 1 : 3
        r_triangle_aggregate(connectivity(triangle, point), 1:3) = r_triangle_aggregate(connectivity(triangle, point), 1:3) + F(triangle, :) * area;
        r_triangle_aggregate(connectivity(triangle, point), 4) = r_triangle_aggregate(connectivity(triangle, point), 4) + area;
    end
end
% Scale it all down to find an average.
for i = 1 : size(radial, 1)
    radial(i, :) = r_triangle_aggregate(i, 1:3) / r_triangle_aggregate(i, 4);
end

% TODO: remove the following-----------------------------------------------
% Synthetic centerline + centerline collar data.
centerline_points = [0.5 0.5 0;
    0.5 0.5 0.1;
    0.5 0.5 0.125;
    0.5 0.5 0.2;
    0.5 0.5 0.312;
    0.5 0.5 0.371;
    0.5 0.5 0.4;
    0.5 0.5 0.5;
    0.5 0.5 0.65;
    0.5 0.5 0.75;
    0.5 0.5 0.9;
    0.5 0.5 0.98;
    0.5 0.5 1];
% collars
collars = cell(size(centerline_points, 1), 1);
for clp = 1 : size(centerline_points, 1)
    for point = 1 : size(points, 1)
        if (abs(points(point, 3) - centerline_points(clp, 3)) < 0.1)
            collars{clp} = cat(1, collars{clp}, point);
        end
    end
end
% unit tangent vectors of centerline points
unit_tangents = zeros(size(centerline_points, 1), 3);
for point = 1 : size(unit_tangents, 1)
    unit_tangents(point, :) = [0, 0, 1];
end
%--------------------------------------------------------------------------


% Compute a tangential (aka axial) vector.
% Reverse mapping of collars centerline points to the collar.
collars_rev = cell(size(points, 1), 1);
for clp = 1 : size(centerline_points, 1)
    for i = 1 : size(collars{clp}, 1)
        collar_point = collars{clp}(i);
        collars_rev{collar_point} = cat(1, collars_rev{collar_point}, clp);
    end
end

% Create an orthonormal frame for each point on the cylinder.
axial = zeros(size(points, 1), 3); % Stores axial vector in [x,y,z]
axial_aggregate = zeros(size(points, 1), 4); % Stores [x,y,z,#clps]
for point = 1 : size(points, 1)
    clps = collars_rev{point}; % Corresponding centerline-points.
    % Sum the centerline points.
     for clp = 1 : size(clps, 1)
        axial_aggregate(point, 1:3) = axial_aggregate(point, 1:3) + unit_tangents(clp, 1:3);
        axial_aggregate(point, 4) = axial_aggregate(point, 4) + 1;
    end
end

% Take average of these centerline-points.
for point = 1 : size(axial, 1)
    axial(point, :) = axial_aggregate(point, 1:3) / axial_aggregate(point, 4);
end

% Get the circumferential vector by taking cross product of radial x axial.
circ = zeros(size(points, 1), 3);
for point = 1 : size(points, 1)
    circ(point, :) = cross(radial(point, :), axial(point, :));
end

% SANITY CHECK THAT DOT PRODS ARE ROUGHLY 0
radial_axial = zeros(size(points, 1), 1);
problematics = [];
radial_circ = zeros(size(points, 1), 1);
axial_circ = zeros(size(points, 1), 1);
for i = 1 : size(points, 1)
    if (dot(radial(i, :), axial(i, :)) < 10^(-5))
        radial_axial(i) = 0;
        radial_axial(i) = dot(radial(i, :), axial(i, :));
    else
        radial_axial(i) = dot(radial(i, :), axial(i, :));
        problematics = [problematics; points(i, :)];
    end
%     if dot(radial(i, :), circ(i, :)) < 10^(-5)
%         radial_circ(i) = 0;
%     else
        radial_circ(i) = dot(radial(i, :), circ(i, :));
%     end
%     if dot(axial(i, :), circ(i, :)) < 10^(-5)
%         axial_circ(i) = 0;
%     else
        axial_circ(i) = dot(axial(i, :), circ(i, :));
%     end
end

% Normallize all the frame vectors.
for i = 1 : size(points, 1)
    axial(i, :) = axial(i, :) / norm(axial(i, :));
    radial(i, :) = radial(i, :) / norm(radial(i, :));
    if sum(circ(i, :)) ~= 0
        circ(i, :) = circ(i, :) / norm(circ(i, :));
    else
        circ(i, :) = [0, 0, -1.5];
    end
end


% Visuallize the components of shear stress at each point.
% This code is the same as the code from computing mag_w.
adotw = zeros(size(points, 1), 1); % The axial component of w; at each point, is a scalar.
rdotw = zeros(size(points, 1), 1); % The radial component of w; scalars.
cdotw = zeros(size(points, 1), 1); % The circumferential component of w; scalars.
all_w = zeros(size(points, 1), 3); % The entirety of w; at each point, is a vector.
for point = 1 : size(points, 1)
    x0 = points(point, 1);
    y0 = points(point, 2);
    z0 = points(point, 3);
    x_idx = binary_search_voxel(x0, x) - 1; % We subtract 1 because these are indices in w, which is smaller than the original.
    y_idx = binary_search_voxel(y0, y) - 1;
    z_idx = binary_search_voxel(z0, z) - 1;
    % Beginning of special cases.
    if (x0 == x(x_idx))
        if (y0 == y(y_idx))
            if (z0 == z(z_idx))
                % Most special case: the corner case.
                adotw(point) = dot(axial(point, :), w{x_idx, y_idx, z_idx});
                rdotw(point) = dot(radial(point, :), w{x_idx, y_idx, z_idx});
                cdotw(point) = dot(circ(point, :), w{x_idx, y_idx, z_idx});
                all_w(point) = w{x_idx, y_idx, z_idx};
            else
                % On a x-y edge
                adotw(point) = dot(axial(point, :), ...
                    1/2 * (w{x_idx, y_idx, z_idx} + ...
                    w{x_idx, y_idx, z_idx+1}));
                rdotw(point) = dot(radial(point, :), ...
                    1/2 * (w{x_idx, y_idx, z_idx} + ...
                    w{x_idx, y_idx, z_idx+1}));
                cdotw(point) = dot(circ(point, :), ...
                    1/2 * (w{x_idx, y_idx, z_idx} + ...
                    w{x_idx, y_idx, z_idx+1}));
                all_w(point) = 1/2 * (w{x_idx, y_idx, z_idx} + ...
                    x{s_idx, y_idx, z_idx+1});
            end
        else
            if (z0 == z(z_idx))
                % On a x-z edge
                adotw(point) = dot(axial(point, :), ...
                    1/2 * (w{x_idx, y_idx, z_idx} + ...
                    w{x_idx, y_idx+1, z_idx}));
                rdotw(point) = dot(radial(point, :), ...
                    1/2 * (w{x_idx, y_idx, z_idx} + ...
                    w{x_idx, y_idx+1, z_idx}));
                cdotw(point) = dot(circ(point, :), ...
                    1/2 * (w{x_idx, y_idx, z_idx} + ...
                    w{x_idx, y_idx+1, z_idx}));
                all_w(point) = 1/2 * (w{x_idx, y_idx, z_idx} + ...
                    w{x_idx, y_idx+1, z_idx});
            else
                % On a x face
                adotw(point) = dot(axial(point, :), ...
                    1/4 * (w{x_idx, y_idx, z_idx} + ...
                    w{x_idx, y_idx+1, z_idx} + ...
                    w{x_idx, y_idx, z_idx+1} + ...
                    w{x_idx, y_idx+1, z_idx+1}));
                rdotw(point) = dot(radial(point, :), ...
                    1/4 * (w{x_idx, y_idx, z_idx} + ...
                    w{x_idx, y_idx+1, z_idx} + ...
                    w{x_idx, y_idx, z_idx+1} + ...
                    w{x_idx, y_idx+1, z_idx+1}));
                cdotw(point) = dot(circ(point, :), ...
                    1/4 * (w{x_idx, y_idx, z_idx} + ...
                    w{x_idx, y_idx+1, z_idx} + ...
                    w{x_idx, y_idx, z_idx+1} + ...
                    w{x_idx, y_idx+1, z_idx+1}));
                all_w(point) = 1/4 * (w{x_idx, y_idx, z_idx} + ...
                    w{x_idx, y_idx+1, z_idx} + ...
                    w{x_idx, y_idx, z_idx+1} + ...
                    w{x_idx, y_idx+1, z_idx+1});
            end
        end
    else
        if (y0 == y(y_idx))
            if (z0 == z(z_idx))
                % On a y-z edge
                adotw(point) = dot(axial(point, :), ...
                    1/2 * (w{x_idx, y_idx, z_idx} + ...
                    w{x_idx+1, y_idx, z_idx}));
                rdotw(point) = dot(radial(point, :), ...
                    1/2 * (w{x_idx, y_idx, z_idx} + ...
                    w{x_idx+1, y_idx, z_idx}));
                cdotw(point) = dot(circ(point, :), ...
                    1/2 * (w{x_idx, y_idx, z_idx} + ...
                    w{x_idx+1, y_idx, z_idx}));
                all_w(point) = 1/2 * (w{x_idx, y_idx, z_idx} + ...
                    w{x_idx+1, y_idx, z_idx});
            else
                % On a y face
                adotw(point) = dot(axial(point, :), ...
                    1/4 * (w{x_idx, y_idx, z_idx} + ...
                    w{x_idx+1, y_idx, z_idx} + ...
                    w{x_idx, y_idx, z_idx+1} + ...
                    w{x_idx+1, y_idx, z_idx+1}));
                rdotw(point) = dot(radial(point, :), ...
                    1/4 * (w{x_idx, y_idx, z_idx} + ...
                    w{x_idx+1, y_idx, z_idx} + ...
                    w{x_idx, y_idx, z_idx+1} + ...
                    w{x_idx+1, y_idx, z_idx+1}));
                cdotw(point) = dot(circ(point, :), ...
                    1/4 * (w{x_idx, y_idx, z_idx} + ...
                    w{x_idx+1, y_idx, z_idx} + ...
                    w{x_idx, y_idx, z_idx+1} + ...
                    w{x_idx+1, y_idx, z_idx+1}));
                all_w(point) = 1/4 * (w{x_idx, y_idx, z_idx} + ...
                    w{x_idx+1, y_idx, z_idx} + ...
                    w{x_idx, y_idx, z_idx+1} + ...
                    w{x_idx+1, y_idx, z_idx+1});
            end
        else
            if (z0 == z(z_idx))
                % On a z face
                adotw(point) = dot(axial(point, :), ...
                    1/4 * (w{x_idx, y_idx, z_idx} + ...
                    w{x_idx+1, y_idx, z_idx} + ...
                    w{x_idx, y_idx+1, z_idx} + ...
                    w{x_idx+1, y_idx+1, z_idx}));
                rdotw(point) = dot(radial(point, :), ...
                    1/4 * (w{x_idx, y_idx, z_idx} + ...
                    w{x_idx+1, y_idx, z_idx} + ...
                    w{x_idx, y_idx+1, z_idx} + ...
                    w{x_idx+1, y_idx+1, z_idx}));
                cdotw(point) = dot(circ(point, :), ...
                    1/4 * (w{x_idx, y_idx, z_idx} + ...
                    w{x_idx+1, y_idx, z_idx} + ...
                    w{x_idx, y_idx+1, z_idx} + ...
                    w{x_idx+1, y_idx+1, z_idx}));
                all_w(point) = 1/4 * (w{x_idx, y_idx, z_idx} + ...
                    w{x_idx+1, y_idx, z_idx} + ...
                    w{x_idx, y_idx+1, z_idx} + ...
                    w{x_idx+1, y_idx+1, z_idx});
            else
                % End of special cases, only generic case left: average all 8 points around this point.
                adotw(point) = dot(axial(point, :), ...
                    1/8 * (w{x_idx, y_idx, z_idx} + ...
                    w{x_idx+1, y_idx, z_idx} + ...
                    w{x_idx, y_idx+1, z_idx} + ...
                    w{x_idx, y_idx, z_idx+1} + ...
                    w{x_idx+1, y_idx+1, z_idx} + ...
                    w{x_idx+1, y_idx, z_idx+1} + ...
                    w{x_idx, y_idx+1, z_idx+1} + ...
                    w{x_idx+1, y_idx+1, z_idx+1}));
                rdotw(point) = dot(radial(point, :), ...
                    1/8 * (w{x_idx, y_idx, z_idx} + ...
                    w{x_idx+1, y_idx, z_idx} + ...
                    w{x_idx, y_idx+1, z_idx} + ...
                    w{x_idx, y_idx, z_idx+1} + ...
                    w{x_idx+1, y_idx+1, z_idx} + ...
                    w{x_idx+1, y_idx, z_idx+1} + ...
                    w{x_idx, y_idx+1, z_idx+1} + ...
                    w{x_idx+1, y_idx+1, z_idx+1}));
                cdotw(point) = dot(circ(point, :), ...
                    1/8 * (w{x_idx, y_idx, z_idx} + ...
                    w{x_idx+1, y_idx, z_idx} + ...
                    w{x_idx, y_idx+1, z_idx} + ...
                    w{x_idx, y_idx, z_idx+1} + ...
                    w{x_idx+1, y_idx+1, z_idx} + ...
                    w{x_idx+1, y_idx, z_idx+1} + ...
                    w{x_idx, y_idx+1, z_idx+1} + ...
                    w{x_idx+1, y_idx+1, z_idx+1}));
                all_w(point) = 1/8 * (w{x_idx, y_idx, z_idx} + ...
                    w{x_idx+1, y_idx, z_idx} + ...
                    w{x_idx, y_idx+1, z_idx} + ...
                    w{x_idx, y_idx, z_idx+1} + ...
                    w{x_idx+1, y_idx+1, z_idx} + ...
                    w{x_idx+1, y_idx, z_idx+1} + ...
                    w{x_idx, y_idx+1, z_idx+1} + ...
                    w{x_idx+1, y_idx+1, z_idx+1});
            end
        end
    end
end


% Write data to a VTK file.
filename_out = "output.vtk";
fid = fopen(filename_out, 'w');
% Write VTK DataFile Version.
fprintf(fid, "# vtk DataFile Version 2.0\n");
% Write title.
fprintf(fid, "Test Output\n");
% Write format data.
fprintf(fid, "ASCII\n");
% Write header.
fprintf(fid, "DATASET POLYDATA\n");

% Write points.
fprintf(fid, "POINTS %d float\n", size(points, 1));
for point = 1 : size(points, 1)
    fprintf(fid, "%f %f %f\n", points(point, 1), points(point, 2), points(point, 3));
end

% Write connectivites.
fprintf(fid, "POLYGONS %d %d\n", size(connectivity, 1), size(connectivity, 1) * 4);
for triangle = 1 : size(connectivity, 1)
    fprintf(fid, "3 %d %d %d\n", connectivity(triangle, 1)-1, connectivity(triangle, 2)-1, connectivity(triangle, 3)-1);
end

% % Write scalars.
fprintf(fid, "POINT_DATA %d\n", size(points, 1));
fprintf(fid, "SCALARS axial_component float 1\n");
fprintf(fid, "LOOKUP_TABLE default\n");
for point = 1 : size(points, 1)
    fprintf(fid, "%f\n", adotw(point));
end
% 
fprintf(fid, "SCALARS radial_component float 1\n");
fprintf(fid, "LOOKUP_TABLE default\n");
for point = 1 : size(points, 1)
    fprintf(fid, "%f\n", rdotw(point));
end
% 
fprintf(fid, "SCALARS circ_component float 1\n");
fprintf(fid, "LOOKUP_TABLE default\n");
for point = 1 : size(points, 1)
    fprintf(fid, "%f\n", cdotw(point));
end


% Write vectors.
fprintf(fid, "VECTORS radial_vector_unit float\n");
for point = 1 : size(points, 1)
    fprintf(fid, "%f %f %f\n", radial(point, 1), radial(point, 2), radial(point, 3));
end
fprintf(fid, "VECTORS axial_vector_unit float\n");
for point = 1 : size(points, 1)
    fprintf(fid, "%f %f %f\n", axial(point, 1), axial(point, 2), axial(point, 3));
end
fprintf(fid, "VECTORS circ_vector_unit float\n");
for point = 1 : size(points, 1)
    fprintf(fid, "%f %f %f\n", circ(point, 1), circ(point, 2), circ(point, 3));
end 
fprintf(fid, "VECTORS w float\n");
for point = 1 : size(points, 1)
    fprintf(fid, "%f %f %f\n", all_w(point, 1), all_w(point, 2), all_w(point, 3));
end 
fclose(fid);









%% Helper functions.

% Write to a VTK file
function vtk_write_tr(filename, title, format, points, connectivity, scalars)
% Writes a triangulation to a VTK file.
% Inputs:
%   - filename: a string representing the name of the file to write to. Be
%   sure to include .vtk at the end
%   - title: a string representing the title of the vtk file.
%   - format: the format of the data, either "ASCII" or "BINARY"
%   - points: the points, in an n x 3 array (x,y,z)
%   - connectivity: the connectivity, in an n x 3 array (p1,p2,p3).

% Write data to a VTK file.
fid = fopen(filename, 'w');
% Write VTK DataFile Version.
fprintf(fid, "# vtk DataFile Version 2.0\n");
% Write title.
fprintf(fid, strcat(title, "\n"));
% Write format data.
fprintf(fid, strcat(format, "\n"));
% Write header.
fprintf(fid, "DATASET POLYDATA\n");

% Write points.
fprintf(fid, "POINTS %d float\n", size(points, 1));
for point = 1 : size(points, 1)
    fprintf(fid, "%f %f %f\n", points(point, 1), points(point, 2), points(point, 3));
end

% Write connectivites.
fprintf(fid, "POLYGONS %d %d\n", size(connectivity, 1), size(connectivity, 1) * 4);
for triangle = 1 : size(connectivity, 1)
    fprintf(fid, "3 %d %d %d\n", connectivity(triangle, 1)-1, connectivity(triangle, 2)-1, connectivity(triangle, 3)-1);
end

% TODO: scalars
% % Write scalars.
% fprintf(fid, "POINT_DATA %d\n", size(points, 1));
% fprintf(fid, "SCALARS w_magnitude float 1\n");
% fprintf(fid, "LOOKUP_TABLE w_mags\n");
% for point = 1 : size(points, 1)
%     fprintf(fid, "%f\n", w_cylinder(point));
% end

fclose(fid);
end

% Compute area of a triangle given by 3 points (x,y,z)
function area = area_of_triangle(point1, point2, point3)
% Inputs:
%   - point1,2,3: each point of the triangle, given by (x,y,z)
vec1 = point2 - point1;
vec2 = point3 - point1;
area = norm(cross(vec1, vec2)) / 2;
end


% Binary search function to find maximum voxel coordinate <= input.
function idx = binary_search_voxel(value, vector)
% Returns the index in vector correpsonding to the coordinate of the
% maximum voxel coordinate less or equal to the input value.
% Inputs:
%   - value: the value we are trying to find the index for
%   - vector: the vector we are trying to find value in. This is a
%     horizontal vector.
% Outputs:
%   - idx: the index of the maximum coordinate in vector that is less than
%   or equal to input.
lo = 1;
hi = size(vector, 2);
while (lo <= hi)
    mid = floor((hi + lo) / 2);
    % Edge case
    if (mid == size(vector, 2))
        break;
        % If we've found the voxel.
    elseif (value >= vector(mid) && value < vector(mid + 1))
        idx = mid;
        return;
        % If we want to search downwards.
    elseif (value < vector(mid))
        hi = mid;
        % If we want to search upwards.
    else
        lo = mid + 1;
    end
end
% If we reach here, value is not in the vector.
idx = -1;
end

function  imshow3D( Img, disprange, initS )
%IMSHOW3D displays 3D grayscale or RGB images in a slice by slice fashion
%with mouse-based slice browsing and window and level adjustment control,
%and auto slice browsing control.
%
% Usage:
% imshow3D ( Image )
% imshow3D ( Image , [] )
% imshow3D ( Image , [LOW HIGH] )
% imshow3D ( Image , [] , initsn )
%
%    Image:      3D image MxNxKxC (K slices of MxN images) C is either 1
%                (for grayscale images) or 3 (for RGB images)
%    [LOW HIGH]: display range that controls the display intensity range of
%                a grayscale image (default: the broadest available range)
%    initsn:     The slice number to be displayed initially (default:
%                mid-slice number)
%
% Use the scroll bar or mouse scroll wheel to switch between slices. To
% adjust window and level values keep the mouse right button pressed, and
% drag the mouse up and down (for level adjustment) or right and left (for
% window adjustment). Window and level adjustment control works only for
% grayscale images.
% "Play" button displays all the slices as a sequence of frames. The time
% interval value can also be adjusted (default time interval is 100 ms).
%
% "Auto W/L" button adjust the window and level automatically for grayscale
% images.
%
% While "Fine Tune" checkbox is checked the window/level adjustment gets 16
% times less sensitive to mouse movement, to make it easier to control
% display intensity rang.
%
% Note: The sensitivity of mouse-based window and level adjustment is set
% based on the user-defined display intensity range; the wider the range,
% the more sensitivity to mouse drag.
%
% Note: IMSHOW3DFULL is a newer version of IMSHOW3D (also available on
% MathWorks) that displays 3D grayscale or RGB images from three
% perpendicular views (i.e., axial, sagittal, and coronal).
%
%   Example
%   --------
%       % To display an image (MRI example)
%       load mri
%       Image = squeeze(D);
%       figure,
%       imshow3D(Image)
%
%       % To display the image, and adjust the display range
%       figure,
%       imshow3D(Image,[20 100]);
%
%       % To define the initial slice number
%       figure,
%       imshow3D(Image,[],5);
%
%   See also IMSHOW.
%
% - Maysam Shahedi (mshahedi@gmail.com)
% - Released: 1.0.0   Date: 2013/04/15
% - Revision: 1.1.0   Date: 2013/04/19
% - Revision: 1.5.0   Date: 2016/09/22
% - Revision: 1.6.0   Date: 2018/06/07
% - Revision: 1.6.1   Date: 2018/10/29
%
sno = size(Img,3);  % number of slices
S = round(sno/2);
PlayFlag = false;   % Play flag, playing when it is 'True'
Tinterv = 100;
global InitialCoord;
MinV = 0;
MaxV = max(Img(:));
LevV = (double( MaxV) + double(MinV)) / 2;
Win = double(MaxV) - double(MinV);
WLAdjCoe = (Win + 1)/1024;
FineTuneC = [1 1/16];    % Regular/Fine-tune mode coefficients
if isa(Img,'uint8')
    MaxV = uint8(Inf);
    MinV = uint8(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'uint16')
    MaxV = uint16(Inf);
    MinV = uint16(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'uint32')
    MaxV = uint32(Inf);
    MinV = uint32(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'uint64')
    MaxV = uint64(Inf);
    MinV = uint64(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'int8')
    MaxV = int8(Inf);
    MinV = int8(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'int16')
    MaxV = int16(Inf);
    MinV = int16(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'int32')
    MaxV = int32(Inf);
    MinV = int32(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'int64')
    MaxV = int64(Inf);
    MinV = int64(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'logical')
    MaxV = 0;
    MinV = 1;
    LevV =0.5;
    Win = 1;
    WLAdjCoe = 0.1;
end
SFntSz = 9;
txtFntSz = 10;
LVFntSz = 9;
WVFntSz = 9;
BtnSz = 10;
if (nargin < 3)
    S = round(sno/2);
else
    S = initS;
    if S > sno
        S = sno;
        warning('Initial slice number out of range');
    elseif S < 1
        S = 1;
        warning('Initial slice number out of range');
    end
end
if (nargin < 2)
    [Rmin Rmax] = WL2R(Win, LevV);
elseif numel(disprange) == 0
    [Rmin Rmax] = WL2R(Win, LevV);
else
    LevV = (double(disprange(2)) + double(disprange(1))) / 2;
    Win = double(disprange(2)) - double(disprange(1));
    WLAdjCoe = (Win + 1)/1024;
    [Rmin Rmax] = WL2R(Win, LevV);
end
clf
axes('position',[0,0.2,1,0.8]), imshow(squeeze(Img(:,:,S,:)), [Rmin Rmax])
FigPos = get(gcf,'Position');
S_Pos = [30 45 uint16(FigPos(3)-100)+1 20];
Stxt_Pos = [30 65 uint16(FigPos(3)-100)+1 15];
Wtxt_Pos = [20 18 60 20];
Wval_Pos = [75 20 50 20];
Ltxt_Pos = [130 18 45 20];
Lval_Pos = [170 20 50 20];
Btn_Pos = [240 20 70 20];
ChBx_Pos = [320 20 80 20];
Play_Pos = [uint16(FigPos(3)-100)+40 45 30 20];
Time_Pos = [uint16(FigPos(3)-100)+35 20 40 20];
Ttxt_Pos = [uint16(FigPos(3)-100)-50 18 90 20];
% W/L Button styles:
WL_BG = ones(Btn_Pos(4),Btn_Pos(3),3)*0.85;
WL_BG(1,:,:) = 1; WL_BG(:,1,:) = 1; WL_BG(:,end-1,:) = 0.6; WL_BG(:,end,:) = 0.4; WL_BG(end,:,:) = 0.4;
% Play Button styles:
Play_BG = ones(Play_Pos(4),Play_Pos(3),3)*0.85;
Play_BG(1,:,:) = 1; Play_BG(:,1,:) = 1; Play_BG(:,end-1,:) = 0.6; Play_BG(:,end,:) = 0.4; Play_BG(end,:,:) = 0.4;
Play_Symb = [0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1; 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1; 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1;...
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1;...
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1;...
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1; 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1; 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1;...
    0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
Play_BG(floor((Play_Pos(4)-13)/2)+1:floor((Play_Pos(4)-13)/2)+13,floor(Play_Pos(3)/2)-7:floor(Play_Pos(3)/2)+6,:) = ...
    repmat(Play_Symb,1,1,3) .* Play_BG(floor((Play_Pos(4)-13)/2)+1:floor((Play_Pos(4)-13)/2)+13,floor(Play_Pos(3)/2)-7:floor(Play_Pos(3)/2)+6,:);
Pause_BG = ones(Play_Pos(4),Play_Pos(3),3)*0.85;
Pause_BG(1,:,:) = 1; Pause_BG(:,1,:) = 1; Pause_BG(:,end-1,:) = 0.6; Pause_BG(:,end,:) = 0.4; Pause_BG(end,:,:) = 0.4;
Pause_Symb = repmat([0, 0, 0, 1, 1, 1, 1, 0, 0, 0],13,1);
Pause_BG(floor((Play_Pos(4)-13)/2)+1:floor((Play_Pos(4)-13)/2)+13,floor(Play_Pos(3)/2)-5:floor(Play_Pos(3)/2)+4,:) = ...
    repmat(Pause_Symb,1,1,3) .* Pause_BG(floor((Play_Pos(4)-13)/2)+1:floor((Play_Pos(4)-13)/2)+13,floor(Play_Pos(3)/2)-5:floor(Play_Pos(3)/2)+4,:);
if sno > 1
    shand = uicontrol('Style', 'slider','Min',1,'Max',sno,'Value',S,'SliderStep',[1/(sno-1) 10/(sno-1)],'Position', S_Pos,'Callback', {@SliceSlider, Img});
    stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String',sprintf('Slice# %d / %d',S, sno), 'FontSize', SFntSz);
    playhand = uicontrol('Style', 'pushbutton','Position', Play_Pos, 'Callback' , @Play);
    set(playhand, 'cdata', Play_BG)
    ttxthand = uicontrol('Style', 'text','Position', Ttxt_Pos,'String','Interval (ms): ',  'FontSize', txtFntSz);
    timehand = uicontrol('Style', 'edit','Position', Time_Pos,'String',sprintf('%d',Tinterv), 'BackgroundColor', [1 1 1], 'FontSize', LVFntSz,'Callback', @TimeChanged);
else
    stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String','2D image', 'FontSize', SFntSz);
end
ltxthand = uicontrol('Style', 'text','Position', Ltxt_Pos,'String','Level: ',  'FontSize', txtFntSz);
wtxthand = uicontrol('Style', 'text','Position', Wtxt_Pos,'String','Window: ',  'FontSize', txtFntSz);
lvalhand = uicontrol('Style', 'edit','Position', Lval_Pos,'String',sprintf('%6.0f',LevV), 'BackgroundColor', [1 1 1], 'FontSize', LVFntSz,'Callback', @WinLevChanged);
wvalhand = uicontrol('Style', 'edit','Position', Wval_Pos,'String',sprintf('%6.0f',Win), 'BackgroundColor', [1 1 1], 'FontSize', WVFntSz,'Callback', @WinLevChanged);
Btnhand = uicontrol('Style', 'pushbutton','Position', Btn_Pos,'String','Auto W/L', 'FontSize', BtnSz, 'Callback' , @AutoAdjust);
set(Btnhand, 'cdata', WL_BG)
ChBxhand = uicontrol('Style', 'checkbox','Position', ChBx_Pos,'String','Fine-tune', 'FontSize', txtFntSz);
set (gcf, 'WindowScrollWheelFcn', @mouseScroll);
set (gcf, 'ButtonDownFcn', @mouseClick);
set(get(gca,'Children'),'ButtonDownFcn', @mouseClick);
set(gcf,'WindowButtonUpFcn', @mouseRelease)
set(gcf,'ResizeFcn', @figureResized)
% -=< Figure resize callback function >=-
    function figureResized(object, eventdata)
        FigPos = get(gcf,'Position');
        S_Pos = [30 45 uint16(FigPos(3)-100)+1 20];
        Stxt_Pos = [30 65 uint16(FigPos(3)-100)+1 15];
        Play_Pos = [uint16(FigPos(3)-100)+40 45 30 20];
        Time_Pos = [uint16(FigPos(3)-100)+35 20 40 20];
        Ttxt_Pos = [uint16(FigPos(3)-100)-50 18 90 20];
        if sno > 1
            set(shand,'Position', S_Pos);
            set(playhand, 'Position', Play_Pos)
            set(ttxthand, 'Position', Ttxt_Pos)
            set(timehand, 'Position', Time_Pos)
        end
        set(stxthand,'Position', Stxt_Pos);
        set(ltxthand,'Position', Ltxt_Pos);
        set(wtxthand,'Position', Wtxt_Pos);
        set(lvalhand,'Position', Lval_Pos);
        set(wvalhand,'Position', Wval_Pos);
        set(Btnhand,'Position', Btn_Pos);
        set(ChBxhand,'Position', ChBx_Pos);
    end
% -=< Slice slider callback function >=-
    function SliceSlider (hObj,event, Img)
        S = round(get(hObj,'Value'));
        set(get(gca,'children'),'cdata',squeeze(Img(:,:,S,:)))
        caxis([Rmin Rmax])
        if sno > 1
            set(stxthand, 'String', sprintf('Slice# %d / %d',S, sno));
        else
            set(stxthand, 'String', '2D image');
        end
    end
% -=< Mouse scroll wheel callback function >=-
    function mouseScroll (object, eventdata)
        UPDN = eventdata.VerticalScrollCount;
        S = S - UPDN;
        if (S < 1)
            S = 1;
        elseif (S > sno)
            S = sno;
        end
        if sno > 1
            set(shand,'Value',S);
            set(stxthand, 'String', sprintf('Slice# %d / %d',S, sno));
        else
            set(stxthand, 'String', '2D image');
        end
        set(get(gca,'children'),'cdata',squeeze(Img(:,:,S,:)))
    end
% -=< Mouse button released callback function >=-
    function mouseRelease (object,eventdata)
        set(gcf, 'WindowButtonMotionFcn', '')
    end
% -=< Mouse click callback function >=-
    function mouseClick (object, eventdata)
        MouseStat = get(gcbf, 'SelectionType');
        if (MouseStat(1) == 'a')        %   RIGHT CLICK
            InitialCoord = get(0,'PointerLocation');
            set(gcf, 'WindowButtonMotionFcn', @WinLevAdj);
        end
    end
% -=< Window and level mouse adjustment >=-
    function WinLevAdj(varargin)
        PosDiff = get(0,'PointerLocation') - InitialCoord;
        Win = Win + PosDiff(1) * WLAdjCoe * FineTuneC(get(ChBxhand,'Value')+1);
        LevV = LevV - PosDiff(2) * WLAdjCoe * FineTuneC(get(ChBxhand,'Value')+1);
        if (Win < 1)
            Win = 1;
        end
        [Rmin, Rmax] = WL2R(Win,LevV);
        caxis([Rmin, Rmax])
        set(lvalhand, 'String', sprintf('%6.0f',LevV));
        set(wvalhand, 'String', sprintf('%6.0f',Win));
        InitialCoord = get(0,'PointerLocation');
    end
% -=< Window and level text adjustment >=-
    function WinLevChanged(varargin)
        LevV = str2double(get(lvalhand, 'string'));
        Win = str2double(get(wvalhand, 'string'));
        if (Win < 1)
            Win = 1;
        end
        [Rmin, Rmax] = WL2R(Win,LevV);
        caxis([Rmin, Rmax])
    end
% -=< Window and level to range conversion >=-
    function [Rmn Rmx] = WL2R(W,L)
        Rmn = L - (W/2);
        Rmx = L + (W/2);
        if (Rmn >= Rmx)
            Rmx = Rmn + 1;
        end
    end
% -=< Window and level auto adjustment callback function >=-
    function AutoAdjust(object,eventdata)
        Win = double(max(Img(:))-min(Img(:)));
        Win (Win < 1) = 1;
        LevV = double(min(Img(:)) + (Win/2));
        [Rmin, Rmax] = WL2R(Win,LevV);
        caxis([Rmin, Rmax])
        set(lvalhand, 'String', sprintf('%6.0f',LevV));
        set(wvalhand, 'String', sprintf('%6.0f',Win));
    end
% -=< Play button callback function >=-
    function Play (hObj,event)
        PlayFlag = ~PlayFlag;
        if PlayFlag
            set(playhand, 'cdata', Pause_BG)
        else
            set(playhand, 'cdata', Play_BG)
        end
        while PlayFlag
            S = S + 1;
            if (S > sno)
                S = 1;
            end
            set(shand,'Value',S);
            set(stxthand, 'String', sprintf('Slice# %d / %d',S, sno));
            set(get(gca,'children'),'cdata',squeeze(Img(:,:,S,:)))
            pause(Tinterv/1000)
        end
    end
% -=< Time interval adjustment callback function>=-
    function TimeChanged(varargin)
        Tinterv = str2double(get(timehand, 'string'));
    end

end
% -=< Maysam Shahedi (mshahedi@gmail.com), October 29, 2018>=-