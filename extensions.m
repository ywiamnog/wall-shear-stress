% Requires:
%   - Input .stl file is composed of approximately evenly distributed
%   triangles with no "holes".
%   - The anuli of the input .stl file are flat along the x, y, and z
%   planes.
%   - That is, it should be a safe assumption that all points on the anuli
%   are within 0.1*(median edge length) of the boundary value.
compute_all_extensions("testvessel.stl", 0, 0, 0, 0, 1, 1);

function compute_all_extensions(filename, num_x_pos, num_x_neg, num_y_pos, num_y_neg, num_z_pos, num_z_neg)
% Computes and adds all extensions to input .stl file.
% Inputs:
%   - filename: a string, the name of the .stl file we are reading
%   - num_x_pos: an integer, the number of extensions expected on the x 
%     axis in the positive direction
%   - num_x_neg: an integer, the number of extensions expected on the x
%     axis in the negative direction.
%   - other inputs: similar to num_x_pos/neg.

%% Setup
% Read in input .stl file.
input = stlread(filename);
original_points = output.Points;
original_connectivity = output.ConnectivityList;

% Compute preliminary triangle neighbors (triangles that share an edge).
triangle_neighbors = map_triangle_neighbors(original_connectivity);
% Compute preliminary point neighbors (points that share an edge).
point_neighbors = map_point_neighbors(original_points);

% Compute tolerance
tol = compute_tol(original_points, original_connectivity);

% Find bounds of vessel.
[x_min, x_max, y_min, y_max, z_min, z_max] = find_bounds(original_points);

% Find all triangles whose centroid are within tolerance of x/y/z_min/max.
[x_min_bin, x_max_bin, y_min_bin, y_max_bin, z_min_bin, z_max_bin] = compute_bins(tol, x_min, x_max, y_min, y_max, z_min, z_max, original_points, original_connectivity);

% Find and remove all anuli from the vessel.
%%%%%%%%% PICK UP HERE

%% Extend in the positive x direction.
for num = 1 : num_x_pos
    % TODO: try to do a try-catch here to catch wrong input number of
    % extensions errors??
end

%% Extend in the negative x direction.
for num = 1 : num_x_neg
end

%% Extend in the positive y direction.
for num = 1 : num_y_pos
end

%% Extend in the negative y direction.
for num = 1 : num_y_neg
end

%% Extend in the positive z direction.
for num = 1 : num_z_pos
end

%% Extend in the negative z direction.
for num = 1 : num_z_neg
end

end




%==========================================================================
%% Helper Functions

% Find x/y/x_min/max bounds for the vessel.
function [x_min, x_max, y_min, y_max, z_min, z_max] = find_bounds(points)
x_max = max(points(:,1));
x_min = min(points(:,1));
y_max = max(points(:,2));
y_min = min(points(:,2));
z_max = max(points(:,3));
z_min = min(points(:,3));
end


% Maps each triangle index to the indices of triangle that the triangle
% shares an edge with.
function [triangle_neighbors] = map_triangle_neighbors(connectivity)
% Creates a mapping of triangle_idx : neighbors_idx
% The key is the row index, neighbors are stored in the rows. 
% If neighbor = 0, that means there is no neighbor
% Inputs:
%   - connectivity: the input connectivity.
triangle_neighbors = cell(size(connectivity, 1), 1);
for i = 1 : size(connectivity, 1)
    for j = i : size(connectivity, 1)
        if (i ~= j)
            triangle1 = connectivity(i, :);
            triangle2 = connectivity(j, :);
            if (share_edge(triangle1, triangle2))
                neighbors{i} = cat(1, neighbors{i}, j);
                neighbors{j} = cat(1, neighbors{j}, i);
            end
        end
    end
end
end


% Maps each point index to the indices of points that the triangle shares
% an edge with.
function [point_neighbors] = map_point_neighbors(points, connectivity)
% Creates a mapping of point : neighbors, where neighbors are other points
% that this point shares an edge with
% The key is the row index, neighbors are stored in the rows. 
% If neighbor = 0, that means there is no neighbor
% Inputs:
%   - points is in the form [x, y, z]
%   - connectivity is in the form [point1, point2, point3]
point_neighbors = cell(size(points, 1), 1);
for i = 1 : size(connectivity, 1)
    point1 = connectivity(i, 1);
    point2 = connectivity(i, 2);
    point3 = connectivity(i, 3);
    point_neighbors{point1} = cat(1, point_neighbors{point1}, [point2, point3]);
    point_neighbors{point2} = cat(1, point_neighbors{point2}, [point1, point3]);
    point_neighbors{point3} = cat(1, point_neighbors{point3}, [point1, point2]);
end
% Remove duplicates
for i = 1 : size(point_neighbors)
    point_neighbors{i} = unique(point_neighbors{i}(:).');
end
end


% Returns a boolean representing whether triangle1 and triangle2 share an
% edge.
function [share_edge_bool] = share_edge(triangle1, triangle2)
% Inputs:
%   - triangle1, triangle2: triangles in a connectivity list. Both are 1x3
%   arrays that represent p1, p2, p3 point indices of a triangle.
num_shared_points = 0;
share_edge_bool = false; % Initialize to false.
% Iterate through points
for point1 = 1 : 3 % Point in triangle1.
    for point2 = 1 : 3 % Point in triangle2.
        if (triangle1(point1) == triangle2(point2))
            num_shared_points = num_shared_points + 1;
        end
    end
end
% If the two triangles share exactly two points, they share an edge.
if (num_shared_points == 2)
    share_edge_bool = true;
end
end


% Computes tolerance based on 0.1 * median length of an edge
function [tol] = compute_tol(points, connectivity)
num_edges = 0;
edge_lengths = zeros(size(connectivity, 1) * 3, 1);
for i = 1 : size(connectivity, 1)
    p1 = connectivity(i, 1);
    p2 = connectivity(i, 2);
    p3 = connectivity(i, 3);

    num_edges = num_edges + 1;
    edge_lengths(num_edges) = abs(norm(points(p1, :) - points(p2, :)));
    num_edges = num_edges + 1;
    edge_lengths(num_edges) = abs(norm(points(p1, :) - points(p3, :)));
    num_edges = num_edges + 1;
    edge_lengths(num_edges) = abs(norm(points(p2, :) - points(p3, :)));
end
tol = 0.1 * median(edge_lengths);
end


% Compute the triangles whose centroid are within tolerance of the bounds
% of the bounding box, and places them into a bin.
function [x_min_bin, x_max_bin, y_min_bin, y_max_bin, z_min_bin, z_max_bin] = compute_bins(tol, x_min, x_max, y_min, y_max, z_min, z_max, original_points, original_connectivity)
% Inputs:
%   - x/y/z_min/max: the boundary values of the bins.

% Initialize.
x_min_bin = [];
x_max_bin = [];
y_min_bin = [];
y_max_bin = [];
z_min_bin = [];
z_max_bin = [];

for t = 1 : size(connectivity, 1)
    % Find centroid of triangle.
    p1 = points(connectivity(t, 1), :);
    p2 = points(connectivity(t, 2), :);
    p3 = points(connectivity(t, 3), :);
    centroid = (p1 + p2 + p3) ./ 3;

    % Determine if centroid is within tol of any side of the bounding box.
    % If within tolerance, add to correpsonding bin.
    if (abs(centroid(1) - x_min) < tol)
        x_min_bin = cat(1, x_min_bin, connectivity(t, :));
    end
    if (abs(centroid(1) - x_max) < tol)
        x_max_bin = cat(1, x_max_bin, connectivity(t, :));
    end
    if (abs(centroid(1) - y_min) < tol)
        y_min_bin = cat(1, y_min_bin, connectivity(t, :));
    end
    if (abs(centroid(1) - y_max) < tol)
        y_max_bin = cat(1, y_max_bin, connectivity(t, :));
    end
    if (abs(centroid(1) - z_min) < tol)
        z_min_bin = cat(1, z_min_bin, connectivity(t, :));
    end
    if (abs(centroid(1) - z__max) < tol)
        z_max_bin = cat(1, z_max_bin, connectivity(t, :));
    end
end
end




    