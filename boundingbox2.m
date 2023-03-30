% % Read in input .stl file
% output = stlread("testvessel.stl");
% points = output.Points;
% connectivity = output.ConnectivityList;

% % boolean flags representing if anuli are expected on corresponding face
% x_min_bool = false;
% x_max_bool = false;
% y_min_bool = false;
% y_max_bool = false;
% z_min_bool = true;
% z_max_bool = true;
% testing = true;

% % Compute tolerance
% [tol, edge_lengths, num_edges] = compute_tol(points, connectivity);

% % make neighbors maps
% triangle_neighbors = maptriangleneighbors(connectivity);
% points_neighbors = mappointneighbors(points, connectivity);

% % Compute bounds and create bounding box;
% [x_min, x_max, y_min, y_max, z_min, z_max, new_connectivity, new_points] = bounding_box(points, connectivity);
% bounds = [x_min, x_max, y_min, y_max, z_min, z_max];
% if (testing)
%     TR = triangulation(new_connectivity,new_points);
%     stlwrite(TR,"testvesselwithbox.stl",'text')
% end

% % Compute bins
% [x_min_bin, x_max_bin, y_min_bin, y_max_bin, z_min_bin, z_max_bin, testing_centroids] ...
%     = create_bins(bounds, points, connectivity, tol);
% if (testing)
%     z_bins = cat(1, z_min_bin, z_max_bin);
%     stlwrite(triangulation(z_bins, points), "anuli.stl", "text");
% end

% % Remove anuli to get the rest of the surfaces
% [removed_anuli_conn] = remove_anuli(connectivity, z_min_bin, z_max_bin);
% if (testing)
%     TR = triangulation(removed_anuli_conn,new_points);
%     stlwrite(TR,"removedanuli.stl",'text')
% end

% % Use flooding algorithm to get other lateral surfaces
% [inside, ouside] = flood(removed_anuli_conn);
% if (testing)
%     TR = triangulation(inside,new_points);
%     stlwrite(TR,"latsurf1.stl",'text')
%     TR = triangulation(ouside,new_points);
%     stlwrite(TR,"latsurf2.stl",'text')
% end

% The edges on the rim are the edges that are both in a triangle on a
% latsurf and a triangle in a zbin.
z_bot_out = find_rim(z_min_bin, ouside, points);
z_bot_in = find_rim(z_min_bin, inside, points);
z_top_out = find_rim(z_max_bin, ouside, points);
z_top_in = find_rim(z_max_bin, inside, points);

% Extend the anuli outwards
[extended_z_min_points, extended_z_min_connectivity] = extend_anulus(z_min_bin, points, connectivity);
[extended_z_max_points, extended_z_max_connectivity] = extend_anulus(z_max_bin, points, connectivity);
if (testing)
    TR = triangulation(extended_z_min_connectivity,extended_z_min_points);
    stlwrite(TR,"z_min_extended.stl",'text')
    TR = triangulation(extended_z_max_connectivity,extended_z_max_points);
    stlwrite(TR,"z_max_extended.stl",'text')
end


%==========================================================================



% compute_tol : Computes tolerance based on 0.1 * median length of an edge.
function [tol, edge_lengths, num_edges] = compute_tol(points, connectivity)
% Computes a rough value for the tolerance
num_edges = 0;
edge_lengths = zeros(size(connectivity, 1) * 3, 1);
for triangle = 1 : size(connectivity, 1)
    i = connectivity(triangle, 1);
    j = connectivity(triangle, 2);
    k = connectivity(triangle, 3);

    num_edges = num_edges + 1;
    edge_lengths(num_edges) = abs(norm(points(i, :) - points(j, :)));
    num_edges = num_edges + 1;
    edge_lengths(num_edges) = abs(norm(points(j, :) - points(k, :)));
    num_edges = num_edges + 1;
    edge_lengths(num_edges) = abs(norm(points(i, :) - points(k, :)));
end
sort(edge_lengths);
tol = 0.1 * median(edge_lengths);
end

function [shares_edge] = share_edge(triangle1, triangle2)
% Determines if the two input triangles share an edge.
num_common_points = 0;
shares_edge = false; % initialize to false
% iterate through points and set to false if no edge shared
for point1 = 1 : 3 % point in triangle1
    for point2 = 1 : 3 % point in triangle2
        if (triangle1(point1)==triangle2(point2))
            num_common_points = num_common_points + 1;
        end
    end
end
% if 2 points in common, they share an edge
if num_common_points==2
    shares_edge = true;
end
end

function [neighbors] = maptriangleneighbors(connectivity)
% Creates a mapping of triangle : neighbors
% The key is the row index, neighbors are stored in the rows. 
% If neighbor = 0, that means there is no neighbor
neighbors = cell(size(connectivity, 1), 1);
for i = 1 : size(connectivity)
    for j = 1 : size(connectivity)
        if (i ~= j)
            triangle1 = connectivity(i, :);
            triangle2 = connectivity(j, :);
            if (share_edge(triangle1, triangle2))
                neighbors{i} = [neighbors{i}; j];
                neighbors{j} = [neighbors{j}; i];
                % TODO: remove double-counting
            end
        end
    end
end
% remove double-counting
for i = 1 : size(neighbors)
    neighbors{i} = unique(neighbors{i}(:).');
end
end

function [neighbors] = mappointneighbors(points, connectivity)
% Creates a mapping of point : neighbors, where neighbors are other points
% that this point shares an edge with
% The key is the row index, neighbors are stored in the rows. 
% If neighbor = 0, that means there is no neighbor
% Inputs:
%   - points is in the form [x, y, z]
%   - connectivity is in the form [point1, point2, point3]--orig_conn. 
neighbors = cell(size(points, 1), 1);
for i = 1 : size(connectivity, 1)
    point1 = connectivity(i, 1);
    point2 = connectivity(i, 2);
    point3 = connectivity(i, 3);
    neighbors{point1} = [neighbors{point1}; point2; point3];
    neighbors{point2} = [neighbors{point2}; point1; point3];
    neighbors{point3} = [neighbors{point3}; point1; point2];
end
% remove double-counting
for i = 1 : size(neighbors)
    neighbors{i} = unique(neighbors{i}(:).');
end
end


% bounding_box : Computes a bounding box and returns x/y/z_min/max as well as new connectivity and points.
function[x_min, x_max, y_min, y_max, z_min, z_max, new_connectivity, new_points] = ...
    bounding_box(points, connectivity)
% Creates a bounding box based on the input .stl file
% Inputs:
%   points - an array of points
%   connectivity - an array representing triangles

% Get x, y, z, max and mins
x_max = max(points(:,1));
x_min = min(points(:,1));
y_max = max(points(:,2));
y_min = min(points(:,2));
z_max = max(points(:,3));
z_min = min(points(:,3));

% Create coordinates and connectivities of bounding box
box_points = [x_min, y_min, z_min; %1
    x_min, y_min, z_max;    %2
    x_min, y_max, z_min;    %3
    x_min, y_max, z_max;    %4
    x_max, y_min, z_min;    %5
    x_max, y_min, z_max;    %6
    x_max, y_max, z_min;    %7
    x_max, y_max, z_max];  %8

box_connectivity = [5, 1, 6;
    2, 1, 6;
    4, 2, 6;
    6, 4, 8;
    6, 5, 7;
    7, 6, 8;
    3, 1, 5;
    5, 3, 7;
    2, 1, 3;
    3, 2, 4;
    4, 3, 8;
    7, 3, 8];
points_size = size(points);
points_rows = points_size(1);
box_connectivity = box_connectivity + points_rows;

% Concatenate bounding box points and connectivities to original
new_points = cat(1, points, box_points);
new_connectivity = cat(1, connectivity, box_connectivity);
end


% create_bins : Computes bins of points at each bound
function [x_min_bin, x_max_bin, y_min_bin, y_max_bin, z_min_bin, ...
    z_max_bin, testing_centroids]= create_bins(bounds, points, connectivity, tol)
% This function creates 6 bins, containing the triangles whose centroids
% are within tolerance of the edge of the bounding box.
% Inputs:
%   bounds - an array [x_min, x_max, y_min, y_max, z_min, z_max]
%   points - an array of point coordinates
%   connectivity - an array of triangle connectivities
% Outputs:
%   bins -  an array [x_min bin, x_max bin, y_min bin, y_max bin, 
%           z_min bin, z_max bin]
%           Each row of the array represents a triangle on the boundary.
x_min = bounds(1); x_min_bin = []; x_min_points = []; x_min_conn = [];
x_max = bounds(2); x_max_bin = []; x_max_points = []; x_max_conn = [];
y_min = bounds(3); y_min_bin = []; y_min_points = []; y_min_conn = [];
y_max = bounds(4); y_max_bin = []; y_max_points = []; y_max_conn = [];
z_min = bounds(5); z_min_bin = []; z_min_points = []; z_min_conn = [];
z_max = bounds(6); z_max_bin = []; z_max_points = []; z_max_conn = [];
testing_centroids = [];

% Find triangles whose centroids are on x, y, z - min, max
for triangle = 1:length(connectivity)
    % t1x t1y t1z
    % t2x t2y t2z
    % t3x t3y t3z
    triangle_coords = zeros(3,3); % each row is a point in the triangle
    % Get coordinates of this triangle
    for node_idx = 1:3
        node_num = connectivity(triangle, node_idx);
        triangle_coords(node_idx, :) = points(node_num, :);
    end

    % Find centroid of triangle
    centroid = sum(triangle_coords, 1) ./ 3;
    testing_centroids = [testing_centroids; centroid];

    % Check if centroid is within tolerance of the bounding box.
    % If within tolerance, add to corresponding bin.
%     diff = abs(centroid(1) - x_min);
    if abs(centroid(1) - x_min) < tol
        x_min_bin = [x_min_bin; connectivity(triangle, :)];
    end

%     diff = abs(centroid(1) - x_max);
    if abs(centroid(1) - x_max) < tol
        x_max_bin = [x_max_bin; connectivity(triangle, :)];
    end

%     diff = abs(centroid(2) - y_min);
    if abs(centroid(2) - y_min) < tol
        y_min_bin = [y_min_bin; connectivity(triangle, :)];
    end

%     diff = abs(centroid(2) - y_max);
    if abs(centroid(2) - y_max) < tol
        y_max_bin = [y_max_bin; connectivity(triangle, :)];
    end

%     diff = abs(centroid(3) - z_min);
    if abs(centroid(3) - z_min) < tol
        z_min_bin = [z_min_bin; connectivity(triangle, :)];
    end

%     diff = abs(centroid(3) - z_max);
    if abs(centroid(3) - z_max) < tol
        z_max_bin = [z_max_bin; connectivity(triangle, :)];
    end
end
end


% remove_anuli : removes anuli from the rest of the surfaces
function [removed_anuli_conn] = remove_anuli(original_connectivity, z_min_bin, z_max_bin)
% Remove min bin triangles from original connectivity
removed_anuli_conn = original_connectivity; % make a copy
for orig_elem = 1 : size(original_connectivity, 1) % TODO: very brute force, can it be better?
    for anuli_elem = 1 : size(z_min_bin)
        if (sum(removed_anuli_conn(orig_elem, :) ~= z_min_bin(anuli_elem, :)) == 0|| ...
                sum(removed_anuli_conn(orig_elem, :) ~= z_max_bin(anuli_elem, :)) == 0)
            removed_anuli_conn(orig_elem, :) = [0, 0, 0]; % Set to be all zeros
        end
    end
end
% remove all zero rows
removed_anuli_conn( ~any(removed_anuli_conn,2), : ) = [];
disp("Completed removing anuli")
end


% flood : flooding algorithm to get the other two lateral surfaces using a FIFO queue.
function [lat_surf1_conn, lat_surf2_conn] = flood(removed_anuli_conn)
% Compute triangle neighbors
neighbors = maptriangleneighbors(removed_anuli_conn);
% make copy of connectivity for flooding 
lat_surf1_conn = [];
lat_surf2_conn = removed_anuli_conn; 
% BFS
visited = zeros(size(removed_anuli_conn, 1), 1);
visited(1) = true;
% queue holds triangle indices from orig_connectivity
queue = []; % initialize queue to be an empty queue. Horizonal array.
queue = cat(2, queue, 1);
% add to lat_surf1_conn and remove from lat_surf2_conn
lat_surf1_conn = cat(1, lat_surf1_conn, removed_anuli_conn(1, :));
lat_surf2_conn(1, :) = [0,0,0];
while (size(queue, 2) > 0) % while queue is not empty
    % dequeue j
    j = queue(1); 
    queue(1, 1) = 0;
    queue( :, all(~queue,1) ) = [];
    % for each neighbor h of j
    nbrs = neighbors{j};
    for nbridx = 1 : size(nbrs, 2)
        h = nbrs(nbridx);
        if (~visited(h))
            visited(h) = true;
            queue = cat(2, queue, h); % enqueue h
            % add to lat_surf1_conn and remove from lat_surf2_conn
            lat_surf1_conn = cat(1, lat_surf1_conn, removed_anuli_conn(h, :));
            lat_surf2_conn(h, :) = [0,0,0];
        end
    end
end
% remove empty rows in lat_surf2_conn
lat_surf2_conn( ~any(lat_surf2_conn,2), : ) = [];

end


% find_rim : Compute points on the rims by determining which edges are shared between
% an anuli and a lateral surface.
% NOTE: latsurf1 is inside, latsurf2 is outside
function [rim] = find_rim(surface1, surface2, points)
% Computes a list of points on the rim shared between two surfaces.
% Inputs:
%   - surface1,2: surfaces, represented as list of n x 3 of connectivities
% Outputs:
%   - rim: the rim between surface1 and surface2, represented as a list of n x 3 of coordinate points.
rim = [];
for i = 1 : size(surface1, 1) % i is a triangle index
    for j = 1 : size(surface2, 1) % j is a triangle index
        for n = 1 : 3 % n is a point index in a triangle
            for m = 1 : 3 % m is a point index in a triangle
                if (surface1(i, n)==surface2(j, m))
                    rim = cat(1, rim, points(surface1(i, n), :));
                end
            end
        end
    end
end
end


% extend_anulus : extend all points on an anuli outward
function [new_points, new_connectivity] = extend_anulus(anulus, points, connectivity)
% Given an input anulus, extends it outwards.
% Inputs:
%   - anulus: in the form of n x 3 of connectivities
%   - points, connectivity: as expected
% Outputs:
%   - new_points, new_connectivity: as expected with extended anulus
%     Note that these have not been connected to the main body yet

% get all points on the anulus
anulus_points = [points(anulus(1, 1), :), anulus(1, 1);
    points(anulus(1, 2), :), anulus(1, 2);
    points(anulus(1, 3), :), anulus(1, 3)]; % formatted x, y, z, pointidx
anulus_conn = [1, 2, 3];
for i = 2 : size(anulus, 1)
    % add all points not already added
    for p = 1 : 3
        pointidx = anulus(i, p);
        if (~ismember(pointidx, anulus_points(:, 4)))
            anulus_points = cat(1, anulus_points, [points(pointidx, :), pointidx]);
        end
    end
    % add to anulus_conn
    p1 = find(anulus_points(:, 4) == anulus(i, 1)); %returns the x-th element of the array that is this value.
    p2 = find(anulus_points(:, 4) == anulus(i, 2));
    p3 = find(anulus_points(:, 4) == anulus(i, 3));
    anulus_conn = cat(1, anulus_conn, [p1, p2, p3]);
end

% SANITY CHECK THAT ANULUS_CONN AND ANULUS_POINTS ARE CORRECT
TR = triangulation(anulus_conn,anulus_points(:, 1:3));
stlwrite(TR,"temp.stl",'text')

% Compute center of vessel
centroid = sum(points(1:3), 1) / size(points, 1);

% Extend all points on anulus by factor around centroid
for i = 1 : size(anulus_points, 1)
    anulus_points(i, 1:3) = 2 *(anulus_points(i, 1:3) - centroid) + centroid;
end

% Connect all points on anulus
anulus_conn = anulus_conn + size(points, 1);

new_points = cat(1, points, anulus_points(:, 1:3));
new_connectivity = cat(1, connectivity, anulus_conn);
end

function [max_diameter] = find_diameter(bin)
% Helper function: Given a set of points in the bin, finds the max diameter.
% The max diameter = the largest distance of (point - point)
% Input: bin of points, formatted as [x, y, z, pointidx]
max_diameter = 0;
for i = 1 : size(bin, 1)
    for j = 1 : size(bin, 1)
        diameter = abs(norm(bin(i, 1:3) - (bin(j, 1:3))));
        max_diameter = max(max_diameter, diameter);
    end
end
disp(max_diameter)
end