% Requires:.
%   - Input .stl file is composed of approximately evenly distributed
%   triangles with no "holes".
%   - The anuli of the input .stl file are flat along the x, y, and z
%   planes.
%   - That is, it should be a safe assumption that all points on the anuli
%   are within 0.1*(median edge length) of the boundary value.

compute_all_extensions("xtest.stl", 1, 1, 0, 0, 0, 0);
compute_all_extensions("ytest.stl", 0, 0, 1, 1, 0, 0);
compute_all_extensions("ztest.stl", 0, 0, 0, 0, 1, 1);

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
original_points = input.Points;
original_connectivity = input.ConnectivityList;

% ------------------------REMOVE---------------
% num_x_pos = 1;
% num_x_neg = 1;
% num_y_pos = 0;
% num_y_neg = 0;
% num_z_pos = 0;
% num_z_neg = 0;
%------------------------------------------------

% Compute preliminary triangle neighbors (triangles that share an edge).
triangle_neighbors = map_triangle_neighbors(original_connectivity);
% Compute preliminary point neighbors (points that share an edge).
point_neighbors = map_point_neighbors(original_points, original_connectivity);

% Compute tolerance
tol = compute_tol(original_points, original_connectivity);

% Find bounds of vessel.
[x_min, x_max, y_min, y_max, z_min, z_max] = find_bounds(original_points, original_connectivity);

% Find all triangles whose centroid are within tolerance of x/y/z_min/max.
% These bins are of the form [x, y, z, connectivity_idx].
[x_min_bin, x_max_bin, y_min_bin, y_max_bin, z_min_bin, z_max_bin] = compute_bins(tol, x_min, x_max, y_min, y_max, z_min, z_max, original_points, original_connectivity); % DEBUG: nothing in y_max+bin, y_min_bin, z_max_bin.

% Find and remove all anuli from the vessel.
temp_connectivity = original_connectivity;
if (num_x_neg > 0)
    removed_anuli_conn = remove_anulus(temp_connectivity, x_min_bin);
    temp_connectivity = removed_anuli_conn;
end
if (num_x_pos > 0)
    removed_anuli_conn = remove_anulus(temp_connectivity, x_max_bin);
    temp_connectivity = removed_anuli_conn;
end
if (num_y_neg > 0)
    removed_anuli_conn = remove_anulus(temp_connectivity, y_min_bin);
    temp_connectivity = removed_anuli_conn;
end
if (num_y_pos > 0)
    removed_anuli_conn = remove_anulus(temp_connectivity, y_max_bin);
    temp_connectivity = removed_anuli_conn;
end
if (num_z_neg > 0)
    removed_anuli_conn = remove_anulus(temp_connectivity, z_min_bin);
    temp_connectivity = removed_anuli_conn;
end
if (num_z_pos > 0)
    removed_anuli_conn = remove_anulus(temp_connectivity, z_max_bin);
    temp_connectivity = removed_anuli_conn;
end
removed_anuli_conn = temp_connectivity;
% Remove all zero rows.
removed_anuli_conn( ~any(removed_anuli_conn,2), : ) = [];
TRtesting = triangulation(removed_anuli_conn, original_points);
stlwrite(TRtesting, "sanitycheck2.stl", 'text');

% Use flooding algorithm to get other lateral surfaces.
[inside, outside] = flood_lateral(removed_anuli_conn);
TR = triangulation(inside,original_points);
stlwrite(TR,"latsurf1.stl",'text')
TR = triangulation(outside,original_points);
stlwrite(TR,"latsurf2.stl",'text')

% Separate the anuli on each face from each other.
% These anuli are stored in cell arrays of size num_anuli, where each cell
% contains all the triangles in that anulus.

x_max_anuli_cell = separate_anuli(x_max_bin, num_x_pos);
x_min_anuli_cell = separate_anuli(x_min_bin, num_x_neg);
y_max_anuli_cell = separate_anuli(y_max_bin, num_y_pos);
y_min_anuli_cell = separate_anuli(y_min_bin, num_y_neg);
z_max_anuli_cell = separate_anuli(z_max_bin, num_z_pos);
z_min_anuli_cell = separate_anuli(z_min_bin, num_z_neg);

% Extend anuli.
new_points = original_points;
new_connectivity = original_connectivity;

%% X positive axis:
for i = 1 : num_x_pos % For each anulus,
    anulus = x_max_anuli_cell{i};
    % Determine points on the rim and store them in a bucket.
    out_rim = find_rim(anulus, outside, original_points);
    in_rim = find_rim(anulus, inside, original_points);
    % Extend the anulus outwards
    % FIX: This one looks weird. Possibly out_rim incorrect?
    [new_points, new_connectivity] = extend_anulus(1, anulus, 1.5, 1, in_rim, out_rim, point_neighbors, new_points, new_connectivity);
end

%% X negative axis:
for i = 1 : num_x_neg % For each anulus,
    anulus = x_min_anuli_cell{i};
    % Determine points on the rim and store them in a bucket.
    out_rim = find_rim(anulus, outside, original_points);
    in_rim = find_rim(anulus, inside, original_points);
    % Extend the anulus outwards
    [new_points, new_connectivity] = extend_anulus(1, anulus, 1.5, -1, in_rim, out_rim, point_neighbors, new_points, new_connectivity);
end

%% Y positive axis:
for i = 1 : num_y_pos % For each anulus,
    anulus = y_max_anuli_cell{i};
    % Determine points on the rim and store them in a bucket.
    out_rim = find_rim(anulus, outside, original_points);
    in_rim = find_rim(anulus, inside, original_points);
    % Extend the anulus outwards
    [new_points, new_connectivity] = extend_anulus(2, anulus, 1.5, 1, in_rim, out_rim, point_neighbors, new_points, new_connectivity);
end

%% Y negative axis:
for i = 1 : num_y_neg % For each anulus,
    anulus = y_min_anuli_cell{i};
    % Determine points on the rim and store them in a bucket.
    out_rim = find_rim(anulus, outside, original_points);
    in_rim = find_rim(anulus, inside, original_points);
    % Extend the anulus outwards
    [new_points, new_connectivity] = extend_anulus(2, anulus, 1.5, -1, in_rim, out_rim, point_neighbors, new_points, new_connectivity);
end

%% Z positive axis:
for i = 1 : num_z_pos % For each anulus,
    anulus = z_max_anuli_cell{i};
    % Determine points on the rim and store them in a bucket.
    out_rim = find_rim(anulus, outside, original_points);
    in_rim = find_rim(anulus, inside, original_points);
    % Extend the anulus outwards
    [new_points, new_connectivity] = extend_anulus(3, anulus, 1.5, 1, in_rim, out_rim, point_neighbors, new_points, new_connectivity);
end

%% Z negative axis:
for i = 1 : num_z_neg % For each anulus,
    anulus = z_min_anuli_cell{i};
    % Determine points on the rim and store them in a bucket.
    out_rim = find_rim(anulus, outside, original_points);
    in_rim = find_rim(anulus, inside, original_points);
    % Extend the anulus outwards
    [new_points, new_connectivity] = extend_anulus(3, anulus, 1.5, -1, in_rim, out_rim, point_neighbors, new_points, new_connectivity);
end

% Remove original anuli.
if (num_x_pos > 0)
new_connectivity = remove_anulus(new_connectivity, x_max_bin);
end
if (num_x_neg > 0)
new_connectivity = remove_anulus(new_connectivity, x_min_bin);
end
if (num_y_pos > 0)
new_connectivity = remove_anulus(new_connectivity, y_max_bin);
end
if (num_y_neg > 0)
new_connectivity = remove_anulus(new_connectivity, y_min_bin);
end
if (num_z_pos > 0)
new_connectivity = remove_anulus(new_connectivity, z_max_bin);
end
if (num_z_neg > 0)
new_connectivity = remove_anulus(new_connectivity, z_min_bin);
end
new_connectivity( ~any(new_connectivity,2), : ) = [];

%% Write output.
TR = triangulation(new_connectivity, new_points);
stlwrite(TR, "xyzoutput.stl", 'text')
%%%%%%%%% PICK UP HERE


% Determine points on the rims and store them in a bucket.
% ======FLOODING ALGORITHM FOR RIMS???======
%   - Store rims in cell arrays, where the n-th index corresponds to the
%   n-th anulus on that face of the box
%   - Have 12 total cell arrays, one for each face of the box, including
%   inside and outside rims.

% x_min_rims_out = find_rim(x_min_bin, outside, points);
% x_min_rims_in = find_rim(x_min_bin, inside, points);
% x_max_rims_out = find_rim(x_max_bin, outside, points);
% x_max_rims_in = find_rim(x_max_bin, inside, points);
% y_min_rims_out = find_rim(y_min_bin, outside, points);
% y_min_rims_in = find_rim(y_min_bin, inside, points);
% y_max_rims_out = find_rim(y_max_bin, outside, points);
% y_max_rims_in = find_rim(y_max_bin, inside, points);
% z_min_rims_out = find_rim(z_min_bin, outside, points);
% z_min_rims_in = find_rim(z_min_bin, inside, points);
% z_max_rims_out = find_rim(z_max_bin, outside, points);
% z_max_rims_in = find_rim(z_max_bin, inside, points);
% all_rim_points = cat(1, all_rim_points, find_rim());

% Split each bin into individual rims, stored in cells.
% x_min_rims_out_cells = separate_rims(x_min_rims_out, point_neighbors);

% Extend each anulus outwards.
end



%==========================================================================
%% Helper Functions

% Find x/y/x_min/max bounds for the vessel.
function [x_min, x_max, y_min, y_max, z_min, z_max] = find_bounds(points, connectivity)
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

TR = triangulation(new_connectivity, new_points);
stlwrite(TR, "sanitycheck1.stl", 'text')
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
    for j = 1 : size(connectivity, 1)
        if (i ~= j)
            triangle1 = connectivity(i, :);
            triangle2 = connectivity(j, :);
            if (share_edge(triangle1, triangle2))
                triangle_neighbors{i} = cat(1, triangle_neighbors{i}, j);
                triangle_neighbors{j} = cat(1, triangle_neighbors{j}, i);
            end
        end
    end
end
% Remove double-counting
for i = 1 : size(triangle_neighbors)
    triangle_neighbors{i} = unique(triangle_neighbors{i}(:).');
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
tol = 0.15 * median(edge_lengths);
end


% Compute the triangles whose centroid are within tolerance of the bounds
% of the bounding box, and places them into a bin.
function [x_min_bin, x_max_bin, y_min_bin, y_max_bin, z_min_bin, z_max_bin] = compute_bins(tol, x_min, x_max, y_min, y_max, z_min, z_max, points, connectivity)
% Inputs:
%   - x/y/z_min/max: the boundary values of the bins.
% Outputs:
%   - bins containing triangles that are on the boundaries of the boxes.

% Initialize.
x_min_bin = [];
x_max_bin = [];
y_min_bin = [];
y_max_bin = [];
z_min_bin = [];
z_max_bin = [];

for t = 1 : size(connectivity, 1)
    % Find centroid of triangle.
    %     p1 = points(connectivity(t, 1), :);
    %     p2 = points(connectivity(t, 2), :);
    %     p3 = points(connectivity(t, 3), :);
    %     centroid = (p1 + p2 + p3) ./ 3;
    triangle_coords = zeros(3,3); % each row is a point in the triangle
    % Get coordinates of this triangle
    for node_idx = 1:3
        node_num = connectivity(t, node_idx);
        triangle_coords(node_idx, :) = points(node_num, :);
    end

    % Find centroid of triangle
    centroid = sum(triangle_coords, 1) ./ 3;

    % Determine if centroid is within tol of any side of the bounding box.
    % If within tolerance, add to correpsonding bin.
    if (abs(centroid(1) - x_min) < tol)
        x_min_bin = cat(1, x_min_bin, [connectivity(t, :), t]);
    end
    if (abs(centroid(1) - x_max) < tol)
        x_max_bin = cat(1, x_max_bin, [connectivity(t, :), t]);
    end
    if (abs(centroid(2) - y_min) < tol)
        y_min_bin = cat(1, y_min_bin, [connectivity(t, :), t]);
    end
    if (abs(centroid(2) - y_max) < tol)
        y_max_bin = cat(1, y_max_bin, [connectivity(t, :), t]);
    end
    if (abs(centroid(3) - z_min) < tol)
        z_min_bin = cat(1, z_min_bin, [connectivity(t, :), t]);
    end
    if (abs(centroid(3) - z_max) < tol)
        z_max_bin = cat(1, z_max_bin, [connectivity(t, :), t]);
    end
end
end


% Remove all anuli from one side of the bounding box.
function [removed_anuli_conn] = remove_anulus(connectivity, remove_bin)
% Inputs:
%   - remove_bin: the bin of triangles to remove
% Outputs:
%   - removed_anuli_conn: connectivity of main body after anulus is removed.
removed_anuli_conn = connectivity; % Make a copy.
for idx = 1 : size(remove_bin, 1)
    removed_anuli_conn(remove_bin(idx, 4), :) = [0, 0, 0];
end
end


% Flooding algorithm to get the other two lateral surfaces (inside,
% outside) using a FIFO queue.
function [lat_surf1_conn, lat_surf2_conn] = flood_lateral(removed_anuli_conn)
% Compute triangle neighbors.
neighbors = map_triangle_neighbors(removed_anuli_conn);
% Make copy of connectivity for flooding.
lat_surf1_conn = [];
lat_surf2_conn = removed_anuli_conn;
% BFS.
visited = zeros(size(removed_anuli_conn, 1), 1);
visited(1) = true;
% Queue holds triangle indices from orig_connectivity.
queue = []; % Initialize queue to be an empty queue. Horizonal array.
queue = cat(2, queue, 1);
% Add to lat_surf1_conn and remove from lat_surf2_conn.
lat_surf1_conn = cat(1, lat_surf1_conn, removed_anuli_conn(1, :));
lat_surf2_conn(1, :) = [0,0,0];
while (size(queue, 2) > 0) % While queue is not empty.
    % Dequeue j.
    j = queue(1);
    queue(1, 1) = 0;
    queue( :, all(~queue,1) ) = [];
    % For each neighbor h of j.
    nbrs = neighbors{j};
    for nbridx = 1 : size(nbrs, 2)
        h = nbrs(nbridx);
        if (~visited(h))
            visited(h) = true;
            queue = cat(2, queue, h); % Enqueue h.
            % Add to lat_surf1_conn and remove from lat_surf2_conn.
            lat_surf1_conn = cat(1, lat_surf1_conn, removed_anuli_conn(h, :));
            lat_surf2_conn(h, :) = [0,0,0];
        end
    end
end
% Remove empty rows in lat_surf2_conn.
lat_surf2_conn( ~any(lat_surf2_conn,2), : ) = [];
end


% Compute points on the rims by determining which points are shared between
% an anulus and a lateral surface.
function [rim] = find_rim(surface1, surface2, points)
% Computes a list of points on the rim shared between two surfaces.
% Inputs:
%   - surface1,2: surfaces, represented as list of n x 3 of connectivities
% Outputs:
%   - rim: the rim between surface1 and surface2, represented as a list of n x 4 of x,y,z,idx.
rim = [];
for t1 = 1 : size(surface1, 1) % t1 is a triangle index.
    for t2 = t1 : size(surface2, 1) % t2 is a triangle index.
        for i = 1 : 3 % i is a point index within triangle t1.
            for j = 1 : 3 % j is a point index within triangle t2.
                if (surface1(t1, i) == surface2(t2, j))
                    rim = cat(1, rim, [points(surface1(t1, i), :), surface1(t1, i)]);
                end
            end
        end
    end
end
rim = unique(rim,'row','stable');
end


% % Compute points on the rim using a flooding algorithm and determining
% % which points are shared between an anulus and a lateral surface.
% function [rims] = find_rims(surface1, surface2, points, point_neighbors, num_anuli)
% % All points that are in common between the two surface. all_common_points
% % is a n x 4 array of x,y,z,points_idx.
% all_common_points = find_common_points(surface1, surface2, point);
% rims = cell(num_anuli); % Separates the anuli from each other.
%
% % BFS
% starting_point = all_common_points(1, 4); % Initialize starting_point.
% all_common_points(1, :) = []; % Remove starting_point from all_common_points.
%
%
% queue = []; % Initialize queue to be an empty queue. Horizonal array.
% while (size(all_common_points, 1) > 0)
%     j = queue(1);
%     queue(1, 1) = 0;
%     queue( : , all(~queue, 1)) = [];
%     % foreach neighbor h of j, do
%         % if v_h == False, then
%             % v_h <- True;
%             % enqueue(h);
%
% end
%
%
%
%
% % Iterate through surfaces to find first point in common.
%
% % Remove the point from points array.
%
% % BFS to find all other points on the array using the point to neighbors
% % mapping, considering only points on both surfaces.
%
% end


% Separate the anuli on one boundary from each other.
function[anuli_cell] = separate_anuli(bin, num_anuli)
% Inputs:
%   - bin: the x/y/z_min/max_bin of all points on boundary of a face of box.
%   - num_anuli: the expected number of anuli on this face.
%     Assumption: there are exactly num_anuli groups of triangles on this
%     face.
% Outputs:
%   - anuli_cell: a cell array of size num_anuli, where each cell contains
%   the triangles in one anulus.
% Compute triangle neighbors of the bin.
nbrs = map_triangle_neighbors(bin); % Map of triangle_idx : idxs of nbrs, indices relative to input bin.

% Do BFS to separate the anuli from each other and store in a cell array.
anuli_cell = cell(num_anuli, 1);
for i = 1 : num_anuli
    % Initialize everything.
    visited = zeros(size(bin, 1), 1); % idxs correspond to row idx of bin.
    anulus = []; % Triangles in the anulus we are currently separating from the rest.
    anulus = cat(1, anulus, bin(1, 1:3));
    queue = []; % Initialize queue to be an empty queue of rowidxs. Vertical array.
    queue = cat(1, queue, 1);
    visited(1) = true;
    bin(1, :) = [0, 0, 0, 0]; % Remove from original.
    
    while (size(queue, 1) > 0)
        % Dequeue j.
        j = queue(1, :);
        queue(1, :) = [];
        % For each neighbor h of j
        j_nbrs = nbrs{j}; % Array of neighbor idxs of j.
        for nbr_idx = 1 : size(j_nbrs, 2)
            h = j_nbrs(nbr_idx); % h is an index of bin.
            if (~visited(h))
                visited(h) = true;
                queue = cat(1, queue, h); % Enqueue h.
                % Add to anulus and remove from bin.
                anulus = cat(1, anulus, bin(h, 1:3));
                % Find row idx of bin that contains triangle number h.
                bin(h, :) = [0, 0, 0, 0];
            end
        end
    end
    % Place into anuli_cell.
    anuli_cell{i} = anulus;
end
end

% Extend all points on an anulus outwards, including surface and rim
% points. Removes interior points of original anulus.
function [new_points, new_connectivity] = extend_anulus(axis, anulus, extension_factor, direction, in , out, points_neighbors, orig_points, orig_connectivity)
% Given an input anulus, extends it outwards.
% Inputs:
%   - axis: x/y/z correspond to 1/2/3, respectively
%   - anulus: The triangles in the anulus, in the form of n x 3 conns.
%   - extension_factor: factor by which the points are extended.
%   - direction: -1 if in negative direction, +1 if in positive direction.
%   - in, out: the points on the inside and outside rims of the
%   corresponding anulus.
%   - orig_points/connectivity: as expected.
% Outputs:
%   - new_points/connectivity: as expected with extended anulus.

% Get all points on the anulus.
%Formatted x, y, z, pointidx.
anulus_points = [orig_points(anulus(1, 1), :), anulus(1, 1); % Point 1, triangle 1.
    orig_points(anulus(1, 2), :), anulus(1, 2); % Point 2, triangle 1.
    orig_points(anulus(1, 3), :), anulus(1, 3)]; % Point 3, triangle 1.
anulus_conn = [1, 2, 3];
for i = 2 : size(anulus, 1)
    % Add all points not already added.
    for p = 1 : 3
        point_idx = anulus(i, p);
        if (~ismember(point_idx, anulus_points(:, 4)))
            anulus_points = cat(1, anulus_points, [orig_points(point_idx, :), point_idx]);
        end
    end
    % Add to anulus_conn.
    p1 = find(anulus_points(:, 4) == anulus(i, 1)); %returns the x-th element of the array that is this value.
    p2 = find(anulus_points(:, 4) == anulus(i, 2));
    p3 = find(anulus_points(:, 4) == anulus(i, 3));
    anulus_conn = cat(1, anulus_conn, [p1, p2, p3]);
end

% TODO: insert sanity check here if needed.

% Compute center of anulus:
centroid = sum(anulus_points(:, 1:3), 1) / size(anulus_points, 1);

new_anulus_points = zeros(size(anulus_points, 1), 4);
% Extend all points on anulus by factor around centroid.
for i = 1 : size(anulus_points, 1)
    new_anulus_points(i, 1:3) = extension_factor * (anulus_points(i, 1:3) - centroid(:, 1:3)) + centroid(:, 1:3);
    new_anulus_points(i, 4) = anulus_points(i, 4);
end

% Translate points outwards.
diameter = find_diameter(new_anulus_points);
new_anulus_points(:, axis) = new_anulus_points(:, axis) + 0.2 * diameter * direction;

% Connect all points on anulus.
anulus_conn = anulus_conn + size(orig_points, 1);

new_points = cat(1, orig_points, new_anulus_points(:, 1:3));
new_connectivity = cat(1, orig_connectivity, anulus_conn);
% Connect points with original points.
% Connect in.
original_centroid = sum(anulus_points(:, 1:3), 1) / size(anulus_points, 1); % centroid of original points on anulus
new_centroid = sum(new_anulus_points(:, 1:3)) / size(new_anulus_points, 1); % centroid of new points on anulus
vector_through_centroids = new_centroid - original_centroid; % vector going through the two centroids.

% Base case: nothing added yet, add first triangle.
numadded = 0;
in_nbrs = cell(size(orig_points, 1), 1);
for i = 1 : size(in, 1)
    pointidx = in(i, 4);
    for nbridx = 1 : size(points_neighbors{pointidx}, 2)
        nbr = points_neighbors{pointidx}(nbridx);
        if ismember(nbr, in(:, 4))
            in_nbrs{pointidx} = [in_nbrs{pointidx}, nbr];
        end
    end
end
point1 = in(1, 4); % Find first element to use.
point2 = find(new_anulus_points(:, 4)==point1) + size(orig_points, 1);
point3 = in_nbrs{point1}(1);
in_added = [point1, point2, point3];

% depending on the dot product of vector through centroids and vector to point_centroid, each new triangle is either ccw or cw.
direction_new_triangle = cross(new_points(point2, 1:3) - new_points(point1, 1:3), new_points(point3, 1:3) - new_points(point1, 1:3));
norm_vec_new_triangle = direction_new_triangle/norm(direction_new_triangle);
if (dot(vector_through_centroids, norm_vec_new_triangle) > 0)
    direction = 1;
else
    direction = 2;
end

% new_connectivity = cat(1, new_connectivity, [point1, point3, point2]);
new_connectivity = concat_triangle(new_connectivity, point1, point2, point3, direction);

point1 = point3;
% find extended neighbor of point 1
point3 = find(new_anulus_points(:, 4)==point1) + size(orig_points, 1);
% for i = 1 : size(new_anulus_points)
%     if new_anulus_points(i, 4)==point1
%         point3 = i + size(orig_points, 1);
%         in_added = cat(2, in_added, point3);
%         break;
%     end
% end
% new_connectivity = cat(1, new_connectivity, [point1, point3, point2]);
new_connectivity = concat_triangle(new_connectivity, point1, point2, point3, direction);

numadded = numadded + 4;
% iterative case
while numadded < size(in, 1) * 2 - 1
    point2 = point3;
    if ismember(in_nbrs{point1}(1), in_added)
        point3 = in_nbrs{point1}(2);
    else
        point3 = in_nbrs{point1}(1);
    end
    in_added = [in_added, point3];
    % new_connectivity = cat(1, new_connectivity, [point1, point3, point2]);
    new_connectivity = concat_triangle(new_connectivity, point1, point2, point3, direction);

    point1 = point3;
    % find extended neighbor of point 1
    for i = 1 : size(new_anulus_points, 1)
        if new_anulus_points(i, 4)==point1
            point3 = i + size(orig_points, 1);
            break;
        end
    end
    in_added = [in_added, point3];
    % new_connectivity = cat(1, new_connectivity, [point1, point3, point2]);
    new_connectivity = concat_triangle(new_connectivity, point1, point2, point3, direction);

    numadded = numadded + 2;
end
point2 = point3;
point3 = in_added(1);
% new_connectivity = cat(1, new_connectivity, [point1, point3, point2]);
new_connectivity = concat_triangle(new_connectivity, point1, point2, point3, direction);

point1 = point3;
point3 = in_added(2);

% new_connectivity = cat(1, new_connectivity, [point1, point3, point2]);
new_connectivity = concat_triangle(new_connectivity, point1, point2, point3, direction);

TR = triangulation(new_connectivity, new_points);
stlwrite(TR, "xyzoutput1.stl", 'text')



% connect out
% Base case: nothing added yet, add first triangle.
numadded = 0;
out_nbrs = cell(size(orig_points, 1), 1);
for i = 1 : size(out, 1)
    pointidx = out(i, 4);
    for nbridx = 1 : size(points_neighbors{pointidx}, 2)
        nbr = points_neighbors{pointidx}(nbridx);
        if ismember(nbr, out(:, 4))
            out_nbrs{pointidx} = [out_nbrs{pointidx}, nbr];
        end
    end
end
point1 = out(1, 4); % find first element to use.
point2 = find(new_anulus_points(:, 4)==point1) + size(orig_points, 1);
point3 = out_nbrs{point1}(1);
out_added = [point1, point2, point3];
% depending on the dot product of vector through centroids and vector to point_centroid, each new triangle is either ccw or cw.
direction_new_triangle = cross(new_points(point2, 1:3) - new_points(point1, 1:3), new_points(point3, 1:3) - new_points(point1, 1:3));
norm_vec_new_triangle = direction_new_triangle/norm(direction_new_triangle);
if (dot(vector_through_centroids, norm_vec_new_triangle) > 0)
    direction = 2;
else
    direction = 1;
end
% new_connectivity = cat(1, new_connectivity, [point1, point3, point2]);
new_connectivity = concat_triangle(new_connectivity, point1, point2, point3, direction);

point1 = point3;
% find extended neighbor of point 1
for i = 1 : size(new_anulus_points)
    if new_anulus_points(i, 4)==point1
        point3 = i + size(orig_points, 1);
        out_added = cat(2, out_added, point3);
        break;
    end
end

% new_connectivity = cat(1, new_connectivity, [point1, point3, point2]);
new_connectivity = concat_triangle(new_connectivity, point1, point2, point3, direction);

numadded = numadded + 4;
% iterative case
while numadded < size(out, 1) * 2 - 1
    point2 = point3;
    if ismember(out_nbrs{point1}(1), out_added)
        point3 = out_nbrs{point1}(2);
    else
        point3 = out_nbrs{point1}(1);
    end
    out_added = [out_added, point3];

    % new_connectivity = cat(1, new_connectivity, [point1, point3, point2]);
    new_connectivity = concat_triangle(new_connectivity, point1, point2, point3, direction);

    point1 = point3;
    % find extended neighbor of point 1
    for i = 1 : size(new_anulus_points, 1)
        if new_anulus_points(i, 4)==point1
            point3 = i + size(orig_points, 1);
            break;
        end
    end
    out_added = [out_added, point3];

    % new_connectivity = cat(1, new_connectivity, [point1, point3, point2]);
    new_connectivity = concat_triangle(new_connectivity, point1, point2, point3, direction);

    numadded = numadded + 2;
end
point2 = point3;
point3 = out_added(1);

% new_connectivity = cat(1, new_connectivity, [point1, point3, point2]);
new_connectivity = concat_triangle(new_connectivity, point1, point2, point3, direction);

point1 = point3;
point3 = out_added(2);

% new_connectivity = cat(1, new_connectivity, [point1, point3, point2]);
new_connectivity = concat_triangle(new_connectivity, point1, point2, point3, direction);

TR = triangulation(new_connectivity, new_points);
stlwrite(TR, "xyzoutput2.stl", 'text')

end

function [connectivity] = concat_triangle(connectivity, point1, point2, point3, direction)
% Concatenates one triangle of three points to connectivity in a direction
% that ensures the right normal vector direction.
% Inputs:
%   - connectivity: as expected
%   - point1, point2, point3: the three points we add to the connectivity list.
%   - direction: 1 if points go clockwise, 2 if points go counter-clockwise.
if (direction == 1)
    connectivity = [connectivity; point1, point2, point3];
else
    connectivity = [connectivity; point1, point3, point2];
end
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







