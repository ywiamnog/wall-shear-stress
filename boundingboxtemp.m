% Read in input .stl file
output = stlread("testvessel.stl");
points = output.Points;
connectivity = output.ConnectivityList;

% TODO: compute avg edge length, tol = 0.1(avg edge length)

% % bounds is an array [x_min, x_max, y_min, y_max, z_min, z_max]
% [x_min, x_max, y_min, y_max, z_min, z_max] = ...
%     bounding_box(points, connectivity);
% bounds = [x_min, x_max, y_min, y_max, z_min, z_max];
% % bounding_box("testvessel.stl", "box.stl", false)
% 
% % Create bins. Contains triangle numbers in the original connectivity list.
% [x_min_bin, x_max_bin, y_min_bin, y_max_bin, z_min_bin, z_max_bin, testing_centroids] = ...
%     create_bins(bounds, points, connectivity);
% testing_centroids = sortrows(testing_centroids, 3);
% 
% % create extensions
% % [bot_in, bot_out, top_in, top_out] = create_extensions(points, connectivity, z_min_bin, z_max_bin, z_min, z_max, 0.07) 
% [conn_copy] = remove_anuli(connectivity, z_min_bin, z_max_bin);
% [top, bot] = split_top_bot(conn_copy, points, z_min, z_max, 0.07);
% [top_out, top_in, bot_out, bot_in] = distribute_bins(conn_copy, top, bot, points);
[top_out_new, top_in_new, bot_out_new, bot_in_new, points, connectivity] = extend_points(top_out, top_in, bot_out, bot_in, points, connectivity);
[top_out_nbrs, top_in_nbrs, bot_out_nbrs, bot_in_nbrs] = connectpoints (top_out, top_in, bot_out, bot_in, top_out_new, top_in_new, bot_out_new, bot_in_new, points, connectivity);


% make neighbors map
% triangleneighbors = maptriangleneighbors(connectivity);
% pointsneighbors = mappointneighbors(points, connectivity);


function[x_min, x_max, y_min, y_max, z_min, z_max] = ...
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

% Write new .stl file
TR = triangulation(new_connectivity,new_points);
% TR = triangulation(box_connectivity,box_points); % For testing purposes
stlwrite(TR,"testvesselwithbox.stl",'text')

end



function [x_min_bin, x_max_bin, y_min_bin, y_max_bin, z_min_bin, ...
    z_max_bin, testing_centroids]= create_bins(bounds, points, connectivity)
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

tol = 0.07;

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
stlwrite(triangulation(z_max_bin, points), "z_max_bin.stl", "text");
stlwrite(triangulation(z_min_bin, points), "z_min_bin.stl", "text");

end


% helper function to calc min/max/mean triangle edge lengths to select tol.


% create extensions
function [conn_copy] = remove_anuli(original_connectivity, z_min_bin, z_max_bin)
% Remove min bin triangles from original connectivity
conn_copy = original_connectivity; % make a copy
for orig_elem = 1 : size(original_connectivity, 1) % very brute force, can it be better?
    for end_elem = 1 : size(z_min_bin)
        if (sum(conn_copy(orig_elem, 1) ~= z_min_bin(end_elem)) == 0 || ...
                sum(conn_copy(orig_elem, 1) ~= z_max_bin(end_elem)) == 0)
            conn_copy(orig_elem, :) = [0, 0, 0]; % Set to be all zeros
        end
    end
end
% remove all zero rows
conn_copy( ~any(conn_copy,2), : ) = [];
disp("Completed removing anuli")
end


function [top, bot] = split_top_bot(conn_copy, original_points, z_min, z_max, tol)
% split the triangles into two bins of points: top vs bot
% both of the form [x, y, z, pointidx, triangleidx]
top = [];
bot = [];
for idx = 1 : size(conn_copy, 1)
    triangle = conn_copy(idx, :);
    % determine if on rim
    num_edges_shared = 0;
    for test_idx = 1 : size(conn_copy, 1)
        test_tri = conn_copy(test_idx, :);
        if idx~=test_idx
            shares_edge = share_edge(conn_copy(triangle, :), conn_copy(test_tri, :));
            if shares_edge==true
                num_edges_shared = num_edges_shared + 1;
            end
        end
    end
    % if on rim, identify points on rim and determine if on top or bottom
    if num_edges_shared < 3
        for point = 1 : 3
            if abs(original_points(triangle(point), 3) - z_min) < tol % TODO: expand to be for more than just z
                bot = [bot; original_points(triangle(point), 1:3), triangle(point), idx];
            elseif abs(original_points(triangle(point), 3) - z_max) < tol
                top = [top; original_points(triangle(point), 1:3), triangle(point), idx];
            end
        end
    end
end
disp("Split into top and bottom");
end


function [top_out, top_in, bot_out, bot_in] = distribute_bins(conn_copy, top, bot, original_points)
% bot/top_in/out are formatted as n x 4 arrays of x, y, z, point idx
% Determine centroids of top and bottom
centroid_top = sum(top(:, 1:3), 1) / size(top(:, 1:3), 1);
centroid_bot = sum(bot(:, 1:3), 1) / size(bot(:, 1:3), 1);
disp("Computed centroids");

% Place each point into its corresponding bin
top_out = [];
top_in = [];
bot_out = [];
bot_in = [];
% compute normal vectors
tr_top = triangulation(conn_copy, original_points);
norm_vecs = faceNormal(tr_top);
disp("Computed normal vectors")

% for top
for idx = 1 : size(top, 1)
    triangle = top(idx, 4);
    % calculate normal vector
    norm_vec = norm_vecs(triangle, :);
    % calculate vector of centroid - pointonrim
    diff_vec = centroid_top - top(idx, 1:3);
    % dot product the two vectors
    % if dot prod < 0, outside
    if (dot(norm_vec, diff_vec) < 0)
        top_out = [top_out; top(idx, 1:4)];
    else
        top_in = [top_in; top(idx, 1:4)];
    end
end

% for bot
for idx = 1 : size(bot, 1)
    triangle = bot(idx, 4);
    % calculate normal vector
    norm_vec = norm_vecs(triangle, :);
    % calculate vector of centroid - pointonrim
    diff_vec = centroid_bot - bot(idx, 1:3);
    % dot product the two vectors
    % if dot prod < 0, outside
    if (dot(norm_vec, diff_vec) < 0)
        bot_out = [bot_out; bot(idx, 1:4)];
    else
        bot_in = [bot_in; bot(idx, 1:4)];
    end
end

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


function [top_out_new, top_in_new, bot_out_new, bot_in_new, points, connectivity] = extend_points(top_out, top_in, bot_out, bot_in, points, connectivity)
% while top/bot_in/out are formatted as n x 1 arrays of point indices,
% top/bot_in/out_new are formatted as n x 3 arrays of point coordinates.
top_out_new = zeros(size(top_out, 1), 3);
top_in_new = zeros(size(top_in, 1), 3);
bot_out_new = zeros(size(bot_out, 1), 3);
bot_in_new = zeros(size(bot_in, 1), 3);

centroid_top_in = sum(top_in(:, 1:3), 1) / size(top_in(:, 1:3), 1);
centroid_top_out = sum(top_out(:, 1:3), 1) / size(top_out(:, 1:3), 1);
centroid_bot_in = sum(bot_in(:, 1:3), 1) / size(bot_in(:, 1:3), 1);
centroid_bot_out = sum(bot_out(:, 1:3), 1) / size(bot_out(:, 1:3), 1);

% Scale by ~10% of diameter
diameter_top_in = finddiameter(top_in, points);
diameter_top_out = finddiameter(top_out, points);
diameter_bot_in = finddiameter(bot_in, points);
diameter_bot_out = finddiameter(bot_out, points);

for i = 1 : size(top_out, 1)
    top_out_new(i, 1:3) = 0.1 * diameter_top_out * (top_out(i, 1:3) - centroid_top_out) + centroid_top_out;
end
top_out_new = [top_out_new, top_out(:, 4)];
for i = 1 : size(top_in, 1)
    top_in_new(i, 1:3) = 0.1 * diameter_top_in * (top_in(i, 1:3) - centroid_top_in) + centroid_top_in;
end
top_in_new = [top_in_new, top_in(:, 4)];
for i = 1 : size(bot_in, 1)
    bot_in_new(i, 1:3) = 0.1 * diameter_bot_in * (bot_in(i, 1:3) - centroid_bot_in) + centroid_bot_in;
end
bot_in_new = [bot_in_new, bot_in(:, 4)];
for i = 1 : size(bot_out, 1)
    bot_out_new(i, 1:3) = 0.1 * diameter_bot_out * (bot_out(i, 1:3) - centroid_bot_out) + centroid_bot_out;
end
bot_out_new = [bot_out_new, bot_out(:, 4)];
end


function [top_out_nbrs, top_in_nbrs, bot_out_nbrs, bot_in_nbrs] = connectpoints (top_out, top_in, bot_out, bot_in, top_out_new, top_in_new, bot_out_new, bot_in_new, points, connectivity)
% get neighbors of the points
pointsneighbors = mappointneighbors(points, connectivity); % all points
top_out_nbrs = cell(size(points, 1), 1);
for i = 1 : size(top_out_new, 1)
    pointidx = top_out_new(i, 4);
    for nbridx = 1 : size(pointsneighbors{pointidx})
        nbr = pointsneighbors{pointidx}(nbridx);
        if ismember(nbr, top_out_new)
            top_out_nbrs{pointidx} = [top_out_nbrs{pointidx}, nbr];
        end
    end
end
top_in_nbrs = cell(size(points, 1), 1);
for i = 1 : size(top_in_new, 1)
    pointidx = top_in_new(i, 4);
    for nbr = 1 : size(pointsneighbors{pointidx})
        if ismember(pointsneighbors{pointidx}(nbr), top_in_new)
            top_in_nbrs{pointidx} = [top_in_nbrs{pointidx}, pointsneighbors{pointidx}(nbr)];
        end
    end
end
bot_out_nbrs = cell(size(points, 1), 1);
for i = 1 : size(bot_out_new, 1)
    pointidx = bot_out_new(i, 4);
    for nbr = 1 : size(pointsneighbors{pointidx})
        if ismember(pointsneighbors{pointidx}(nbr), bot_out_new)
            bot_out_nbrs{pointidx} = [bot_out_nbrs{pointidx}, pointsneighbors{pointidx}(nbr)];
        end
    end
end
bot_in_nbrs = cell(size(points, 1), 1);
for i = 1 : size(bot_in_new, 1)
    pointidx = bot_in_new(i, 4);
    for nbr = 1 : size(pointsneighbors{pointidx})
        if ismember(pointsneighbors{pointidx}(nbr), bot_in_new)
            bot_in_nbrs{pointidx} = [bot_in_nbrs{pointidx}, pointsneighbors{pointidx}(nbr)];
        end
    end
end

% CONNECT TOP OUT
% add points to points list and extend connectivities
points = [points; top_out_new(:, 1:3)];

% base case
numadded = 0; % add points to this when used
point1 = top_out(1, 4); % find first element to use. TODO: change
point2 = top_out_new(1, 4) + size(points, 1); % point corresponding to point1
point3 = top_out_nbrs{point1}(1);
top_out_added = [top_out_added, point1, point2, point3];
connectivity = [connectivity; point1, point2, point3];
point1 = point3;
% find extended neighbor of point 1
for i = 1 : size(top_out_new)
    if top_out_new(i, 4)==point1
        point3 = i + size(points, 1);
        top_out_added = [top_out_added, point3];
        break;
    end
end
connectivity = [connectivity; point1, point2, point3];
numadded = numadded + 4;
% iterative case
while numadded < size(top_out_new, 1) * 2 - 1
    point2 = point3;
    if ismember(top_out_nbrs{point1}(1), top_out_added)
        point3 = top_out_nbrs{point1}(2);
    else
        point3 = top_out_nbrs{point1}(1);
    end
    top_out_added = [top_out_added, point3];
    connectivity = [connectivity; point1, point2, point3];
    point1 = point3;
    % find extended neighbor of point 1
    for i = 1 : size(top_out_new)
        if top_out_new(i, 4)==point1
            point3 = i + size(points, 1);
            break;
        end
    end
    top_out_added = [top_out_added, point3];
    connectivity = [connectivity; point1, point2, point3];
    numadded = numadded + 2;
end
point2 = point3;
point3 = top_out_added(1);
connectivity = [connectivity; point1, point2, point3];
point1 = point3;
point3 = top_out_added(2);
connectivity = [connectivity; point1, point2, point3];


% CONNECT TOP IN
% add points to points list and extend connectivities
points = [points; top_in_new(:, 1:3)];

% base case
numadded = 0; % add points to this when used
point1 = top_in(1, 4); % find first element to use. TODO: change
point2 = top_in_new(1, 4) + size(points, 1); % point corresponding to point1
point3 = top_in_nbrs{point1}(1);
top_in_added = [top_in_added, point1, point2, point3];
connectivity = [connectivity; point1, point2, point3];
point1 = point3;
% find extended neighbor of point 1
for i = 1 : size(top_in_new)
    if top_in_new(i, 4)==point1
        point3 = i + size(points, 1);
        top_in_added = [top_in_added, point3];
        break;
    end
end
connectivity = [connectivity; point1, point2, point3];
numadded = numadded + 4;
% iterative case
while numadded < size(top_in_new, 1) * 2 - 1
    point2 = point3;
    if ismember(top_in_nbrs{point1}(1), top_in_added)
        point3 = top_in_nbrs{point1}(2);
    else
        point3 = top_in_nbrs{point1}(1);
    end
    top_in_added = [top_in_added, point3];
    connectivity = [connectivity; point1, point2, point3];
    point1 = point3;
    % find extended neighbor of point 1
    for i = 1 : size(top_in_new)
        if top_in_new(i, 4)==point1
            point3 = i + size(points, 1);
            break;
        end
    end
    top_in_added = [top_in_added, point3];
    connectivity = [connectivity; point1, point2, point3];
    numadded = numadded + 2;
end
point2 = point3;
point3 = top_in_added(1);
connectivity = [connectivity; point1, point2, point3];
point1 = point3;
point3 = top_in_added(2);
connectivity = [connectivity; point1, point2, point3];


% CONNECT BOT IN
% add points to points list and extend connectivities
points = [points; bot_in_new(:, 1:3)];

% base case
numadded = 0; % add points to this when used
point1 = bot_in(1, 4); % find first element to use. TODO: change
point2 = bot_in_new(1, 4) + size(points, 1); % point corresponding to point1
point3 = bot_in_nbrs{point1}(1);
bot_in_added = [bot_in_added, point1, point2, point3];
connectivity = [connectivity; point1, point2, point3];
point1 = point3;
% find extended neighbor of point 1
for i = 1 : size(bot_in_new)
    if bot_in_new(i, 4)==point1
        point3 = i + size(points, 1);
        bot_in_added = [bot_in_added, point3];
        break;
    end
end
connectivity = [connectivity; point1, point2, point3];
numadded = numadded + 4;
% iterative case
while numadded < size(bot_in_new, 1) * 2 - 1
    point2 = point3;
    if ismember(bot_in_nbrs{point1}(1), bot_in_added)
        point3 = bot_in_nbrs{point1}(2);
    else
        point3 = bot_in_nbrs{point1}(1);
    end
    bot_in_added = [bot_in_added, point3];
    connectivity = [connectivity; point1, point2, point3];
    point1 = point3;
    % find extended neighbor of point 1
    for i = 1 : size(bot_in_new)
        if bot_in_new(i, 4)==point1
            point3 = i + size(points, 1);
            break;
        end
    end
    bot_in_added = [bot_in_added, point3];
    connectivity = [connectivity; point1, point2, point3];
    numadded = numadded + 2;
end
point2 = point3;
point3 = bot_in_added(1);
connectivity = [connectivity; point1, point2, point3];
point1 = point3;
point3 = bot_in_added(2);
connectivity = [connectivity; point1, point2, point3];


% CONNECT BOT OUT
% add points to points list and extend connectivities
points = [points; bot_out_new(:, 1:3)];

% base case
numadded = 0; % add points to this when used
point1 = bot_out(1, 4); % find first element to use. TODO: change
point2 = bot_out_new(1, 4) + size(points, 1); % point corresponding to point1
point3 = bot_out_nbrs{point1}(1);
bot_out_added = [bot_out_added, point1, point2, point3];
connectivity = [connectivity; point1, point2, point3];
point1 = point3;
% find extended neighbor of point 1
for i = 1 : size(bot_out_new)
    if bot_out_new(i, 4)==point1
        point3 = i + size(points, 1);
        bot_out_added = [bot_out_added, point3];
        break;
    end
end
connectivity = [connectivity; point1, point2, point3];
numadded = numadded + 4;
% iterative case
while numadded < size(bot_out_new, 1) * 2 - 1
    point2 = point3;
    if ismember(bot_out_nbrs{point1}(1), bot_out_added)
        point3 = bot_out_nbrs{point1}(2);
    else
        point3 = bot_out_nbrs{point1}(1);
    end
    bot_out_added = [bot_out_added, point3];
    connectivity = [connectivity; point1, point2, point3];
    point1 = point3;
    % find extended neighbor of point 1
    for i = 1 : size(bot_out_new)
        if bot_out_new(i, 4)==point1
            point3 = i + size(points, 1);
            break;
        end
    end
    bot_out_added = [bot_out_added, point3];
    connectivity = [connectivity; point1, point2, point3];
    numadded = numadded + 2;
end
point2 = point3;
point3 = bot_out_added(1);
connectivity = [connectivity; point1, point2, point3];
point1 = point3;
point3 = bot_out_added(2);
connectivity = [connectivity; point1, point2, point3];

    % filter through points neighbors to get just the points
    % connect each point to its neighbors next one:
    % . ._._._._._._._._._.  new points
    %  /|/|/|/|/|/|/|/|/|/|  connections to be added
    % ' ' ' ' ' ' ' ' ' '    original points--connections alr exist

    % TODO: change top/bot_in/out to be array of point indices
end

function [maxdiameter] = finddiameter(bin, points)
% Helper function: Given a set of points in the bin, finds the max diameter.
% The max diameter = the largest distance of (point - point)
% Input: bin of points, formatted as [x, y, z, pointidx]
maxdiameter = 0;
pointidxs = bin(:, 4);
for i = 1 : size(bin, 1)
    for j = 1 : size(bin, 1)
        diameter = norm(points(pointidxs(i), 1:3) - points(pointidxs(j), 1:3));
        maxdiameter = max(maxdiameter, diameter);
    end
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

function [tol] = computetol(points, connectivity)
% Computes a rough value for the tolerance
end
