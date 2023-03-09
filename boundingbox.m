% Read in input .stl file
output = stlread("testvessel.stl");
points = output.Points;
connectivity = output.ConnectivityList;

% bounds is an array [x_min, x_max, y_min, y_max, z_min, z_max]
[x_min, x_max, y_min, y_max, z_min, z_max] = ...
    bounding_box(points, connectivity);
bounds = [x_min, x_max, y_min, y_max, z_min, z_max];
% bounding_box("testvessel.stl", "box.stl", false)

% Create bins. Contains triangle numbers in the original connectivity list.
[x_min_bin, x_max_bin, y_min_bin, y_max_bin, z_min_bin, z_max_bin, testing_centroids] = ...
    create_bins(bounds, points, connectivity);
testing_centroids = sortrows(testing_centroids, 3);

% create extensions
[bot_in, bot_out, top_in, top_out] = create_extensions(points, connectivity, z_min_bin, z_max_bin, z_min, z_max, 0.07) 


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


function [bot_in, bot_out, top_in, top_out] = create_extensions(original_points, original_connectivity, z_min_bin, z_max_bin, z_min, z_max, tol) % This is specific to a vessel that is along the Z axis.
% create extensions
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

% split the triangles into two bins of points: top vs bot
% both of the form [x, y, z, triangleidx]
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
%         str = ["compared triangles ", idx, " and ", test_idx];
%         disp(str)
    end
    % if on rim, identify points on rim and determine if on top or bottom
    if num_edges_shared < 3
        for point = 1 : 3
            if abs(original_points(triangle(point), 3) - z_min) < tol % TODO: expand to be for more than just z
                bot = [bot; original_points(triangle(point)), triangle];
            elseif abs(original_points(triangle(point), 3) - z_max) < tol
                top = [top; original_points(triangle(point)), triangle];
            end
        end
    end
end
disp("Split into top and bottom");

% Determine centroids of top and bottom
centroid_top = sum(top, 1) / size(top, 1);
centroid_bot = sum(bot, 1) / size(bot, 1);
disp("Computed centroids");

% Place each point into its corresponding bin
top_out = [];
top_in = [];
bot_out = [];
bot_in = [];
% compute normal vectors
tr_top = triangulation(top, original_points);
norm_vecs = faceNormal(tr_top);
disp("Computed normal vectors")

% for top
for idx = 1 : size(top, 1)
    triangle = top(idx, 4);
    % calculate normal vector
    norm_vec = norm_vecs(triangle);
    % calculate vector of centroid - pointonrim
    diff_vec = centroid_top - top(idx, 1:3);
    % dot product the two vectors
    % if dot prod < 0, outside
    if (dot(norm_vec, diff_vec) < 0)
        top_out = [top_out; top(idx)];
    else
        top_in = [top_in; top(idx)];
    end
end

% for bot
for idx = 1 : size(top, 1)
    triangle = top(idx, 4);
    % calculate normal vector
    norm_vec = norm_vecs(triangle);
    % calculate vector of centroid - pointonrim
    diff_vec = centroid_bot - bot(idx, 1:3);
    % dot product the two vectors
    % if dot prod < 0, outside
    if (dot(norm_vec, diff_vec) < 0)
        bot_out = [bot_out; top(idx)];
    else
        bot_in = [bot_in; top(idx)];
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

