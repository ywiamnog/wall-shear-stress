% clear; clc; close all;

% Given an input stl and an input centerline file, compute an orthonormal
% frame at each point in the stl.
% Requires:
%   - Input stl and input centerline vtk are aligned properly.
% Inputs:
%   - stl_filename: the name of the input slt file.
%   - centerline_filename: the name of the input centerline vtk file.
% Outputs:
%   - axial, radial, circ: each are nx3 matrices, where matrix(point, :) is
%   the vector for that point.

% centerline_filename = "centerlinepoints.vtk";
% collars_foldername = "collars";

function [all_points, axial, radial, circ] = computeframe(centerline_filename, collars_foldername)

% Read in vtk file.
centerline = vtk_read_centerline(centerline_filename);
points_centerline = centerline.points;
tangents = centerline.tangents; 

% Maps centerline idx : [x,y,z] of collar points, [p1,p2,p3] of connectivity, [radials of collar points], [tangents of collar points].
collars = cell(size(points_centerline, 1), 3); 
% Initialize collars.
for clp = 1 : size(points_centerline, 1)
    collars_filename = strcat(collars_foldername, "/", num2str(clp - 1), ".vtk");
    collar_TR = read_collar_from_vtk(collars_filename);
    collars{clp, 1} = collar_TR.points; % Points.
    collars{clp, 2} = collar_TR.connectivity; % Connectivity.
    collars{clp, 3} = zeros(size(collars{clp, 1}, 1), 4); % Radial vector (normal vector).
    collars{clp, 4} = zeros(size(collars{clp, 1}, 1), 4); % Axial vector (tangent vector).
    collars{clp, 5} = zeros(size(collars{clp, 1}, 1), 3); % Circumferential vector.
end

% Add normal vectors.
for clp = 1 : size(points_centerline, 1)
    collar_points = collars{clp, 1};
    collar_connectivity = collars{clp, 2};
    TR = triangulation(collar_connectivity, collar_points);
    N = faceNormal(TR);
    radial_aggregate_temp = zeros(size(collar_points, 1), 4);
    for triangle = 1 : size(collar_connectivity, 1)
        point1 = collar_points(collar_connectivity(triangle, 1), 1:3);
        point2 = collar_points(collar_connectivity(triangle, 2), 1:3);
        point3 = collar_points(collar_connectivity(triangle, 3), 1:3);
        area = area_of_triangle(point1, point2, point3);
        % Weighting process
        for point = 1 : 3
            radial_aggregate_temp(collar_connectivity(triangle, point), 1:3) = radial_aggregate_temp(collar_connectivity(triangle, point), 1:3) + N(triangle, :) * area;
            radial_aggregate_temp(collar_connectivity(triangle, point), 4) = radial_aggregate_temp(collar_connectivity(triangle, point), 4) + area;
        end
    end
    collars{clp, 3} = collars{clp, 3} + radial_aggregate_temp;
end

% Normalize normal vectors
for clp = 1 : size(points_centerline, 1)
    for point = 1 : size(collars{clp, 1}, 1)
        collars{clp, 3}(point, 1:3) = collars{clp, 3}(point, 1:3) / collars{clp, 3}(point, 4);
        collars{clp, 3}(point, 1:3) = collars{clp, 3}(point, 1:3) / norm(collars{clp, 3}(point, 1:3));
    end
end


% Add tangent vectors.
for clp = 1 : size(points_centerline, 1)
    collar_points = collars{clp, 1};
    tangent = tangents(clp, :);
    % Add tangent to each point in this collar.
    for point = 1 : size(collar_points, 1)
        collars{clp, 4}(point, 1:3) = collars{clp, 4}(point, 1:3) + tangent;
        collars{clp, 4}(point, 4) = collars{clp, 4}(point, 4) + 1;
        % Go through other collars and try to find these points.
        for other_clp = 1 : size(points_centerline, 1)
            % Skip if is current centerline point.
            if (other_clp == clp) 
                continue;
            end
            % Search for point in other collar.
            for other_point = 1 : size(collars{other_clp, 1}, 1)
                if (collars{other_clp, 1}(other_point, :) ~= collar_points(point, :))
                    continue;
                end
                % Add tangent to point in other collar.
                collars{other_clp, 4}(other_point, 1:3) = collars{other_clp, 4}(other_point, 1:3) + tangent;
                collars{other_clp, 4}(other_point, 4) = collars{other_clp, 4}(other_point, 4) + 1;
            end
        end
    end
end

% Normalize tangent vectors.
for clp = 1 : size(points_centerline, 1)
    for point = 1 : size(collars{clp, 1}, 1)
        collars{clp, 4}(point, 1:3) = collars{clp, 4}(point, 1:3) / collars{clp, 4}(point, 4);
        collars{clp, 4}(point, 1:3) = collars{clp, 4}(point, 1:3) / norm(collars{clp, 4}(point, 1:3));
        % Compute circumferential vectors.
    end
end 

% % Compute circumferential vectors.
% for clp = 1 : size(points_centerline, 1)
%     for point = 1 : size(collars{clp, 1}, 1)
%         collars{clp, 5}(point, :) = cross(collars{clp, 3}(point, 1:3), collars{clp, 4}(point, 1:3));
%         % Normalize.
%         collars{clp, 5}(point, :) = collars{clp, 5}(point, :) / norm(collars{clp, 5}(point, :));
%     end
% end

% Consolidate all data into individual matrices.
% All of their row indices correspond with each other.
all_points = [];
radial = [];
axial = [];
circ = [];
for clp = 1 : size(points_centerline, 1)
    all_points = [all_points; collars{clp, 1}(:, 1:3)];
%     radial = [radial; collars{clp, 3}(:, 1:3)];
    axial = [axial; collars{clp, 4}(:, 1:3)];
%     circ = [circ; collars{clp, 5}(:, 1:3)];
end


%% Write out to VTK.
% Write data to a VTK file.
filename_out = "frame.vtk";
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
fprintf(fid, "POINTS %d float\n", size(all_points, 1));
for point = 1 : size(all_points, 1)
    fprintf(fid, "%f %f %f\n", all_points(point, 1), all_points(point, 2), all_points(point, 3));
end

% % Write connectivites.
% fprintf(fid, "POLYGONS %d %d\n", size(connectivity, 1), size(connectivity, 1) * 4);
% for triangle = 1 : size(connectivity, 1)
%     fprintf(fid, "3 %d %d %d\n", connectivity(triangle, 1)-1, connectivity(triangle, 2)-1, connectivity(triangle, 3)-1);
% end

fprintf(fid, "POINT_DATA %d\n", size(all_points, 1));

% % Write vectors.
% fprintf(fid, "VECTORS radial_vector_unit float\n");
% for point = 1 : size(all_points, 1)
%     fprintf(fid, "%f %f %f\n", radial(point, 1), radial(point, 2), radial(point, 3));
% end
fprintf(fid, "VECTORS axial_vector_unit float\n");
for point = 1 : size(all_points, 1)
    fprintf(fid, "%f %f %f\n", axial(point, 1), axial(point, 2), axial(point, 3));
end
% fprintf(fid, "VECTORS circ_vector_unit float\n");
% for point = 1 : size(all_points, 1)
%     fprintf(fid, "%f %f %f\n", circ(point, 1), circ(point, 2), circ(point, 3));
% end

fclose(fid);



%% Helper functions
% Compute area of a triangle given by 3 points (x,y,z)
    function area = area_of_triangle(point1, point2, point3)
        % Inputs:
        %   - point1,2,3: each point of the triangle, given by (x,y,z)
        vec1 = point2 - point1;
        vec2 = point3 - point1;
        area = norm(cross(vec1, vec2)) / 2;
    end


    function centerline = vtk_read_centerline(filename)
        % Reads from a VTK file.
        % Requires the VTK file to have only one POINTS section detailing the
        % points that the centerline consists of. Everything else will be
        % discarded.
        %% Read multi-data VTK file
        fid = fopen(filename,'r');
        % initialize variables
        centerline.points = [];
        centerline.tangents = [];
        % flag for searching a title line to find where the numbers are
        flag1 = 0;
        while ~feof(fid)
            % read header
            str = fgets(fid);
            str = strip(str);
            if (strlength(str) >= 6 && strcmp(str(1:6), "POINTS"))
                separate = split(str);
                num_points = str2double(separate{2});
                centerline.points = zeros(num_points, 3);
                for i = 1 : num_points
                    str = fgets(fid);
                    str = strip(str);
                    separate = split(str);
                    centerline.points(i, 1) = str2double(separate{1});
                    centerline.points(i, 2) = str2double(separate{2});
                    centerline.points(i, 3) = str2double(separate{3});
                end
            end
        end
        fclose(fid);

        % Compute centerline tangents.
        centerline.tangents = zeros(size(centerline.points, 1), 3);
        for i = 2 : size(centerline.points, 1) - 1
            centerline.tangents(i, :) = ((centerline.points(i+1, :) - centerline.points(i-1, :))) / norm(((centerline.points(i+1, :) - centerline.points(i-1, :))));
        end
        centerline.tangents(1, :) = centerline.tangents(2, :);
        centerline.tangents(size(centerline.tangents, 1), :) = centerline.tangents(size(centerline.tangents, 1)-1, :);

    end


    function collar = read_collar_from_vtk(filename)
        % The number of collars should be equal to the number of centerline points.
        % Inputs:
        %   - filename: the name of the vtk file we are reading, the file contains
        %   one collar
        %   - centerline: the centerline structure that this collar belongs to.
        %   - centerline_idx: the index of the centerline that this collar
        %   belongs to.

        %% Read multi-data VTK file
        fid = fopen(filename,'r');

        % initialize variables
        collar.points = [];
        collar.connectivity = [];
        while ~feof(fid)
            str = fgets(fid);
            if (strlength(str) >= 6 && strcmp(str(1:6), "POINTS"))
                separate = split(str);
                num_points = str2double(separate{2});
                collar.points = zeros(num_points, 3);
                for i = 1 : num_points
                    str = fgets(fid);
                    separate = split(str);
                    collar.points(i, 1) = str2double(separate{1});
                    collar.points(i, 2) = str2double(separate{2});
                    collar.points(i, 3) = str2double(separate{3});
                end
            elseif (strlength(str) >= 8 && strcmp(str(1:8), "POLYGONS"))
                separate = split(str);
                num_triangles = str2double(separate{2});
                collar.connectivity = zeros(num_triangles, 3);
                for i = 1 : num_triangles
                    str = fgets(fid);
                    separate = split(str);
                    collar.connectivity(i, 1) = str2double(separate{2});
                    collar.connectivity(i, 2) = str2double(separate{3});
                    collar.connectivity(i, 3) = str2double(separate{4});
                end
            end
        end
        fclose(fid);

    end

end

