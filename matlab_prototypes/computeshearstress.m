% clear; % Clear workspace.
clc; % Clear command window.
close all;

% Given an input VTK file of velocity data, computes shear stress at each
% voxel
% Inputs:
%   - velocity_filename: the name of the VTK file that is to be read.
% Outputs:
%   - E: a 3D matrix, where each element stores a 3x3 matrix of gradients.
%   These gradient matrices will later be used to interpolate onto the
%   actual data. 

% Inputs:
velocity_filename = "0.vtk";
% stl_filename = "test_cylinder_long.stl";
centerline_filename = "centerlinepoints_ascending_aorta.vtk";
collars_foldername= "collars_ascending_aorta";
mu = 1; % Viscosity.

% Read in VTK file.
vtk_data = vtk_read(velocity_filename);
translation = vtk_data.translation;
transformation = vtk_data.transformation;
num_points = vtk_data.num_points;
% points = vtk_data.points;
dimensions = vtk_data.dimensions;
velocity = vtk_data.velocity;

% Number of nodes in each direction.
num_x = vtk_data.dimensions(1);
num_y = vtk_data.dimensions(2);
num_z = vtk_data.dimensions(3);
% Step between nodes in each direction.
hx = 1 / (num_x - 1);
hy = 1 / (num_y - 1);
hz = 1 / (num_z - 1);
% Vectors representing axis ticks along each direction.
x = linspace(0, num_x-1, num_x);
y = linspace(0, num_y-1, num_y);
z = linspace(0, num_z-1, num_z);


%% Calculate gradients.
% Calculate partial derivatives.
du1dx = zeros(num_x, num_y, num_z);
du1dy = zeros(num_x, num_y, num_z);
du1dz = zeros(num_x, num_y, num_z);
du2dx = zeros(num_x, num_y, num_z);
du2dy = zeros(num_x, num_y, num_z);
du2dz = zeros(num_x, num_y, num_z);
du3dx = zeros(num_x, num_y, num_z);
du3dy = zeros(num_x, num_y, num_z);
du3dz = zeros(num_x, num_y, num_z);
% For general points
for i = 2 : num_x - 1
    for j = 2 : num_y - 1
        for k = 2 : num_z - 1
            du1dx(i, j, k) = (velocity{i+1, j, k}(1) - velocity{i-1, j, k}(1)) / (2 * hx);
            du1dy(i, j, k) = (velocity{i, j+1, k}(1) - velocity{i, j-1, k}(1)) / (2 * hy);
            du1dz(i, j, k) = (velocity{i, j, k+1}(1) - velocity{i, j, k-1}(1)) / (2 * hz);
            du2dx(i, j, k) = (velocity{i+1, j, k}(2) - velocity{i-1, j, k}(2)) / (2 * hx);
            du2dy(i, j, k) = (velocity{i, j+1, k}(2) - velocity{i, j-1, k}(2)) / (2 * hy);
            du2dz(i, j, k) = (velocity{i, j, k+1}(2) - velocity{i, j, k-1}(2)) / (2 * hz);
            du3dx(i, j, k) = (velocity{i+1, j, k}(3) - velocity{i-1, j, k}(3)) / (2 * hx);
            du3dy(i, j, k) = (velocity{i, j+1, k}(3) - velocity{i, j-1, k}(3)) / (2 * hy);
            du3dz(i, j, k) = (velocity{i, j, k+1}(3) - velocity{i, j, k-1}(3)) / (2 * hz);
        end
    end
end

% % todo...TEMPORARILY: ALL EDGE CASES ARE JUST 0.
% % Edge cases. 
% % For each face voxel, I chose to use the value from the voxel next to it.
% % For x-faces.
% temp_x = [1, num_x];
% temp_y = [1, num_y];
% temp_z = [1, num_z];
% for j = 2 : num_y - 1
%     for k = 2 : num_z - 1
%         du1dx(1, j, k) = du1dx(2, j, k);
%         du1dy(1, j, k) = du1dy(2, j, k);
%         du1dz(1, j, k) = du1dz(2, j, k);
%         du2dx(1, j, k) = du2dx(2, j, k);
%         du2dy(1, j, k) = du2dy(2, j, k);
%         du2dz(1, j, k) = du2dz(2, j, k);
%         du3dx(1, j, k) = du3dx(2, j, k);
%         du3dy(1, j, k) = du3dy(2, j, k);
%         du3dz(1, j, k) = du3dz(2, j, k);
%         du1dx(num_x, j, k) = du1dx(num_x-1, j, k);
%         du1dy(num_x, j, k) = du1dy(num_x-1, j, k);
%         du1dz(num_x, j, k) = du1dz(num_x-1, j, k);
%         du2dx(num_x, j, k) = du2dx(num_x-1, j, k);
%         du2dy(num_x, j, k) = du2dy(num_x-1, j, k);
%         du2dz(num_x, j, k) = du2dz(num_x-1, j, k);
%         du3dx(num_x, j, k) = du3dx(num_x-1, j, k);
%         du3dy(num_x, j, k) = du3dy(num_x-1, j, k);
%         du3dz(num_x, j, k) = du3dz(num_x-1, j, k);
%     end
% end
% 
% % For y-faces.
% for i = 2 : num_x - 1
%     for k = 2 : num_z - 1
%         du1dx(i, 1, k) = du1dx(i, 2, k);
%         du1dy(i, 1, k) = du1dy(i, 2, k);
%         du1dz(i, 1, k) = du1dz(i, 2, k);
%         du2dx(i, 1, k) = du2dx(i, 2, k);
%         du2dy(i, 1, k) = du2dy(i, 2, k);
%         du2dz(i, 1, k) = du2dz(i, 2, k);
%         du3dx(i, 1, k) = du3dx(i, 2, k);
%         du3dy(i, 1, k) = du3dy(i, 2, k);
%         du3dz(i, 1, k) = du3dz(i, 2, k);
%         du1dx(i, num_y, k) = du1dx(i, num_y-1, k);
%         du1dy(i, num_y, k) = du1dy(i, num_y-1, k);
%         du1dz(i, num_y, k) = du1dz(i, num_y-1, k);
%         du2dx(i, num_y, k) = du2dx(i, num_y-1, k);
%         du2dy(i, num_y, k) = du2dy(i, num_y-1, k);
%         du2dz(i, num_y, k) = du2dz(i, num_y-1, k);
%         du3dx(i, num_y, k) = du3dx(i, num_y-1, k);
%         du3dy(i, num_y, k) = du3dy(i, num_y-1, k);
%         du3dz(i, num_y, k) = du3dz(i, num_y-1, k);
%     end
% end
% 
% % For z-faces.
% for i = 2 : num_x - 1
%     for j = 2 : num_y - 1
%         du1dx(i, j, 1) = du1dx(i, j, 2);
%         du1dy(i, j, 1) = du1dy(i, j, 2);
%         du1dz(i, j, 1) = du1dz(i, j, 2);
%         du2dx(i, j, 1) = du2dx(i, j, 2);
%         du2dy(i, j, 1) = du2dy(i, j, 2);
%         du2dz(i, j, 1) = du2dz(i, j, 2);
%         du3dx(i, j, 1) = du3dx(i, j, 2);
%         du3dy(i, j, 1) = du3dy(i, j, 2);
%         du3dz(i, j, 1) = du3dz(i, j, 2);
%         du1dx(i, j, num_z) = du1dx(i, j, num_z-1);
%         du1dy(i, j, num_z) = du1dy(i, j, num_z-1);
%         du1dz(i, j, num_z) = du1dz(i, j, num_z-1);
%         du2dx(i, j, num_z) = du2dx(i, j, num_z-1);
%         du2dy(i, j, num_z) = du2dy(i, j, num_z-1);
%         du2dz(i, j, num_z) = du2dz(i, j, num_z-1);
%         du3dx(i, j, num_z) = du3dx(i, j, num_z-1);
%         du3dy(i, j, num_z) = du3dy(i, j, num_z-1);
%         du3dz(i, j, num_z) = du3dz(i, j, num_z-1);
%     end
% end
% 
% % For each edge voxel, I chose to use the values from the two face voxels
% % next to it.
% % For xy-edges.
% for k = 2 : num_z - 1
%     du1dx(1, 1, k) = (du1dx(1, 2, k) + du1dx(2, 1, k)) / 2;
%     du1dy(1, 1, k) = (du1dx(1, 2, k) + du1dx(2, 1, k)) / 2;
%     du1dz(1, 1, k) = (du1dx(1, 2, k) + du1dx(2, 1, k)) / 2;
%     du2dx(1, 1, k) = (du1dx(1, 2, k) + du1dx(2, 1, k)) / 2;
%     du2dy(1, 1, k) = (du1dx(1, 2, k) + du1dx(2, 1, k)) / 2;
%     du2dz(1, 1, k) = (du1dx(1, 2, k) + du1dx(2, 1, k)) / 2;
%     du3dx(1, 1, k) = (du1dx(1, 2, k) + du1dx(2, 1, k)) / 2;
%     du3dy(1, 1, k) = (du1dx(1, 2, k) + du1dx(2, 1, k)) / 2;
%     du3dz(1, 1, k) = (du1dx(1, 2, k) + du1dx(2, 1, k)) / 2;
%     du1dx(1, num_y, k) = (du1dx(1, num_y-1, k) + du1dx(2, num_y, k)) / 2;
%     du1dy(1, num_y, k) = (du1dx(1, num_y-1, k) + du1dx(2, num_y, k)) / 2;
%     du1dz(1, num_y, k) = (du1dx(1, num_y-1, k) + du1dx(2, num_y, k)) / 2;
%     du2dx(1, num_y, k) = (du1dx(1, num_y-1, k) + du1dx(2, num_y, k)) / 2;
%     du2dy(1, num_y, k) = (du1dx(1, num_y-1, k) + du1dx(2, num_y, k)) / 2;
%     du2dz(1, num_y, k) = (du1dx(1, num_y-1, k) + du1dx(2, num_y, k)) / 2;
%     du3dx(1, num_y, k) = (du1dx(1, num_y-1, k) + du1dx(2, num_y, k)) / 2;
%     du3dy(1, num_y, k) = (du1dx(1, num_y-1, k) + du1dx(2, num_y, k)) / 2;
%     du3dz(1, num_y, k) = (du1dx(1, num_y-1, k) + du1dx(2, num_y, k)) / 2;
%     du1dx(num_x, 1, k) = (du1dx(num_x, 2, k) + du1dx(num_x-1, 1, k)) / 2;
%     du1dy(num_x, 1, k) = (du1dx(num_x, 2, k) + du1dx(num_x-1, 1, k)) / 2;
%     du1dz(num_x, 1, k) = (du1dx(num_x, 2, k) + du1dx(num_x-1, 1, k)) / 2;
%     du2dx(num_x, 1, k) = (du1dx(num_x, 2, k) + du1dx(num_x-1, 1, k)) / 2;
%     du2dy(num_x, 1, k) = (du1dx(num_x, 2, k) + du1dx(num_x-1, 1, k)) / 2;
%     du2dz(num_x, 1, k) = (du1dx(num_x, 2, k) + du1dx(num_x-1, 1, k)) / 2;
%     du3dx(num_x, 1, k) = (du1dx(num_x, 2, k) + du1dx(num_x-1, 1, k)) / 2;
%     du3dy(num_x, 1, k) = (du1dx(num_x, 2, k) + du1dx(num_x-1, 1, k)) / 2;
%     du3dz(num_x, 1, k) = (du1dx(num_x, 2, k) + du1dx(num_x-1, 1, k)) / 2;
%     du1dx(num_x, num_y, k) = (du1dx(num_x, num_y-1, k) + du1dx(num_x-1, num_y, k)) / 2;
%     du1dy(num_x, num_y, k) = (du1dx(num_x, num_y-1, k) + du1dx(num_x-1, num_y, k)) / 2;
%     du1dz(num_x, num_y, k) = (du1dx(num_x, num_y-1, k) + du1dx(num_x-1, num_y, k)) / 2;
%     du2dx(num_x, num_y, k) = (du1dx(num_x, num_y-1, k) + du1dx(num_x-1, num_y, k)) / 2;
%     du2dy(num_x, num_y, k) = (du1dx(num_x, num_y-1, k) + du1dx(num_x-1, num_y, k)) / 2;
%     du2dz(num_x, num_y, k) = (du1dx(num_x, num_y-1, k) + du1dx(num_x-1, num_y, k)) / 2;
%     du3dx(num_x, num_y, k) = (du1dx(num_x, num_y-1, k) + du1dx(num_x-1, num_y, k)) / 2;
%     du3dy(num_x, num_y, k) = (du1dx(num_x, num_y-1, k) + du1dx(num_x-1, num_y, k)) / 2;
%     du3dz(num_x, num_y, k) = (du1dx(num_x, num_y-1, k) + du1dx(num_x-1, num_y, k)) / 2;
% end
% 
% % For xz-edges.
% for j = 2 : num_y - 1
%     du1dx(1, j, 1) = (du1dx(1, j, 2) + du1dx(2, j, 1)) / 2;
%     du1dy(1, j, 1) = (du1dx(1, j, 2) + du1dx(2, j, 1)) / 2;
%     du1dz(1, j, 1) = (du1dx(1, j, 2) + du1dx(2, j, 1)) / 2;
%     du2dx(1, j, 1) = (du1dx(1, j, 2) + du1dx(2, j, 1)) / 2;
%     du2dy(1, j, 1) = (du1dx(1, j, 2) + du1dx(2, j, 1)) / 2;
%     du2dz(1, j, 1) = (du1dx(1, j, 2) + du1dx(2, j, 1)) / 2;
%     du3dx(1, j, 1) = (du1dx(1, j, 2) + du1dx(2, j, 1)) / 2;
%     du3dy(1, j, 1) = (du1dx(1, j, 2) + du1dx(2, j, 1)) / 2;
%     du3dz(1, j, 1) = (du1dx(1, j, 2) + du1dx(2, j, 1)) / 2;
%     du1dx(1, j, num_z) = (du1dx(1, j, num_z-1) + du1dx(2, j, num_z)) / 2;
%     du1dy(1, j, num_z) = (du1dx(1, j, num_z-1) + du1dx(2, j, num_z)) / 2;
%     du1dz(1, j, num_z) = (du1dx(1, j, num_z-1) + du1dx(2, j, num_z)) / 2;
%     du2dx(1, j, num_z) = (du1dx(1, j, num_z-1) + du1dx(2, j, num_z)) / 2;
%     du2dy(1, j, num_z) = (du1dx(1, j, num_z-1) + du1dx(2, j, num_z)) / 2;
%     du2dz(1, j, num_z) = (du1dx(1, j, num_z-1) + du1dx(2, j, num_z)) / 2;
%     du3dx(1, j, num_z) = (du1dx(1, j, num_z-1) + du1dx(2, j, num_z)) / 2;
%     du3dy(1, j, num_z) = (du1dx(1, j, num_z-1) + du1dx(2, j, num_z)) / 2;
%     du3dz(1, j, num_z) = (du1dx(1, j, num_z-1) + du1dx(2, j, num_z)) / 2;
%     du1dx(num_x, j, 1) = (du1dx(num_x, j, 2) + du1dx(num_x-1, j, 1)) / 2;
%     du1dy(num_x, j, 1) = (du1dx(num_x, j, 2) + du1dx(num_x-1, j, 1)) / 2;
%     du1dz(num_x, j, 1) = (du1dx(num_x, j, 2) + du1dx(num_x-1, j, 1)) / 2;
%     du2dx(num_x, j, 1) = (du1dx(num_x, j, 2) + du1dx(num_x-1, j, 1)) / 2;
%     du2dy(num_x, j, 1) = (du1dx(num_x, j, 2) + du1dx(num_x-1, j, 1)) / 2;
%     du2dz(num_x, j, 1) = (du1dx(num_x, j, 2) + du1dx(num_x-1, j, 1)) / 2;
%     du3dx(num_x, j, 1) = (du1dx(num_x, j, 2) + du1dx(num_x-1, j, 1)) / 2;
%     du3dy(num_x, j, 1) = (du1dx(num_x, j, 2) + du1dx(num_x-1, j, 1)) / 2;
%     du3dz(num_x, j, 1) = (du1dx(num_x, j, 2) + du1dx(num_x-1, j, 1)) / 2;
%     du1dx(num_x, j, num_z) = (du1dx(num_x, j, num_z-1) + du1dx(num_x-1, j, num_z)) / 2;
%     du1dy(num_x, j, num_z) = (du1dx(num_x, j, num_z-1) + du1dx(num_x-1, j, num_z)) / 2;
%     du1dz(num_x, j, num_z) = (du1dx(num_x, j, num_z-1) + du1dx(num_x-1, j, num_z)) / 2;
%     du2dx(num_x, j, num_z) = (du1dx(num_x, j, num_z-1) + du1dx(num_x-1, j, num_z)) / 2;
%     du2dy(num_x, j, num_z) = (du1dx(num_x, j, num_z-1) + du1dx(num_x-1, j, num_z)) / 2;
%     du2dz(num_x, j, num_z) = (du1dx(num_x, j, num_z-1) + du1dx(num_x-1, j, num_z)) / 2;
%     du3dx(num_x, j, num_z) = (du1dx(num_x, j, num_z-1) + du1dx(num_x-1, j, num_z)) / 2;
%     du3dy(num_x, j, num_z) = (du1dx(num_x, j, num_z-1) + du1dx(num_x-1, j, num_z)) / 2;
%     du3dz(num_x, j, num_z) = (du1dx(num_x, j, num_z-1) + du1dx(num_x-1, j, num_z)) / 2;
% end
% 
% % For yz-edges.
% for i = 2 : num_x - 1
%     du1dx(i, 1, 1) = (du1dx(i, 1, 2) + du1dx(i, 2, 1)) / 2;
%     du1dy(i, 1, 1) = (du1dx(i, 1, 2) + du1dx(i, 2, 1)) / 2;
%     du1dz(i, 1, 1) = (du1dx(i, 1, 2) + du1dx(i, 2, 1)) / 2;
%     du2dx(i, 1, 1) = (du1dx(i, 1, 2) + du1dx(i, 2, 1)) / 2;
%     du2dy(i, 1, 1) = (du1dx(i, 1, 2) + du1dx(i, 2, 1)) / 2;
%     du2dz(i, 1, 1) = (du1dx(i, 1, 2) + du1dx(i, 2, 1)) / 2;
%     du3dx(i, 1, 1) = (du1dx(i, 1, 2) + du1dx(i, 2, 1)) / 2;
%     du3dy(i, 1, 1) = (du1dx(i, 1, 2) + du1dx(i, 2, 1)) / 2;
%     du3dz(i, 1, 1) = (du1dx(i, 1, 2) + du1dx(i, 2, 1)) / 2;
%     du1dx(i, 1, num_z) = (du1dx(i, 1, num_z-1) + du1dx(i, 2, num_z)) / 2;
%     du1dy(i, 1, num_z) = (du1dx(i, 1, num_z-1) + du1dx(i, 2, num_z)) / 2;
%     du1dz(i, 1, num_z) = (du1dx(i, 1, num_z-1) + du1dx(i, 2, num_z)) / 2;
%     du2dx(i, 1, num_z) = (du1dx(i, 1, num_z-1) + du1dx(i, 2, num_z)) / 2;
%     du2dy(i, 1, num_z) = (du1dx(i, 1, num_z-1) + du1dx(i, 2, num_z)) / 2;
%     du2dz(i, 1, num_z) = (du1dx(i, 1, num_z-1) + du1dx(i, 2, num_z)) / 2;
%     du3dx(i, 1, num_z) = (du1dx(i, 1, num_z-1) + du1dx(i, 2, num_z)) / 2;
%     du3dy(i, 1, num_z) = (du1dx(i, 1, num_z-1) + du1dx(i, 2, num_z)) / 2;
%     du3dz(i, 1, num_z) = (du1dx(i, 1, num_z-1) + du1dx(i, 2, num_z)) / 2;
%     du1dx(i, num_y, 1) = (du1dx(i, num_y, 2) + du1dx(i, num_y-1, 1)) / 2;
%     du1dy(i, num_y, 1) = (du1dx(i, num_y, 2) + du1dx(i, num_y-1, 1)) / 2;
%     du1dz(i, num_y, 1) = (du1dx(i, num_y, 2) + du1dx(i, num_y-1, 1)) / 2;
%     du2dx(i, num_y, 1) = (du1dx(i, num_y, 2) + du1dx(i, num_y-1, 1)) / 2;
%     du2dy(i, num_y, 1) = (du1dx(i, num_y, 2) + du1dx(i, num_y-1, 1)) / 2;
%     du2dz(i, num_y, 1) = (du1dx(i, num_y, 2) + du1dx(i, num_y-1, 1)) / 2;
%     du3dx(i, num_y, 1) = (du1dx(i, num_y, 2) + du1dx(i, num_y-1, 1)) / 2;
%     du3dy(i, num_y, 1) = (du1dx(i, num_y, 2) + du1dx(i, num_y-1, 1)) / 2;
%     du3dz(i, num_y, 1) = (du1dx(i, num_y, 2) + du1dx(i, num_y-1, 1)) / 2;
%     du1dx(i, num_y, num_z) = (du1dx(i, num_y, num_z-1) + du1dx(i, num_y-1, num_z)) / 2;
%     du1dy(i, num_y, num_z) = (du1dx(i, num_y, num_z-1) + du1dx(i, num_y-1, num_z)) / 2;
%     du1dz(i, num_y, num_z) = (du1dx(i, num_y, num_z-1) + du1dx(i, num_y-1, num_z)) / 2;
%     du2dx(i, num_y, num_z) = (du1dx(i, num_y, num_z-1) + du1dx(i, num_y-1, num_z)) / 2;
%     du2dy(i, num_y, num_z) = (du1dx(i, num_y, num_z-1) + du1dx(i, num_y-1, num_z)) / 2;
%     du2dz(i, num_y, num_z) = (du1dx(i, num_y, num_z-1) + du1dx(i, num_y-1, num_z)) / 2;
%     du3dx(i, num_y, num_z) = (du1dx(i, num_y, num_z-1) + du1dx(i, num_y-1, num_z)) / 2;
%     du3dy(i, num_y, num_z) = (du1dx(i, num_y, num_z-1) + du1dx(i, num_y-1, num_z)) / 2;
%     du3dz(i, num_y, num_z) = (du1dx(i, num_y, num_z-1) + du1dx(i, num_y-1, num_z)) / 2;
% end


% For each corner voxel, I chose to use the values from the three edge
% voxels next to it.


% Put all of the individual gradients together.
gradients = cell(num_x, num_y, num_z);
for i = 1 : num_x
    for j = 1 : num_y
        for k = 1 : num_z
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
        end
    end
end

% Compute E on the Cartesian grid.
E = cell(num_x, num_y, num_z); % E is nxnxn, each elem storing a 3x3 -- 5-D.
for i = 2 : num_x-1 % NOTE: same bounds as gradient.
    for j = 2 : num_y-1
        for k = 2 : num_z-1
            E{i, j, k} = mu * gradients{i, j, k} + (gradients{i, j, k})'; % 3x3.
        end
    end
end
% Edge cases. 
% For each face voxel, I chose to use the value from the voxel next to it.
% For x-faces.
% temp_x = [1, num_x];
% temp_y = [1, num_y];
% temp_z = [1, num_z];
for j = 2 : num_y - 1
    for k = 2 : num_z - 1
        E{1, j, k} = E{2, j, k};
        E{num_x, j, k} = E{num_x - 1, j, k};
    end
end

% For y-faces.
for i = 2 : num_x - 1
    for k = 2 : num_z - 1
        E{i, 1, k} = E{i, 2, k};
        E{i, num_y, k} = E{i, num_y-1, k};
    end
end

% For z-faces.
for i = 2 : num_x - 1
    for j = 2 : num_y - 1
        E{i, j, 1} = E{i, j, 2};
        E{i, j, num_z} = E{i, j, num_z-1};
    end
end

% For each edge voxel, I chose to use the values from the two face voxels
% next to it.
% For xy-edges.
for k = 2 : num_z - 1
    E{1, 1, k} = (E{1, 2, k} + E{2, 1, k}) / 2;
    E{1, num_y, k} = (E{1, num_y-1, k} + E{2, num_y, k}) / 2;
    E{num_x, 1, k} = (E{num_x, 2, k} + E{num_x-1, 1, k}) / 2;
    E{num_x, num_y, k} = (E{num_x, num_y-1, k} + E{num_x-1, num_y, k}) / 2;
end

% For xz-edges.
for j = 2 : num_y - 1
    E{1, j, 1} = (E{1, j, 2} + E{2, j, 1}) / 2;
    E{1, j, num_z} = (E{1, j, num_z-1} + E{2, j, num_z}) / 2;
    E{num_x, j, 1} = (E{num_x, j, 2} + E{num_x-1, j, 1}) / 2;
    E{num_x, j, num_z} = (E{num_x, j, num_z-1} + E{num_x-1, j, num_z}) / 2;
end

% For yz-edges.
for i = 2 : num_x - 1
    E{i, 1, 1} = (E{i, 1, 2} + E{i, 2, 1}) / 2;
    E{i, 1, num_z} = (E{i, 1, num_z-1} + E{i, 2, num_z}) / 2;
    E{i, num_y, 1} = (E{i, num_y, 2} + E{i, num_y-1, 1}) / 2;
    E{i, num_y, num_z} = (E{i, num_y, num_z-1} + E{i, num_y-1, num_z}) / 2;
end

% For corners.
E{1, 1, 1} = (E{1, 1, 2} + E{1, 2, 1} + E{2, 1, 1}) / 3;
E{1, 1, num_z} = (E{1, 1, num_z-1} + E{1, 2, num_z} + E{2, 1, num_z}) / 3;
E{1, num_y, 1} = (E{1, num_y, 2} + E{1, num_y-1, 1} + E{2, num_y, 1}) / 3;
E{num_x, 1, 1} = (E{num_x, 1, 2} + E{num_x, 2, 1} + E{num_x-1, 1, 1}) / 3;
E{num_x, num_y, 1} = (E{num_x, num_y, 2} + E{num_x, num_y-1, 1} + E{num_x-1, num_y, 1}) / 3;
E{num_x, 1, num_z} = (E{num_x, 1, num_z-1} + E{num_x, 2, num_z} + E{num_x-1, 1, num_z}) / 3;
E{1, num_y, num_z} = (E{1, num_y, num_z-1} + E{1, num_y-1, num_z} + E{2, num_y, num_z}) / 3;
E{num_x, num_y, num_z} = (E{num_x, num_y, num_z-1} + E{num_x, num_y-1, num_z} + E{num_x-1, num_y, num_z}) / 3;

%% Compute an orthonormal frame at each point of the mesh data.
[all_points, axial, radial, circ] = computeframe(centerline_filename, collars_foldername);

% Modify points_mesh so that all points are also axis aligned.
%   NOTE: we previously used the regular Cartesian grid (axis aligned, 
%   based at [0,0,0]). The actual data is not this regular, and a 
%   translation and transformation are necessary to successfully binary
%   search for it.


%% Interpolate E and calculate w = E*r at each point on the mesh.
w_interpolated = zeros(size(all_points, 1), 3); % w = E * radial_unit. This is interpolated.
E_int = [0,0,0; 0,0,0; 0,0,0]; % The interpolated value of E at each point.
for point = 1 : size(all_points, 1)
    coords_orig = all_points(point, 1:3);
    coords_axis_aligned = (transformation \ (coords_orig' - translation))' + [1,1,1];
    x_val = coords_axis_aligned(1);
    y_val = coords_axis_aligned(2);
    z_val = coords_axis_aligned(3);
    % NOTE: REPLACED BINARY_SEARCH_VOXEL WITH ROUND DOWN + 1.
    x_idx = floor(x_val) + 1;
    y_idx = floor(y_val) + 1;
    z_idx = floor(z_val) + 1;
%     x_idx = binary_search_voxel(x_val, x);
%     y_idx = binary_search_voxel(y_val, y);
%     z_idx = binary_search_voxel(z_val, z);
    % Interpolation starts here.
    if (x_val == x(x_idx))
        if (y_val == y(y_idx))
            if (z_val == z(z_idx))
                % Case 1: The corner case.
                E_int = E{x_idx, y_idx, z_idx};
                else
                % Case 2: On an xy edge.
                E_temp = 1/2 * (E{x_idx, y_idx, z_idx} + ...
                    E{x_idx, y_idx, z_idx+1});
            end
        else
            if (z_val == z(z_idx))
                % Case 3: On an xz edge.
                E_temp = 1/2 * (E{x_idx, y_idx, z_idx} + ...
                    E{x_idx, y_idx+1, z_idx});
            else
                % case 4: On an x face.
                E_temp = 1/4 * (E{x_idx, y_idx, z_idx} + ...
                    E{x_idx, y_idx+1, z_idx} + ...
                    E{x_idx, y_idx, z_idx+1} + ...
                    E{x_idx, y_idx+1, z_idx+1});
            end
        end
    else
        if (y_val == y(y_idx))
            if(z_val == z(z_idx))
                % case 5: On a yz edge.
                E_temp = 1/2 * (E{x_idx, y_idx, z_idx} + ...
                    E{x_idx+1, y_idx, z_idx});
            else
                % Case 6: On a y face.
                E_temp = 1/4 * (E{x_idx, y_idx, z_idx} + ...
                    E{x_idx+1, y_idx, z_idx} + ...
                    E{x_idx, y_idx, z_idx+1} + ...
                    E{x_idx+1, y_idx, z_idx+1});
            end
        else 
            if (z_val == z(z_idx))
                % Case 7: On a z face.
                E_temp =  1/4 * (E{x_idx, y_idx, z_idx} + ...
                    E{x_idx+1, y_idx, z_idx} + ...
                    E{x_idx, y_idx+1, z_idx} + ...
                    E{x_idx+1, y_idx+1, z_idx});
            else 
                % Case 8: Generic case.
                E_temp = 1/8 * (E{x_idx, y_idx, z_idx} + ...
                    E{x_idx+1, y_idx, z_idx} + ...
                    E{x_idx, y_idx+1, z_idx} + ...
                    E{x_idx, y_idx, z_idx+1} + ...
                    E{x_idx+1, y_idx+1, z_idx} + ...
                    E{x_idx+1, y_idx, z_idx+1} + ...
                    E{x_idx, y_idx+1, z_idx+1} + ...
                    E{x_idx+1, y_idx+1, z_idx+1});
            end
        end
    end
    w_interpolated(point, :) = E_temp * radial(point, :)';
end





%% Helper functions.
function vtk_data = vtk_read(filename)
% Reads from a VTK file. 
% Requires the VTK file to have:
%   - one POINTS section detailing the points of interest. These points
%   should be a regular Cartesian grid.
%   - a VECTORS section detailing the velocity at each point.
%   - Everything else will be discarded.
%% Read multi-data VTK file
fid = fopen(filename,'r');
% initialize variables
vtk_data.translation = [];
vtk_data.transformation = [];
vtk_data.num_points = 0;
vtk_data.dimensions = [];
vtk_data.velocity = [];

while ~feof(fid)
    % read header
    str = fgets(fid);
    str = strip(str);
    if (strlength(str) >= 10 && strcmp(str(1:10), "AFFINE_MAP"))
        vtk_data.translation = zeros(3, 1);
        vtk_data.transformation = zeros(3, 3);
        separate = split(str);
        vtk_data.translation(1) = str2double(separate{2});
        vtk_data.translation(2) = str2double(separate{3});
        vtk_data.translation(3) = str2double(separate{4});
        vtk_data.transformation(1, 1) = str2double(separate{5});
        vtk_data.transformation(2, 1) = str2double(separate{6});
        vtk_data.transformation(3, 1) = str2double(separate{7});
        vtk_data.transformation(1, 2) = str2double(separate{8});
        vtk_data.transformation(2, 2) = str2double(separate{9});
        vtk_data.transformation(3, 2) = str2double(separate{10});
        vtk_data.transformation(1, 3) = str2double(separate{11});
        vtk_data.transformation(2, 3) = str2double(separate{12});
        vtk_data.transformation(3, 3) = str2double(separate{13});
    elseif (strlength(str) >= 10 && strcmp(str(1:10), "DIMENSIONS"))
        vtk_data.dimensions = zeros(3, 1);
        separate = split(str);
        vtk_data.dimensions(1) = str2double(separate{2});
        vtk_data.dimensions(2) = str2double(separate{3});
        vtk_data.dimensions(3) = str2double(separate{4});
    elseif (strlength(str) >= 6 && strcmp(str(1:6), "POINTS"))
        separate = split(str);
        vtk_data.num_points = str2double(separate{2});
%         for i = 1 : vtk_data.num_points
%             str = fgets(fid);
%             str = strip(str);
%             separate = split(str);
%             vtk_data.points(i, 1) = str2double(separate{1});
%             vtk_data.points(i, 2) = str2double(separate{2});
%             vtk_data.points(i, 3) = str2double(separate{3});
%         end
    elseif (strlength(str) >= 7 && strcmp(str(1:7), "VECTORS"))
        vtk_data.velocity = cell(vtk_data.dimensions(1), vtk_data.dimensions(2), vtk_data.dimensions(3));
        for k = 1 : vtk_data.dimensions(3)
            for j = 1 : vtk_data.dimensions(2)
                for i = 1 : vtk_data.dimensions(1)
                    str = fgets(fid);
                    str = strip(str);
                    separate = split(str);
                    vtk_data.velocity{i, j, k} = [0, 0, 0];
                    vtk_data.velocity{i, j, k}(1) = str2double(separate{1});
                    vtk_data.velocity{i, j, k}(2) = str2double(separate{2});
                    vtk_data.velocity{i, j, k}(3) = str2double(separate{3});
                end
            end
        end
    else
        continue;
    end
end
fclose(fid);
end

% Binary search function to find maximum voxel coordinate <= input.
function idx = binary_search_voxel(value, vector)
% Returns the index in vector correpsonding to the coordinate of the
% maximum voxel coordinate less or equal to the input value.
% Inputs:
%   - value: the value we are trying to find the index for
%   - vector: the vector we are trying to find value in. This is a
%     horizontal vector.
%   - translation: the vector by which we translate the vector
%   - transformation: the matrix by which we transform the vector
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
    elseif ((value < vector(lo)) || (value > vector(hi)))
        disp("value not found");
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


