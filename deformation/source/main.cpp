#include "common_defs.h"
#include "misc_io.h"
#include "vtk_io.h"
#include "misc_math.h"

Triangulation deform(Triangulation stl_input, Index<3> n, Vector<3> lower, Vector<3> upper);
Vector<3> deformation_formula(Vector<3> coords);


int main(int argc, char** argv)
{
//     // assert(argc == 7);
//     //
//     // const std::string filename_of_vtk_with_axis_alignment(argv[1]);
//     // const std::string axis_aligned_surface_filename(argv[2]);
//     // const std::string axis_aligned_aorta_skeleton_filename(argv[3]);
//     // const std::string axis_aligned_t_skeleton_filename(argv[4]);
//     // const std::string surface_filename(argv[5]);
//     // const std::string skeleton_filename(argv[6]);
//     //
//     // // extract affinity from filename_of_vtk_with_axis_alignment
//     // Grid<3, 3> grid;
//     // VTK_IO::read(grid, filename_of_vtk_with_axis_alignment);
//     // const Affinity<3,3> affinity = grid.affinity;
//     //
//     // // apply affinity to the surface  
//     // // Miscellaneous_IO::transform_obj(affinity, axis_aligned_surface_filename, surface_filename);    
//     // Triangulation triangulation_axis_aligned;
//     // Miscellaneous_IO::read_from_obj(triangulation_axis_aligned, axis_aligned_surface_filename);
//     // Triangulation triangulation = triangulation_axis_aligned.transform(affinity);
//     // Miscellaneous_IO::write_to_obj(triangulation, surface_filename);
//     //
//     // // read skeletons into polyline trees and combine the trees
//     // MyPolylineTree<3> aorta;
//     // VTK_IO::read(aorta, axis_aligned_aorta_skeleton_filename);
//     //
//     // // MyPolylineTree<3> t;
//     // // VTK_IO::read(t, axis_aligned_t_skeleton_filename);
//     //
//     // std::vector<MyPolyline<3>> polylines_combined;
//     // polylines_combined.insert(polylines_combined.end(), aorta.polylines.begin(), aorta.polylines.end());
//     // // polylines_combined.insert(polylines_combined.end(), t.polylines.begin(), t.polylines.end());
//     // MyPolylineTree<3> combined_axis_aligned(polylines_combined);
//     //
//     // // apply affinity to the polyline tree  
//     // MyPolylineTree<3> combined = combined_axis_aligned.transform(affinity);
//     // VTK_IO::write(combined, skeleton_filename);

    // std::cout << "hello world" << std::endl;
    
    Triangulation mesh_input;
    std::string filename = "/home/yimo/Downloads/bunny.obj";
    Index<3> n = {10, 10, 10};
    Vector<3> lower = {-50, -50, -50};
    Vector<3> upper = {100, 100, 100};
    Miscellaneous_IO::read_from_obj(mesh_input, filename);
    Triangulation stl_output = deform(mesh_input, n, lower, upper);

    return 0;
}


Triangulation deform(Triangulation stl_input, Index<3> n, Vector<3> lower, Vector<3> upper)
{
    // Make grid.
    uint64_t num_x = n[0];
    uint64_t num_y = n[1];
    uint64_t num_z = n[2];
    float x_min = lower[0];
    float y_min = lower[1];
    float z_min = lower[2];
    float x_max = upper[0];
    float y_max = upper[1];
    float z_max = upper[2];
    LabelledGrid<3, 3, Vector<3>> deformation_grid(Grid<3, 3>(n, lower, upper), Vector<3>(0, 0, 0));
    
    // Fill in deformation grid with respective deformations, given by the deformation formula.
    for (int x = 0; x < num_x; x++) {
        for (int y = 0; y < num_y; y++) {
            for (int z = 0; z < num_z; z++) {
                // Get the rank.
                uint64_t rank = deformation_grid.indexer.to_rank(Index<3>(x, y, z));
                // Get the actual value of coordinates at this place.
                Vector<3> coords = {0, 0, 0};
                coords[0] = lower[0] + x * (upper[0] - lower[0]) / (n[0] - 1);
                coords[1] = lower[1] + y * (upper[1] - lower[1]) / (n[1] - 1);
                coords[2] = lower[2] + z * (upper[2] - lower[2]) / (n[2] - 1);
                // Set label to be the deformation vector.
                Vector<3> translation = deformation_formula(coords);
                deformation_grid.set_label(rank, translation);
            }
        }
    }

    // Iterate through stl file and do the actual deformation.
    std::vector<Vector<3>> new_vertices;
    for (int i = 0; i < stl_input.vertices.size(); i++) {
        std::cout << "===Starting " << i << ": " << std::endl; // REMOVE:
        std::cout << "Original point: " << stl_input.vertices[i].transpose() <<std::endl; // REMOVE:

        float x_orig = stl_input.vertices[i][0];
        float y_orig = stl_input.vertices[i][1];
        float z_orig = stl_input.vertices[i][2];

        // Check if this point is outside the grid.
        if (x_orig < x_min || x_orig > x_max) {
            std::cout << "Point " << i << ": " << stl_input.vertices[i].transpose() << " out of x bound for grid." << std::endl;
            exit(1);
        } else if (y_orig < y_min || y_orig > y_max) {
            std::cout << "Point " << i << ": " << stl_input.vertices[i].transpose() << " out of y bound for grid." << std::endl;
            exit(1);
        } else if (z_orig < z_min || z_orig > z_max) {
            std::cout << "Point " << i << ": " << stl_input.vertices[i].transpose() << " out of z bound for grid." << std::endl;
            exit(1);
        }

        Vector<3> deformation = deformation_grid.interpolate_label(stl_input.vertices[i]);
        std::cout << "Deformation:    " << deformation.transpose() << std::endl; // REMOVE:

        // Deform the point as needed.
        Vector<3> new_vertex = {0, 0, 0};
        new_vertex[0] = x_orig + deformation[0];
        new_vertex[1] = y_orig + deformation[1];
        new_vertex[2] = z_orig + deformation[2];
        std::cout << "New point:      " << new_vertex.transpose() << std::endl; // REMOVE:

        // Check if new point is within bounds.
        if (new_vertex[0] < x_min || new_vertex[0] > x_max) {
            std::cout << "Point " << i << ": " << stl_input.vertices[i].transpose() << " out of x bound after deformation: " << new_vertex.transpose() << std::endl;
            exit(1);
        } else if (new_vertex[1] < y_min || new_vertex[1] > y_max) {
            std::cout << "Point " << i << ": " << stl_input.vertices[i].transpose() << " out of y bound after deformation: " << new_vertex.transpose() <<std::endl;
            exit(1);
        } else if (new_vertex[2] < z_min || new_vertex[2] > z_max) {
            std::cout << "Point " << i << ": " << stl_input.vertices[i].transpose() << " out of z bound after deformation: " << new_vertex.transpose() <<std::endl;
            exit(1);
        }

        // Add new point to new vertex list.
        new_vertices.push_back(new_vertex);
    }

//     return Triangulation(new_vertices, stl_input.triangles);
    Triangulation output(new_vertices, stl_input.triangles);
    Miscellaneous_IO::write_to_obj(output, "bunny_rotation_transformation.obj");
    return output;

}

/* Deformation formula. 
 * Deforms a given (x,y,z) into an output (dx,dy,dz) representative of the translation needed to get the deformed point.
 * Requires:
 *  - Input ijk is an Index<3> corresponding to the i, j, k voxel numbers in the 3d grid.
 *  - Output new_ijk is a Vector<3> correpsonding to the translation that the numbers go through. 
 */
Vector<3> deformation_formula(Vector<3> coords) 
{
    std::cout << coords.transpose() << std::endl; //REMOVE:
    Vector<3> translation = {0, 0, 0};
    translation[0] = -coords[1] - coords[0];
    translation[1] = coords[0] - coords[1];
    // translation[0] = ((coords[0] - 0.5) * (4 * pow(coords[2] - 1, 2)) * 0.4 + 0.5) - coords[0];
    // translation[1] = ((coords[1] - 0.5) * (4 * pow(coords[2] - 1, 2)) * 0.4 + 0.5) - coords[1];
    std::cout << translation.transpose() << std::endl; //REMOVE:
    std::cout << "" << std::endl; // REMOVE:

    return translation;
}

/* Helper function to determine which voxel a point lies in.
 * The voxel that a point lies in will determine the translation the point undergoes in its deformation.
 * Requires:
 *  - The dimension that we are looking at should be axis-aligned. If it is not axis aligned, transform the point first via an affine mapping.
 *    The axis that we are mapping onto should also be axis-aligned.
 *  - Input value is the number that we are trying to find in the voxels.
 *  - Input n is the number of voxels along the dimension we are looking at.
 *  - Input lower is the lower bound of the dimension we are looking at.
 *  - Input upper is the upper bound of the dimension we are looking at.
 *  - Output is the voxel number at which value resides.
 */
uint64_t determine_voxel(float value, uint64_t n, float lower, float upper)
{
    float interval_size = (upper - lower + 1) / n;
    uint64_t voxel_less_than_or_equal_to = (value - lower) / interval_size;
    return voxel_less_than_or_equal_to;
}



