#include "common_defs.h"
#include "misc_io.h"
#include "polyline.h"
#include "vtk_io.h"
#include "triangulation.h"
#include "labelled_grid.h"
#include "div_grad_curl.h"


class Point;
struct compVectorLex;
void compute_and_save_wall_shear_stress(const MyPolyline<3>& centerline, 
    const std::vector<Triangulation>& collars, 
    const std::vector<LabelledGrid<3,3,Vector<3>>>& velocity_timeslices, 
    const std::vector<LabelledGrid<3,3,MatrixRect<3,3>>>& velocity_derivative_timeslices, 
    const Triangulation& surface,
    const std::string& output_directory);
std::map<Vector<3>, Point, compVectorLex>  compute_frame(const MyPolyline<3>& centerline, const std::vector<Triangulation>& collars, const Triangulation& surface);
std::map<Vector<3>, Point, compVectorLex>  compute_and_write_wall_shear_stress(const std::vector<LabelledGrid<3,3,Vector<3>>>& velocity_timeslices, const std::vector<LabelledGrid<3,3,MatrixRect<3,3>>>& velocity_derivative_timeslices, const Triangulation& surface, 
    std::map<Vector<3>, Point, compVectorLex> points_map, std::vector<Vector<3>> all_points, std::vector<Vector<3>> axial_vectors, std::vector<Vector<3>> radial_vectors, std::vector<Vector<3>> circ_vectors, const std::string& output_directory, double mu, int t);
double norm(Vector<3> arr);

void compute_and_save_residuals_along_centerline(
    const MyPolyline<3>& centerline, 
    const std::vector<Triangulation>& collars, 
    const std::vector<LabelledGrid<3,3,MatrixRect<3,3>>>& velocity_derivative_timeslices, 
    const Triangulation& surface,
    std::map<Vector<3>, Point, compVectorLex> points_map,
    const std::string& output_directory, int t, double mu);
void compute_and_write_time_averaged_shear_stress(std::map<Vector<3>, Point, compVectorLex> points_map, const std::string& output_directory);
void compute_and_save_average_wall_shear_stress_along_centerline2(
    const MyPolyline<3>& centerline, 
    const std::vector<Triangulation>& collars, 
    const std::vector<LabelledGrid<3,3,MatrixRect<3,3>>>& velocity_derivative_timeslices, 
    const Triangulation& surface,
    std::map<Vector<3>, Point, compVectorLex> points_map,
    const std::string& output_directory, int t);
void compute_and_save_average_wall_shear_stress_along_centerline(
    const MyPolyline<3>& centerline, 
    const std::vector<Triangulation>& collars, 
    const std::vector<LabelledGrid<3,3,MatrixRect<3,3>>>& velocity_derivative_timeslices, 
    const std::string& output_directory);


/* Define a custom comparator to support ordering in points_map. */
struct compVectorLex {
    bool operator()(const Vector<3>& lhs, const Vector<3>& rhs) const {
        return (lhs[0] < rhs[0]) ||
        ((lhs[0] == rhs[0]) && (lhs[1] < rhs[1])) ||
        ((lhs[0] == rhs[0]) && (lhs[1] == rhs[1]) && (lhs[2] < rhs[2]));
    }
};


class Point
{
    private:
        Vector<3> coords;
        Vector<3> radial;
        Vector<3> tangent;
        Vector<3> circ;
        Vector<4> radial_aggregate;
        Vector<4> tangent_aggregate;
        std::vector<Vector<3>> shear_stress_over_timeslices;
        std::vector<Vector<3>> my_axial_component_over_timeslices;
        std::vector<Vector<3>> my_radial_component_over_timeslices;
        std::vector<Vector<3>> my_circ_component_over_timeslices;
        std::vector<Vector<3>> dans_axial_component_over_timeslices;
        std::vector<Vector<3>> dans_radial_component_over_timeslices;
        std::vector<Vector<3>> dans_circ_component_over_timeslices;


    public:
        Point() {
            this->coords = {0, 0, 0};
            this->radial = {0, 0, 0};
            this->tangent = {0, 0, 0};
            this->circ = {0, 0, 0};
            this->radial_aggregate = {0, 0, 0, 0};
            this->tangent_aggregate = {0, 0, 0, 0};
        }

        Point(Vector<3> coords) {
            this->coords = coords;
            this->radial = {0, 0, 0};
            this->tangent = {0, 0, 0};
            this->circ = {0, 0, 0};
            this->radial_aggregate = {0, 0, 0, 0};
            this->tangent_aggregate = {0, 0, 0, 0};
        }

        void add_radial(Vector<3> radial) {
            this->radial_aggregate[0] += radial[0];
            this->radial_aggregate[1] += radial[1];
            this->radial_aggregate[2] += radial[2];
            this->radial_aggregate[3] += 1;
        }

        void add_tangent(Vector<3> tangent) {
            this->tangent_aggregate[0] += tangent[0];
            this->tangent_aggregate[1] += tangent[1];
            this->tangent_aggregate[2] += tangent[2];
            this->tangent_aggregate[3] += 1;
        }

        void compute_total_radial() {
            // Take average.
            this->radial[0] = this->radial_aggregate[0] / this->radial_aggregate[3];
            this->radial[1] = this->radial_aggregate[1] / this->radial_aggregate[3];
            this->radial[2] = this->radial_aggregate[2] / this->radial_aggregate[3];
            // Normalize.
            double magnitude = norm(this->radial);
            this->radial[0] /= magnitude;
            this->radial[1] /= magnitude;
            this->radial[2] /= magnitude;
        }

        void compute_total_tangent() {
            // Take average.
            this->tangent[0] = this->tangent_aggregate[0] / this->tangent_aggregate[3];
            this->tangent[1] = this->tangent_aggregate[1] / this->tangent_aggregate[3];
            this->tangent[2] = this->tangent_aggregate[2] / this->tangent_aggregate[3];
            // Normalize.
            double magnitude = norm(this->tangent);
            this->tangent[0] /= magnitude;
            this->tangent[1] /= magnitude;
            this->tangent[2] /= magnitude;
        }

        void compute_and_set_circ() {
            // Compute cross product of radial and tangent to get circ.
            this->circ[0] = this->radial[1] * this->tangent[2] - this->tangent[1] * this->radial[2];
            this->circ[1] = this->radial[2] * this->tangent[0] - this->tangent[2] * this->radial[0];
            this->circ[2] = this->radial[0] * this->tangent[1] - this->tangent[0] * this->radial[1];
            // Normalize.
            double magnitude = norm(this->circ);
            this->circ[0] /= magnitude;
            this->circ[1] /= magnitude;
            this->circ[2] /= magnitude;
        }

        void set_tangent(Vector<3> tangent) {
            this->tangent = tangent;
        }

        void set_radial(Vector<3> radial) {
            this->radial = radial;
        }

        void set_circ(Vector<3> circ) {
            this->circ = circ;
        }

        void add_shear_stress(Vector<3> shear_stress) {
            this->shear_stress_over_timeslices.push_back(shear_stress);
        }

        void add_my_components_timeslices(Vector<3> axial, Vector<3> radial, Vector<3> circ) {
            this->my_axial_component_over_timeslices.push_back(axial);
            this->my_radial_component_over_timeslices.push_back(radial);
            this->my_circ_component_over_timeslices.push_back(circ);
        }

        void add_dans_components_timeslices(Vector<3> axial, Vector<3> radial, Vector<3> circ) {
            this->dans_axial_component_over_timeslices.push_back(axial);
            this->dans_radial_component_over_timeslices.push_back(radial);
            this->dans_circ_component_over_timeslices.push_back(circ);
        }

        Vector<3> get_tangent() {
            return this->tangent;
        }

        Vector<3> get_radial() {
            return this->radial;
        }

        Vector<3> get_circ() {
            return this-> circ;
        }

        std::vector<Vector<3>> get_shear_stresses() {
            return this->shear_stress_over_timeslices;
        }

        Vector<3> get_shear_stress(int time) {
            return this->shear_stress_over_timeslices[time];
        }

        Vector<3> get_my_axial_component(int time) {
            return this->my_axial_component_over_timeslices[time];
        }

        Vector<3> get_my_radial_component(int time) {
            return this->my_radial_component_over_timeslices[time];
        }

        Vector<3> get_my_circ_component(int time) {
            return this->my_circ_component_over_timeslices[time];
        }

        Vector<3> get_dans_axial_component(int time) {
            return this->dans_axial_component_over_timeslices[time];
        }

        Vector<3> get_dans_radial_component(int time) {
            return this->dans_radial_component_over_timeslices[time];
        }

        Vector<3> get_dans_circ_component(int time) {
            return this->dans_circ_component_over_timeslices[time];
        }


};


void compute_and_save_wall_shear_stress(
    const MyPolyline<3>& centerline, 
    const std::vector<Triangulation>& collars, 
    const std::vector<LabelledGrid<3,3,Vector<3>>>& velocity_timeslices, 
    const std::vector<LabelledGrid<3,3,MatrixRect<3,3>>>& velocity_derivative_timeslices, 
    const Triangulation& surface,
    const std::string& output_directory)
{
    /* Throughout the entire code, we use a viscosity constant mu = 4 * 10^(-3) Pa-s. */

    {
        // // CREATION OF SYNTHETIC DATA: CYLINDER ------------------------------------------
        // // LabelledGrid<3,3,Vector<3>> velocity_timeslices_debug_pre_use(Grid<3,3>(velocity_derivative_timeslices[6].indexer.n, Affinity<3,3>(MatrixRect<3,3>::Identity(), velocity_derivative_timeslices[6].affinity.translation_part)), Vector<3>(0,0,1));
        // LabelledGrid<3,3,Vector<3>> velocity_timeslices_debug_pre_use(Grid<3,3>(Index<3>(20,20,20), Affinity<3,3>(MatrixRect<3,3>::Identity() / 20)), Vector<3>(0,0,5));
        // for (int rank = 0; rank < velocity_timeslices_debug_pre_use.num_cells(); rank++) {
        //     Index<3> idx = velocity_timeslices_debug_pre_use.indexer.to_index(rank);
        //     Vector<3> pos = velocity_timeslices_debug_pre_use.position(rank);
        //     if (pow(pos[0] - 0.5, 2) + pow(pos[1] - 0.5, 2) < pow(0.25, 2) && pos[2] > 0.1 && pos[1] < 0.9) {
        //         velocity_timeslices_debug_pre_use.set_label(rank, Vector<3>(0,0,-8));
        //     } else {
        //         velocity_timeslices_debug_pre_use.set_label(rank, Vector<3>(0, 0, 1));
        //     }
        // }
        // VTK_IO::write(velocity_timeslices_debug_pre_use, "test_cylinder.vtk", false, "mag"); // Visuallize vector field.
        // const std::vector<LabelledGrid<3,3,Vector<3>>> velocity_timeslices_debug(10, velocity_timeslices_debug_pre_use);
        // // Synethic surface.
        // Triangulation surface_debug;
        // Miscellaneous_IO::read_from_obj(surface_debug, false,"/home/yimo/Documents/DEBUG/test_cylinder2.obj");
        // // Transform surface to fit inside where I want it to be.
        // double x_min = 0;
        // double x_max = 1;
        // double y_min = 0;
        // double y_max = 1;
        // double z_min = 0;
        // double z_max = 1;
        // double z_lower_lim = 0.11;
        // double z_upper_lim = 0.89;
        // double radius = 0.25;
        // for (int i = 0; i < surface_debug.get_vertices().size(); i++) {
        //     Vector<3> point = surface_debug.get_vertices()[i];
        //     if (point[0] < x_min) {
        //         x_min = point[0];
        //     }
        //     if (point[0] > x_max) {
        //         x_max = point[0];
        //     }
        //     if (point[1] < y_min) {
        //         y_min = point[1];
        //     }
        //     if (point[1] > y_max) {
        //         y_max = point[1];
        //     }
        //     if (point[2] < z_min) {
        //         z_min = point[2];
        //     }
        //     if (point[2] > z_max) {
        //         z_max = point[2];
        //     }
        // }
        // std::vector<Vector<3>> new_vertices;
        // for (int i = 0; i < surface_debug.get_vertices().size(); i++) {
        //     Vector<3> new_vertex = surface_debug.get_vertices()[i];
        //     new_vertex[0] = radius * 2 * ((surface_debug.get_vertices()[i][0] - x_min) / (x_max - x_min) - 0.5) + 0.5;
        //     new_vertex[1] = radius * 2 * ((surface_debug.get_vertices()[i][1] - y_min) / (y_max - y_min) - 0.5) + 0.5;
        //     new_vertex[2] = (z_upper_lim - z_lower_lim) * ((surface_debug.get_vertices()[i][2] - z_min) / (z_max - z_min) - 0.5) + 0.5;
        //     new_vertices.push_back(new_vertex);
        // }
        // //
        // // Compute vertex normals.
        // std::vector<Vector<3>> vertex_normals;
        // for (int i = 0; i < new_vertices.size(); i++) {
        //     Vector<3> normal = new_vertices[i] - Vector<3>(0.5, 0.5, new_vertices[i][2]);
        //     normal = normal / norm(normal); // Normalize.
        //     vertex_normals.push_back(normal);
        // }
        // Triangulation new_surface_debug(new_vertices, surface_debug.get_triangles(), vertex_normals, std::vector<uint64_t>(surface_debug.get_triangles().size(), 0));
        // // Visualize new surface.
        // std::string filename = "new_surface_debug.vtk";
        //     std::ofstream file(filename);
        //     assert(file.is_open());
        //     assert(file);
        //     file << "# vtk DataFile Version 4.2" << std::endl;
        //     file << "definition precedes existence" << std::endl;
        //     file << "ASCII" << std::endl;
        //     file << "DATASET POLYDATA" << std::endl;
        //     // Write points.
        //     file << "POINTS " << new_surface_debug.get_vertices().size() << " " << "double" << std::endl;
        //     for (int point_idx = 0; point_idx < new_surface_debug.get_vertices().size(); point_idx++) {
        //         file << new_surface_debug.get_vertices()[point_idx][0] << " " << new_surface_debug.get_vertices()[point_idx][1] << " " << new_surface_debug.get_vertices()[point_idx][2] << std::endl;
        //     }
        //     // Write connectivity.
        //     file << "POLYGONS " << new_surface_debug.get_triangles().size() << " " << 4 * new_surface_debug.get_triangles().size() << std::endl;
        //     for (int triangle_idx = 0; triangle_idx < new_surface_debug.get_triangles().size(); triangle_idx++) {
        //         file << 3 << " " << new_surface_debug.get_triangles()[triangle_idx][0] << " " << new_surface_debug.get_triangles()[triangle_idx][1] << " " << new_surface_debug.get_triangles()[triangle_idx][2] << std::endl;
        //     }
        // // Synethic velocity derivative field.
        // std::vector<LabelledGrid<3,3,MatrixRect<3,3>>> velocity_derivative_timeslices_debug;
        // for (auto velocity_timeslice : velocity_timeslices_debug)
        //     velocity_derivative_timeslices_debug.push_back(DIV_GRAD_CURL::derivative(velocity_timeslice));
        // //
        // // Synthetic centerline.
        // std::vector<Vector<3>> centerline_points;
        // for (int i = 0; i < 100; i++) {
        //     Vector<3> centerline_point = {0.5, 0.5, 0.01*i};
        //     centerline_points.push_back(centerline_point);
        // }
        // MyPolyline<3> centerline_debug(centerline_points);
        // //
        // // Synethic collars.
        // std::vector<Triangulation> collars_debug;
        // for (int i = 0; i < 100; i++) {
        //     std::vector<Vector<3>> vertices;
        //     std::vector<Index<3>> triangles;
        //     std::vector<Vector<3>> normals;
        //     for (int idx = 0; idx < new_surface_debug.get_vertices().size(); idx++) {
        //         Vector<3> coords = new_surface_debug.get_vertices()[idx];
        //         if (abs(coords[2] - 0.01*i) < 0.01) { // Add to collars.
        //             vertices.push_back(coords);
        //             normals.push_back(new_surface_debug.get_vertex_normals()[idx]);
        //         }
        //     }
        //     Triangulation one_collar_for_debug(vertices, triangles, normals, std::vector<u_int64_t>());
        //     collars_debug.push_back(one_collar_for_debug);
        // }
        // //
        // // Get and save the affinity.
        // Affinity<3, 3> affinity = velocity_derivative_timeslices_debug[0].affinity;
        // //
        // /* 
        // * Compute an orthonormal frame.
        // */
        // // std::map<Vector<3>, Point, compVectorLex> points_map = compute_frame(centerline_debug, collars_debug, new_surface_debug); 
        // std::map<Vector<3>, Point, compVectorLex> points_map = compute_frame(centerline_debug, collars_debug, new_surface_debug);
        // //
        // // Post-processing of compute_frame(): save components of points_map to a vectors of Vectors. 
        // // The order of these points should be in the same order as they are in surface.get_vertices().
        // std::vector<Vector<3>> all_points;
        // std::vector<Vector<3>> axial_vectors;
        // std::vector<Vector<3>> radial_vectors;
        // std::vector<Vector<3>> circ_vectors;
        // for (Vector<3> coords : new_surface_debug.get_vertices()) {
        //     Point point = points_map[coords];
        //     all_points.push_back(coords);
        //     axial_vectors.push_back(point.get_tangent());
        //     radial_vectors.push_back(point.get_radial());
        //     circ_vectors.push_back(point.get_circ());
        // }    
        // //
        // /* 
        // * Compute wall shear stress on each individual timeslice.
        // */
        // double mu = 4 * pow(10, -3);
        // points_map = compute_and_write_wall_shear_stress(velocity_timeslices_debug, velocity_derivative_timeslices_debug, new_surface_debug, points_map, all_points, axial_vectors, radial_vectors, circ_vectors, output_directory, mu, 6); 
        // compute_and_save_average_wall_shear_stress_along_centerline2(centerline_debug, collars_debug, velocity_derivative_timeslices_debug, surface_debug, points_map, output_directory, 6);
        // compute_and_save_average_wall_shear_stress_along_centerline(centerline_debug, collars_debug, velocity_derivative_timeslices_debug, output_directory);
        // // ------------------------------------
    }

    {
        // // CREATION OF SYNTHETIC DATA: DESCENDING AORTA --------------------------------------
        // // Synethic surface.
        // // Miscellaneous_IO::read_from_obj(surface, false,"/home/yimo/Documents/DEBUG/descending_aorta.obj");
        // //
        // // Create synthetic velocity field. 
        // // Rewrite format of surface to use winding number. 
        // std::vector<std::array<Vector<3>,3>> triangle_list;
        // // LabelledGrid<3,3,Vector<3>> velocity_timeslices_debug_pre_use(Grid<3,3>(velocity_derivative_timeslices[6].indexer.n, velocity_derivative_timeslices[6].affinity), Vector<3>(0,0,1));
        // LabelledGrid<3,3,Vector<3>> velocity_timeslices_debug_pre_use(Grid<3,3>(Index<3>(20,20,100), Affinity<3,3>(MatrixRect<3,3>::Identity()*2, Vector<3>(0,30,-120))), Vector<3>(0,0,1));
        // for (int triangle_idx = 0; triangle_idx < surface.get_triangles().size(); triangle_idx++) {
        //     Index<3> triangle = surface.get_triangles()[triangle_idx];
        //     Vector<3> point1 = surface.get_vertices()[triangle[0]];
        //     Vector<3> point2 = surface.get_vertices()[triangle[1]];
        //     Vector<3> point3 = surface.get_vertices()[triangle[2]];
        //     std::array<Vector<3>, 3> triangle_with_points = {point1, point2, point3};
        //     triangle_list.push_back(triangle_with_points);
        // }
        // for (int rank = 0; rank < velocity_timeslices_debug_pre_use.num_cells(); rank++) {
        //     // Compute winding number.
        //     if (MiscellaneousMath::winding_number_3d(velocity_timeslices_debug_pre_use.position(rank), triangle_list) > 0.5) {
        //         velocity_timeslices_debug_pre_use.set_label(rank, Vector<3>(0,0,-11));
        //     }
        // }
        // VTK_IO::write(velocity_timeslices_debug_pre_use, "test_descending_aorta.vtk", false, "mag"); // Visuallize vector field.
        // const std::vector<LabelledGrid<3,3,Vector<3>>> velocity_timeslices_debug(10, velocity_timeslices_debug_pre_use);
        // //
        // // Synethic velocity derivative field.
        // std::vector<LabelledGrid<3,3,MatrixRect<3,3>>> velocity_derivative_timeslices_debug;
        // for (auto velocity_timeslice : velocity_timeslices_debug)
        //     velocity_derivative_timeslices_debug.push_back(DIV_GRAD_CURL::derivative(velocity_timeslice));
        // //
        // // Get and save the affinity.
        // Affinity<3, 3> affinity = velocity_derivative_timeslices_debug[0].affinity;
        // //
        // /* 
        // * Compute an orthonormal frame.
        // */
        // // std::map<Vector<3>, Point, compVectorLex> points_map = compute_frame(centerline_debug, collars_debug, new_surface_debug); 
        // std::map<Vector<3>, Point, compVectorLex> points_map = compute_frame(centerline, collars, surface);
        // //
        // // Post-processing of compute_frame(): save components of points_map to a vectors of Vectors. 
        // // The order of these points should be in the same order as they are in surface.get_vertices().
        // std::vector<Vector<3>> all_points;
        // std::vector<Vector<3>> axial_vectors;
        // std::vector<Vector<3>> radial_vectors;
        // std::vector<Vector<3>> circ_vectors;
        // for (Vector<3> coords : surface.get_vertices()) {
        //     Point point = points_map[coords];
        //     all_points.push_back(coords);
        //     axial_vectors.push_back(point.get_tangent());
        //     radial_vectors.push_back(point.get_radial());
        //     circ_vectors.push_back(point.get_circ());
        // }    
        // //
        // /* 
        // * Compute wall shear stress on each individual timeslice.
        // */
        // double mu = 4 * pow(10, -3);
        // points_map = compute_and_write_wall_shear_stress(velocity_timeslices_debug, velocity_derivative_timeslices_debug, surface, points_map, all_points, axial_vectors, radial_vectors, circ_vectors, output_directory, mu, 6); 
        // compute_and_save_average_wall_shear_stress_along_centerline2(centerline, collars, velocity_derivative_timeslices_debug, surface, points_map, output_directory, 6);
        // compute_and_save_average_wall_shear_stress_along_centerline(centerline, collars, velocity_derivative_timeslices_debug, output_directory);
        // //==================================================================================================
    }


    {
        // USE WITH REAL DATA

        /* 
        * Compute an orthonormal frame.
        */ 
        std::map<Vector<3>, Point, compVectorLex> points_map = compute_frame(centerline, collars, surface);
    
        // Post-processing of compute_frame(): save components of points_map to a vectors of Vectors. 
        // The order of these points should be in the same order as they are in surface.get_vertices().
        std::vector<Vector<3>> all_points;
        std::vector<Vector<3>> axial_vectors;
        std::vector<Vector<3>> radial_vectors;
        std::vector<Vector<3>> circ_vectors;
        for (Vector<3> coords : surface.get_vertices()) {
            Point point = points_map[coords];
            all_points.push_back(coords);
            axial_vectors.push_back(point.get_tangent());
            radial_vectors.push_back(point.get_radial());
            circ_vectors.push_back(point.get_circ());
        }    
    
        /* 
        * Compute wall shear stress on each individual timeslice.
        */
        double mu = 4 * pow(10, -3); // Viscosity constant.
        for (int time = 0; time < velocity_timeslices.size(); time++) {
            points_map = compute_and_write_wall_shear_stress(velocity_timeslices, velocity_derivative_timeslices, surface, points_map, all_points, axial_vectors, radial_vectors, circ_vectors, output_directory, mu, time); 
        }

        for (int time = 0; time < velocity_timeslices.size(); time++) {
            compute_and_save_residuals_along_centerline(centerline, collars, velocity_derivative_timeslices, surface, points_map, output_directory, time, mu);
            compute_and_save_average_wall_shear_stress_along_centerline2(centerline, collars, velocity_derivative_timeslices, surface, points_map, output_directory, time);
        }

        compute_and_write_time_averaged_shear_stress(points_map, output_directory);

    }

    std::cout << "DONE" << std::endl;
}


/* Function that computes an orthonormal frame at each point. 
 * Output points_map maps a point's coordinates (Vector<3>) to a Point object. This helps bundle all the information about this point together. */
std::map<Vector<3>, Point, compVectorLex> compute_frame(const MyPolyline<3>& centerline, const std::vector<Triangulation>& collars, const Triangulation& surface) {
    // Initialize and populate points map to map coordinates of points : Point.
    std::map<Vector<3>, Point, compVectorLex> points_map;

    // Read centerline and correspond the centerline points with their respective collars.
    std::vector<Vector<3>> centerline_tangents;
    for (int t = 0; t < centerline.vertices.size(); t++) {
        centerline_tangents.push_back(centerline.tangent(t));
    }
    
    // Populate the points_map and the points' corresponding collars.
    for (int collar_i = 0; collar_i < collars.size(); collar_i++) {
        for (int v = 0; v < collars[collar_i].get_vertices().size(); v++) {
            Vector<3> coords = {0, 0, 0};
            coords[0] = collars[collar_i].get_vertices()[v][0];
            coords[1] = collars[collar_i].get_vertices()[v][1];
            coords[2] = collars[collar_i].get_vertices()[v][2];
            // Add point to points map if not already in there.
            if (points_map.find(coords) == points_map.end()) {
                points_map[coords] = Point(coords);
            }
        }
    }

    // Calculate tangents at each point.
    for (int clp = 0; clp < centerline.vertices.size(); clp++) { // Iterate over clps / collars.
        // Get tangent.
        Vector<3> tangent = centerline.tangent(clp);
        for (int v = 0; v < collars[clp].get_vertices().size(); v++) { // Iterate through points within a collar.
            // Get the point.
            Vector<3> coords = {0, 0, 0};
            coords[0] = collars[clp].get_vertices()[v][0];
            coords[1] = collars[clp].get_vertices()[v][1];
            coords[2] = collars[clp].get_vertices()[v][2];
            // Add tangent.
            points_map[coords].add_tangent(tangent);
        }
    }

    // Calculate normals at each point.
    for (int clp = 0; clp < centerline.vertices.size(); clp++) { // Iterate over clps / collars.
        for (int v = 0; v < collars[clp].get_vertices().size(); v++) { // Iterate through points within a collar.
            // Get the point.
            Vector<3> coords = {0, 0, 0};
            coords[0] = collars[clp].get_vertices()[v][0];
            coords[1] = collars[clp].get_vertices()[v][1];
            coords[2] = collars[clp].get_vertices()[v][2];
            // Add normal.
            points_map[coords].add_radial(collars[clp].get_vertex_normals()[v]);
        }
    }

    // Iterate over points_map to average and normallize all vectors.
    for (auto it : points_map) {
        points_map[it.first].compute_total_radial(); // Average and normalize radials.
        points_map[it.first].compute_total_tangent(); // Average and normalize tangents.
    }

    // Shift radial vector to be normal to axial vector.
    for (auto it : points_map) {
        // TODO: Figure out if we need to force the frame to be perfectly orthonormal. Uncomment the next three lines if forcing orthonormality.
        // // Compute projection of radial onto tangential, and subtract that from tangential to get a tangential vector that is actually normal to the radial ((I - a*a^T)*r).
        // Vector<3> new_radial = (MatrixRect<3, 3>::Identity() - points_map[it.first].get_tangent() * points_map[it.first].get_tangent().transpose()) * points_map[it.first].get_radial();
        // points_map[it.first].set_radial(new_radial / norm(new_radial));
        points_map[it.first].compute_and_set_circ(); // Compute circs from cross prod of radial and tangent.
    }

    // Iterate over surface to include all not included points.
    for (Vector<3> point : surface.get_vertices()) { // REVIEW:
        // Skip if already included in map.
        if (points_map.find(point) != points_map.end()) {
            continue;
        }
        points_map[point] = Point(point);
        points_map[point].set_tangent(Vector<3>(0, 0, 0));
        points_map[point].set_radial(Vector<3>(0, 0, 0));
        points_map[point].set_circ(Vector<3>(0, 0, 0));
    }

    return points_map;
}


/* Function that calculates real actual wall shear stress. */
std::map<Vector<3>, Point, compVectorLex> compute_and_write_wall_shear_stress(const std::vector<LabelledGrid<3,3,Vector<3>>>& velocity_timeslices, const std::vector<LabelledGrid<3,3,MatrixRect<3,3>>>& velocity_derivative_timeslices, const Triangulation& surface, 
    std::map<Vector<3>, Point, compVectorLex> points_map, std::vector<Vector<3>> all_points, std::vector<Vector<3>> axial_vectors, std::vector<Vector<3>> radial_vectors, std::vector<Vector<3>> circ_vectors, 
    const std::string& output_directory, double mu, int t) 
{
    // NOTE: E on a LabelledGrid is already computed (in velocity_derivative_timeslices).

    // Interpolate E with respect to the points in the input surface mesh.
    std::vector<Vector<3>> w_interpolated; // This will be ordered in the same way that all_points, axial_vectors, radial_vectors, and circ_vectors are.

    for (int point_idx = 0; point_idx < all_points.size(); point_idx++) {
        
        MatrixRect<3,3> interpolated_E = velocity_derivative_timeslices[t].interpolate_label(all_points[point_idx]);
        Vector<3> shear_stress = mu * (interpolated_E + interpolated_E.transpose()) / 2 * radial_vectors[point_idx];

        // Add shear stress to corresponding Point object.
        points_map[all_points[point_idx]].add_shear_stress(shear_stress);

        w_interpolated.push_back(shear_stress);
    }

    // Compute components of wall shear stress
    std::vector<double> axial_components;
    std::vector<double> radial_components;
    std::vector<double> circ_components;
    for (int point_idx = 0; point_idx < all_points.size(); point_idx++) {
        axial_components.push_back(axial_vectors[point_idx].dot(w_interpolated[point_idx]));
        radial_components.push_back(radial_vectors[point_idx].dot(w_interpolated[point_idx]));
        circ_components.push_back(circ_vectors[point_idx].dot(w_interpolated[point_idx]));
    }

    /* 
    * Write this timeslice to VTK file.
    */ 
    std::string filename = output_directory + "/" + std::to_string(t) + ".vtk";
    std::ofstream file(filename);
    assert(file.is_open());
    assert(file);
    file << "# vtk DataFile Version 4.2" << std::endl;
    file << "definition precedes existence" << std::endl;
    file << "ASCII" << std::endl;
    file << "DATASET POLYDATA" << std::endl;

    // Write points.
    file << "POINTS " << all_points.size() << " " << "double" << std::endl;
    for (int point_idx = 0; point_idx < all_points.size(); point_idx++) {
        file << all_points[point_idx][0] << " " << all_points[point_idx][1] << " " << all_points[point_idx][2] << std::endl;
    }

    // Write connectivity.
    file << "POLYGONS " << surface.get_triangles().size() << " " << 4 * surface.get_triangles().size() << std::endl;
    for (int triangle_idx = 0; triangle_idx < surface.get_triangles().size(); triangle_idx++) {
        file << 3 << " " << surface.get_triangles()[triangle_idx][0] << " " << surface.get_triangles()[triangle_idx][1] << " " << surface.get_triangles()[triangle_idx][2] << std::endl;
    }
    
    // Write shear stress components.
    file << "POINT_DATA " << all_points.size() << std::endl;
    file << "SCALARS " << "axial_component_of_shear_stress " << "double " << "1" << std::endl;
    file << "LOOKUP_TABLE default" << std::endl;
    for (int point_idx = 0; point_idx < axial_components.size(); point_idx++) {
        file << axial_components[point_idx] << std::endl;
    }
    file << "SCALARS " << "radial_component_of_shear_stress " << "double " << "1" << std::endl;
    file << "LOOKUP_TABLE default" << std::endl;
    for (int point_idx = 0; point_idx < radial_components.size(); point_idx++) {
        file << radial_components[point_idx] << std::endl;
    }
    file << "SCALARS " << "circ_component_of_shear_stress " << "double " << "1" << std::endl;
    file << "LOOKUP_TABLE default" << std::endl;
    for (int point_idx = 0; point_idx < circ_components.size(); point_idx++) {
        file << circ_components[point_idx] << std::endl;
    }

    // Write frame vectors.for (int v = 0; v < collars[clp].get_vertices().size(); v++) { // Iterate through points within a collar.
    file << "VECTORS " << "axial_vector_unit " << "double" << std::endl;
    for (int point_idx = 0; point_idx < axial_vectors.size(); point_idx++) {
        file << axial_vectors[point_idx][0] << " " << axial_vectors[point_idx][1] << " " << axial_vectors[point_idx][2] << std::endl;
    }
    file << "VECTORS " << "radial_vector_unit " << "double" << std::endl;
    for (int point_idx = 0; point_idx < radial_vectors.size(); point_idx++) {
        file << radial_vectors[point_idx][0] << " " << radial_vectors[point_idx][1] << " " << radial_vectors[point_idx][2] << std::endl;
    }
    file << "VECTORS " << "circ_vector_unit " << "double" << std::endl;
    for (int point_idx = 0; point_idx < circ_vectors.size(); point_idx++) {
        file << circ_vectors[point_idx][0] << " " << circ_vectors[point_idx][1] << " " << circ_vectors[point_idx][2] << std::endl;
    }

    // Write shear stress vector.
    file << "VECTORS " << "shear_stress " << "double" << std::endl;
    for (int point_idx = 0; point_idx < w_interpolated.size(); point_idx++) {
        file << w_interpolated[point_idx][0] << " " << w_interpolated[point_idx][1] << " " << w_interpolated[point_idx][2] << std::endl;
    }

    file.close();
    return points_map;
}


/* Helper function that finds the norm of an input 1x3 array.
 * Requires:
 *  - Input vector is of type Vector<3>. Input is 1-D in one direction, 3-D in the other direction.
 */
double norm(Vector<3> arr) 
{
    double sum_sq = 0;
    for (int i = 0; i < 3; i++) {
        sum_sq += pow(arr[i], 2);
    }
    return sqrt(sum_sq);
}


void compute_and_save_average_wall_shear_stress_along_centerline2(
    const MyPolyline<3>& centerline, 
    const std::vector<Triangulation>& collars, 
    const std::vector<LabelledGrid<3,3,MatrixRect<3,3>>>& velocity_derivative_timeslices, 
    const Triangulation& surface,
    std::map<Vector<3>, Point, compVectorLex> points_map,
    const std::string& output_directory, int t)
{
    // Calculate an average shear stress vector at each centerline point.
    std::vector<double> avg_axial_along_centerline;
    std::vector<double> avg_radial_along_centerline;
    std::vector<double> avg_circ_along_centerline;
    const std::vector<double> distances_from_start = MiscellaneousMath::distances_from_start<3>(centerline.vertices);
    // Iterate over collars.
    for (int clp = 0; clp < centerline.vertices.size(); clp++) {
        double total_shear_stress_axial = 0;
        double total_shear_stress_radial = 0;
        double total_shear_stress_circ = 0;
        // Compute averages within a collar, iterating through points inside the collar.
        for (Vector<3> coords : collars[clp].get_vertices()) {
            Vector<3> temp = points_map[coords].get_shear_stress(t);
            Vector<3> a = points_map[coords].get_tangent();
            Vector<3> r = points_map[coords].get_radial();
            Vector<3> c = points_map[coords].get_circ();
            total_shear_stress_axial += temp.dot(a);
            total_shear_stress_radial += temp.dot(r);
            total_shear_stress_circ += temp.dot(c);
        }
        // The y values we are outputting for the plot. 
        avg_axial_along_centerline.push_back(total_shear_stress_axial / collars[clp].get_vertices().size());
        avg_radial_along_centerline.push_back(total_shear_stress_radial / collars[clp].get_vertices().size());
        avg_circ_along_centerline.push_back(total_shear_stress_circ / collars[clp].get_vertices().size());

    }

    // Plot magnitude of shear stress over the centerline (distance from start)
    const std::vector<Vector<2>> zipped_axial = Miscellaneous_IO::zip(distances_from_start, avg_axial_along_centerline);
    VTK_IO::write<2>(zipped_axial, output_directory + std::string("avg_axial") + std::to_string(t) + (".vtk"));

    const std::vector<Vector<2>> zipped_radial = Miscellaneous_IO::zip(distances_from_start, avg_radial_along_centerline);
    VTK_IO::write<2>(zipped_radial, output_directory + std::string("avg_radial") + std::to_string(t) + (".vtk"));

    const std::vector<Vector<2>> zipped_circ = Miscellaneous_IO::zip(distances_from_start, avg_circ_along_centerline);
    VTK_IO::write<2>(zipped_circ, output_directory + std::string("avg_circ") + std::to_string(t) + (".vtk"));
}

void compute_and_write_time_averaged_shear_stress(std::map<Vector<3>, Point, compVectorLex> points_map, const std::string& output_directory)
{
    std::string filename = output_directory + std::string("time_averaged_shear_stresses.txt");
    std::ofstream file(filename);
    assert(file.is_open());
    assert(file);

    for (auto it : points_map){ 
        double total_shear_stress;
        for (int t = 0; t < points_map[it.first].get_shear_stresses().size(); t++) {
            total_shear_stress += norm(points_map[it.first].get_shear_stress(t));
        }
        total_shear_stress /= points_map[it.first].get_shear_stresses().size();
        file << total_shear_stress << std::endl;
    }

    file.close();
}

void compute_and_save_residuals_along_centerline(
    const MyPolyline<3>& centerline, 
    const std::vector<Triangulation>& collars, 
    const std::vector<LabelledGrid<3,3,MatrixRect<3,3>>>& velocity_derivative_timeslices, 
    const Triangulation& surface,
    std::map<Vector<3>, Point, compVectorLex> points_map,
    const std::string& output_directory, int it, double mu) 
{
    const std::vector<double> distances_from_start = MiscellaneousMath::distances_from_start<3>(centerline.vertices);

    /* CALCULATE DAN'S COMPUTATIONS */
    std::vector<Vector<3>> centerline_tangents;
    for (uint64_t is=0; is<centerline.vertices.size(); ++is)
        centerline_tangents.push_back(centerline.tangent(is));
    const LabelledGrid<3,3,MatrixRect<3,3>> velocity_derivatives = velocity_derivative_timeslices[it];

    std::vector<double> wss_axial;
    std::vector<double> wss_radial;
    std::vector<double> wss_tangential; // for lack of a better word

    for (uint64_t is=0; is<centerline.vertices.size(); ++is)
    {
        const Vector<3> t = centerline_tangents[is];                    // axial direction
        const uint64_t num_collar_vertices = collars[is].get_vertices().size();
        
        double wss_axial_collar_sum = 0;
        double wss_radial_collar_sum = 0;
        double wss_tangential_collar_sum = 0;
        for (uint64_t ic=0; ic < num_collar_vertices; ++ic)
        {
            const Vector<3> n = collars[is].get_vertex_normals()[ic];         // radial direction

            // TODO: Figure out if we need to force the frame to be perfectly orthonormal. Uncomment the next three lines if forcing orthonormality
            // // Compute projection of radial onto tangential, and subtract that from tangential to get a tangential vector that is actually normal to the radial ((I - a*a^T)*r).
            // Vector<3> new_n = (MatrixRect<3, 3>::Identity() - t * t.transpose()) * n;
            // new_n /= norm(new_n);
            Vector<3> new_n = n;
            const Vector<3> b = new_n.cross(t).normalized();                // "tangential" direction
            assert(fabs(b.norm() - 1) < 1e-6);

            const MatrixRect<3,3> D = velocity_derivatives.interpolate_label(collars[is].get_vertices()[ic]);
            const MatrixRect<3,3> S = mu * (D + D.transpose())/2;    // symmetric quadratic form

            // compute S(n,t), S(n,n), S(n,b):
            const Vector<3> temp = S*new_n;
            wss_axial_collar_sum += temp.dot(t);
            wss_radial_collar_sum += temp.dot(new_n);
            wss_tangential_collar_sum += temp.dot(b);
        }

        wss_axial.push_back(wss_axial_collar_sum / num_collar_vertices);
        wss_radial.push_back(wss_radial_collar_sum / num_collar_vertices);
        wss_tangential.push_back(wss_tangential_collar_sum / num_collar_vertices);        
    }


    /* CALCULATE MY COMPUTATIONS */
    // Calculate an average shear stress vector at each centerline point.
    std::vector<double> avg_axial_along_centerline;
    std::vector<double> avg_radial_along_centerline;
    std::vector<double> avg_circ_along_centerline;
    // Iterate over collars.
    for (int clp = 0; clp < centerline.vertices.size(); clp++) {
        double total_shear_stress_axial = 0;
        double total_shear_stress_radial = 0;
        double total_shear_stress_circ = 0;
        // Compute averages within a collar, iterating through points inside the collar.
        for (Vector<3> coords : collars[clp].get_vertices()) {
            Vector<3> temp = points_map[coords].get_shear_stress(it);
            Vector<3> a = points_map[coords].get_tangent();
            Vector<3> r = points_map[coords].get_radial();
            Vector<3> c = points_map[coords].get_circ();
            total_shear_stress_axial += temp.dot(a);
            total_shear_stress_radial += temp.dot(r);
            total_shear_stress_circ += temp.dot(c);
        }
        // The y values we are outputting for the plot. 
        avg_axial_along_centerline.push_back(total_shear_stress_axial / collars[clp].get_vertices().size());
        avg_radial_along_centerline.push_back(total_shear_stress_radial / collars[clp].get_vertices().size());
        avg_circ_along_centerline.push_back(total_shear_stress_circ / collars[clp].get_vertices().size());
    }


    /* CALCULATE RESIDUALS */
    assert(wss_axial.size() == avg_axial_along_centerline.size());
    assert(wss_radial.size() == avg_radial_along_centerline.size());
    assert(wss_tangential.size() == avg_circ_along_centerline.size());
    std::vector<double> axial_residuals; //= wss_axial - avg_axial_along_centerline;
    std::vector<double> radial_residuals; // = wss_radial - avg_radial_along_centerline;
    std::vector<double> circ_residuals; //= wss_tangential - avg_circ_along_centerline;
    for (int i = 0; i < wss_axial.size(); i++) {
        axial_residuals.push_back(wss_axial[i] - avg_axial_along_centerline[i]);
        radial_residuals.push_back(wss_radial[i] - avg_radial_along_centerline[i]);
        circ_residuals.push_back(wss_tangential[i] - avg_circ_along_centerline[i]);
    }

    {
    const std::vector<Vector<2>> zipped = Miscellaneous_IO::zip(distances_from_start, axial_residuals);
    VTK_IO::write<2>(zipped, output_directory + std::string("axial_residuals") + std::to_string(it) + (".vtk"));
    }
    {
    const std::vector<Vector<2>> zipped = Miscellaneous_IO::zip(distances_from_start, radial_residuals);
    VTK_IO::write<2>(zipped, output_directory + std::string("radial_residuals") + std::to_string(it) + (".vtk"));
    }
    {
    const std::vector<Vector<2>> zipped = Miscellaneous_IO::zip(distances_from_start, circ_residuals);
    VTK_IO::write<2>(zipped, output_directory + std::string("circ_residuals") + std::to_string(it) + (".vtk"));
    }
    
}

void compute_and_save_average_wall_shear_stress_along_centerline(
    const MyPolyline<3>& centerline, 
    const std::vector<Triangulation>& collars, 
    const std::vector<LabelledGrid<3,3,MatrixRect<3,3>>>& velocity_derivative_timeslices, 
    const std::string& output_directory)
{
    assert(collars.size() == centerline.vertices.size());

    // gather some prelimenary data (constant over timeslices):

    double mu = 4 * pow(10, -3);
    const std::vector<double> distances_from_start = MiscellaneousMath::distances_from_start<3>(centerline.vertices);

    std::vector<Vector<3>> centerline_tangents;
    for (uint64_t is=0; is<centerline.vertices.size(); ++is)
        centerline_tangents.push_back(centerline.tangent(is));

    for (uint64_t it=0; it<velocity_derivative_timeslices.size(); ++it)
    // for (uint64_t it=6; it < 7; it++)
    {
        // compute wall shear stresses for the current timeslice:

        const LabelledGrid<3,3,MatrixRect<3,3>> velocity_derivatives = velocity_derivative_timeslices[it];

        std::vector<double> wss_axial;
        std::vector<double> wss_radial;
        std::vector<double> wss_tangential; // for lack of a better word

        for (uint64_t is=0; is<centerline.vertices.size(); ++is)
        {
            const Vector<3> t = centerline_tangents[is];                    // axial direction
            assert(fabs(t.norm() - 1) < 1e-6);

            const uint64_t num_collar_vertices = collars[is].get_vertices().size();
            
            double wss_axial_collar_sum = 0;
            double wss_radial_collar_sum = 0;
            double wss_tangential_collar_sum = 0;
            for (uint64_t ic=0; ic < num_collar_vertices; ++ic)
            {
                const Vector<3> n = collars[is].get_vertex_normals()[ic];         // radial direction
                assert(fabs(n.norm() - 1) < 1e-6);


                // TODO: Figure out if we need to force the frame to be perfectly orthonormal. Uncomment the next three lines if forcing orthonormality
                // // Compute projection of radial onto tangential, and subtract that from tangential to get a tangential vector that is actually normal to the radial ((I - a*a^T)*r).
                // Vector<3> new_n = (MatrixRect<3, 3>::Identity() - t * t.transpose()) * n;
                // new_n /= norm(new_n);
                Vector<3> new_n = n;



                const Vector<3> b = new_n.cross(t).normalized();                // "tangential" direction
                assert(fabs(b.norm() - 1) < 1e-6);

                const MatrixRect<3,3> D = velocity_derivatives.interpolate_label(collars[is].get_vertices()[ic]);
                const MatrixRect<3,3> S = mu * (D + D.transpose())/2;    // symmetric quadratic form

                // compute S(n,t), S(n,n), S(n,b):
                const Vector<3> temp = S*new_n;
                wss_axial_collar_sum += temp.dot(t);
                wss_radial_collar_sum += temp.dot(new_n);
                wss_tangential_collar_sum += temp.dot(b);
            }

            wss_axial.push_back(wss_axial_collar_sum / num_collar_vertices);
            wss_radial.push_back(wss_radial_collar_sum / num_collar_vertices);
            wss_tangential.push_back(wss_tangential_collar_sum / num_collar_vertices);        
        }

        // save each component of the wall shear stress for the current timeslice:

        // const std::string timeslice_directory(output_directory + std::to_string(it));
        // Miscellaneous_IO::initialize_directory(timeslice_directory);

        {
            const std::vector<Vector<2>> zipped = Miscellaneous_IO::zip(distances_from_start, wss_tangential);
            VTK_IO::write<2>(zipped, output_directory + std::string("wss_tangential") + std::to_string(it) + (".vtk"));
        }
        {
            const std::vector<Vector<2>> zipped = Miscellaneous_IO::zip(distances_from_start, wss_axial);
            VTK_IO::write<2>(zipped, output_directory + std::string("wss_axial") + std::to_string(it) + (".vtk"));
        }
        {
            const std::vector<Vector<2>> zipped = Miscellaneous_IO::zip(distances_from_start, wss_radial);
            VTK_IO::write<2>(zipped, output_directory + std::string("wss_radial") + std::to_string(it) + (".vtk"));
        }
    }
}
