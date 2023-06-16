#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Constrained_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>
#include <chrono>
#include <fstream>
#include <iostream>
typedef CGAL::Simple_cartesian<double>               Kernel;
typedef Kernel::Point_3                              Point_3;
typedef CGAL::Surface_mesh<Point_3>                  Surface_mesh;
typedef boost::graph_traits<Surface_mesh>::halfedge_descriptor  halfedge_descriptor;
typedef boost::graph_traits<Surface_mesh>::edge_descriptor      edge_descriptor;

namespace SMS = CGAL::Surface_mesh_simplification;

// BGL property map which indicates whether an edge is marked as non-removable
struct Border_is_constrained_edge_map
{
  const Surface_mesh* sm_ptr;
  typedef edge_descriptor                                       key_type;
  typedef bool                                                  value_type;
  typedef value_type                                            reference;
  typedef boost::readable_property_map_tag                      category;
  Border_is_constrained_edge_map(const Surface_mesh& sm) : sm_ptr(&sm) {}
  friend value_type get(const Border_is_constrained_edge_map& m, const key_type& edge) {
    return CGAL::is_border(edge, *m.sm_ptr);
  }
};
typedef SMS::Constrained_placement<SMS::Midpoint_placement<Surface_mesh>,
                                  Border_is_constrained_edge_map > Placement;

int main(int argc, char** argv)
{
  Surface_mesh surface_mesh;
  const std::string filename = (argc > 1) ? argv[1] : "bacb.off";
  std::ifstream is(filename);
  if(!is || !(is >> surface_mesh))
  {
    std::cerr << "Failed to read input mesh: " << filename << std::endl;
    return EXIT_FAILURE;
  }
  if(!CGAL::is_triangle_mesh(surface_mesh))
  {
    std::cerr << "Input geometry is not triangulated." << std::endl;
    return EXIT_FAILURE;
  }

      
  Surface_mesh::Property_map<halfedge_descriptor, std::pair<Point_3, Point_3> > constrained_halfedges;
  constrained_halfedges = surface_mesh.add_property_map<halfedge_descriptor,std::pair<Point_3, Point_3> >("h:vertices").first;
  size_t n_border=0;
  for(halfedge_descriptor hd : halfedges(surface_mesh))
  {
    if(CGAL::is_border(hd, surface_mesh))
    {
      constrained_halfedges[hd] = std::make_pair(surface_mesh.point(source(hd, surface_mesh)),
                                                surface_mesh.point(target(hd, surface_mesh)));
      ++n_border;
    }
  }

  std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();
  // In this example, the simplification stops when the number of undirected edges
  // drops below 10% of the initial count
  double stop_ratio = (argc > 2) ? std::stod(argv[2]) : 0.1;
  SMS::Count_ratio_stop_predicate<Surface_mesh> stop(stop_ratio);
  Border_is_constrained_edge_map bem(surface_mesh);
  int r = SMS::edge_collapse(surface_mesh, stop, CGAL::parameters::edge_is_constrained_map(bem)
                                                  .get_placement(Placement(bem)));
  std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();
  std::cout << "\nFinished!\n" << r << " edges removed.\n" << surface_mesh.number_of_edges() << " final edges.\n";
  std::cout << "Time elapsed: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() << "ms" << std::endl;
  CGAL::IO::write_polygon_mesh((argc > 3) ? argv[3] : "out.off", surface_mesh, CGAL::parameters::stream_precision(17));
  return EXIT_SUCCESS;
}