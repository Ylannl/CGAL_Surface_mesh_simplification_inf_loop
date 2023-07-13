# CGAL_Surface_mesh_simplification_inf_loop

See issue: https://github.com/CGAL/cgal/issues/7529

- bacb.off: initial problematic mesh. It was fixed based on [this suggestion](https://github.com/CGAL/cgal/issues/7529#issuecomment-1594731945).
- b9ef.off: mesh that still gives an error even after the above fix.
- CMakeLists.txt generated with cgal_create_CMakeLists utility.
- simplify.cpp source code to reproduce issue. Need to give path bacb.off file as first argument.
