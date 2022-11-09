// All rights reserved.
//

#ifndef included_MeshGen
#define included_MeshGen

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <sstream>
#include <math.h>

// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/perf_log.h"
#include "libmesh/boundary_info.h"
#include "libmesh/utility.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/zero_function.h"
#include "libmesh/const_function.h"
#include "libmesh/parsed_function.h"

// For systems of equations the DenseSubMatrix
// and DenseSubVector provide convenient ways for
// assembling the element matrix and vector on a
// component-by-component basis.
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"

#include "libmesh/nonlinear_solver.h"
#include "libmesh/nonlinear_implicit_system.h"

// The definition of a geometric element
#include "libmesh/elem.h"

#include <libmesh/mesh_modification.h>
#include <libmesh/mesh_tools.h>

#include "libmesh/face_tri3.h"

#include "libmesh/cell_hex27.h"
#include "libmesh/cell_tet4.h"

#include <map>

#include <iostream>
#include <fstream>

#include "libmesh/string_to_enum.h"
#include "libmesh/getpot.h"

using namespace std;
using namespace libMesh;

// DarcySolver is a static class that provides data and functions
// required to sovle the darcy flow.
class MeshGen
{
//private:
  //  EquationSystems & es;

public:
	MeshGen ();
    ~MeshGen();
	
	static int idij(int nx, int i, int j);
	static int idijk(int nx, int ny, int i, int j, int k);
	static void create_mesh(double xmin, double xmax, double ymin, double ymax, int nx, int ny, Mesh & mesh);
	static void create_mesh_Cooks2D(double xmin, double xmax, double ymin, double ymax, int nx, int ny, Mesh & mesh);
	static void create_mesh_cube2D(double xmin, double xmax, double ymin, double ymax, int nx, int ny, Mesh & mesh);
	static void create_mesh_cube_swing_2D(double xmin, double xmax, double ymin, double ymax, int nx, int ny, Mesh & mesh);
	static void create_mesh_Cooks3D(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, int nx, int ny, int nz, Mesh & mesh);
//---------------------------------------------------------

};


#endif
