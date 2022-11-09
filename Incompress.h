// All rights reserved.
//

#ifndef included_Incompress
#define included_Incompress

// C++ include files that we need
#include <algorithm>
#include <iostream>
#include <math.h>
#include <sstream>

// Basic include file needed for the mesh functionality.
#include "libmesh/boundary_info.h"
#include "libmesh/const_function.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"
#include "libmesh/equation_systems.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/fe.h"
#include "libmesh/libmesh.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/parsed_function.h"
#include "libmesh/perf_log.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/transient_system.h"
#include "libmesh/utility.h"
#include "libmesh/zero_function.h"

// For systems of equations the DenseSubMatrix
// and DenseSubVector provide convenient ways for
// assembling the element matrix and vector on a
// component-by-component basis.
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"

#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/nonlinear_solver.h"

// The definition of a geometric element
#include "libmesh/elem.h"

#include <libmesh/mesh_modification.h>
#include <libmesh/mesh_tools.h>

#include <map>

#include <fstream>
#include <iostream>

#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"

#include "GeomPar.h"
#include "MatVecOper.h"
#include "TensorDer.h"

using namespace std;
using namespace libMesh;

// DarcySolver is a static class that provides data and functions
// required to sovle the darcy flow.
class Incompress {

public:
  Incompress();
  ~Incompress();

  static int time_itr;
  static double beta_s, dt, incomp_res, incomp_cur, pressure_old;
  static DenseVector<double> dincompdu;
  static double dincompdp;

  static int steady_stab;

  static void incomp_disp();
  static void incomp_vel(EquationSystems &es, const Elem *elem, unsigned int qp,
                         const std::vector<std::vector<Real>> &phi);
  static void compute_incomp_res(EquationSystems &es, const Elem *elem,
                                 unsigned int qp,
                                 const std::vector<std::vector<Real>> &phi);
  static void
  compute_incomp_res_der_dis(const std::vector<std::vector<Real>> &phi,
                             unsigned int qp, unsigned int dof_j);
  static void
  compute_incomp_res_der_vel(const std::vector<std::vector<Real>> &phi,
                             unsigned int qp, unsigned int dof_j);
  static void compute_incomp_res_der(const std::vector<std::vector<Real>> &phi,
                                     unsigned int qp, unsigned int dof_j);
  static void define_systems(EquationSystems &es);
  static void solve_incomp(EquationSystems &es);
  static void assemble_incomp(EquationSystems &es,
                              const std::string &libmesh_dbg_var(system_name));
  //---------------------------------------------------------
};

#endif
