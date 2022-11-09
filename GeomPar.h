// All rights reserved.
//

#ifndef included_GeomPar
#define included_GeomPar

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

#include "MatVecOper.h"

using namespace std;
using namespace libMesh;

// DarcySolver is a static class that provides data and functions
// required to sovle the darcy flow.
class GeomPar {

public:
  GeomPar();
  ~GeomPar();

  static DenseMatrix<double> I, grad_u, grad_v, F, FT, C, FInv, FInvTra, CInv;
  static DenseVector<double> grad_pre, acc_vec;

  static double detF, pressure;

  static void init_geoPar();
  static void
  compute_geoPar(EquationSystems &es, const Elem *elem, unsigned int qp,
                 const std::vector<std::vector<Real>> &phi,
                 const std::vector<std::vector<RealGradient>> &dphi);
  static void compute_FT();
  static void compute_C();
  static void compute_FInv();
  static void compute_FInvTra();
  static void compute_CInv();
  //---------------------------------------------------------
};

#endif
