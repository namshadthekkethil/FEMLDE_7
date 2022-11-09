// All rights reserved.
//

#ifndef included_TensorDer
#define included_TensorDer

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
#include "HOModel.h"
#include "MatVecOper.h"

typedef struct tensorDer {
  double dJFinvTdF[MESH_DIMENSION][MESH_DIMENSION][MESH_DIMENSION]
                  [MESH_DIMENSION];
} tensorDer;

using namespace std;
using namespace libMesh;

// DarcySolver is a static class that provides data and functions
// required to sovle the darcy flow.
class TensorDer {

public:
  TensorDer();
  ~TensorDer();

  static tensorDer tensorder;

  static DenseVector<double> gradNA, gradNB;

  static DenseMatrix<double> gradNAx, gradNAy;
#if (MESH_DIMENSION == 3)
  static DenseMatrix<double> gradNAz;
#endif

  static DenseMatrix<double> dJdF, dFdux, dFduy, dFduz, dJFinvTdux, dJFinvTduy,
      dJFinvTduz;

  static void compute_gradNA(const std::vector<std::vector<RealGradient>> &dphi,
                             unsigned int qp, unsigned int dof_i);
  static void
  compute_gradNAxyz(const std::vector<std::vector<RealGradient>> &dphi,
                    unsigned int qp, unsigned int dof_i);
  static void compute_gradNB(const std::vector<std::vector<RealGradient>> &dphi,
                             unsigned int qp, unsigned int dof_j);
  static void compute_dJdF();
  static void compute_dJFinvTdF();
  static void compute_dJFinvTdu();
  //---------------------------------------------------------
};

#endif
