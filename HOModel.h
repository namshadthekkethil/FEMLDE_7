// All rights reserved.
//

#ifndef included_HOModel
#define included_HOModel

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

typedef struct structHO {
  double dCdE[MESH_DIMENSION][MESH_DIMENSION][MESH_DIMENSION][MESH_DIMENSION];
  double d2Jm23dE2[MESH_DIMENSION][MESH_DIMENSION][MESH_DIMENSION]
                  [MESH_DIMENSION];
  double d2I1bardE2[MESH_DIMENSION][MESH_DIMENSION][MESH_DIMENSION]
                   [MESH_DIMENSION];
  double d2I4fbardE2[MESH_DIMENSION][MESH_DIMENSION][MESH_DIMENSION]
                    [MESH_DIMENSION];
  double d2I4sbardE2[MESH_DIMENSION][MESH_DIMENSION][MESH_DIMENSION]
                    [MESH_DIMENSION];
  double d2I8fsbardE2[MESH_DIMENSION][MESH_DIMENSION][MESH_DIMENSION]
                     [MESH_DIMENSION];
  double Diso[MESH_DIMENSION][MESH_DIMENSION][MESH_DIMENSION][MESH_DIMENSION];
  double Df[MESH_DIMENSION][MESH_DIMENSION][MESH_DIMENSION][MESH_DIMENSION];
  double Ds[MESH_DIMENSION][MESH_DIMENSION][MESH_DIMENSION][MESH_DIMENSION];
  double Dfs[MESH_DIMENSION][MESH_DIMENSION][MESH_DIMENSION][MESH_DIMENSION];
  double D[MESH_DIMENSION][MESH_DIMENSION][MESH_DIMENSION][MESH_DIMENSION];
  double Dp[MESH_DIMENSION][MESH_DIMENSION][MESH_DIMENSION][MESH_DIMENSION];
} structHO;

using namespace std;
using namespace libMesh;

// DarcySolver is a static class that provides data and functions
// required to sovle the darcy flow.
class HOModel {

public:
  HOModel();
  ~HOModel();
  static double aHO, bHO, afHO, bfHO, asHO, bsHO, afsHO, bfsHO;

  static double Jm23, I1, I1bar, dWdI1bar, d2WdI1bar2;
  static DenseMatrix<double> dJm23dE, dI1dE, dI1bardE, SIso;

  static double I4f, I4fbar, dWdI4fbar, d2WdI4fbar2;
  static DenseMatrix<double> dI4fdE, dI4fbardE, Sf;

  static double I4s, I4sbar, dWdI4sbar, d2WdI4sbar2;
  static DenseMatrix<double> dI4sdE, dI4sbardE, Ss;

  static double I8fs, I8fsbar, dWdI8fsbar, d2WdI8fsbar2;
  static DenseMatrix<double> dI8fsdE, dI8fsbardE, Sfs;

  static DenseMatrix<double> S;

  static double dWdI1, dWdI4f, dWdI4s, dWdI8fs;

  static double mu;

  static DenseMatrix<double> Sp;

  static structHO Dmat;

  static void compute_HO_PK2(EquationSystems &es, const Elem *elem);
  static void compute_HO_PK2_iso();
  static void compute_HO_D_iso();
  static void compute_HO_D();
  static void compute_HO_PK2_aniso(EquationSystems &es, const Elem *elem);
  static void compute_HO_D_aniso();
  static void compute_PK2_iso();
  static void compute_PK2_aniso();
  static void compute_D_iso();
  static void compute_D_aniso();
  static void compute_D_p();
  static void compute_mu(EquationSystems &es, const Elem *elem);

  //---------------------------------------------------------
};

#endif
