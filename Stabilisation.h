// All rights reserved.
//

#ifndef included_Stabilisation
#define included_Stabilisation

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
#include "TensorDer.h"

typedef void (*COMPUTEMU)(EquationSystems &es, const Elem *elem);

using namespace std;
using namespace libMesh;

// DarcySolver is a static class that provides data and functions
// required to sovle the darcy flow.
class Stabilisation {

public:
  Stabilisation();
  ~Stabilisation();

  static COMPUTEMU compute_mu;

  static double alpha_stab, tau_stab, rho_s, res_stab;
  static double mu;
  static double deltat_mu, dt;

  static DenseVector<double> JFinvTragradPre, prime_vec, u_prime_old_vec,
      v_prime_old_vec, u_prime_vec;
  static DenseVector<double> JFinvTragradNA;

  static DenseVector<double> dJFinvTduxgradNA, dJFinvTduygradNA,
      dJFinvTduzgradNA;
  static DenseVector<double> dJFinvTduxgradP, dJFinvTduygradP, dJFinvTduzgradP;

  static DenseVector<double> NBex, NBey, NBez;

  static DenseVector<double> duprimedux, duprimeduy, duprimeduz;

  static DenseVector<double> dstabdu;

  static DenseVector<double> JFinvTragradNB;

  static DenseVector<double> duprimedp;

  static DenseVector<double> dofi_cur;
  static DenseMatrix<double> dofij_cur;

  static int num_nodes;

  static double dstabdp;

  static int time_itr;

  static int steady_stab;

  static void init_stab(EquationSystems &es);
  static void compute_tau(EquationSystems &es);
  static void compute_prime();
  static void compute_stab();
  static void compute_stab_der(const std::vector<std::vector<Real>> &phi,
                               unsigned int qp, unsigned int dof_j);
  static void define_systems(EquationSystems &es);
  static void assemble_prime(EquationSystems &es,
                             const std::string &libmesh_dbg_var(system_name));
  static void solve_prime(EquationSystems &es);
  static void compute_u_prime(EquationSystems &es, const Elem *elem,
                              unsigned int qp,
                              const std::vector<std::vector<Real>> &phi);
  static void compute_count_dofi(EquationSystems &es);
  static void compute_count_dofij(EquationSystems &es);
  static void compute_dofi_cur(EquationSystems &es, const Elem *elem);
  static void compute_dofij_cur(EquationSystems &es, const Elem *elem);
  //---------------------------------------------------------
};

#endif
