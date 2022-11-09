
#include "Stabilisation.h"

double Stabilisation::alpha_stab, Stabilisation::tau_stab, Stabilisation::rho_s,
    Stabilisation::res_stab;
double Stabilisation::mu;

double Stabilisation::deltat_mu, Stabilisation::dt;

DenseVector<double> Stabilisation::JFinvTragradPre, Stabilisation::prime_vec,
    Stabilisation::u_prime_old_vec, Stabilisation::v_prime_old_vec,
    Stabilisation::u_prime_vec;
DenseVector<double> Stabilisation::JFinvTragradNA;

DenseVector<double> Stabilisation::dJFinvTduxgradNA,
    Stabilisation::dJFinvTduygradNA, Stabilisation::dJFinvTduzgradNA;
DenseVector<double> Stabilisation::dJFinvTduxgradP,
    Stabilisation::dJFinvTduygradP, Stabilisation::dJFinvTduzgradP;

DenseVector<double> Stabilisation::NBex, Stabilisation::NBey,
    Stabilisation::NBez;

DenseVector<double> Stabilisation::duprimedux, Stabilisation::duprimeduy,
    Stabilisation::duprimeduz;

DenseVector<double> Stabilisation::dstabdu;

DenseVector<double> Stabilisation::JFinvTragradNB;

DenseVector<double> Stabilisation::duprimedp;

DenseVector<double> Stabilisation::dofi_cur;
DenseMatrix<double> Stabilisation::dofij_cur;

int Stabilisation::num_nodes;

double Stabilisation::dstabdp;

int Stabilisation::time_itr;

int Stabilisation::steady_stab;

using namespace std;

// Bring in everything from the libMesh namespace
using namespace libMesh;

Stabilisation::Stabilisation() {}
Stabilisation::~Stabilisation() {}

void Stabilisation::init_stab(EquationSystems &es) {

  compute_count_dofi(es);
  compute_count_dofij(es);
}

void Stabilisation::compute_tau(EquationSystems &es) {
  const MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  NonlinearImplicitSystem &system =
      es.get_system<NonlinearImplicitSystem>("NonlinearElasticity");

  const unsigned int u_var = system.variable_number("u");
  const DofMap &dof_map = system.get_dof_map();

  FEType fe_type = dof_map.variable_type(u_var);
  UniquePtr<FEBase> fe(FEBase::build(dim, fe_type));
  QGauss qrule(dim, CONSTANT);
  fe->attach_quadrature_rule(&qrule);

  const std::vector<Real> &JxW = fe->get_JxW();
  const std::vector<std::vector<Real>> &phi = fe->get_phi();
  const std::vector<std::vector<RealGradient>> &dphi = fe->get_dphi();

  double hbyck_min = 1.0e6;

  double h2bymu_min = 1.0e6;
  double h2bymu_glob = 0.0;

  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
      mesh.active_local_elements_end();
  for (; el != end_el; ++el) {
    const Elem *elem = *el;

    fe->reinit(elem);

    GeomPar::compute_geoPar(es, elem, 0, phi, dphi);
    compute_mu(es, elem);

    double h_elem = min(elem->length(0, 1), elem->length(1, 2));
    h_elem = min(h_elem, elem->length(0, 2));

#if (MESH_DIMENSION == 3)
    h_elem = min(h_elem, elem->length(0, 3));
    h_elem = min(h_elem, elem->length(1, 3));
    h_elem = min(h_elem, elem->length(2, 3));
#endif

    if (steady_stab == 0) {
      double hbyck_elem = h_elem / sqrt(mu / rho_s);

      if (hbyck_elem < hbyck_min)
        hbyck_min = hbyck_elem;

      // cout << "mu=" << mu << " hbyck_elem=" << hbyck_elem
      //      << " hbyck_min=" << h2bymu_min << endl;
    }

    else if (steady_stab == 1) {
      double h2bymu_elem = (h_elem * h_elem) / mu;

      if (h2bymu_elem < h2bymu_min)
        h2bymu_min = h2bymu_elem;

      // cout << "mu=" << mu << " h2bymu_elem=" << h2bymu_elem
      //      << " h2bymu_min=" << h2bymu_min << endl;
    }
  }
  // cout << "h2bymu_min=" << h2bymu_min << endl;
  if (steady_stab == 0) {
    MPI_Allreduce(&hbyck_min, &deltat_mu, 1, MPI_DOUBLE, MPI_MIN,
                  MPI_COMM_WORLD);
    tau_stab = 0.5 * alpha_stab * max(deltat_mu / 100.0, min(deltat_mu, dt));
  } else if (steady_stab == 1) {
    // MPI_Allreduce(&h2bymu_min, &h2bymu_glob, 1, MPI_DOUBLE, MPI_MIN,
    //               MPI_COMM_WORLD);
    h2bymu_glob = h2bymu_min;
    tau_stab = 0.5 * alpha_stab * h2bymu_glob;
  }
}

void Stabilisation::compute_prime() {
  JFinvTragradPre.resize(MESH_DIMENSION);
  GeomPar::FInvTra.vector_mult(JFinvTragradPre, GeomPar::grad_pre);
  JFinvTragradPre.scale(GeomPar::detF);

  prime_vec.resize(MESH_DIMENSION);
  if (steady_stab == 0) {
    prime_vec.add(rho_s, GeomPar::acc_vec);
    prime_vec.add(-1.0, JFinvTragradPre);
    prime_vec.scale(-tau_stab / rho_s);
  } else if (steady_stab == 1) {
    prime_vec.add(tau_stab, JFinvTragradPre);
  }
}

void Stabilisation::compute_dofi_cur(EquationSystems &es, const Elem *elem) {
  System &dofi_system = es.get_system<System>("dofiCount");
  unsigned int dofi_var = dofi_system.variable_number("dofiVar");
  const DofMap &dof_map_dofi = dofi_system.get_dof_map();
  std::vector<dof_id_type> dof_indices_dofi;

  dof_map_dofi.dof_indices(elem, dof_indices_dofi, dofi_var);

  num_nodes = dof_indices_dofi.size();

  dofi_cur.resize(num_nodes);
  for (unsigned int j = 0; j < num_nodes; j++) {
    dofi_cur(j) = dofi_system.current_solution(dof_indices_dofi[j]);
  }
}

void Stabilisation::compute_dofij_cur(EquationSystems &es, const Elem *elem) {
  System &dofij_system = es.get_system<System>("dofijSystem");
  unsigned int dofij01_var = dofij_system.variable_number("dofij_01");
  unsigned int dofij02_var = dofij_system.variable_number("dofij_02");
#if (MESH_DIMENSION == 3)
  unsigned int dofij03_var = dofij_system.variable_number("dofij_03");
#endif
  unsigned int dofij12_var = dofij_system.variable_number("dofij_12");
#if (MESH_DIMENSION == 3)
  unsigned int dofij13_var = dofij_system.variable_number("dofij_13");
  unsigned int dofij23_var = dofij_system.variable_number("dofij_23");
#endif
  const DofMap &dof_map_dofij = dofij_system.get_dof_map();

  const int dof_indices_dofij_01 =
      elem->dof_number(dofij_system.number(), dofij01_var, 0);
  const int dof_indices_dofij_02 =
      elem->dof_number(dofij_system.number(), dofij02_var, 0);
#if (MESH_DIMENSION == 3)
  const int dof_indices_dofij_03 =
      elem->dof_number(dofij_system.number(), dofij03_var, 0);
#endif
  const int dof_indices_dofij_12 =
      elem->dof_number(dofij_system.number(), dofij12_var, 0);
#if (MESH_DIMENSION == 3)
  const int dof_indices_dofij_13 =
      elem->dof_number(dofij_system.number(), dofij13_var, 0);
  const int dof_indices_dofij_23 =
      elem->dof_number(dofij_system.number(), dofij23_var, 0);
#endif

  dofij_cur.resize(num_nodes, num_nodes);
  for (unsigned int j = 0; j < num_nodes; j++) {
    dofij_cur(j, j) = dofi_cur(j);
  }

  dofij_cur(0, 1) = dofij_system.current_solution(dof_indices_dofij_01);
  dofij_cur(0, 2) = dofij_system.current_solution(dof_indices_dofij_02);
#if (MESH_DIMENSION == 3)
  dofij_cur(0, 3) = dofij_system.current_solution(dof_indices_dofij_03);
#endif
  dofij_cur(1, 2) = dofij_system.current_solution(dof_indices_dofij_12);
#if (MESH_DIMENSION == 3)
  dofij_cur(1, 3) = dofij_system.current_solution(dof_indices_dofij_13);
  dofij_cur(2, 3) = dofij_system.current_solution(dof_indices_dofij_23);
#endif

  dofij_cur(1, 0) = dofij_cur(0, 1);
  dofij_cur(2, 0) = dofij_cur(0, 2);
  dofij_cur(2, 1) = dofij_cur(1, 2);
#if (MESH_DIMENSION == 3)
  dofij_cur(3, 0) = dofij_cur(0, 3);
  dofij_cur(3, 1) = dofij_cur(1, 3);
  dofij_cur(3, 2) = dofij_cur(2, 3);
#endif
}

void Stabilisation::compute_stab() {

  JFinvTragradNA.resize(MESH_DIMENSION);
  GeomPar::FInvTra.vector_mult(JFinvTragradNA, TensorDer::gradNA);
  JFinvTragradNA.scale(GeomPar::detF);
  if (steady_stab == 0) {
    if (time_itr < 2)
      res_stab = -MatVecOper::contractVec(u_prime_vec, JFinvTragradNA);
    else
      res_stab = -MatVecOper::contractVec(prime_vec, JFinvTragradNA);
  } else if (steady_stab == 1) {
    res_stab = -MatVecOper::contractVec(prime_vec, JFinvTragradNA);
  }
}

void Stabilisation::compute_stab_der(const std::vector<std::vector<Real>> &phi,
                                     unsigned int qp, unsigned int dof_j) {

  dJFinvTduxgradNA.resize(MESH_DIMENSION);
  dJFinvTduygradNA.resize(MESH_DIMENSION);
  dJFinvTduzgradNA.resize(MESH_DIMENSION);

  TensorDer::dJFinvTdux.vector_mult(dJFinvTduxgradNA, TensorDer::gradNA);
  TensorDer::dJFinvTduy.vector_mult(dJFinvTduygradNA, TensorDer::gradNA);
#if (MESHDIMENSION == 3)
  TensorDer::dJFinvTduz.vector_mult(dJFinvTduzgradNA, TensorDer::gradNA);
#endif

  dJFinvTduxgradP.resize(MESH_DIMENSION);
  dJFinvTduygradP.resize(MESH_DIMENSION);
  dJFinvTduzgradP.resize(MESH_DIMENSION);
  TensorDer::dJFinvTdux.vector_mult(dJFinvTduxgradP, GeomPar::grad_pre);
  TensorDer::dJFinvTduy.vector_mult(dJFinvTduygradP, GeomPar::grad_pre);
#if (MESHDIMENSION == 3)
  TensorDer::dJFinvTduz.vector_mult(dJFinvTduzgradP, GeomPar::grad_pre);
#endif

  DenseVector<double> NBex(MESH_DIMENSION);
  DenseVector<double> NBey(MESH_DIMENSION);
#if (MESHDIMENSION == 3)
  DenseVector<double> NBez(MESH_DIMENSION);
#endif
  NBex(0) = phi[dof_j][qp];
  NBex(1) = phi[dof_j][qp];
#if (MESHDIMENSION == 3)
  NBex(2) = phi[dof_j][qp];
#endif

  duprimedux.resize(MESH_DIMENSION);
  duprimeduy.resize(MESH_DIMENSION);
  duprimeduz.resize(MESH_DIMENSION);

  if (steady_stab == 0) {
    duprimedux.add(-tau_stab * (1.0 / (dt * dt)), NBex);
    duprimeduy.add(-tau_stab * (1.0 / (dt * dt)), NBey);
#if (MESHDIMENSION == 3)
    duprimeduz.add(-tau_stab * (1.0 / (dt * dt)), NBez);
#endif

    duprimedux.add(tau_stab / rho_s, dJFinvTduxgradP);
    duprimeduy.add(tau_stab / rho_s, dJFinvTduygradP);
#if (MESHDIMENSION == 3)
    duprimeduz.add(tau_stab / rho_s, dJFinvTduzgradP);
#endif
  } else if (steady_stab == 1) {

    duprimedux.add(tau_stab, dJFinvTduxgradP);
    duprimeduy.add(tau_stab, dJFinvTduygradP);
#if (MESHDIMENSION == 3)
    duprimeduz.add(tau_stab, dJFinvTduzgradP);
#endif
  }

  dstabdu.resize(MESH_DIMENSION);

  dstabdu(0) = -MatVecOper::contractVec(prime_vec, dJFinvTduxgradNA);
  dstabdu(1) = -MatVecOper::contractVec(prime_vec, dJFinvTduygradNA);
#if (MESHDIMENSION == 3)
  dstabdu(2) = -MatVecOper::contractVec(prime_vec, dJFinvTduzgradNA);
#endif

  dstabdu(0) += -MatVecOper::contractVec(duprimedux, JFinvTragradNA);
  dstabdu(1) += -MatVecOper::contractVec(duprimeduy, JFinvTragradNA);
#if (MESHDIMENSION == 3)
  dstabdu(2) += -MatVecOper::contractVec(duprimeduz, JFinvTragradNA);
#endif

  JFinvTragradNB.resize(MESH_DIMENSION);
  GeomPar::FInvTra.vector_mult(JFinvTragradNB, TensorDer::gradNB);
  JFinvTragradNB.scale(GeomPar::detF);

  duprimedp.resize(MESH_DIMENSION);
  if (steady_stab == 0)
    duprimedp.add(tau_stab / rho_s, JFinvTragradNB);
  else if (steady_stab == 1)
    duprimedp.add(tau_stab, JFinvTragradNB);

  dstabdp = -MatVecOper::contractVec(duprimedp, JFinvTragradNA);

  if (steady_stab == 0) {
    if (time_itr < 2) {
      dstabdu.scale(0.5 * dt);
      dstabdp *= 0.5 * dt;
    }
  }
}

void Stabilisation::define_systems(EquationSystems &es) {
  LinearImplicitSystem &system_prime =
      es.add_system<LinearImplicitSystem>("fineDisSystem");
  system_prime.add_variable("uxPrime", FIRST, LAGRANGE);
  system_prime.add_variable("uyPrime", FIRST, LAGRANGE);
#if (MESH_DIMENSION == 3)
  system_prime.add_variable("uzPrime", FIRST, LAGRANGE);
#endif
  system_prime.add_variable("vxPrime", FIRST, LAGRANGE);
  system_prime.add_variable("vyPrime", FIRST, LAGRANGE);
#if (MESH_DIMENSION == 3)
  system_prime.add_variable("vzPrime", FIRST, LAGRANGE);
#endif

  System &system_dofi = es.add_system<System>("dofiCount");
  system_dofi.add_variable("dofiVar", FIRST, LAGRANGE);

  ExplicitSystem &dofij_system = es.add_system<ExplicitSystem>("dofijSystem");
  dofij_system.add_variable("dofij_01", CONSTANT, MONOMIAL);
  dofij_system.add_variable("dofij_02", CONSTANT, MONOMIAL);
#if (MESH_DIMENSION == 3)
  dofij_system.add_variable("dofij_03", CONSTANT, MONOMIAL);
#endif
  dofij_system.add_variable("dofij_12", CONSTANT, MONOMIAL);
#if (MESH_DIMENSION == 3)
  {
    dofij_system.add_variable("dofij_13", CONSTANT, MONOMIAL);
    dofij_system.add_variable("dofij_23", CONSTANT, MONOMIAL);
  }
#endif

  system_prime.attach_assemble_function(assemble_prime);
}

void Stabilisation::assemble_prime(
    EquationSystems &es, const std::string &libmesh_dbg_var(system_name)) {
  const MeshBase &mesh = es.get_mesh();

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the Stokes system object.
  LinearImplicitSystem &system =
      es.get_system<LinearImplicitSystem>("fineDisSystem");

  unsigned int u_var;

  // Numeric ids corresponding to each variable in the system
  u_var = system.variable_number("uxPrime");

  // Get the Finite Element type for "u".  Note this will be
  // the same as the type for "v".
  FEType fe_vel_type = system.variable_type(u_var);

  // Build a Finite Element object of the specified type for
  // the velocity variables.
  UniquePtr<FEBase> fe_vel(FEBase::build(dim, fe_vel_type));

  // A Gauss quadrature rule for numerical integration.
  // Let the FEType object decide what order rule is appropriate.
  QGauss qrule(dim, fe_vel_type.default_quadrature_order());

  // Tell the finite element objects to use our quadrature rule.
  fe_vel->attach_quadrature_rule(&qrule);

  UniquePtr<FEBase> fe_face(FEBase::build(dim, fe_vel_type));
  QGauss qface(dim - 1,
               fe_vel_type.default_quadrature_order()); // Not sure what the
                                                        // best accuracy is here

  fe_face->attach_quadrature_rule(&qface);

  const std::vector<double> &JxW = fe_vel->get_JxW();
  const std::vector<std::vector<double>> &phi = fe_vel->get_phi();
  const std::vector<std::vector<RealGradient>> &dphi = fe_vel->get_dphi();

  const DofMap &dof_map = system.get_dof_map();

  // Define data structures to contain the element matrix
  // and right-hand-side vector contribution.  Following
  // basic finite element terminology we will denote these
  // "Ke" and "Fe".
  DenseMatrix<double> Ke;

#if (MESH_DIMENSION == 2)
  DenseSubMatrix<Number> Ke_var[4][4] = {
      {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke),
       DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)},
      {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke),
       DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)},
      {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke),
       DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)},
      {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke),
       DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)}};
#elif (MESH_DIMENSION == 3)
  DenseSubMatrix<Number> Ke_var[6][6] = {
      {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke),
       DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke),
       DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)},
      {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke),
       DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke),
       DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)},
      {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke),
       DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke),
       DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)},
      {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke),
       DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke),
       DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)},
      {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke),
       DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke),
       DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)},
      {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke),
       DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke),
       DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)},
  };
#endif

#if (MESH_DIMENSION == 2)
  DenseVector<Number> Fe;
  DenseSubVector<Number> Fe_var[4] = {
      DenseSubVector<Number>(Fe), DenseSubVector<Number>(Fe),
      DenseSubVector<Number>(Fe), DenseSubVector<Number>(Fe)};
#elif (MESH_DIMENSION == 3)
  DenseVector<Number> Fe;
  DenseSubVector<Number> Fe_var[6] = {
      DenseSubVector<Number>(Fe), DenseSubVector<Number>(Fe),
      DenseSubVector<Number>(Fe), DenseSubVector<Number>(Fe),
      DenseSubVector<Number>(Fe), DenseSubVector<Number>(Fe)};
#endif

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;
  std::vector<std::vector<dof_id_type>> dof_indices_uvw(2 * dim);

  // Now we will loop over all the elements in the mesh that
  // live on the local processor. We will compute the element
  // matrix and right-hand-side contribution.  Since the mesh
  // will be refined we want to only consider the ACTIVE elements,
  // hence we use a variant of the active_elem_iterator.
  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
      mesh.active_local_elements_end();

  for (; el != end_el; ++el) {
    // Store a pointer to the element we are currently
    // working on.  This allows for nicer syntax later.
    const Elem *elem = *el;

    // Get the degree of freedom indices for the
    // current element.  These define where in the global
    // matrix and right-hand-side this element will
    // contribute to.
    dof_map.dof_indices(elem, dof_indices);

    for (unsigned int var = 0; var < 2 * dim; var++) {
      dof_map.dof_indices(elem, dof_indices_uvw[var], var);
    }

    const unsigned int n_dofs = dof_indices.size();
    const unsigned int n_var_dofs = dof_indices_uvw[0].size();

    fe_vel->reinit(elem);

    Ke.resize(n_dofs, n_dofs);
    for (unsigned int var_i = 0; var_i < 2 * dim; var_i++)
      for (unsigned int var_j = 0; var_j < 2 * dim; var_j++)
        Ke_var[var_i][var_j].reposition(var_i * n_var_dofs, var_j * n_var_dofs,
                                        n_var_dofs, n_var_dofs);
    Fe.resize(n_dofs);
    for (unsigned int var = 0; var < 2 * dim; var++)
      Fe_var[var].reposition(var * n_var_dofs, n_var_dofs);

    for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {

      GeomPar::compute_geoPar(es, elem, qp, phi, dphi);
      Stabilisation::compute_prime();
      Stabilisation::compute_u_prime(es, elem, qp, phi);

      for (unsigned int dof_i = 0; dof_i < n_var_dofs; dof_i++) {
        for (unsigned int var_i = 0; var_i < dim; var_i++) {
          Fe_var[var_i](dof_i) += u_prime_vec(var_i) * phi[dof_i][qp] * JxW[qp];
          Fe_var[var_i + dim](dof_i) +=
              Stabilisation::prime_vec(var_i) * phi[dof_i][qp] * JxW[qp];
        }

        for (unsigned int dof_j = 0; dof_j < n_var_dofs; dof_j++) {
          for (unsigned int var_i = 0; var_i < 2 * dim; var_i++) {
            Ke_var[var_i][var_i](dof_i, dof_j) +=
                phi[dof_i][qp] * phi[dof_j][qp] * JxW[qp];
          }
        }
      }
    } // end of the quadrature point qp-loop

    system.matrix->add_matrix(Ke, dof_indices);
    system.rhs->add_vector(Fe, dof_indices);
  }
}

void Stabilisation::solve_prime(EquationSystems &es) {
  LinearImplicitSystem &system =
      es.get_system<LinearImplicitSystem>("fineDisSystem");

  system.solve();
}

void Stabilisation::compute_u_prime(EquationSystems &es, const Elem *elem,
                                    unsigned int qp,
                                    const std::vector<std::vector<Real>> &phi) {
  System &system_prime = es.get_system<System>("fineDisSystem");
  std::vector<std::vector<dof_id_type>> dof_indices_prime(2 * MESH_DIMENSION);
  const DofMap &dof_map_prime = system_prime.get_dof_map();

  for (unsigned int var = 0; var < 2 * MESH_DIMENSION; var++) {
    dof_map_prime.dof_indices(elem, dof_indices_prime[var], var);
  }

  const unsigned int n_var_dofs = dof_indices_prime[0].size();

  u_prime_old_vec.resize(MESH_DIMENSION);
  v_prime_old_vec.resize(MESH_DIMENSION);
  for (unsigned int var_i = 0; var_i < MESH_DIMENSION; var_i++) {
    for (unsigned int j = 0; j < n_var_dofs; j++) {
      u_prime_old_vec(var_i) += phi[j][qp] * system_prime.current_solution(
                                                 dof_indices_prime[var_i][j]);
      v_prime_old_vec(var_i) +=
          phi[j][qp] * system_prime.current_solution(
                           dof_indices_prime[var_i + MESH_DIMENSION][j]);
    }
  }

  u_prime_vec.resize(MESH_DIMENSION);
  u_prime_vec.add(1.0, u_prime_old_vec);
  u_prime_vec.add(0.5 * dt, v_prime_old_vec);
  u_prime_vec.add(0.5 * dt, prime_vec);
}

void Stabilisation::compute_count_dofi(EquationSystems &es) {
  const MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  System &system = es.get_system<System>("dofiCount");
  unsigned int dofi_var = system.variable_number("dofiVar");

  const DofMap &dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(dofi_var);
  UniquePtr<FEBase> fe(FEBase::build(dim, fe_type));
  QGauss qrule(dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule(&qrule);

  std::vector<dof_id_type> dof_indices_var;

  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
      mesh.active_local_elements_end();

  for (; el != end_el; ++el) {
    const Elem *elem = *el;
    dof_map.dof_indices(elem, dof_indices_var, dofi_var);

    const unsigned int n_var_dofs = dof_indices_var.size();

    fe->reinit(elem);

    for (unsigned int j = 0; j < n_var_dofs; j++) {

      Point pj;
      pj = elem->point(j);
      set<const Elem *> neighborSet_j;

      elem->find_point_neighbors(pj, neighborSet_j);

      system.solution->set(dof_indices_var[j], neighborSet_j.size());
    }
  }
  system.solution->close();
  system.solution->localize(*system.current_local_solution);
}

void Stabilisation::compute_count_dofij(EquationSystems &es) {
  const MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  ExplicitSystem &system = es.get_system<ExplicitSystem>("dofijSystem");
  unsigned int dofij01_var = system.variable_number("dofij_01");
  unsigned int dofij02_var = system.variable_number("dofij_02");
#if (MESH_DIMENSION == 3)
  unsigned int dofij03_var = system.variable_number("dofij_03");
#endif
  unsigned int dofij12_var = system.variable_number("dofij_12");
#if (MESH_DIMENSION == 3)
  unsigned int dofij13_var = system.variable_number("dofij_13");
  unsigned int dofij23_var = system.variable_number("dofij_23");
#endif

  const DofMap &dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(dofij01_var);
  UniquePtr<FEBase> fe(FEBase::build(dim, fe_type));
  QGauss qrule(dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule(&qrule);

  std::vector<dof_id_type> dof_indices_var_01;
  std::vector<dof_id_type> dof_indices_var_02;
#if (MESH_DIMENSION == 3)
  std::vector<dof_id_type> dof_indices_var_03;
#endif
  std::vector<dof_id_type> dof_indices_var_12;
#if (MESH_DIMENSION == 3)
  std::vector<dof_id_type> dof_indices_var_13;
  std::vector<dof_id_type> dof_indices_var_23;
#endif

  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
      mesh.active_local_elements_end();

  for (; el != end_el; ++el) {
    const Elem *elem = *el;
    dof_map.dof_indices(elem, dof_indices_var_01, dofij01_var);
    dof_map.dof_indices(elem, dof_indices_var_02, dofij02_var);
#if (MESH_DIMENSION == 3)
    dof_map.dof_indices(elem, dof_indices_var_03, dofij03_var);
#endif
    dof_map.dof_indices(elem, dof_indices_var_12, dofij12_var);
#if (MESH_DIMENSION == 3)
    dof_map.dof_indices(elem, dof_indices_var_13, dofij13_var);
    dof_map.dof_indices(elem, dof_indices_var_23, dofij23_var);
#endif

    fe->reinit(elem);

    Point p0, p1, p2;
#if (MESH_DIMENSION == 3)
    Point p3;
#endif

    p0 = elem->point(0);
    p1 = elem->point(1);
    p2 = elem->point(2);
#if (MESH_DIMENSION == 3)
    p3 = elem->point(3);
#endif

    set<const Elem *> neighborSet_01, neighborSet_02, neighborSet_12;
#if (MESH_DIMENSION == 3)
    set<const Elem *> neighborSet_03, neighborSet_13, neighborSet_23;
#endif

    elem->find_edge_neighbors(p0, p1, neighborSet_01);
    elem->find_edge_neighbors(p0, p2, neighborSet_02);
    elem->find_edge_neighbors(p1, p2, neighborSet_12);
#if (MESH_DIMENSION == 3)
    elem->find_edge_neighbors(p0, p3, neighborSet_03);
    elem->find_edge_neighbors(p1, p3, neighborSet_13);
    elem->find_edge_neighbors(p2, p3, neighborSet_23);
#endif

    system.solution->set(dof_indices_var_01[0], neighborSet_01.size());
    system.solution->set(dof_indices_var_02[0], neighborSet_02.size());
#if (MESH_DIMENSION == 3)
    system.solution->set(dof_indices_var_03[0], neighborSet_03.size());
#endif
    system.solution->set(dof_indices_var_12[0], neighborSet_12.size());
#if (MESH_DIMENSION == 3)
    system.solution->set(dof_indices_var_13[0], neighborSet_13.size());
    system.solution->set(dof_indices_var_23[0], neighborSet_23.size());
#endif
  }

  // Should call close and update when we set vector entries directly
  system.solution->close();
  system.update();
}

COMPUTEMU Stabilisation::compute_mu = NULL;
