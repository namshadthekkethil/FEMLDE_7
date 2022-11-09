
#include "Incompress.h"

int Incompress::time_itr;
double Incompress::beta_s, Incompress::dt, Incompress::incomp_res,
    Incompress::incomp_cur, Incompress::pressure_old;

DenseVector<double> Incompress::dincompdu;
double Incompress::dincompdp;

int Incompress::steady_stab;

using namespace std;

// Bring in everything from the libMesh namespace
using namespace libMesh;

Incompress::Incompress() {}
Incompress::~Incompress() {}

void Incompress::incomp_disp() {
  incomp_res = (GeomPar::detF - 1.0) - ((1.0 / beta_s) * GeomPar::pressure);
}

void Incompress::incomp_vel(EquationSystems &es, const Elem *elem,
                            unsigned int qp,
                            const std::vector<std::vector<Real>> &phi) {
  System &incomp_system = es.get_system<System>("incomp_system");
  unsigned int incomp_var = incomp_system.variable_number("incomp_var");
  const DofMap &dof_map_incomp = incomp_system.get_dof_map();
  std::vector<dof_id_type> dof_indices_incomp;
  dof_map_incomp.dof_indices(elem, dof_indices_incomp, incomp_var);

  System &pressure_old_system = es.get_system<System>("displacement_old");
  unsigned int pressure_old_var =
      pressure_old_system.variable_number("pressure_old");
  const DofMap &dof_map_pressure_old = pressure_old_system.get_dof_map();
  std::vector<dof_id_type> dof_indices_pressure_old;
  dof_map_pressure_old.dof_indices(elem, dof_indices_pressure_old,
                                   pressure_old_var);

  const unsigned int n_var_dofs = dof_indices_incomp.size();

  incomp_cur = 0.0;
  pressure_old = 0.0;
  for (unsigned int j = 0; j < n_var_dofs; j++) {
    incomp_cur +=
        phi[j][qp] * (incomp_system.current_solution(dof_indices_incomp[j]));
    pressure_old +=
        phi[j][qp] *
        (pressure_old_system.current_solution(dof_indices_pressure_old[j]));
  }

  incomp_res = GeomPar::detF *
                   MatVecOper::contractMat(GeomPar::FInvTra, GeomPar::grad_v) -
               ((1.0 / beta_s) * ((GeomPar::pressure - pressure_old) / dt)) +
               incomp_cur;
}

void Incompress::compute_incomp_res(EquationSystems &es, const Elem *elem,
                                    unsigned int qp,
                                    const std::vector<std::vector<Real>> &phi) {
  if (steady_stab == 0) {
    if (time_itr < 2)
      incomp_disp();
    else
      incomp_vel(es, elem, qp, phi);
  } else if (steady_stab == 1) {
    incomp_disp();
  }
}

void Incompress::compute_incomp_res_der_dis(
    const std::vector<std::vector<Real>> &phi, unsigned int qp,
    unsigned int dof_j) {
  dincompdu.resize(MESH_DIMENSION);
  dincompdu(0) = MatVecOper::contractMat(TensorDer::dJdF, TensorDer::dFdux);
  dincompdu(1) = MatVecOper::contractMat(TensorDer::dJdF, TensorDer::dFduy);
#if (MESH_DIMENSION == 3)
  dincompdu(2) = MatVecOper::contractMat(TensorDer::dJdF, TensorDer::dFduz);
#endif
  dincompdp = -(1.0 / beta_s) * phi[dof_j][qp];
}

void Incompress::compute_incomp_res_der_vel(
    const std::vector<std::vector<Real>> &phi, unsigned int qp,
    unsigned int dof_j) {
  dincompdu.resize(MESH_DIMENSION);
  dincompdu(0) =
      MatVecOper::contractMat(TensorDer::dJFinvTdux, GeomPar::grad_v);
  dincompdu(1) =
      MatVecOper::contractMat(TensorDer::dJFinvTduy, GeomPar::grad_v);
#if (MESH_DIMENSION == 3)
  dincompdu(2) =
      MatVecOper::contractMat(TensorDer::dJFinvTduz, GeomPar::grad_v);
#endif

  dincompdu(0) += GeomPar::detF *
                  MatVecOper::contractMat(GeomPar::FInvTra, TensorDer::dFdux) *
                  (1.0 / dt);
  dincompdu(1) += GeomPar::detF *
                  MatVecOper::contractMat(GeomPar::FInvTra, TensorDer::dFduy) *
                  (1.0 / dt);
#if (MESH_DIMENSION == 3)
  dincompdu(2) += GeomPar::detF *
                  MatVecOper::contractMat(GeomPar::FInvTra, TensorDer::dFduz) *
                  (1.0 / dt);
#endif

  dincompdp = -((1.0 / beta_s) / dt) * phi[dof_j][qp];
}

void Incompress::compute_incomp_res_der(
    const std::vector<std::vector<Real>> &phi, unsigned int qp,
    unsigned int dof_j) {
  if (steady_stab == 0) {
    if (time_itr < 2)
      compute_incomp_res_der_dis(phi, qp, dof_j);
    else
      compute_incomp_res_der_vel(phi, qp, dof_j);
  } else if (steady_stab == 1) {
    compute_incomp_res_der_dis(phi, qp, dof_j);
  }
}

void Incompress::define_systems(EquationSystems &es) {
  LinearImplicitSystem &system_incomp =
      es.add_system<LinearImplicitSystem>("incomp_system");
  system_incomp.add_variable("incomp_var", FIRST, LAGRANGE);

  system_incomp.attach_assemble_function(assemble_incomp);
}

void Incompress::solve_incomp(EquationSystems &es) {
  LinearImplicitSystem &system_incomp =
      es.get_system<LinearImplicitSystem>("incomp_system");

  system_incomp.solve();
}

void Incompress::assemble_incomp(
    EquationSystems &es, const std::string &libmesh_dbg_var(system_name)) {
  const MeshBase &mesh = es.get_mesh();

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the Stokes system object.
  LinearImplicitSystem &system =
      es.get_system<LinearImplicitSystem>("incomp_system");

  unsigned int u_var;

  // Numeric ids corresponding to each variable in the system
  u_var = system.variable_number("incomp_var");

  FEType fe_vel_type = system.variable_type(u_var);
  UniquePtr<FEBase> fe_vel(FEBase::build(dim, fe_vel_type));
  QGauss qrule(dim, fe_vel_type.default_quadrature_order());
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

  DenseMatrix<double> Ke;
  DenseVector<Number> Fe;

  std::vector<dof_id_type> dof_indices;

  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
      mesh.active_local_elements_end();

  for (; el != end_el; ++el) {
    // Store a pointer to the element we are currently
    // working on.  This allows for nicer syntax later.
    const Elem *elem = *el;

    dof_map.dof_indices(elem, dof_indices);

    const unsigned int n_dofs = dof_indices.size();

    fe_vel->reinit(elem);

    Ke.resize(n_dofs, n_dofs);
    Fe.resize(n_dofs);

    for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {

      GeomPar::compute_geoPar(es, elem, qp, phi, dphi);
      incomp_vel(es, elem, qp, phi);

      for (unsigned int dof_i = 0; dof_i < n_dofs; dof_i++) {

        Fe(dof_i) += incomp_res * phi[dof_i][qp] * JxW[qp];

        for (unsigned int dof_j = 0; dof_j < n_dofs; dof_j++) {

          Ke(dof_i, dof_j) += phi[dof_i][qp] * phi[dof_j][qp] * JxW[qp];
        }
      }
    }
    system.matrix->add_matrix(Ke, dof_indices);
    system.rhs->add_vector(Fe, dof_indices);
  }
}
