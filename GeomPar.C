
#include "GeomPar.h"

DenseMatrix<double> GeomPar::I, GeomPar::grad_u, GeomPar::grad_v, GeomPar::F,
    GeomPar::FT, GeomPar::C, GeomPar::FInv, GeomPar::FInvTra, GeomPar::CInv;
DenseVector<double> GeomPar::grad_pre, GeomPar::acc_vec;
double GeomPar::detF, GeomPar::pressure;

using namespace std;

// Bring in everything from the libMesh namespace
using namespace libMesh;

GeomPar::GeomPar() {}
GeomPar::~GeomPar() {}

void GeomPar::init_geoPar() {
  I.resize(MESH_DIMENSION, MESH_DIMENSION);
  I(0, 0) = 1.0;
  I(1, 1) = 1.0;
#if (MESH_DIMENSION == 3)
  I(2, 2) = 1.0;
#endif
}

void GeomPar::compute_geoPar(
    EquationSystems &es, const Elem *elem, unsigned int qp,
    const std::vector<std::vector<Real>> &phi,
    const std::vector<std::vector<RealGradient>> &dphi) {

  NonlinearImplicitSystem &system =
      es.get_system<NonlinearImplicitSystem>("NonlinearElasticity");
  std::vector<std::vector<dof_id_type>> dof_indices_var(MESH_DIMENSION + 1);
  const DofMap &dof_map = system.get_dof_map();

  System &system_vel = es.get_system<System>("velocity");
  std::vector<std::vector<dof_id_type>> dof_indices_vel(MESH_DIMENSION + 1);
  const DofMap &dof_map_vel = system_vel.get_dof_map();

  for (unsigned int var = 0; var < MESH_DIMENSION + 1; var++) {
    dof_map.dof_indices(elem, dof_indices_var[var], var);
    dof_map_vel.dof_indices(elem, dof_indices_vel[var], var);
  }

  const unsigned int n_var_dofs = dof_indices_var[0].size();

  grad_u.resize(MESH_DIMENSION, MESH_DIMENSION);
  grad_v.resize(MESH_DIMENSION, MESH_DIMENSION);
  for (unsigned int var_i = 0; var_i < MESH_DIMENSION; var_i++) {
    for (unsigned int var_j = 0; var_j < MESH_DIMENSION; var_j++)
      for (unsigned int j = 0; j < n_var_dofs; j++) {
        grad_u(var_i, var_j) +=
            dphi[j][qp](var_j) *
            system.current_solution(dof_indices_var[var_i][j]);
        grad_v(var_i, var_j) +=
            dphi[j][qp](var_j) *
            system_vel.current_solution(dof_indices_vel[var_i][j]);
      }
  }

  pressure = 0.0;
  for (unsigned int j = 0; j < n_var_dofs; j++) {
    pressure += phi[j][qp] *
                system.current_solution(dof_indices_var[MESH_DIMENSION][j]);
  }

  grad_pre.resize(MESH_DIMENSION);
  for (unsigned int var_j = 0; var_j < MESH_DIMENSION; var_j++)
    for (unsigned int j = 0; j < n_var_dofs; j++) {
      grad_pre(var_j) +=
          dphi[j][qp](var_j) *
          system.current_solution(dof_indices_var[MESH_DIMENSION][j]);
    }

  System &system_acc = es.get_system<System>("acceleration");
  const DofMap &dof_map_acc = system_acc.get_dof_map();
  std::vector<std::vector<dof_id_type>> dof_indices_acc(MESH_DIMENSION);

  for (unsigned int var = 0; var < MESH_DIMENSION; var++) {
    dof_map_acc.dof_indices(elem, dof_indices_acc[var], var);
  }

  acc_vec.resize(MESH_DIMENSION);
  for (unsigned int var_i = 0; var_i < MESH_DIMENSION; var_i++) {
    for (unsigned int j = 0; j < n_var_dofs; j++) {
      acc_vec(var_i) +=
          phi[j][qp] * system_acc.current_solution(dof_indices_acc[var_i][j]);
    }
  }

  F.resize(MESH_DIMENSION, MESH_DIMENSION);
  F.add(1.0, grad_u);
  F.add(1.0, I);
  detF = MatVecOper::detMat(F);
  compute_FT();
  compute_C();
  compute_FInv();
  compute_FInvTra();
  compute_CInv();
}

void GeomPar::compute_FT() {
  FT.resize(MESH_DIMENSION, MESH_DIMENSION);
  MatVecOper::transposeMat(F, FT);
}

void GeomPar::compute_C() {
  C.resize(MESH_DIMENSION, MESH_DIMENSION);
  C.add(1.0, FT);
  C.right_multiply(F);
}

void GeomPar::compute_FInv() {
  FInv.resize(MESH_DIMENSION, MESH_DIMENSION);
  MatVecOper::inverseMat(F, FInv);
}

void GeomPar::compute_FInvTra() {
  FInvTra.resize(MESH_DIMENSION, MESH_DIMENSION);
  MatVecOper::transposeMat(FInv, FInvTra);
}

void GeomPar::compute_CInv() {
  CInv.resize(MESH_DIMENSION, MESH_DIMENSION);
  CInv.add(1.0, FInv);
  CInv.right_multiply(FInvTra);
}
