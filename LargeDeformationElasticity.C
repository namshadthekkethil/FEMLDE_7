#include "LargeDeformationElasticity.h"

using namespace libMesh;
using namespace std;

LargeDeformationElasticity::~LargeDeformationElasticity() {}

void LargeDeformationElasticity::define_systems() {
  NonlinearImplicitSystem &system =
      es.add_system<NonlinearImplicitSystem>("NonlinearElasticity");

  unsigned int u_var, v_var, w_var;

  cout << "before defining systems" << endl;

  u_var = system.add_variable("u", FIRST, LAGRANGE);

  v_var = system.add_variable("v", FIRST, LAGRANGE);

#if (MESH_DIMENSION == 3)
  { w_var = system.add_variable("w", FIRST, LAGRANGE); }
#endif

  system.add_variable("pressure", FIRST, LAGRANGE);

  System &system_dis_old = es.add_system<System>("displacement_old");
  system_dis_old.add_variable("u_old", FIRST, LAGRANGE);
  system_dis_old.add_variable("v_old", FIRST, LAGRANGE);
#if (MESH_DIMENSION == 3)
  system_dis_old.add_variable("w_old", FIRST, LAGRANGE);
#endif
  system_dis_old.add_variable("pressure_old", FIRST, LAGRANGE);

  System &system_dis_n1 = es.add_system<System>("displacement_n1");
  system_dis_n1.add_variable("u_n1", FIRST, LAGRANGE);
  system_dis_n1.add_variable("v_n1", FIRST, LAGRANGE);
#if (MESH_DIMENSION == 3)
  system_dis_n1.add_variable("w_n1", FIRST, LAGRANGE);
#endif
  system_dis_n1.add_variable("pressure_n1", FIRST, LAGRANGE);

  System &system_vel = es.add_system<System>("velocity");
  system_vel.add_variable("vel_x", FIRST, LAGRANGE);
  system_vel.add_variable("vel_y", FIRST, LAGRANGE);
#if (MESH_DIMENSION == 3)
  system_vel.add_variable("vel_z", FIRST, LAGRANGE);
#endif
  system_vel.add_variable("pressure_vel", FIRST, LAGRANGE);
#if (POROUS == 1)
  system_vel.add_variable("m_vel", FIRST, LAGRANGE);
#endif

  System &system_acc = es.add_system<System>("acceleration");
  system_acc.add_variable("acc_x", FIRST, LAGRANGE);
  system_acc.add_variable("acc_y", FIRST, LAGRANGE);
#if (MESH_DIMENSION == 3)
  system_acc.add_variable("acc_z", FIRST, LAGRANGE);
#endif
  system_acc.add_variable("pressure_acc", FIRST, LAGRANGE);
#if (POROUS == 1)
  system_acc.add_variable("m_acc", FIRST, LAGRANGE);
#endif

  System &system_vel_old = es.add_system<System>("velocity_old");
  system_vel_old.add_variable("vel_x_old", FIRST, LAGRANGE);
  system_vel_old.add_variable("vel_y_old", FIRST, LAGRANGE);
#if (MESH_DIMENSION == 3)
  system_vel_old.add_variable("vel_z_old", FIRST, LAGRANGE);
#endif
  system_vel_old.add_variable("pressure_vel_old", FIRST, LAGRANGE);
#if (POROUS == 1)
  system_vel_old.add_variable("m_vel_old", FIRST, LAGRANGE);
#endif

  System &system_vel_n1 = es.add_system<System>("velocity_n1");
  system_vel_n1.add_variable("vel_x_n1", FIRST, LAGRANGE);
  system_vel_n1.add_variable("vel_y_n1", FIRST, LAGRANGE);
#if (MESH_DIMENSION == 3)
  system_vel_n1.add_variable("vel_z_n1", FIRST, LAGRANGE);
#endif
  system_vel_n1.add_variable("pressure_vel_n1", FIRST, LAGRANGE);
#if (POROUS == 1)
  system_vel_n1.add_variable("m_vel_n1", FIRST, LAGRANGE);
#endif

  System &system_bcid = es.add_system<System>("bcidSystem");
  system_bcid.add_variable("bcidVar", FIRST, LAGRANGE);

  // Also, initialize an ExplicitSystem to store stresses
  ExplicitSystem &stress_system = es.add_system<ExplicitSystem>("StressSystem");
  stress_system.add_variable("sigma_00", CONSTANT, MONOMIAL);
  stress_system.add_variable("sigma_01", CONSTANT, MONOMIAL);
#if (MESH_DIMENSION == 3)
  stress_system.add_variable("sigma_02", CONSTANT, MONOMIAL);
#endif
  stress_system.add_variable("sigma_11", CONSTANT, MONOMIAL);
#if (MESH_DIMENSION == 3)
  {
    stress_system.add_variable("sigma_12", CONSTANT, MONOMIAL);
    stress_system.add_variable("sigma_22", CONSTANT, MONOMIAL);
  }
#endif

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

  ExplicitSystem &pmono_system = es.add_system<ExplicitSystem>("pMonoSystem");
  pmono_system.add_variable("pressureMono", CONSTANT, MONOMIAL);

  ExplicitSystem &J_system = es.add_system<ExplicitSystem>("JSystem");
  J_system.add_variable("JVar", CONSTANT, MONOMIAL);

  System &system_nodeid = es_cur.add_system<System>("NodeIdSystem");
  system_nodeid.add_variable("NodeIdVar", FIRST, LAGRANGE);
}

void LargeDeformationElasticity::pre_solve() { update_nodeid(); }

void LargeDeformationElasticity::solve_lde() {
  NonlinearImplicitSystem &system =
      es.get_system<NonlinearImplicitSystem>("NonlinearElasticity");
  System &system_dis_old = es.get_system<System>("displacement_old");
  System &system_dis_n1 = es.get_system<System>("displacement_n1");

  System &system_vel = es.get_system<System>("velocity");
  System &system_acc = es.get_system<System>("acceleration");
  System &system_vel_old = es.get_system<System>("velocity_old");
  System &system_vel_n1 = es.get_system<System>("velocity_n1");

  old_new(system_vel_n1, system_vel_old);
  old_new(system_vel_old, system_vel);

  old_new(system_dis_n1, system_dis_old);
  old_new(system_dis_old, system);

  libMesh::out << "time=" << ttime << ", Performing solve " << std::endl;

  system.solve();

  libMesh::out << "System solved at nonlinear iteration "
               << system.n_nonlinear_iterations()
               << " , final nonlinear residual norm: "
               << system.final_nonlinear_residual() << std::endl
               << std::endl;
}

void LargeDeformationElasticity::old_new(System &system_old,
                                         System &system_new) {
  system_old.solution->zero();
  system_old.solution->add(1.0, *system_new.current_local_solution);
  system_old.solution->close();
  system_old.solution->localize(*system_old.current_local_solution);
}

void LargeDeformationElasticity::update_nodeid() {
  MeshBase &mesh = es_cur.get_mesh();
  System &system_nodeid = es_cur.get_system<System>("NodeIdSystem");
  NumericVector<double> &nodeid_solution = *(system_nodeid.solution);

  libMesh::MeshBase::const_node_iterator node_it = mesh.local_nodes_begin();
  libMesh::MeshBase::const_node_iterator node_end = mesh.local_nodes_end();

  for (; node_it != node_end; ++node_it) {
    const libMesh::Node *nd = *node_it;
    const unsigned int dof_num_nodeid =
        nd->dof_number(system_nodeid.number(), 0, 0);
    nodeid_solution.set(dof_num_nodeid, nd->id());
  }
  nodeid_solution.close();
  nodeid_solution.localize(*system_nodeid.current_local_solution);
}

void LargeDeformationElasticity::move_mesh() {
  MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  MeshBase &mesh_cur = es_cur.get_mesh();

  NonlinearImplicitSystem &system =
      es.get_system<NonlinearImplicitSystem>("NonlinearElasticity");
  System &system_nodeid = es_cur.get_system<System>("NodeIdSystem");

  NumericVector<Number> &dis_solution = *(system.solution);
  dis_solution.close();
  // dis_solution.localize(*system.current_local_solution);

  NumericVector<double> &nodeid_solution = *(system_nodeid.solution);
  nodeid_solution.close();
  nodeid_solution.localize(*system_nodeid.current_local_solution);

  vector<double> dis_sol_vec;
  dis_solution.localize(dis_sol_vec);

  vector<double> nodeid_sol_vec;
  nodeid_solution.localize(nodeid_sol_vec);

  for (unsigned int n = 0; n < mesh.n_nodes(); n++) {
    int nodeid_n = nodeid_sol_vec[n];
    Point p_node = mesh.node_ref(nodeid_n);
    if (dim == 3) {
#if (POROUS == 0)
      Point p_dis(dis_sol_vec[(dim + 1) * n], dis_sol_vec[(dim + 1) * n + 1],
                  dis_sol_vec[(dim + 1) * n + 2]);
#elif (POROUS == 1)
      Point p_dis(dis_sol_vec[(dim + 2) * n], dis_sol_vec[(dim + 2) * n + 1],
                  dis_sol_vec[(dim + 2) * n + 2]);
#endif
      mesh_cur.node_ref(nodeid_n) = p_node + p_dis;
    }
    if (dim == 2) {
#if (POROUS == 0)
      Point p_dis(dis_sol_vec[(dim + 1) * n], dis_sol_vec[(dim + 1) * n + 1],
                  0.0);
#elif (POROUS == 1)
      Point p_dis(dis_sol_vec[(dim + 2) * n], dis_sol_vec[(dim + 2) * n + 1],
                  0.0);
#endif
      mesh_cur.node_ref(nodeid_n) = p_node + p_dis;
    }
  }
}

void LargeDeformationElasticity::compute_stresses() {
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

  // Also, get a reference to the ExplicitSystem
  ExplicitSystem &stress_system = es.get_system<ExplicitSystem>("StressSystem");
  const DofMap &stress_dof_map = stress_system.get_dof_map();
  unsigned int sigma_vars[3 * dim - 3];
  if (dim == 3) {
    sigma_vars[0] = stress_system.variable_number("sigma_00");
    sigma_vars[1] = stress_system.variable_number("sigma_01");
    sigma_vars[2] = stress_system.variable_number("sigma_02");
    sigma_vars[3] = stress_system.variable_number("sigma_11");
    sigma_vars[4] = stress_system.variable_number("sigma_12");
    sigma_vars[5] = stress_system.variable_number("sigma_22");
  }

  if (dim == 2) {
    sigma_vars[0] = stress_system.variable_number("sigma_00");
    sigma_vars[1] = stress_system.variable_number("sigma_01");
    sigma_vars[2] = stress_system.variable_number("sigma_11");
  }

  // Storage for the stress dof indices on each element
  std::vector<dof_id_type> stress_dof_indices_var;

  // To store the stress tensor on each element
  DenseMatrix<Number> elem_avg_stress_tensor(dim, dim);

  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
      mesh.active_local_elements_end();

  for (; el != end_el; ++el) {
    const Elem *elem = *el;

    fe->reinit(elem);

    GeomPar::compute_geoPar(es, elem, 0, phi, dphi);
    HOModel::compute_HO_PK2(es, elem);

    // clear the stress tensor
    elem_avg_stress_tensor.resize(dim, dim);

    elem_avg_stress_tensor.add(1.0, HOModel::S);

    elem_avg_stress_tensor.scale(1. / GeomPar::detF);
    elem_avg_stress_tensor.left_multiply(GeomPar::F);
    elem_avg_stress_tensor.right_multiply_transpose(GeomPar::F);

#if (MESH_DIMENSION == 2)
    stress_dof_map.dof_indices(elem, stress_dof_indices_var, 0);
    stress_system.solution->set(stress_dof_indices_var[0],
                                elem_avg_stress_tensor(0, 0));
    stress_dof_map.dof_indices(elem, stress_dof_indices_var, 1);
    stress_system.solution->set(stress_dof_indices_var[0],
                                elem_avg_stress_tensor(0, 1));
    stress_dof_map.dof_indices(elem, stress_dof_indices_var, 2);
    stress_system.solution->set(stress_dof_indices_var[0],
                                elem_avg_stress_tensor(1, 1));
#endif

#if (MESH_DIMENSION == 3)
    stress_dof_map.dof_indices(elem, stress_dof_indices_var, 0);
    stress_system.solution->set(stress_dof_indices_var[0],
                                elem_avg_stress_tensor(0, 0));
    stress_dof_map.dof_indices(elem, stress_dof_indices_var, 1);
    stress_system.solution->set(stress_dof_indices_var[0],
                                elem_avg_stress_tensor(0, 1));
    stress_dof_map.dof_indices(elem, stress_dof_indices_var, 2);
    stress_system.solution->set(stress_dof_indices_var[0],
                                elem_avg_stress_tensor(0, 2));
    stress_dof_map.dof_indices(elem, stress_dof_indices_var, 3);
    stress_system.solution->set(stress_dof_indices_var[0],
                                elem_avg_stress_tensor(1, 1));
    stress_dof_map.dof_indices(elem, stress_dof_indices_var, 4);
    stress_system.solution->set(stress_dof_indices_var[0],
                                elem_avg_stress_tensor(1, 2));
    stress_dof_map.dof_indices(elem, stress_dof_indices_var, 5);
    stress_system.solution->set(stress_dof_indices_var[0],
                                elem_avg_stress_tensor(2, 2));
#endif
  }

  // Should call close and update when we set vector entries directly
  stress_system.solution->close();
  stress_system.update();
}

void LargeDeformationElasticity::compute_J() {
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

  // Also, get a reference to the ExplicitSystem
  ExplicitSystem &J_system = es.get_system<ExplicitSystem>("JSystem");
  const DofMap &J_dof_map = J_system.get_dof_map();
  unsigned int J_var;
  J_var = J_system.variable_number("JVar");

  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
      mesh.active_local_elements_end();

  // int count_Jn = 0;

  for (; el != end_el; ++el) {
    const Elem *elem = *el;

    fe->reinit(elem);

    GeomPar::compute_geoPar(es, elem, 0, phi, dphi);

    const int dof_index = elem->dof_number(J_system.number(), 0, 0);

    J_system.solution->set(dof_index, GeomPar::detF);
  }
  // cout<<count_Jn<<endl;

  // Should call close and update when we set vector entries directly
  J_system.solution->close();
  J_system.solution->localize(*J_system.current_local_solution);
}

void LargeDeformationElasticity::compute_pmono() {
  const MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  NonlinearImplicitSystem &system =
      es.get_system<NonlinearImplicitSystem>("NonlinearElasticity");

  unsigned int displacement_vars[dim];
  displacement_vars[0] = system.variable_number("u");
  displacement_vars[1] = system.variable_number("v");
  if (dim == 3)
    displacement_vars[2] = system.variable_number("w");
  const unsigned int u_var = system.variable_number("u");

  const DofMap &dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(u_var);
  UniquePtr<FEBase> fe(FEBase::build(dim, fe_type));
  QGauss qrule(dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule(&qrule);

  const std::vector<Real> &JxW = fe->get_JxW();
  const std::vector<std::vector<Real>> &phi = fe->get_phi();
  const std::vector<std::vector<RealGradient>> &dphi = fe->get_dphi();

  // Also, get a reference to the ExplicitSystem
  ExplicitSystem &pmono_system = es.get_system<ExplicitSystem>("pMonoSystem");
  const DofMap &pmono_dof_map = pmono_system.get_dof_map();
  unsigned int pmono_var;
  pmono_var = pmono_system.variable_number("pressureMono");

  // Storage for the stress dof indices on each element
  std::vector<std::vector<dof_id_type>> dof_indices_var(dim + 1);
  std::vector<dof_id_type> pmono_dof_indices_var;

  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
      mesh.active_local_elements_end();

  for (; el != end_el; ++el) {
    const Elem *elem = *el;

    for (unsigned int var = 0; var < dim + 1; var++)
      dof_map.dof_indices(elem, dof_indices_var[var], var);

    const unsigned int n_var_dofs = dof_indices_var[0].size();

    fe->reinit(elem);

    Number pmono_cur = 0.0;

    for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {
      double pmono_ = 0.0;

      for (unsigned int j = 0; j < n_var_dofs; j++)
        pmono_ += phi[j][qp] * system.current_solution(dof_indices_var[dim][j]);

      pmono_cur += -pmono_ * JxW[qp];
      // cout<<"J="<<F.det()<<endl;
    }

    pmono_cur *= (1. / elem->volume());

    const int dof_index = elem->dof_number(pmono_system.number(), 0, 0);

    pmono_system.solution->set(dof_index, pmono_cur);
  }

  // Should call close and update when we set vector entries directly
  pmono_system.solution->close();
  pmono_system.update();
}

PREITR LargeDeformationElasticity::pre_itr = NULL;
PARAMQPRESID LargeDeformationElasticity::param_qp_resid = NULL;
PARAMQPJACOB LargeDeformationElasticity::param_qp_jacob = NULL;
PARAMDOFIRESID LargeDeformationElasticity::param_dofi_resid = NULL;
PARAMDOFIJACOB LargeDeformationElasticity::param_dofi_jacob = NULL;
PARAMDOFIDOFJJACOB LargeDeformationElasticity::param_dofi_dofj_jacob = NULL;
BOUNDARYFORCE LargeDeformationElasticity::boundary_force = NULL;
TORSIONRESID LargeDeformationElasticity::torsion_resid = NULL;
TORSIONJACOB LargeDeformationElasticity::torsion_jacob = NULL;
