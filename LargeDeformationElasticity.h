// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

// <h1> Systems Example 7 - Large deformation elasticity (St. Venant-Kirchoff
// material) </h1> \author Lorenzo Zanon \author David Knezevic \date 2014
//
// In this example, we consider an elastic cantilever beam modeled as a St.
// Venant-Kirchoff material (which is an extension of the linear elastic
// material model to the nonlinear regime). The implementation presented here
// uses NonlinearImplicitSystem.
//
// We formulate the PDE on the reference geometry (\Omega) as opposed to the
// deformed geometry (\Omega^deformed). As a result (e.g. see Ciarlet's 3D
// elasticity book, Theorem 2.6-2) the PDE is given as follows:
//
//     \int_\Omega F_im Sigma_mj v_i,j = \int_\Omega f_i v_i + \int_\Gamma g_i
//     v_i ds
//
// where:
//  * F is the deformation gradient, F = I + du/dx (x here refers to reference
//  coordinates).
//  * Sigma is the second Piola-Kirchoff stress, which for the St. Venant
//  Kirchoff model is
//    given by Sigma_ij = C_ijkl E_kl, where E_kl is the strain,
//    E_kl = 0.5 * (u_k,l + u_l,k + u_m,k u_m,l).
//  * f is a body load.
//  * g is a surface traction on the surface \Gamma.
//
// In this example we only consider a body load (e.g. gravity), hence we set g =
// 0.

#ifndef included_LargeDeformationElasticity
#define included_LargeDeformationElasticity

// C++ include files that we need
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

// Various include files needed for the mesh & solver functionality.
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/equation_systems.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/fe.h"
#include "libmesh/getpot.h"
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/zero_function.h"

// The nonlinear solver and system we will be using
#include "libmesh/linear_implicit_system.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/nonlinear_solver.h"

#include <libmesh/mesh_modification.h>
#include <libmesh/mesh_tools.h>

// #include "DarcySolver.h"
#include "MatVecOper.h"
#include "Stabilisation.h"

// class LargeDeformationElasticity;
typedef void (*PREITR)(EquationSystems &es);
typedef void (*PARAMQPRESID)(
    EquationSystems &es, const Elem *elem, unsigned int qp,
    const std::vector<std::vector<Real>> &phi,
    const std::vector<std::vector<RealGradient>> &dphi);
typedef void (*PARAMQPJACOB)(
    EquationSystems &es, const Elem *elem, unsigned int qp,
    const std::vector<std::vector<Real>> &phi,
    const std::vector<std::vector<RealGradient>> &dphi);
typedef void (*PARAMDOFIRESID)(
    const std::vector<std::vector<Real>> &phi,
    const std::vector<std::vector<RealGradient>> &dphi, unsigned int qp,
    unsigned int dof_i);
typedef void (*PARAMDOFIJACOB)(
    const std::vector<std::vector<Real>> &phi,
    const std::vector<std::vector<RealGradient>> &dphi, unsigned int qp,
    unsigned int dof_i);
typedef void (*PARAMDOFIDOFJJACOB)(
    const std::vector<std::vector<Real>> &phi,
    const std::vector<std::vector<RealGradient>> &dphi, unsigned int qp,
    unsigned int dof_i, unsigned int dof_j);

typedef void (*BOUNDARYFORCE)(
    EquationSystems &es, const Elem *elem, unsigned int side,
    const std::vector<std::vector<Real>> &phi_face,
    const std::vector<Real> &JxW_face, const std::vector<Point> &normal_face,
    const std::vector<std::vector<Real>> &phi_face_cur,
    const std::vector<Real> &JxW_face_cur,
    const std::vector<Point> &normal_face_cur, DenseSubVector<double> Re_var[]);
typedef void (*TORSIONRESID)(EquationSystems &es, const Elem *elem,
                             DenseSubVector<double> Re_var[]);
typedef void (*TORSIONJACOB)(
    EquationSystems &es, const Elem *elem,
    DenseSubMatrix<Number> Ke_var[][MESH_DIMENSION + 1]);

// #define MESH_DIMENSION 3

using namespace libMesh;
using namespace std;

class LargeDeformationElasticity
    : public NonlinearImplicitSystem::ComputeResidual,
      public NonlinearImplicitSystem::ComputeJacobian {
  // class LargeDeformationElasticity
  //     : public NonlinearImplicitSystem::ComputeResidualandJacobian {
private:
  EquationSystems &es;
  EquationSystems &es_cur;
  double &ttime;
  double &dt;
  DenseVector<double> &Resid;
  DenseMatrix<double> &Jacob;

public:
  LargeDeformationElasticity(EquationSystems &es_in, EquationSystems &es_in_cur,
                             double &ttime_in, double &dt_in,
                             DenseVector<double> &Resid_in,
                             DenseMatrix<double> &Jacob_in)
      : es(es_in), es_cur(es_in_cur), ttime(ttime_in), dt(dt_in),
        Resid(Resid_in), Jacob(Jacob_in) {}
  ~LargeDeformationElasticity();
  void define_systems();
  void update_nodeid();
  void move_mesh();
  void pre_solve();
  void old_new(System &system_old, System &system_new);
  void solve_lde();
  void solve_lde_stepwise();

  static PREITR pre_itr;
  static PARAMQPRESID param_qp_resid;
  static PARAMQPJACOB param_qp_jacob;
  static PARAMDOFIRESID param_dofi_resid;
  static PARAMDOFIJACOB param_dofi_jacob;
  static PARAMDOFIDOFJJACOB param_dofi_dofj_jacob;
  static BOUNDARYFORCE boundary_force;
  static TORSIONRESID torsion_resid;
  static TORSIONJACOB torsion_jacob;

  virtual void jacobian(const NumericVector<Number> &soln,
                        SparseMatrix<Number> &jacobian,
                        NonlinearImplicitSystem & /*sys*/) {
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

    DenseMatrix<Number> Ke;
#if (MESH_DIMENSION == 2)
    DenseSubMatrix<Number> Ke_var[3][3] = {
        {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke),
         DenseSubMatrix<Number>(Ke)},
        {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke),
         DenseSubMatrix<Number>(Ke)},
        { DenseSubMatrix<Number>(Ke),
          DenseSubMatrix<Number>(Ke),
          DenseSubMatrix<Number>(Ke) }};
#elif (MESH_DIMENSION == 3)
    DenseSubMatrix<Number> Ke_var[4][4] = {
        {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke),
         DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)},
        {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke),
         DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)},
        {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke),
         DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)},
        {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke),
         DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)}};
#endif

    std::vector<dof_id_type> dof_indices;
    std::vector<std::vector<dof_id_type>> dof_indices_var(dim + 1);

    jacobian.zero();

    MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el =
        mesh.active_local_elements_end();

    // std::cout<<"before jacobian"<<std::endl;

    for (; el != end_el; ++el) {
      const Elem *elem = *el;

      dof_map.dof_indices(elem, dof_indices);
      for (unsigned int var = 0; var < dim + 1; var++) {
        dof_map.dof_indices(elem, dof_indices_var[var], var);
      }

      const unsigned int n_dofs = dof_indices.size();
      const unsigned int n_var_dofs = dof_indices_var[0].size();

      fe->reinit(elem);

      Ke.resize(n_dofs, n_dofs);
      for (unsigned int var_i = 0; var_i < dim + 1; var_i++)
        for (unsigned int var_j = 0; var_j < dim + 1; var_j++)
          Ke_var[var_i][var_j].reposition(
              var_i * n_var_dofs, var_j * n_var_dofs, n_var_dofs, n_var_dofs);

      for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {

        param_qp_jacob(es, elem, qp, phi, dphi);

        for (unsigned int dof_i = 0; dof_i < n_var_dofs; dof_i++) {

          param_dofi_jacob(phi, dphi, qp, dof_i);

          for (unsigned int dof_j = 0; dof_j < n_var_dofs; dof_j++) {

            param_dofi_dofj_jacob(phi, dphi, qp, dof_i, dof_j);

            for (unsigned int var_i = 0; var_i < dim + 1; var_i++)
              for (unsigned int var_j = 0; var_j < dim + 1; var_j++) {
                Ke_var[var_i][var_j](dof_i, dof_j) +=
                    Jacob(var_i, var_j) * JxW[qp];
              }
          }
        }
      }

      torsion_jacob(es, elem, Ke_var);

      dof_map.constrain_element_matrix(Ke, dof_indices);
      jacobian.add_matrix(Ke, dof_indices);
    }
    // std::cout<<"after jacobian"<<std::endl;
  }

  /**
   * Evaluate the residual of the nonlinear system.
   */
  virtual void residual(const NumericVector<Number> &soln,
                        NumericVector<Number> &residual,
                        NonlinearImplicitSystem & /*sys*/) {
    const MeshBase &mesh = es.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    pre_itr(es);

    cout << "t_stab=" << Stabilisation::tau_stab << endl;

    NonlinearImplicitSystem &system =
        es.get_system<NonlinearImplicitSystem>("NonlinearElasticity");

    const unsigned int u_var = system.variable_number("u");

    const DofMap &dof_map = system.get_dof_map();

    FEType fe_type = dof_map.variable_type(u_var);
    UniquePtr<FEBase> fe(FEBase::build(dim, fe_type));
    QGauss qrule(dim, CONSTANT);
    fe->attach_quadrature_rule(&qrule);

    UniquePtr<FEBase> fe_face(FEBase::build(dim, fe_type));
    QGauss qface(dim - 1, CONSTANT);
    fe_face->attach_quadrature_rule(&qface);

    /************************************ for current
     * ****************************/
    const MeshBase &mesh_cur = es_cur.get_mesh();
    System &system_cur = es_cur.get_system<System>("NodeIdSystem");

    const unsigned int u_var_cur = system_cur.variable_number("NodeIdVar");

    const DofMap &dof_map_cur = system_cur.get_dof_map();

    FEType fe_type_cur = dof_map_cur.variable_type(u_var_cur);
    UniquePtr<FEBase> fe_cur(FEBase::build(dim, fe_type_cur));
    QGauss qrule_cur(dim, CONSTANT);
    fe_cur->attach_quadrature_rule(&qrule_cur);

    UniquePtr<FEBase> fe_face_cur(FEBase::build(dim, fe_type_cur));
    QGauss qface_cur(dim - 1, CONSTANT);
    fe_face_cur->attach_quadrature_rule(&qface_cur);
    /****************************************************************************/

    const std::vector<Real> &JxW = fe->get_JxW();
    const std::vector<std::vector<Real>> &phi = fe->get_phi();
    const std::vector<std::vector<RealGradient>> &dphi = fe->get_dphi();

    const std::vector<std::vector<Real>> &phi_face = fe_face->get_phi();
    const std::vector<Real> &JxW_face = fe_face->get_JxW();

    const std::vector<Point> &normal_face = fe_face->get_normals();

    const std::vector<std::vector<Real>> &phi_face_cur = fe_face_cur->get_phi();
    const std::vector<Real> &JxW_face_cur = fe_face_cur->get_JxW();
    const std::vector<Point> &normal_face_cur = fe_face_cur->get_normals();

#if (MESH_DIMENSION == 2)
    DenseVector<Number> Re;
    DenseSubVector<Number> Re_var[3] = {DenseSubVector<Number>(Re),
                                        DenseSubVector<Number>(Re),
                                        DenseSubVector<Number>(Re)};
#elif (MESH_DIMENSION == 3)
    DenseVector<Number> Re;
    DenseSubVector<Number> Re_var[4] = {
        DenseSubVector<Number>(Re), DenseSubVector<Number>(Re),
        DenseSubVector<Number>(Re), DenseSubVector<Number>(Re)};
#endif

    std::vector<dof_id_type> dof_indices;
    std::vector<std::vector<dof_id_type>> dof_indices_var(dim + 1);

    residual.zero();

    MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el =
        mesh.active_local_elements_end();

    MeshBase::const_element_iterator el_cur =
        mesh_cur.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el_cur =
        mesh_cur.active_local_elements_end();

    for (; el != end_el; ++el) {
      const Elem *elem = *el;

      const Elem *elem_cur = *el_cur;

      dof_map.dof_indices(elem, dof_indices);
      for (unsigned int var = 0; var < dim + 1; var++) {
        dof_map.dof_indices(elem, dof_indices_var[var], var);
      }

      const unsigned int n_dofs = dof_indices.size();
      const unsigned int n_var_dofs = dof_indices_var[0].size();

      fe->reinit(elem);
      fe_cur->reinit(elem_cur);

      Re.resize(n_dofs);
      for (unsigned int var = 0; var < dim + 1; var++)
        Re_var[var].reposition(var * n_var_dofs, n_var_dofs);

      for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {

        param_qp_resid(es, elem, qp, phi, dphi);

        for (unsigned int dof_i = 0; dof_i < n_var_dofs; dof_i++) {
          param_dofi_resid(phi, dphi, qp, dof_i);

          for (unsigned int i = 0; i < dim + 1; i++) {
            Re_var[i](dof_i) += Resid(i) * JxW[qp];
          }
        }
      }

      for (unsigned int side = 0; side < elem->n_sides(); side++)
        if (elem->neighbor_ptr(side) == libmesh_nullptr) {

          fe_face->reinit(elem, side);
          fe_face_cur->reinit(elem_cur, side);

          boundary_force(es, elem, side, phi_face, JxW_face, normal_face,
                         phi_face_cur, JxW_face_cur, normal_face_cur, Re_var);
        }

      torsion_resid(es, elem, Re_var);

      dof_map.constrain_element_vector(Re, dof_indices);
      residual.add_vector(Re, dof_indices);

      ++el_cur;
    }
  }

  void compute_stresses();
  void compute_J();
  void compute_pmono();

  // unsigned int prssr_loading,act_loading;
};

#endif
