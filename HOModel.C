
#include "HOModel.h"

double HOModel::aHO, HOModel::bHO, HOModel::afHO, HOModel::bfHO, HOModel::asHO,
    HOModel::bsHO, HOModel::afsHO, HOModel::bfsHO;
double HOModel::Jm23, HOModel::I1, HOModel::I1bar, HOModel::dWdI1bar,
    HOModel::d2WdI1bar2;

DenseMatrix<double> HOModel::dJm23dE, HOModel::dI1dE, HOModel::dI1bardE,
    HOModel::SIso;

double HOModel::I4f, HOModel::I4fbar, HOModel::dWdI4fbar, HOModel::d2WdI4fbar2;
DenseMatrix<double> HOModel::dI4fdE, HOModel::dI4fbardE, HOModel::Sf;

double HOModel::I4s, HOModel::I4sbar, HOModel::dWdI4sbar, HOModel::d2WdI4sbar2;
DenseMatrix<double> HOModel::dI4sdE, HOModel::dI4sbardE, HOModel::Ss;

DenseMatrix<double> HOModel::S;

double HOModel::I8fs, HOModel::I8fsbar, HOModel::dWdI8fsbar,
    HOModel::d2WdI8fsbar2;
DenseMatrix<double> HOModel::dI8fsdE, HOModel::dI8fsbardE, HOModel::Sfs;

DenseMatrix<double> HOModel::Sp;

structHO HOModel::Dmat;

double HOModel::dWdI1, HOModel::dWdI4f, HOModel::dWdI4s, HOModel::dWdI8fs;

double HOModel::mu;

using namespace std;

// Bring in everything from the libMesh namespace
using namespace libMesh;

HOModel::HOModel() {}
HOModel::~HOModel() {}

void HOModel::compute_HO_PK2(EquationSystems &es, const Elem *elem) {

  Jm23 = 0.0;
  dJm23dE.resize(MESH_DIMENSION, MESH_DIMENSION);

#if (MESH_DIMENSION == 2)
  Jm23 = pow(GeomPar::detF, -1);
  dJm23dE.add(-Jm23, GeomPar::CInv);
#elif (MESH_DIMENSION == 3)
  Jm23 = pow(GeomPar::detF, -2 / 3.0);
  dJm23dE.add(-(2.0 / 3.0) * Jm23, GeomPar::CInv);
#endif

  for (unsigned int i = 0; i < MESH_DIMENSION; i++)
    for (unsigned int j = 0; j < MESH_DIMENSION; j++)
      for (unsigned int k = 0; k < MESH_DIMENSION; k++)
        for (unsigned int l = 0; l < MESH_DIMENSION; l++) {
          Dmat.dCdE[i][j][k][l] = -(GeomPar::CInv(i, k) * GeomPar::CInv(j, l) +
                                    GeomPar::CInv(i, l) * GeomPar::CInv(j, k));

#if (MESH_DIMENSION == 2)
          Dmat.d2Jm23dE2[i][j][k][l] = -(GeomPar::CInv(i, j) * dJm23dE(k, l) +
                                         Jm23 * Dmat.dCdE[i][j][k][l]);
#elif (MESH_DIMENSION == 3)
          Dmat.d2Jm23dE2[i][j][k][l] =
              -(2.0 / 3.0) * (GeomPar::CInv(i, j) * dJm23dE(k, l) +
                              Jm23 * Dmat.dCdE[i][j][k][l]);
#endif
        }

  compute_HO_PK2_iso();
  compute_HO_PK2_aniso(es, elem);

  S.resize(MESH_DIMENSION, MESH_DIMENSION);
  S.add(1.0, SIso);
  S.add(1.0, Sf);
  S.add(1.0, Ss);
  S.add(1.0, Sfs);

  Sp.resize(MESH_DIMENSION, MESH_DIMENSION);
  Sp.add(GeomPar::pressure * GeomPar::detF, GeomPar::CInv);
}

void HOModel::compute_HO_D() {
  compute_HO_D_iso();
  compute_HO_D_aniso();
  for (unsigned int i = 0; i < MESH_DIMENSION; i++)
    for (unsigned int j = 0; j < MESH_DIMENSION; j++)
      for (unsigned int k = 0; k < MESH_DIMENSION; k++)
        for (unsigned int l = 0; l < MESH_DIMENSION; l++) {
          Dmat.D[i][j][k][l] = Dmat.Diso[i][j][k][l] + Dmat.Df[i][j][k][l] +
                               Dmat.Ds[i][j][k][l] + Dmat.Dfs[i][j][k][l];
        }
  compute_D_p();
}

void HOModel::compute_HO_PK2_iso() {
  I1 = 0.0;
  I1 = GeomPar::C(0, 0) + GeomPar::C(1, 1);
  if (MESH_DIMENSION == 3)
    I1 += GeomPar::C(2, 2);

  I1bar = I1 * Jm23;

  dI1dE.resize(MESH_DIMENSION, MESH_DIMENSION);
  dI1bardE.resize(MESH_DIMENSION, MESH_DIMENSION);

  dI1dE.add(2.0, GeomPar::I);

  dI1bardE.add(I1, dJm23dE);
  dI1bardE.add(Jm23, dI1dE);

#if (MESH_DIMENSION == 2)
  dWdI1bar = (aHO / 2.0) * exp(bHO * (I1bar - 2.0));
#elif (MESH_DIMENSION == 3)
  dWdI1bar = (aHO / 2.0) * exp(bHO * (I1bar - 3.0));
#endif

  compute_PK2_iso();
}

void HOModel::compute_HO_D_iso() {
#if (MESH_DIMENSION == 2)
  d2WdI1bar2 = (aHO / 2.0) * exp(bHO * (I1bar - 2.0)) * bHO;
#elif (MESH_DIMENSION == 3)
  d2WdI1bar2 = (aHO / 2.0) * exp(bHO * (I1bar - 3.0)) * bHO;
#endif

  for (unsigned int i = 0; i < MESH_DIMENSION; i++)
    for (unsigned int j = 0; j < MESH_DIMENSION; j++)
      for (unsigned int k = 0; k < MESH_DIMENSION; k++)
        for (unsigned int l = 0; l < MESH_DIMENSION; l++) {
          Dmat.d2I1bardE2[i][j][k][l] = Dmat.d2Jm23dE2[i][j][k][l] * I1 +
                                        dJm23dE(i, j) * dI1dE(k, l) +
                                        dI1dE(i, j) * dJm23dE(k, l);
        }

  compute_D_iso();
}

void HOModel::compute_HO_PK2_aniso(EquationSystems &es, const Elem *elem) {
  System &f0_system = es.get_system<System>("fiber direction");
  NumericVector<double> &f0_vec = *f0_system.solution;
  int f0_system_num = f0_system.number();

  DenseVector<Number> f0_e(MESH_DIMENSION);
  for (unsigned int d = 0; d < MESH_DIMENSION; ++d) {
    const int dof_index = elem->dof_number(f0_system_num, d, 0);
    f0_e(d) = f0_vec(dof_index);
  }

  System &s0_system = es.get_system<System>("sheet direction");
  NumericVector<double> &s0_vec = *s0_system.solution;
  int s0_system_num = s0_system.number();

  DenseVector<Number> s0_e(MESH_DIMENSION);
  for (unsigned int d = 0; d < MESH_DIMENSION; ++d) {
    const int dof_index = elem->dof_number(s0_system_num, d, 0);
    s0_e(d) = s0_vec(dof_index);
  }

  I4f = MatVecOper::compute_fibre_invariants(f0_e, GeomPar::C, f0_e);
  I4s = MatVecOper::compute_fibre_invariants(s0_e, GeomPar::C, s0_e);
  I8fs = MatVecOper::compute_fibre_invariants(f0_e, GeomPar::C, s0_e);

  I4fbar = Jm23 * I4f;
  I4sbar = Jm23 * I4s;
  I8fsbar = Jm23 * I8fs;

  MatVecOper::compute_diadic(f0_e, f0_e, dI4fdE);
  dI4fdE.scale(2.0);
  dI4fbardE.resize(MESH_DIMENSION, MESH_DIMENSION);
  dI4fbardE.add(I4f, dJm23dE);
  dI4fbardE.add(Jm23, dI4fdE);

  MatVecOper::compute_diadic(s0_e, s0_e, dI4sdE);
  dI4sdE.scale(2.0);
  dI4sbardE.resize(MESH_DIMENSION, MESH_DIMENSION);
  dI4sbardE.add(I4s, dJm23dE);
  dI4sbardE.add(Jm23, dI4sdE);

  MatVecOper::compute_diadic(f0_e, s0_e, dI8fsdE);
  dI8fsdE.scale(2.0);
  dI8fsbardE.resize(MESH_DIMENSION, MESH_DIMENSION);
  dI8fsbardE.add(I8fs, dJm23dE);
  dI8fsbardE.add(Jm23, dI8fsdE);

  dWdI4fbar = afHO *
              exp(bfHO * (max(I4fbar, 1.0) - 1.0) * (max(I4fbar, 1.0) - 1.0)) *
              (max(I4fbar, 1.0) - 1.0);

  dWdI4sbar = asHO *
              exp(bsHO * (max(I4sbar, 1.0) - 1.0) * (max(I4sbar, 1.0) - 1.0)) *
              (max(I4sbar, 1.0) - 1.0);

  dWdI8fsbar = afsHO * exp(bfsHO * (I8fsbar) * (I8fsbar)) * (I8fsbar);

  compute_PK2_aniso();
}

void HOModel::compute_HO_D_aniso() {
  d2WdI4fbar2 =
      afHO * exp(bfHO * (max(I4fbar, 1.0) - 1.0) * (max(I4fbar, 1.0) - 1.0)) *
      (1.0 + 2.0 * bfHO * (max(I4fbar, 1.0) - 1.0) * (max(I4fbar, 1.0) - 1.0));

  d2WdI4sbar2 =
      asHO * exp(bsHO * (max(I4sbar, 1.0) - 1.0) * (max(I4sbar, 1.0) - 1.0)) *
      (1.0 + 2.0 * bsHO * (max(I4sbar, 1.0) - 1.0) * (max(I4sbar, 1.0) - 1.0));

  d2WdI8fsbar2 = afsHO * exp(bfsHO * (I8fsbar) * (I8fsbar)) *
                 (1.0 + 2.0 * bfsHO * (I8fsbar) * (I8fsbar));

  for (unsigned int i = 0; i < MESH_DIMENSION; i++)
    for (unsigned int j = 0; j < MESH_DIMENSION; j++)
      for (unsigned int k = 0; k < MESH_DIMENSION; k++)
        for (unsigned int l = 0; l < MESH_DIMENSION; l++) {
          Dmat.d2I4fbardE2[i][j][k][l] = Dmat.d2Jm23dE2[i][j][k][l] * I4f +
                                         dJm23dE(i, j) * dI4fdE(k, l) +
                                         dI4fdE(i, j) * dJm23dE(k, l);

          Dmat.d2I4sbardE2[i][j][k][l] = Dmat.d2Jm23dE2[i][j][k][l] * I4s +
                                         dJm23dE(i, j) * dI4sdE(k, l) +
                                         dI4sdE(i, j) * dJm23dE(k, l);

          Dmat.d2I8fsbardE2[i][j][k][l] = Dmat.d2Jm23dE2[i][j][k][l] * I8fs +
                                          dJm23dE(i, j) * dI8fsdE(k, l) +
                                          dI8fsdE(i, j) * dJm23dE(k, l);
        }

  compute_D_aniso();
}

void HOModel::compute_PK2_iso() {
  SIso.resize(MESH_DIMENSION, MESH_DIMENSION);
  SIso.add(dWdI1bar, dI1bardE);
}

void HOModel::compute_PK2_aniso() {
  Sf.resize(MESH_DIMENSION, MESH_DIMENSION);
  Sf.add(dWdI4fbar, dI4fbardE);

  Ss.resize(MESH_DIMENSION, MESH_DIMENSION);
  Ss.add(dWdI4sbar, dI4sbardE);

  Sfs.resize(MESH_DIMENSION, MESH_DIMENSION);
  Sfs.add(dWdI8fsbar, dI8fsbardE);
}

void HOModel::compute_D_iso() {
  for (unsigned int i = 0; i < MESH_DIMENSION; i++)
    for (unsigned int j = 0; j < MESH_DIMENSION; j++)
      for (unsigned int k = 0; k < MESH_DIMENSION; k++)
        for (unsigned int l = 0; l < MESH_DIMENSION; l++) {
          Dmat.Diso[i][j][k][l] = d2WdI1bar2 * dI1bardE(i, j) * dI1bardE(k, l) +
                                  dWdI1bar * Dmat.d2I1bardE2[i][j][k][l];
        }
}

void HOModel::compute_D_aniso() {
  for (unsigned int i = 0; i < MESH_DIMENSION; i++)
    for (unsigned int j = 0; j < MESH_DIMENSION; j++)
      for (unsigned int k = 0; k < MESH_DIMENSION; k++)
        for (unsigned int l = 0; l < MESH_DIMENSION; l++) {
          Dmat.Df[i][j][k][l] =
              d2WdI4fbar2 * dI4fbardE(i, j) * dI4fbardE(k, l) +
              dWdI4fbar * Dmat.d2I4fbardE2[i][j][k][l];

          Dmat.Ds[i][j][k][l] =
              d2WdI4sbar2 * dI4sbardE(i, j) * dI4sbardE(k, l) +
              dWdI4sbar * Dmat.d2I4sbardE2[i][j][k][l];

          Dmat.Dfs[i][j][k][l] =
              d2WdI8fsbar2 * dI8fsbardE(i, j) * dI8fsbardE(k, l) +
              dWdI8fsbar * Dmat.d2I8fsbardE2[i][j][k][l];
        }
}

void HOModel::compute_D_p() {
  for (unsigned int i = 0; i < MESH_DIMENSION; i++)
    for (unsigned int j = 0; j < MESH_DIMENSION; j++)
      for (unsigned int k = 0; k < MESH_DIMENSION; k++)
        for (unsigned int l = 0; l < MESH_DIMENSION; l++) {
          Dmat.Dp[i][j][k][l] =
              GeomPar::pressure *
              (GeomPar::CInv(i, j) * GeomPar::detF * GeomPar::CInv(k, l) +
               GeomPar::detF * Dmat.dCdE[i][j][k][l]);
        }
}

void HOModel::compute_mu(EquationSystems &es, const Elem *elem) {
  I1 = 0.0;
  I1 = GeomPar::C(0, 0) + GeomPar::C(1, 1);
  if (MESH_DIMENSION == 3)
    I1 += GeomPar::C(2, 2);

  System &f0_system = es.get_system<System>("fiber direction");
  NumericVector<double> &f0_vec = *f0_system.solution;
  int f0_system_num = f0_system.number();

  DenseVector<Number> f0_e(MESH_DIMENSION);
  for (unsigned int d = 0; d < MESH_DIMENSION; ++d) {
    const int dof_index = elem->dof_number(f0_system_num, d, 0);
    f0_e(d) = f0_vec(dof_index);
  }

  System &s0_system = es.get_system<System>("sheet direction");
  NumericVector<double> &s0_vec = *s0_system.solution;
  int s0_system_num = s0_system.number();

  DenseVector<Number> s0_e(MESH_DIMENSION);
  for (unsigned int d = 0; d < MESH_DIMENSION; ++d) {
    const int dof_index = elem->dof_number(s0_system_num, d, 0);
    s0_e(d) = s0_vec(dof_index);
  }

  I4f = MatVecOper::compute_fibre_invariants(f0_e, GeomPar::C, f0_e);
  I4s = MatVecOper::compute_fibre_invariants(s0_e, GeomPar::C, s0_e);
  I8fs = MatVecOper::compute_fibre_invariants(f0_e, GeomPar::C, s0_e);

#if (MESH_DIMENSION == 2)
  dWdI1 = (aHO / 2.0) * exp(bHO * (I1 - 2.0));
#elif (MESH_DIMENSION == 3)
  dWdI1 = (aHO / 2.0) * exp(bHO * (I1 - 3.0));
#endif

  dWdI4f = afHO * exp(bfHO * (max(I4f, 1.0) - 1.0) * (max(I4f, 1.0) - 1.0)) *
           (max(I4fbar, 1.0) - 1.0);

  dWdI4s = asHO * exp(bsHO * (max(I4s, 1.0) - 1.0) * (max(I4s, 1.0) - 1.0)) *
           (max(I4sbar, 1.0) - 1.0);

  dWdI8fs = afsHO * exp(bfsHO * (I8fs) * (I8fs)) * (I8fs);

  mu = dWdI1 + dWdI4f + dWdI4s + dWdI8fs;
}
