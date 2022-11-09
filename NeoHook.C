
#include "NeoHook.h"

double NeoHook::G;
double NeoHook::Jm23, NeoHook::I1, NeoHook::I1bar, NeoHook::dWdI1bar;

DenseMatrix<double> NeoHook::dJm23dE, NeoHook::dI1dE, NeoHook::dI1bardE,
    NeoHook::SIso;

DenseMatrix<double> NeoHook::S;

DenseMatrix<double> NeoHook::Sp;

structNH NeoHook::Dmat;

double NeoHook::dWdI1;

double NeoHook::mu;

using namespace std;

// Bring in everything from the libMesh namespace
using namespace libMesh;

NeoHook::NeoHook() {}
NeoHook::~NeoHook() {}

void NeoHook::compute_NH_PK2(EquationSystems &es, const Elem *elem) {

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

  compute_NH_PK2_iso();

  S.resize(MESH_DIMENSION, MESH_DIMENSION);
  S.add(1.0, SIso);

  Sp.resize(MESH_DIMENSION, MESH_DIMENSION);
  Sp.add(GeomPar::pressure * GeomPar::detF, GeomPar::CInv);
}

void NeoHook::compute_NH_D() {
  compute_NH_D_iso();
  for (unsigned int i = 0; i < MESH_DIMENSION; i++)
    for (unsigned int j = 0; j < MESH_DIMENSION; j++)
      for (unsigned int k = 0; k < MESH_DIMENSION; k++)
        for (unsigned int l = 0; l < MESH_DIMENSION; l++) {
          Dmat.D[i][j][k][l] = Dmat.Diso[i][j][k][l];
        }
  compute_D_p();
}

void NeoHook::compute_NH_PK2_iso() {
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

  dWdI1bar = (G / 2.0);

  compute_PK2_iso();
}

void NeoHook::compute_NH_D_iso() {

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

void NeoHook::compute_PK2_iso() {
  SIso.resize(MESH_DIMENSION, MESH_DIMENSION);
  SIso.add(dWdI1bar, dI1bardE);
}

void NeoHook::compute_D_iso() {
  for (unsigned int i = 0; i < MESH_DIMENSION; i++)
    for (unsigned int j = 0; j < MESH_DIMENSION; j++)
      for (unsigned int k = 0; k < MESH_DIMENSION; k++)
        for (unsigned int l = 0; l < MESH_DIMENSION; l++) {
          Dmat.Diso[i][j][k][l] = dWdI1bar * Dmat.d2I1bardE2[i][j][k][l];
        }
}

void NeoHook::compute_D_p() {
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

void NeoHook::compute_mu(EquationSystems &es, const Elem *elem) {

  dWdI1 = G / 2.0;

  mu = dWdI1;
}
