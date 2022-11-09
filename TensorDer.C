
#include "TensorDer.h"

tensorDer TensorDer::tensorder;

DenseVector<double> TensorDer::gradNA, TensorDer::gradNB;

DenseMatrix<double> TensorDer::gradNAx, TensorDer::gradNAy;
#if (MESH_DIMENSION == 3)
DenseMatrix<double> TensorDer::gradNAz;
#endif

DenseMatrix<double> TensorDer::dJdF, TensorDer::dFdux, TensorDer::dFduy,
    TensorDer::dFduz;
DenseMatrix<double> TensorDer::dJFinvTdux, TensorDer::dJFinvTduy,
    TensorDer::dJFinvTduz;

using namespace std;

// Bring in everything from the libMesh namespace
using namespace libMesh;

TensorDer::TensorDer() {}
TensorDer::~TensorDer() {}

void TensorDer::compute_gradNA(
    const std::vector<std::vector<RealGradient>> &dphi, unsigned int qp,
    unsigned int dof_i) {

  gradNA.resize(MESH_DIMENSION);
  for (unsigned int i = 0; i < MESH_DIMENSION; i++) {
    gradNA(i) = dphi[dof_i][qp](i);
  }
}

void TensorDer::compute_gradNAxyz(
    const std::vector<std::vector<RealGradient>> &dphi, unsigned int qp,
    unsigned int dof_i) {

  gradNAx.resize(MESH_DIMENSION, MESH_DIMENSION);
  gradNAy.resize(MESH_DIMENSION, MESH_DIMENSION);
#if (MESH_DIMENSION == 3)
  gradNAz.resize(MESH_DIMENSION, MESH_DIMENSION);
#endif
  for (unsigned int i = 0; i < MESH_DIMENSION; i++) {
    gradNAx(0, i) = dphi[dof_i][qp](i);
    gradNAy(1, i) = dphi[dof_i][qp](i);
#if (MESH_DIMENSION == 3)
    gradNAz(2, i) = dphi[dof_i][qp](i);
#endif
  }
}

void TensorDer::compute_gradNB(
    const std::vector<std::vector<RealGradient>> &dphi, unsigned int qp,
    unsigned int dof_j) {

  gradNB.resize(MESH_DIMENSION);
  for (unsigned int i = 0; i < MESH_DIMENSION; i++) {
    gradNB(i) = dphi[dof_j][qp](i);
  }
}

void TensorDer::compute_dJdF() {
  dJdF.resize(MESH_DIMENSION, MESH_DIMENSION);
  dJdF.add(GeomPar::detF, GeomPar::FInvTra);
}

void TensorDer::compute_dJFinvTdF() {
  for (unsigned int i = 0; i < MESH_DIMENSION; i++)
    for (unsigned int j = 0; j < MESH_DIMENSION; j++)
      for (unsigned int k = 0; k < MESH_DIMENSION; k++)
        for (unsigned int l = 0; l < MESH_DIMENSION; l++) {
          tensorder.dJFinvTdF[i][j][k][l] = GeomPar::FInvTra(i, j) * dJdF(k, l);
          tensorder.dJFinvTdF[i][j][k][l] +=
              -0.5 * GeomPar::detF *
              (GeomPar::FInvTra(k, i) * GeomPar::FInvTra(l, j) +
               GeomPar::FInvTra(l, i) * GeomPar::FInvTra(k, j));
        }
}

void TensorDer::compute_dJFinvTdu() {
  dFdux.resize(MESH_DIMENSION, MESH_DIMENSION);
  dFduy.resize(MESH_DIMENSION, MESH_DIMENSION);
  dFduz.resize(MESH_DIMENSION, MESH_DIMENSION);

  for (unsigned int i = 0; i < MESH_DIMENSION; i++) {
    dFdux(0, i) = gradNB(i);
    dFduy(1, i) = gradNB(i);
    if (MESH_DIMENSION == 3)
      dFduz(2, i) = gradNB(i);
  }

  dJFinvTdux.resize(MESH_DIMENSION, MESH_DIMENSION);
  dJFinvTduy.resize(MESH_DIMENSION, MESH_DIMENSION);
  dJFinvTduz.resize(MESH_DIMENSION, MESH_DIMENSION);

  for (unsigned int i = 0; i < MESH_DIMENSION; i++) {
    for (unsigned int j = 0; j < MESH_DIMENSION; j++) {
      for (unsigned int k = 0; k < MESH_DIMENSION; k++) {
        for (unsigned int l = 0; l < MESH_DIMENSION; l++) {
          dJFinvTdux(i, j) += tensorder.dJFinvTdF[i][j][k][l] * dFdux(k, l);
          dJFinvTduy(i, j) += tensorder.dJFinvTdF[i][j][k][l] * dFduy(k, l);
#if (MESHDIMENSION == 3)
          dJFinvTduz(i, j) += tensorder.dJFinvTdF[i][j][k][l] * dFduz(k, l);
#endif
        }
      }
    }
  }
}
