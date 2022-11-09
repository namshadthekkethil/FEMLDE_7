
#include "MatVecOper.h"

using namespace std;

// Bring in everything from the libMesh namespace
using namespace libMesh;

MatVecOper::MatVecOper() {}

MatVecOper::~MatVecOper() {}

void MatVecOper::compute_diadic(DenseVector<double> &a, DenseVector<double> &b,
                                DenseMatrix<double> &c) {
  c.resize(MESH_DIMENSION, MESH_DIMENSION);
  for (unsigned int i = 0; i < MESH_DIMENSION; i++) {
    for (unsigned int j = 0; j < MESH_DIMENSION; j++) {
      c(i, j) = a(i) * b(j);
    }
  }
}

double MatVecOper::compute_fibre_invariants(DenseVector<double> &a0,
                                            DenseMatrix<double> &C,
                                            DenseVector<double> &b0) {
  double invariant = 0.0;
  for (unsigned int i = 0; i < MESH_DIMENSION; i++) {
    for (unsigned int j = 0; j < MESH_DIMENSION; j++) {
      invariant += b0(i) * a0(j) * C(j, i);
    }
  }
  return invariant;
}

double MatVecOper::contractMat(DenseMatrix<double> &A, DenseMatrix<double> &B) {
  double contract_sum = 0.0;
  for (unsigned int i = 0; i < MESH_DIMENSION; i++) {
    for (unsigned int j = 0; j < MESH_DIMENSION; j++) {
      contract_sum += A(i, j) * B(i, j);
    }
  }
  return contract_sum;
}

double MatVecOper::contractVec(DenseVector<double> &A, DenseVector<double> &B) {
  double contract_sum = 0.0;
  for (unsigned int i = 0; i < MESH_DIMENSION; i++) {
    contract_sum += A(i) * B(i);
  }
  return contract_sum;
}

double MatVecOper::detMat(DenseMatrix<double> &A) {
  DenseMatrix<double> Atemp;
  Atemp.resize(MESH_DIMENSION, MESH_DIMENSION);
  Atemp.add(1.0, A);
  double detA = Atemp.det();
  return detA;
}

void MatVecOper::inverseMat(DenseMatrix<double> &A, DenseMatrix<double> &invA) {
  DenseMatrix<double> Atemp;
  Atemp.resize(MESH_DIMENSION, MESH_DIMENSION);
  Atemp.add(1.0, A);
  // cout<<"THAT IS OK THAT IS OK THAT IS OK THAT IS OK"<<endl;

  double detA = Atemp.det();
  // cout<<detA<<" THAT IS OK THAT IS OK THAT IS OK THAT IS OK"<<endl;

  if (MESH_DIMENSION == 2) {
    invA(0, 0) = A(1, 1);
    invA(0, 1) = -A(0, 1);
    invA(1, 0) = -A(1, 0);
    invA(1, 1) = A(0, 0);
  }

  if (MESH_DIMENSION == 3) {
    invA(0, 0) = (A(1, 1) * A(2, 2) - A(2, 1) * A(1, 2));
    invA(0, 1) = (A(0, 2) * A(2, 1) - A(0, 1) * A(2, 2));
    invA(0, 2) = (A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1));
    invA(1, 0) = (A(1, 2) * A(2, 0) - A(1, 0) * A(2, 2));
    invA(1, 1) = (A(0, 0) * A(2, 2) - A(0, 2) * A(2, 0));
    invA(1, 2) = (A(1, 0) * A(0, 2) - A(0, 0) * A(1, 2));
    invA(2, 0) = (A(1, 0) * A(2, 1) - A(2, 0) * A(1, 1));
    invA(2, 1) = (A(2, 0) * A(0, 1) - A(0, 0) * A(2, 1));
    invA(2, 2) = (A(0, 0) * A(1, 1) - A(1, 0) * A(0, 1));
  }
  invA.scale(1.0 / detA);
}

void MatVecOper::transposeMat(DenseMatrix<double> &A, DenseMatrix<double> &At) {
  if (MESH_DIMENSION == 2) {
    At(0, 0) = A(0, 0);
    At(0, 1) = A(1, 0);
    At(1, 0) = A(0, 1);
    At(1, 1) = A(1, 1);
  }

  if (MESH_DIMENSION == 3) {
    At(0, 0) = A(0, 0);
    At(0, 1) = A(1, 0);
    At(0, 2) = A(2, 0);
    At(1, 0) = A(0, 1);
    At(1, 1) = A(1, 1);
    At(1, 2) = A(2, 1);
    At(2, 0) = A(0, 2);
    At(2, 1) = A(1, 2);
    At(2, 2) = A(2, 2);
  }
}

void MatVecOper::crossMat(DenseMatrix<double> &A, DenseMatrix<double> &B,
                          DenseMatrix<double> &C) {
  C(0, 0) = A(1, 1) * B(2, 2) - A(1, 2) * B(2, 1) + A(2, 2) * B(1, 1) -
            A(2, 1) * B(1, 2);
  C(0, 1) = A(1, 2) * B(2, 0) - A(1, 0) * B(2, 2) + A(2, 0) * B(1, 2) -
            A(2, 2) * B(1, 0);
  C(0, 2) = A(1, 0) * B(2, 1) - A(1, 1) * B(2, 0) + A(2, 1) * B(1, 0) -
            A(2, 0) * B(1, 1);

  C(1, 0) = A(0, 2) * B(2, 1) - A(0, 1) * B(2, 2) + A(2, 1) * B(0, 2) -
            A(2, 2) * B(0, 1);
  C(1, 1) = A(2, 2) * B(0, 0) - A(2, 0) * B(0, 2) + A(0, 0) * B(2, 2) -
            A(0, 2) * B(2, 0);
  C(1, 2) = A(2, 0) * B(0, 1) - A(2, 1) * B(0, 0) + A(0, 1) * B(2, 0) -
            A(0, 0) * B(2, 1);

  C(2, 0) = A(0, 1) * B(1, 2) - A(0, 2) * B(1, 1) + A(1, 2) * B(0, 1) -
            A(1, 1) * B(0, 2);
  C(2, 1) = A(0, 2) * B(1, 0) - A(0, 0) * B(1, 2) + A(1, 0) * B(0, 2) -
            A(1, 2) * B(0, 0);
  C(2, 2) = A(0, 0) * B(1, 1) - A(0, 1) * B(1, 0) + A(1, 1) * B(0, 0) -
            A(1, 0) * B(0, 1);
}
