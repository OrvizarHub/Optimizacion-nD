#pragma once

#include "matrix.h"

using namespace std;

void ImprM(Matrix x);
double Norma(Matrix x);
Matrix InverM(Matrix x);
Matrix TransM(Matrix x);
Matrix MultEM(double x, Matrix y);
Matrix MultMM(Matrix x, Matrix y);
Matrix DivME(Matrix x, double y);
Matrix SumarMM(Matrix x, Matrix y);
Matrix RestarMM(Matrix x, Matrix y);
Matrix Iden(int x);