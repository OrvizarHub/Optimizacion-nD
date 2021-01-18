#include "..\linalg\linalg.h"
#include "..\linalg\matrix.h"

using namespace std;

double CostFunction(Matrix x);
Matrix CostFunctionGradient(Matrix x, double deltaX);
Matrix CostFunctionHessian(Matrix x, double deltaX);