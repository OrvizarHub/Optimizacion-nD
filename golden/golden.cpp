#include "..\linalg\linalg.h"
#include "..\linalg\matrix.h"
#include "..\problem\problem.h"

#include "golden.h"

using namespace std;

// Golden-section search
double Minimizar(Matrix x, Matrix sb)
{
    const double tau = 0.381967;
    const double epsilon = 1e-3;

    double a = -5, b = 5;

    double alpha1 = a*(1-tau) + b*tau;
    double alpha2 = a*tau + b*(1-tau);

    double U_alpha1 = CostFunction(SumarMM(x, MultEM(alpha1, sb)));
    double U_alpha2 = CostFunction(SumarMM(x, MultEM(alpha2, sb)));
 
    while (true)
    {
        if(U_alpha1 > U_alpha2)
        {
            a = alpha1;
            alpha1 = alpha2;
            U_alpha1 = U_alpha2;
            alpha2 = tau*a + (1-tau)*b;
            U_alpha2 = CostFunction(SumarMM(x, MultEM(alpha2, sb)));
        }
        else
        {
            b = alpha2;
            alpha2 = alpha1;
            U_alpha2 = U_alpha1;
            alpha1 = tau*b + (1-tau)*a;
            U_alpha1 = CostFunction(SumarMM(x, MultEM(alpha1, sb)));
        }
        
        if(abs(U_alpha1 - U_alpha2) < epsilon)
        {
            break;
        }
    }
        
    return alpha1;
}