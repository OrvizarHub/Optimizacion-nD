#include <math.h>

#include "..\linalg\linalg.h"
#include "..\linalg\matrix.h"

#include "problem.h"

using namespace std;

double CostFunction(Matrix x)
{
    double x1 = x.GetValue(0,0);
    double x2 = x.GetValue(0,1);

    return -20*x1 - 40*x2 + 90*pow(pow(pow(x1,2)+pow(x2-1,2),0.5)-1,2) + 100*pow(pow(pow(x1,2)+pow(x2+1,2),0.5)-1,2);
}

Matrix CostFunctionGradient(Matrix x, double deltaX)
{
    Matrix xvec1(1,2);
    Matrix xvec2(1,2);
    Matrix deriv(1,2);

    for (int i = 0; i < x.Cols(); i++)
    {
        xvec1 = x;
        xvec2 = x;

        xvec1.SetValue(0,i,x.GetValue(0,i) + deltaX);
        xvec2.SetValue(0,i,x.GetValue(0,i) - deltaX);

        deriv.SetValue(0,i,(CostFunction(xvec1)-CostFunction(xvec2))/(2*deltaX));
    }

    return deriv;
}

Matrix CostFunctionHessian(Matrix x, double deltaX)
{
    Matrix temp(x.Rows(), x.Cols());
    double term1, term2, term3, term4;
    Matrix hess(x.Cols(), x.Cols());

    for (int i = 0; i < x.Cols(); i++)
    {
        for (int j = 0; j < x.Cols(); j++)
        {
            if (i==j)
            {
                temp = x;
                temp.SetValue(0,i,x.GetValue(0,i) + deltaX);
                term1 = CostFunction(temp);
                temp.SetValue(0,i,x.GetValue(0,i) - deltaX);
                term2 = CostFunction(temp);
                term3 = CostFunction(x);

                hess.SetValue(i,j, (term1-2*term3+term2)/pow(deltaX, 2));
            }
            else
            {
                temp = x;
                temp.SetValue(0,i,x.GetValue(0,i) + deltaX);
                temp.SetValue(0,j,x.GetValue(0,j) + deltaX);
                term1 = CostFunction(temp);
                
                temp = x;
                temp.SetValue(0,i,x.GetValue(0,i) + deltaX);
                temp.SetValue(0,j,x.GetValue(0,j) - deltaX);
                term2 = CostFunction(temp);

                temp = x;
                temp.SetValue(0,i,x.GetValue(0,i) - deltaX);
                temp.SetValue(0,j,x.GetValue(0,j) + deltaX);
                term3 = CostFunction(temp);

                temp = x;
                temp.SetValue(0,i,x.GetValue(0,i) - deltaX);
                temp.SetValue(0,j,x.GetValue(0,j) - deltaX);
                term4 = CostFunction(temp);

                hess.SetValue(i,j,(term1-term2-term3+term4)/(4*pow(deltaX, 2)));
            }
        }   
    }

    return hess;
}