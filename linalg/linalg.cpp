#include <iostream>
#include <vector>

#include "matrix.h"
#include "linalg.h"

using namespace std;

void ImprM(Matrix x)
{
    for (int i = 0; i < x.Rows(); i++)
    {
        for (int j = 0; j < x.Cols(); j++)
        {
            cout << x.GetValue(i,j) << '\t';
        }
        cout << endl;
    }

    cout << endl;
}

double Norma(Matrix x)
{
    double x1 = x.GetValue(0,0);
    double x2 = x.GetValue(0,1);

    return sqrt(pow(x1, 2) + pow(x2, 2));
}

Matrix InverM(Matrix x)
{
    double a = x.GetValue(0,0);
    double b = x.GetValue(0,1);
    double c = x.GetValue(1,0);
    double d = x.GetValue(1,1);

    Matrix y({{d,-b},{-c,a}});
    double deter = 1/(a*d - b*c);

    for (int i = 0; i < y.Rows(); i++)
    {
        for (int j = 0; j < y.Cols(); j++)
        {
            y.SetValue(i,j, deter*y.GetValue(i,j));
        }
    }
    
    return y;
}

Matrix TransM(Matrix x)
{
    Matrix y(x.Cols(), x.Rows());

    for (int i = 0; i < y.Rows(); i++)
    {
        for (int j = 0; j < y.Cols(); j++)
        {
            y.SetValue(i,j,x.GetValue(j,i));
        }        
    }
    
    return y;
}

Matrix MultEM(double x, Matrix y)
{
    Matrix z(y.Rows(), y.Cols());

    for (int i = 0; i < z.Rows(); i++)
    {
        for (int j = 0; j < z.Cols(); j++)
        {
            z.SetValue(i, j, x*y.GetValue(i,j));
        }
    }
    
    return z;
}

Matrix MultMM(Matrix x, Matrix y)
{
    Matrix z(x.Rows(), y.Cols());

    for (int i = 0; i < x.Rows(); i++)
    {
        for (int j = 0; j < y.Cols(); j++)
        {
            for (int k = 0; k < x.Cols(); k++)
            {
                z.SetValue(i,j,z.GetValue(i,j)+x.GetValue(i,k)*y.GetValue(k,j));
            }   
        }
    }
    
    return z;
}

Matrix DivME(Matrix x, double y)
{
    Matrix z(x.Rows(), x.Cols());

    for (int i = 0; i < z.Rows(); i++)
    {
        for (int j = 0; j < z.Cols(); j++)
        {
            z.SetValue(i,j,x.GetValue(i,j)/y);
        }
        
    }
    
    return z;
}

Matrix SumarMM(Matrix x, Matrix y)
{
    Matrix z(x.Rows(), x.Cols());

    for (int i = 0; i < z.Rows(); i++)
    {
        for (int j = 0; j < z.Cols(); j++)
        {
            z.SetValue(i, j, x.GetValue(i,j) + y.GetValue(i,j));
        }
    }
    
    return z;
}

Matrix RestarMM(Matrix x, Matrix y)
{
    Matrix z(x.Rows(), x.Cols());

    for (int i = 0; i < z.Rows(); i++)
    {
        for (int j = 0; j < z.Cols(); j++)
        {
            z.SetValue(i, j, x.GetValue(i,j) - y.GetValue(i,j));
        }
    }
    
    return z;
}

Matrix Iden(int x)
{
    Matrix y(x, x);

    for (int i = 0; i < y.Rows(); i++)
    {
        for (int j = 0; j < y.Cols(); j++)
        {
            if(i == j)
            {
                y.SetValue(i,j,1);
            }
        }    
    }
    
    return y;
}