#include <iostream>
#include <math.h>

#include "..\utilities\utilities.h"
#include "..\linalg\linalg.h"
#include "..\linalg\matrix.h"
#include "..\problem\problem.h"
#include "..\golden\golden.h"
#include "powell.h"

using namespace std;

Matrix DirInit(int n)
{
    Matrix dir(n+1, n);

    for (int i = 0; i < dir.Rows(); i++)
    {
        for (int j = 0; j < dir.Cols(); j++)
        {
            if(i==j)
            {
                dir.SetValue(i,j,1);
            }
            else
            {
                dir.SetValue(i,j,0);
            }
            
        }
    }
    
    return dir;
}

void Powell()
{
    cout << "#######################################################################" << endl;
    cout << "#                                                                     #" << endl;
    cout << "#                                Powell                               #" << endl;
    cout << "#                                                                     #" << endl;
    cout << "#######################################################################" << endl;

    // Constantes de entorno
    ios::fmtflags oldSet = cout.flags();

    // Constantes
    const double epsilon = 1e-4;
    const int nVar = 2;

    // Variables iniciales
    Matrix dir = DirInit(nVar);
    Matrix x({{-3,2}});
    double Ux_prev = CostFunction(x);

    // Variables del proceso
    double cont = 0;
    Matrix reg;
    double alpha;
    double Ux;
    Matrix y;
    Matrix dir_temp; //sb
    int k;

    cout << "It.\tx1\t\tx2\t\tf(x)\t\tDer" << endl;
    cout << "-----------------------------------------------------------------------" << endl;

    while (true)
    {
        y = x;
        k = 0;

        while (k<nVar)
        {
            alpha = Minimizar(x, dir.GetRow(k));
            x = SumarMM(x, MultEM(alpha, dir.GetRow(k)));
            Ux = CostFunction(x);
            k++;
        }

        dir.SetRow(k, RestarMM(x,y));
        alpha = Minimizar(x, dir.GetRow(k));
        x = SumarMM(x, MultEM(alpha, dir.GetRow(k)));
        Ux = CostFunction(x);

        dir_temp = dir;
        for (int j = 0; j < nVar; j++)
        {
            dir.SetRow(j, dir.GetRow(j+1));
        }
        
        cont++;

        reg.AddValue({cont, x.GetValue(0,0), x.GetValue(0,1), Ux});
        MostrarEstado(cont, x, Ux, Matrix(1,nVar), oldSet);

        if(abs(Ux_prev - Ux) < epsilon)
        {
            cout << "-----------------------------------------------------------------------" << endl;
            MostrarEstado(cont, x, Ux, Matrix(1,nVar), oldSet);
            break;
        }

        Ux_prev = Ux;
    }

    cout.flags(oldSet);
    GuardarRegistro(reg);
}