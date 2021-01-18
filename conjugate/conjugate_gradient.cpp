#include <iostream>
#include <math.h>

#include "..\utilities\utilities.h"
#include "..\linalg\linalg.h"
#include "..\linalg\matrix.h"
#include "..\problem\problem.h"
#include "..\golden\golden.h"
#include "conjugate_gradient.h"

using namespace std;

void ConjugateGradient()
{
    cout << "#######################################################################" << endl;
    cout << "#                                                                     #" << endl;
    cout << "#                          Conjugate Gradient                         #" << endl;
    cout << "#                                                                     #" << endl;
    cout << "#######################################################################" << endl;

     // Constantes de entorno
    ios::fmtflags oldSet = cout.flags();

    // Constantes 
    const double epsilon1 = 1e-4;
    const double epsilon2 = 1e-4;
    const double deltaX = 1e-3;

    // Variables iniciales
    Matrix x({{-3, 2}});
    double Ux_prev = CostFunction(x);

    // Variables de proceso
    double cont = 0;
    Matrix reg;
    Matrix gradient;
    Matrix gradient_prev;
    Matrix sb;
    Matrix sb_prev;
    Matrix x_prev;
    double Ux;
    double alpha;

    cout << "It.\tx1\t\tx2\t\tf(x)\t\tDer" << endl;
    cout << "-----------------------------------------------------------------------" << endl;

    while (true)
    {
        if(cont==0)
        {
            gradient_prev = CostFunctionGradient(x, deltaX);
            sb_prev = MultEM(-1, gradient_prev);
            alpha = Minimizar(x, sb_prev);

            x_prev = SumarMM(x, MultEM(alpha, sb_prev));
            Ux_prev = CostFunction(x_prev);

            cont++;
        }
        else
        {
            gradient = CostFunctionGradient(x_prev, deltaX);
            sb = SumarMM(MultEM(-1, gradient), MultEM(pow(Norma(gradient), 2)/pow(Norma(gradient_prev), 2), sb_prev));
            alpha = Minimizar(x_prev, sb);

            x = SumarMM(x_prev, MultEM(alpha, sb));
            Ux = CostFunction(x);

            cont++;

            reg.AddValue({cont, x.GetValue(0,0), x.GetValue(0,1), Ux});
            MostrarEstado(cont, x, Ux, sb, oldSet);

            if(abs(Ux_prev - Ux) < epsilon1 || Norma(sb) < epsilon2)
            {
                cout << "-----------------------------------------------------------------------" << endl;
                MostrarEstado(cont, x, Ux, sb, oldSet);
                break;
            }

            Ux_prev = Ux;
            gradient_prev = gradient;
            sb_prev = sb;
            x_prev = x;
        }
    }

    cout.flags(oldSet);
    GuardarRegistro(reg);
}