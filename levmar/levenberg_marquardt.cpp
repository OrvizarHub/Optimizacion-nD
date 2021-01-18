#include <iostream>
#include <math.h>

#include "..\utilities\utilities.h"
#include "..\linalg\linalg.h"
#include "..\linalg\matrix.h"
#include "..\problem\problem.h"
#include "levenberg_marquardt.h"

using namespace std;

void LevenbergMarquardt()
{
    cout << "#######################################################################" << endl;
    cout << "#                                                                     #" << endl;
    cout << "#                        Levenberg Marquardt                          #" << endl;
    cout << "#                                                                     #" << endl;
    cout << "#######################################################################" << endl;

    // Constantes de entorno
    ios::fmtflags oldSet = cout.flags();

    // Constantes 
    const double epsilon1 = 1e-3;
    const double epsilon2 = 1e-3;
    const double deltaX = 1e-3;

    // Variables iniciales
    Matrix x({{-3, 2}});
    double Ux_prev = CostFunction(x);

    // Variables de proceso
    double cont = 0;
    Matrix reg;
    Matrix gradient;
    Matrix hessian;
    Matrix sb;
    double Ux;
    double lambda = 1e3;

    cout << "It.\tx1\t\tx2\t\tf(x)\t\tDer" << endl;
    cout << "-----------------------------------------------------------------------" << endl;

    while (true)
    {
        gradient = CostFunctionGradient(x, deltaX);
        hessian = CostFunctionHessian(x, deltaX);
        sb = MultEM(-1, MultMM(InverM(SumarMM(hessian, MultEM(lambda, Iden(x.Cols())))), TransM(gradient)));
        x = SumarMM(x, TransM(sb));
        Ux = CostFunction(x);

        if (Ux < Ux_prev)
        {
            lambda /= 2;
        }
        else
        {
            lambda *=2;
        }

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
    }

    cout.flags(oldSet);
    GuardarRegistro(reg);
}