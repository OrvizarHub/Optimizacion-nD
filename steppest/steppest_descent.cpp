#include <iostream>

#include "..\utilities\utilities.h"
#include "..\linalg\linalg.h"
#include "..\linalg\matrix.h"
#include "..\problem\problem.h"
#include "..\golden\golden.h"
#include "steppest_descent.h"

using namespace std;

void SteppestDescent()
{
    cout << "#######################################################################" << endl;
    cout << "#                                                                     #" << endl;
    cout << "#                    Steppest Descent Algorithm                       #" << endl;
    cout << "#                                                                     #" << endl;
    cout << "#######################################################################" << endl;

    //Constantes ambiente
    ios::fmtflags oldSet = cout.flags();

    //Constantes
    const double epsilon1 = 1e-3;
    const double epsilon2 = 1e-3;
    const double deltaX = 1e-3;

    //Variables iniciales
    Matrix x({{-3, 2}});
    double Ux_prev = CostFunction(x);

    //Variables de proceso
    Matrix sb;
    double alpha;
    double Ux;
    double cont = 0;
    Matrix reg;

    cout << "It.\tx1\t\tx2\t\tf(x)\t\tDer" << endl;
    cout << "-----------------------------------------------------------------------" << endl;

    while (true)
    {
        sb = MultEM(-1, CostFunctionGradient(x, deltaX));
        alpha = Minimizar(x, sb);
        x = SumarMM(x, MultEM(alpha, sb));
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
    }
    
    cout.flags(oldSet);
    GuardarRegistro(reg);
}