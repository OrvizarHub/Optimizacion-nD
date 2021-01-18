#include <iostream>

#include "steppest\steppest_descent.h"
#include "newton\newton.h"
#include "levmar\levenberg_marquardt.h"
#include "conjugate\conjugate_gradient.h"
#include "powell\powell.h"

using namespace std;

void main()
{
    // SteppestDescent();
    //  Newton();
    // LevenbergMarquardt();
    // ConjugateGradient();
    Powell();

    system(R"(python ..\utilities\plot_results.py)");
}