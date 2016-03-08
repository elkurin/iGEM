#include <iostream>
#include <math.h>
#include "function.h"
#include "variable.h"

using namespace std;

double bibun_equation(double t, double x)
{
	//この中にdx/dt = f(x)のf(x)式にあたるものを書いてください
	double r = 0.1;
	double k = 0.01;
	double f = r * x * (1 + k * x);

	return f;
}

double get(void)
{
	//xの初期値を与えてください
	double x = 0.1;

	return x;
}
