#include <iostream>
#include "function.h"
#include "variable.h"

using namespace std;

double equation_eular(double x, double time_bunkai)
{
	x += time_bunkai * bibun_equation(t, x); 
	return x;
}

double equation_runge_kutta(double x, double time_bunkai)
{
	double k1 = bibun_equation(t, x);
	double k2 = bibun_equation(t + time_bunkai / 2, x + k1 * time_bunkai / 2);
	double k3 = bibun_equation(t + time_bunkai / 2, x + k2 * time_bunkai / 2);
	double k4 = bibun_equation(t + time_bunkai, x + k3 * time_bunkai);
	x += time_bunkai / 6 * (k1 + k2 * 2 + k3 * 2 + k4);
	return x;
}

double equation_runge_kutta_5(double x, double time_bunkai)
{
	double k1 = bibun_equation(t, x);
	double k2 = bibun_equation(t + time_bunkai / 2, x + k1 * time_bunkai / 2);
	double k3 = bibun_equation(t + time_bunkai / 4, x + k2 * time_bunkai / 4);
	double k4 = bibun_equation(t + time_bunkai / 2, x + k3 * time_bunkai / 2);
	double k5 = bibun_equation(t + time_bunkai, x + k4 * time_bunkai);
	x += time_bunkai / 10 * (k1 + k2 * 2 + k3 * 4 + k4 * 2 + k5);
	return x;
}

