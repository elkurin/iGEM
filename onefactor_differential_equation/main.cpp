#include <iostream>
#include <fstream>
#include <ostream>
#include "function.h"
#include "variable.h"
#include <stdio.h>

using namespace std;

namespace {
	ofstream take_log("data_equation.log");
}

double time_bunkai = 0.01;
int time_end = 5000;
double t;
double x;

int main(void)
{
	x = get();
	t = 0;
	cout << t << " " << x << endl;
	while (t < time_end) {
		t++;
		x = equation_runge_kutta(x, time_bunkai);
		cout << t * time_bunkai << " " << x << endl;
		take_log << t * time_bunkai << " " << x << endl;
	}

	FILE *gp;
	gp = popen("gnuplot -persist", "w");
	fprintf(gp, "p  \"data_equation.log\" with lines\n");
	pclose(gp);
	
	return 0;
}
