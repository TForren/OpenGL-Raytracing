#include <math.h>
#include <cstddef>
#include <iostream>
#include "QFHeader.h"
using namespace::std;



QuadFormula::QuadFormula() {
	x = 0;
	x1 = 0;
	x2 = 0;
}

double QuadFormula::filterNegatives(int x) {
	if (x == 0) {
		//nothing. our value is already null.
	}
	else if (x <= 0) {
		x = 0; //Our value is negative meaning it is behind us. set to 0
	}
	else if (x > 0) {
		// value is positive so keep it. 
	}
	return x;
}
/*
quadFormula
@params: a, b, c from a quadratic equation
this method will determine the number of results and the result values from the given quadratic equation
if we find 2 values, resultcount will be 2 and x1, x2 will have values. x will be 0
if we find 1 value, resultcount will be 1 and x will have a value. x1,x2 will be 0
0 values means all x's are 0;
*/
void QuadFormula::quadFormula(double a, double b, double c)
	{
		discriminant = ( (b* b) - 4 * a * c);

		if (discriminant > 0) { //we have 2 solutions
			resultcount = 2;
			x = 0;
			x1 = (-b + sqrt(b * b - 4 * a * c)) / (2 * a);
			x2 = (-b - sqrt(b * b - 4 * a * c)) / (2 * a);
		}

		else if (discriminant == 0) {//we have 1 solution
			resultcount = 1;
			x = ((-b + sqrt(b * b - 4 * a * c)) / (2 * a));
		}

		else if (discriminant < 0) {//we have 0 solutuons
			resultcount = 0;
			x = 0;
			x1 = 0;
			x2 = 0;
		}

	}