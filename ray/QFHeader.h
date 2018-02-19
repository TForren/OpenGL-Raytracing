class QuadFormula{

	//QuadFormula 1Constructor

public:
	QuadFormula();
	double filterNegatives(int x);
	void quadFormula(double a, double b, double c);

	double x, x1, x2, discriminant;
	int resultcount;
};