#include<fstream>
#include<math.h>
#include<iostream>
using namespace std;
const int N = 100;
const double A = 0.034;

double m2 (double t){
	return 2 - 3*pow(t,3) + 3*pow(t,2) - 2*M_E;
}
double m1 (double t){
	return -3*pow(t,3)-2;
}
double m (double x){
	return 2*pow(x,4) - 2*exp(x);
}
double courant_number(double A, double tao, double h){
	return 2*A*tao/pow(h,2);
}
double F(double t, double x){
	return -9*t*t + 6*t*x - A*(24*x*x - 2*exp(x));
}
double u(double x, double t){
 return 2 * pow(x,4) - 3 * pow(t, 3) +3 * t * t * x - 2 * exp(x);	
}
	/**
	 * n - число уравнений (строк матрицы)
	 * b - диагональ, лежащая над главной (нумеруется: [0;n-2])
	 * c - главная диагональ матрицы A (нумеруется: [0;n-1])
	 * a - диагональ, лежащая под главной (нумеруется: [1;n-1])
	 * f - правая часть (столбец)
	 * x - решение, массив x будет содержать ответ
	 */
void solveMatrix (int n, double *a, double *c, double *b, double *f, double *x)
{
	double m;
	for (int i = 1; i < n; i++)
	{
		m = a[i]/c[i-1];
		c[i] = c[i] - m*b[i-1];
		f[i] = f[i] - m*f[i-1];
	}

	x[n-1] = f[n-1]/c[n-1];

	for (int i = n - 2; i >= 0; i--)
		x[i]=(f[i]-b[i]*x[i+1])/c[i];

}

int main (){
	ofstream approx, func;
	approx.open("approx.txt");
	func.open("func.txt");
	cout.precision(4);
   	double x[N], t[N], tao, h;
	cin >> tao >> h;
	int i = 0;
	while ((i-1)*h < 1){
		x[i] = i*h;
		i++;
	}	
	
	int x_length = i;
	i = 0;
	while ((i-1)*tao < 1){
		t[i] = i*tao;
		i++;
	}
	int tao_length = i;
	
	double U[tao_length][x_length];
	
	for (int i = 0; i < tao_length; i++)
		for (int j = 0; j < x_length; j++)
			U[i][j] = 0;
	
	
	for (i = 0; i < x_length; i++) 
		U[0][i] = m(x[i]);
	for (i = 0; i < tao_length; i++){
		U[i][0] = m1(t[i]);
		U[i][x_length - 1]  =m2(t[i]);
	} 
	
	
	cout << endl << endl;
	double cour = courant_number(A, tao, h);
	double a[N], b[N], c[N],f[N];
	for (int k = 1; k < tao_length; k++){
		for (i = 0; i < x_length-2; i++){
			a[i]=(cour+1); 
			b[i]=-cour/2; 
			c[i]=-cour/2; 
			f[i] = U[k-1][i+1] + tao * F(t[k], x[i+1]);
		}
		f[0] -= b[0] * U[k][0];
		f[x_length - 3] -= c[x_length - 3] * U[k][x_length - 1] ; 
		b[0] = 0; c[x_length - 3] = 0;
		double result[N];
		solveMatrix(x_length - 2, b, a, c,f, result); 
		
		for (int i = 1; i < x_length - 1; i++){
			U[k][i] = result[i - 1];
		}
	}
	
	cout << endl;
	for (i = 0; i < tao_length; i++){
		for (int j = 0; j < x_length; j++)
			approx << U[i][j] << " ";
		approx << endl;
	}
	cout << endl;
	for (int i = 0; i < tao_length; i++){
		for (int j = 0; j < x_length; j++)
			func << u(x[j], t[i]) << " ";
		func << endl;
	}
	return 0;
}

