#include <iostream>
#include <chrono>
#include "Macierz.h"
using std::cin;
using std::cout;
using std::endl;
#define N 985 
#define e 1
#define f 0

double* pomnozMacierzPrzezWektor(Macierz* M, double* V);
double* residuum(Macierz* A, double* b, double* x);
double norma(double* V, int rozmiar);
int Jacobi(Macierz* M, double* x, double* b, double oczekiwanaNorma, double* czas);
int GaussSeidel(Macierz* M, double* x, double* b, double oczekiwanaNorma, double* czas);
double LU(Macierz* A, double* x, double* b, double* czas);

int main()
{	
	//A
	int a1 = 5 + e, a2 = -1, a3 = a2;
	Macierz* m1 = new Macierz(N, a1, a2, a3);
	double b[N];
	double x[N];
	for (int i = 0; i < N; i++) {
		b[i] = sin(i * (f + 1));
	}
	//B
	cout << "B: " << endl;
	double czasJB = 0;
	double czasGSB = 0;
	cout << "Jacobi potrzebowal " << Jacobi(m1, x, b, 1e-9, &czasJB) << " iteracji. ";
	cout<< "Czas wykonania wyniosl " << czasJB << " sekund." << endl;
	cout << "Gauss-Seidel potrzebowal "<<GaussSeidel(m1, x, b, 1e-9, &czasGSB)<<" iteracji. ";
	cout<< "Czas wykonania wyniosl " << czasGSB << " sekund." << endl;
	cout << "Stosunek czasu wykonywania metody Jacobiego do czasu wykonywania metody Gaussa-Seidla: " << czasJB / czasGSB << endl;
	
	//C
	cout << endl << "C: " << endl;
	double czasJC = 0;
	double czasGSC = 0;
	int a1c = 3, a2c = -1, a3c = a2c;
	Macierz* mc = new Macierz(N, a1c, a2c, a3c);
	for (int i = 0; i < N; i++) {
		b[i] = sin(i * (f + 1));
	}
	cout << "Jacobi, iteracje: " << Jacobi(mc, x, b, 1e-9, &czasJC) << endl;
	cout << "Gauss-Seidel, iteracje: " << GaussSeidel(mc, x, b, 1e-9, &czasGSC) << endl;
	

	//D
	cout << endl << "D: " << endl;
	double czasLU = 0;
	int a1d = 3, a2d = -1, a3d = a2d;
	Macierz* md = new Macierz(N, a1d, a2d, a3d);
	for (int i = 0; i < N; i++) {
		b[i] = sin(i * (f + 1));
	}
	cout << "Metoda faktoryzacji LU, norma z residuum: " << LU(md, x, b, &czasLU) << ". ";
	cout << "Czas wykonania wyniosl " << czasLU << " sekund." << endl;

	//E
	cout << endl << "E: " << endl;
	int rozmiary[5] = { 100, 500, 1000, 2000, 3000 };
	double czasE = 0;
	double czasyTrwaniaJacobi[5];
	double czasyTrwaniaGaussSeidel[5];
	double czasyTrwaniaLU[5];

	for (int i = 0; i < 5; i++) {

		double* bE= new double[rozmiary[i]];
		double* xE = new double[rozmiary[i]];
		for (int j = 0; j < rozmiary[i]; j++) {
			bE[j] = sin(j * (f + 1));
		}
		Macierz* mE = new Macierz(rozmiary[i], a1, a2, a3);

		LU(mE, xE, bE, &czasE);
		czasyTrwaniaLU[i] = czasE;

		Jacobi(mE, xE, bE, 1e-9, &czasE);
		czasyTrwaniaJacobi[i] = czasE;

		GaussSeidel(mE, xE, bE, 1e-9, &czasE);
		czasyTrwaniaGaussSeidel[i] = czasE;

		delete[] bE;
		delete[] xE;
		delete mE;
	}
	cout << "Czasy Jacobi: ";
	for (int i = 0; i < 5; i++) {
		cout << czasyTrwaniaJacobi[i]<<"  ";
	}
	cout << endl;

	cout << "Czasy Gauss-Seidel: ";
	for (int i = 0; i < 5; i++) {
		 cout<< czasyTrwaniaGaussSeidel[i]<<"  ";
	}
	cout << endl;

	cout << "Czasy LU: ";
	for (int i = 0; i < 5; i++) {
		 cout << czasyTrwaniaLU[i]<<"  ";
	}
	cout << endl;
}

double* pomnozMacierzPrzezWektor(Macierz* M, double* V) {
	double* wynik = new double[M->rozmiar]();
	for (int i = 0; i < M->rozmiar; i++) {

		for (int j = 0; j < M->rozmiar; j++) {

			wynik[i] += M->tab[i][j] * V[j];

		}
	}
	return wynik;
}

double* residuum(Macierz* A, double* b, double* x) {
	double* r = new double[A->rozmiar];
	r = pomnozMacierzPrzezWektor(A, x);
	for (int i = 0; i < A->rozmiar; i++) {
		r[i] -= b[i];
	}
	return r;
}


double norma(double* V, int rozmiar) {
	double wynik = 0;
	for (int i = 0; i < rozmiar; i++) {
		wynik += pow(V[i], 2);
	}
	return sqrt(wynik);
}

int Jacobi(Macierz* M, double* x, double* b, double oczekiwanaNorma, double* czas) {
	auto start = std::chrono::high_resolution_clock::now();
	int iteracje = 0;
	double* nastepnyX = new double[M->rozmiar];
	for (int i = 0; i < M->rozmiar; i++) {
		x[i] = 1;
	}

	do {
		for (int i = 0; i < M->rozmiar; i++) {
			double sumaIloczynow = 0;
			for (int j = 0; j < M->rozmiar; j++) {
				if (i != j) {
					sumaIloczynow += M->tab[i][j] * x[j];
				}
			}
			nastepnyX[i] = (b[i] - sumaIloczynow) / M->tab[i][i];
		}
		for (int i = 0; i < M->rozmiar; i++) {
			x[i] = nastepnyX[i];
		}

		iteracje++;
		//cout << norma(residuum(M, b, x), M->rozmiar)<<endl;
	} while (norma(residuum(M, b, x), M->rozmiar) > oczekiwanaNorma);
	if (isnan(norma(residuum(M, b, x), M->rozmiar))) cout << "Jacobi nie zbiega sie!"<< endl;
	double xd = norma(residuum(M, b, x), M->rozmiar);
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	*czas = elapsed.count();
	return iteracje;
}


int GaussSeidel(Macierz* M, double* x, double* b, double oczekiwanaNorma, double* czas) {
	auto start = std::chrono::high_resolution_clock::now();
	int iteracje = 0;
	double* nastepnyX = new double[M->rozmiar];
	for (int i = 0; i < M->rozmiar; i++) {
		x[i] = 1;
	}
	do {
		for (int i = 0; i < M->rozmiar; i++) {
			double sumaIloczynow = 0;

			for (int j = 0; j < i; j++) {
				sumaIloczynow += M->tab[i][j] * nastepnyX[j];
			}

			for (int j = i + 1; j < M->rozmiar; j++) {
				sumaIloczynow += M->tab[i][j] * x[j];
			}

			nastepnyX[i] = (b[i] - sumaIloczynow) / M->tab[i][i];
		}
		for (int i = 0; i < M->rozmiar; i++) {
			x[i] = nastepnyX[i];
		}

		iteracje++;
		//cout << norma(residuum(M, b, x), M->rozmiar) << endl;
	} while (norma(residuum(M, b, x), M->rozmiar) > oczekiwanaNorma);
	if (isnan(norma(residuum(M, b, x), M->rozmiar))) cout << "Gauss-Seidel nie zbiega sie!" << endl;
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	*czas = elapsed.count();
	return iteracje;
}

double LU(Macierz* A, double* x, double* b, double* czas) {
	auto start = std::chrono::high_resolution_clock::now();
	Macierz* U = A->kopia();
	Macierz* L = new Macierz(A->rozmiar, 1, 0, 0);
	for (int i = 0; i < A->rozmiar - 1; i++) {

		for (int j = i + 1; j < A->rozmiar; j++) {
			L->tab[j][i] = U->tab[j][i] / U->tab[i][i];

			for (int k = i; k < A->rozmiar; k++) {
				U->tab[j][k] -= L->tab[j][i] * U->tab[i][k];
			}
		}
	}

	double* y = new double[A->rozmiar];

	for (int i = 0; i < A->rozmiar; i++) {
		double suma = 0;
		for (int j = 0; j < i; j++) {
			suma += L->tab[i][j] * y[j];
		}
		y[i] = (b[i] - suma) / L->tab[i][i];
	}

	for (int i = A->rozmiar - 1; i >= 0; i--) {
		double suma = 0;
		for (int j = i + 1; j < A->rozmiar; j++) {
			suma += U->tab[i][j] * x[j];
		}
		x[i] = (y[i] - suma) / U->tab[i][i];
	}
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	*czas = elapsed.count();
	return norma(residuum(A, b, x), A->rozmiar);
}