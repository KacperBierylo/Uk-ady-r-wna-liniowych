#include "Macierz.h"
#include  <iostream>
using std::cout;
using std::endl;

Macierz::Macierz(int n) {

	this->rozmiar = n;
	this->tab = new double* [n];
	for (int i = 0; i < rozmiar; i++) {
		tab[i] = new double[n]();
	}

}

Macierz::Macierz(int N, int a1, int a2, int a3) {

	this->rozmiar = N;
	this->tab = new double* [N];
	for (int i = 0; i < rozmiar; i++) {
		tab[i] = new double[N];
	}

	for (int i = 0; i < N; i++) {

		for (int j = 0; j < N; j++) {

			if (i == j) tab[i][j] = a1;
			else if ((j - i == 1) || (i - j) == 1) tab[i][j] = a2;
			else if ((j - i == 2) || (i - j) == 2) tab[i][j] = a3;
			else tab[i][j] = 0;
		}
	}
}

void Macierz::wyswietl() {


	for (int i = 0; i < rozmiar; i++) {
		for (int j = 0; j < rozmiar; j++) {
			if (tab[i][j] >= 0)cout << " " << tab[i][j];
			else cout << tab[i][j];
		}
		cout << endl;
	}

}

Macierz* Macierz::kopia() {
	Macierz* nowaMacierz = new Macierz(rozmiar);
	for (int i = 0; i < rozmiar; i++) {
		for (int j = 0; j < rozmiar; j++) {
			nowaMacierz->tab[i][j] = this->tab[i][j];
		}
	}
	return nowaMacierz;
}