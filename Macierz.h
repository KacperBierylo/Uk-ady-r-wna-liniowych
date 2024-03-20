#pragma once

class Macierz {
public:
	int rozmiar;
	double** tab;
	Macierz(int N);
	Macierz(int N, int a1, int a2, int a3);
	void wyswietl();
	Macierz* kopia();
};
