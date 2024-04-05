#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
#include <string>
#include <iostream>

#define DATA double

void euler(DATA y_at_beg, DATA beg, DATA end, int steps, DATA yprim(DATA x, DATA y));
void rk2(DATA y_at_beg, DATA beg, DATA end, int steps, DATA yprim(DATA x, DATA y));
void rk4(DATA y_at_beg, DATA beg, DATA end, int steps, DATA yprim(DATA x, DATA y));
void rk4_2(DATA y_at_beg, DATA beg, DATA end, int steps, DATA yprim(DATA x, DATA y, DATA y_t_tcr));
void euler_2(DATA y_at_beg, DATA beg, DATA end, int steps, DATA yprim(DATA x, DATA y, DATA y_t_tcr));

DATA f(DATA x, DATA y) {
	return pow(x, 2) + y;
}

DATA f2(DATA x, DATA y) {
	return x + y;
}

DATA a1 = 3.5;
DATA a2 = 5;
DATA a3 = 1;
DATA e_dot = 1;
DATA a8 = 0;
DATA p_cr = 0.4;
DATA p0 = 0;


DATA fun_ivm(DATA x, DATA y) {
	DATA t_cr = 1 / a2 * log((p0 - (a1 / a2)) / (p_cr - (a1 / a2)));
	if (x >= t_cr)
		return a1 * e_dot - a2 * y * e_dot - a3 * pow(y, a8) * 1;
	else
		return a1 * e_dot - a2 * y * e_dot;
}

DATA fun_ivm(DATA x, DATA y, DATA y_t_tcr) {
	DATA t_cr = 1 / a2 * log((p0 - (a1 / a2)) / (p_cr - (a1 / a2)));
	if (x >= t_cr)
		return a1 * e_dot - a2 * y * e_dot - a3 * pow(y, a8) * y_t_tcr;
	else
		return a1 * e_dot - a2 * y * e_dot;
}

int main() {
	//for (int i = 1; i <= 50; i++) {
		std::cout << "\nfunkcja pierwsza";
		euler_2(0, 0, 1, 100, fun_ivm);
		// rk2(0.1, 0, 1, i, f);
		// rk4(0.1, 0, 1, i, f);
		// std::cout << "\n\nfunkcja druga";
		// euler(0.1, 0, 1, i, f2);
		// rk2(0.1, 0, 1, i, f2);
		// rk4(0.1, 0, 1, i, f2);
	//}
	//rk4(0.1, 0, 1, 10, f);
		//rk4(0.1, 0, 1, i, f2);
}
DATA calc_val(DATA t, DATA tcr, int steps, DATA beg, DATA end, std::vector<DATA> y) {
	DATA td = t - tcr;
	if (td <= beg)return y[0];
	td -= beg;
	DATA h = (end - beg) / steps;
	int ind = int(floor(td / h));
	DATA rest = td - ind * h;
	DATA proc = rest / h;
	std::cout << "V:" << y[ind] << " " << y[ind + 1] << " " << y[ind] * proc + (1 - proc) * y[ind + 1];
	return y[ind] * proc + (1 - proc) * y[ind + 1];
}

void euler(DATA y_at_beg, DATA beg, DATA end, int steps, DATA yprim(DATA x, DATA y)) {

	DATA h = (end - beg) / steps;

	std::cout << "\n\nPrzedzial [" << beg << "," << end << "] warunek pocz y0=" << y_at_beg << " krok (h) " << h;
	std::cout << std::endl;
	std::vector<DATA> ys;
	ys.resize(steps + 1);
	ys[0] = y_at_beg;

	for (int i = 1; i <= steps; i++) {
		ys[i] = ys[i - 1] + h * yprim(beg + (i - 1) * h, ys[i - 1]);
		std::cout << ys[i] << " " << std::endl;
	}


	std::cout << "\nEuler Wartosc funkcji w punkcie " << ys[ys.size() - 1];

	std::ofstream plik("data.txt");
	for (int i = 0; i < ys.size(); i++)
		plik << ys[i] << std::endl;

	plik.close();

}

void euler_2(DATA y_at_beg, DATA beg, DATA end, int steps, DATA yprim(DATA x, DATA y, DATA y_t_tcr)) {

	DATA h = (end - beg) / steps;

	std::cout << "\n\nPrzedzial [" << beg << "," << end << "] warunek pocz y0=" << y_at_beg << " krok (h) " << h;
	std::cout << std::endl;
	std::vector<DATA> ys;
	ys.resize(steps + 1);
	ys[0] = y_at_beg;

	DATA t_cr = 1 / a2 * log((p0 - (a1 / a2)) / (p_cr - (a1 / a2)));
	std::cout << "\n\nt_cr [" << t_cr << "] \n";

	for (int i = 1; i <= steps; i++) {
		ys[i] = ys[i - 1] + h * yprim(beg + (i - 1) * h, ys[i - 1], calc_val(beg + (i - 1) * h, t_cr, steps, beg, end, ys));
		std::cout << ys[i] << " " << std::endl;
	}


	std::cout << "\nEuler Wartosc funkcji w punkcie " << ys[ys.size() - 1];

	std::ofstream plik("data.txt");
	for (int i = 0; i < ys.size(); i++)
		plik << ys[i] << std::endl;

	plik.close();

}

void rk2(DATA y_at_beg, DATA beg, DATA end, int steps, DATA yprim(DATA x, DATA y)) {

	DATA h = (end - beg) / steps;


	std::cout << "\n\nPrzedzial [" << beg << "," << end << "] warunek pocz y0=" << y_at_beg << " krok (h) " << h;
	std::vector<DATA> ys;
	ys.resize(steps + 1);
	ys[0] = y_at_beg;

	DATA k, k2;

	for (int i = 1; i <= steps; i++) {

		k = yprim(beg + (i - 1) * h, ys[i - 1]);
		k2 = yprim(beg + i * h, ys[i - 1] + h * k);

		ys[i] = ys[i - 1] + h * 0.5 * (k + k2);

		//std::cout << ys[i] << " ";
	}
	std::cout << "\nRK2 Wartosc funkcji w punkcie " << ys[ys.size() - 1];
}

void rk4(DATA y_at_beg, DATA beg, DATA end, int steps, DATA yprim(DATA x, DATA y)) {

	DATA h = (end - beg) / steps;

	std::cout << "\n\nPrzedzial [" << beg << "," << end << "] warunek pocz y0=" << y_at_beg << " krok (h) " << h << std::endl;
	std::vector<DATA> ys;
	ys.resize(steps + 1);
	ys[0] = y_at_beg;

	DATA k, k2, k3, k4;

	for (int i = 1; i <= steps; i++) {

		k = yprim(beg + (i - 1) * h, ys[i - 1]);
		k2 = yprim(beg + (i - 0.5) * h, ys[i - 1] + h * 0.5 * k);
		k3 = yprim(beg + (i - 0.5) * h, ys[i - 1] + h * 0.5 * k2);
		k4 = yprim(beg + i * h, ys[i - 1] + h * k3);

		ys[i] = ys[i - 1] + h * (k + 2 * k2 + 2 * k3 + k4) / 6;
		std::cout << ys[i] << " ";
	}


	std::cout << "\nRK4 Wartosc funkcji w punkcie " << ys[ys.size() - 1];
}

void rk4_2(DATA y_at_beg, DATA beg, DATA end, int steps, DATA yprim(DATA x, DATA y, DATA y_t_tcr)) {

	DATA h = (end - beg) / steps;

	std::cout << "\n\nPrzedzial [" << beg << "," << end << "] warunek pocz y0=" << y_at_beg << " krok (h) " << h << std::endl;
	std::vector<DATA> ys;
	ys.resize(steps + 1);
	ys[0] = y_at_beg;

	DATA k, k2, k3, k4;
	DATA t_cr = 1 / a2 * log((p0 - (a1 / a2)) / (p_cr - (a1 / a2)));

	for (int i = 1; i <= steps; i++) {
		 
		k = yprim(beg + (i - 1) * h, ys[i - 1], calc_val(beg + (i - 1) * h, t_cr, steps, beg, end, ys));
		k2 = yprim(beg + (i - 0.5) * h, ys[i - 1] + h * 0.5 * k, calc_val(beg + (i - 0.5) * h, t_cr, steps, beg, end, ys));
		k3 = yprim(beg + (i - 0.5) * h, ys[i - 1] + h * 0.5 * k2, calc_val(beg + (i - 0.5) * h, t_cr, steps, beg, end, ys));
		k4 = yprim(beg + i * h, ys[i - 1] + h * k3, calc_val(beg + i * h, t_cr, steps, beg, end, ys));

		ys[i] = ys[i - 1] + h * (k + 2 * k2 + 2 * k3 + k4) / 6;
		std::cout << ys[i] << " ";
	}


	std::cout << "\nRK4 Wartosc funkcji w punkcie " << ys[ys.size() - 1];
}