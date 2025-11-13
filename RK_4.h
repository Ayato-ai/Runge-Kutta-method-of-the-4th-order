#pragma once
#include <iostream>
#include <vector>


std::size_t numbers_steps, a, b;
const int p = 4;
double e;
//std::cin << e;

double true_trajectory(double x, double u0) {
	return u0 * exp(-2.5 * x);
}
double test_function(double x, double v) {
	return -2.5 * v;
}
double function_1(double x, double v) {
	return (std::log(x + 1) / (pow(x, 2) + 1)) * v + v - pow(v, 3) * sin(10 * x);
}

std::pair<double, double> Runge_Kytta_4(double(*f)(double, double), double h_n, double x_n, double v_n) {

	double k1 = f(x_n, v_n);
	double k2 = f(x_n + h_n / 2.0, v_n + h_n / 2.0 * k1);
	double k3 = f(x_n + h_n / 2.0, v_n + h_n / 2.0 * k2);
	double k4 = f(x_n + h_n, v_n + h_n * k3);

	x_n = x_n + h_n;
	v_n = v_n + h_n * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
	return { x_n, v_n };
}
double S(double v_n, double v) {
	return abs((v_n - v)) / (pow(2, p) - 1);
}
//void next_point(double v_n, double v, double* h_n, double e) {
//	double S_ = S(v_n, v);
//	if (S > e) {  // Результаты грубые
//		*h_n = *h_n / 2.0;
//	}
//	else if (S < e / pow(2, p + 1)) {
//		*h_n = *h_n * 2.0;
//	}
//	else {
//
//	}
//}

std::pair<double, double> RK_4_OLP(double(*f)(double, double), double x0, double u0, double h, double e) {
	std::pair<double, double> new_point_h;
	std::pair<double, double> new_point_2h;
	double x_n = x0;
	double v_n = u0;
	bool flag = true;
	double S_new = 0;
	
	while (flag) {
		new_point_h = Runge_Kytta_4(f, h, x0, u0);
		new_point_2h = Runge_Kytta_4(test_function, h / 2.0, x0, Runge_Kytta_4(test_function, h / 2.0, x0, u0).second);
		S_new = S(new_point_h.second, new_point_2h.second);

		if (S_new > e)
			h /= 2;
		else if (S_new < e / pow(2, p + 1)) {
			flag = false;
			h *= 2;
		}
		else {
			flag = false;
		}
	}
	return new_point_h;
}