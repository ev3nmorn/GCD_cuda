#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "GF2m.cuh"
#include <iostream>
#include <list>
#include <string>
#include <vector>
#include <algorithm>

class Polynomial {
private:
	int* terms_array_;
	int array_size_;
	GF2m field;

	int get_power_by_value(int coefficient);
	std::string term_to_string(int coefficient, int power);
	int* copy_terms();
	int is_normalized();
	void normalize(int size);
	int dirive_elements(int el1, int el2);

public:
	Polynomial() {};
	Polynomial(GF2m field);
	Polynomial(int* terms, int size, GF2m field);
	int get_array_size() { return array_size_; };
	int* get_terms_array() { return terms_array_; };
	friend Polynomial operator +(Polynomial p1, Polynomial p2);
	friend Polynomial operator -(Polynomial p1, Polynomial p2);
	friend Polynomial operator *(Polynomial p1, Polynomial p2);
	friend Polynomial operator /(Polynomial p1, Polynomial p2);
	friend std::ostream& operator<<(std::ostream &os, Polynomial& poly);
};
