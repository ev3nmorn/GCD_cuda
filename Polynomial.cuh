#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "Term.cuh"
#include <iostream>
#include <list>
#include <string>
#include <vector>
#include <algorithm>

class Polynomial {
private:
	std::vector<Term> terms_;
	GF2m field;

	void sort_by_power();
	void delete_zero_elements();
	bool has_equal_power(Term term);
	void add_to_equal_power(Term term);
	void add_unique(Term term);

public:
	Polynomial(GF2m field);
	Polynomial(std::vector<Term> terms, GF2m field);
	std::vector<Term> get_terms() { return terms_; }
	friend Polynomial operator *(Polynomial p1, Polynomial p2);
	friend Polynomial operator /(Polynomial p1, Polynomial p2);
	friend Polynomial operator +(Polynomial p1, Polynomial p2);
	friend Polynomial operator -(Polynomial p1, Polynomial p2);
	friend std::ostream& operator<<(std::ostream &os, Polynomial& poly);
};
