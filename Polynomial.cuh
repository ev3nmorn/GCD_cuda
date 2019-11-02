#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "Term.cuh"
#include <iostream>
#include <list>
#include <string>
#include <vector>

class Polynomial {
private:
	std::list<Term> terms_;

	void sortByPower();
	void deleteZeroElements();
	bool hasEqualPower(Term term);
	void addToEqualPower(Term term);
	void addUnique(Term term);

public:
	Polynomial();
	Polynomial(std::vector<std::string> terms);
	Polynomial(std::list<Term> terms);
	std::list<Term> getTerms() { return terms_; }
	//static Polynomial GCD(Polynomial p1, Polynomial p2);
	//static Polynomial multiGCD(std::vector<Polynomial> polynomials);
	friend Polynomial operator *(Polynomial p1, Polynomial p2);
	friend Polynomial operator /(Polynomial p1, Polynomial p2);
	friend Polynomial operator +(Polynomial p1, Polynomial p2);
	friend Polynomial operator -(Polynomial p1, Polynomial p2);
	friend std::ostream& operator <<(std::ostream& os, const Polynomial& p);
};
