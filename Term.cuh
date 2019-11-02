#pragma once

#include "GF2m.cuh"
#include <iostream>
#include <string>

class Term {
private:
	int power_;
	std::string coefficient_;

public:
	Term() {};
	Term(std::string coefficient, int power);
	~Term() { };
	int getPower() { return power_; }
	std::string getCoefficient() { return coefficient_; }
	std::string to_string();
	friend Term operator *(Term t1, Term t2);
	friend Term operator /(Term t1, Term t2);
	friend Term operator +(Term t1, Term t2);
	friend Term operator -(Term t1, Term t2);
};