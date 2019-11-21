#pragma once

#include "GF2m.cuh"
#include <iostream>
#include <string>

class Term {
private:
	int coefficient_, power_;
	GF2m field;

public:
	Term(GF2m field) { this->field = field; };
	Term(int coefficient, int power, GF2m field);
	~Term() { };
	int get_power() { return power_; }
	int get_coefficient() { return coefficient_; }
	std::string to_string();
	friend Term operator *(Term t1, Term t2);
	friend Term operator /(Term t1, Term t2);
	friend Term operator +(Term t1, Term t2);
	friend Term operator -(Term t1, Term t2);
};