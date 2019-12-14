#pragma once

#include "Shift_functor.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "thrust/device_vector.h"
#include <iostream>
#include <list>
#include <string>

class GF2m {
private:
	uint32_t* elements_array_;
	int array_size_;
	std::list<std::pair<int, uint32_t>> elementsToPrint_;
	uint32_t fieldMod_, bitsPerElement_, highBitPoly_;

	uint32_t get_high_bit(uint32_t value);
	uint32_t elements_mult(uint32_t op1, uint32_t op2);
	std::string element_to_bit_str(uint32_t value);
	bool is_element_by_value(uint32_t value);

public:
	GF2m(uint32_t irreduciblePolynomial);
	GF2m() {}
	~GF2m() {}
	int get_array_size() { return array_size_; }
	uint32_t* get_field_elements() { return elements_array_; }
	friend std::ostream& operator<<(std::ostream &os, GF2m& field);
};