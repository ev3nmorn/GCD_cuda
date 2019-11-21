#pragma once

#include "Shift_functor.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "thrust/device_vector.h"
#include <iostream>
#include <list>
#include <string>
#include <unordered_map>


class GF2m {
private:
	std::unordered_map<int, uint32_t> elements_;
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
	uint32_t get_elements_count() { return highBitPoly_; }
	std::unordered_map<int, uint32_t> get_elements() { return elements_; }
	friend std::ostream& operator<<(std::ostream &os, GF2m& field);
};