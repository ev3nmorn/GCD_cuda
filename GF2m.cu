#include "GF2m.cuh"
#include <cmath>
#include <iostream>
#include <list>
#include <string>
#include <unordered_map>

uint32_t GF2m::get_high_bit(uint32_t value) {
	uint32_t result = 2147483648;

	while (value < result) {
		result >>= 1;
	}

	return result;
}

uint32_t GF2m::elements_mult(uint32_t op1, uint32_t op2) {
	/*uint32_t result = 0,
		hiBitSet;

	for (uint32_t i = 0; i < GF2m::bitsPerElement_; i++) {
		if ((op2 & 1) != 0) {
			result ^= op1;
		}
		hiBitSet = op1 & (GF2m::highBitPoly_ >> 1);
		op1 <<= 1;
		if (hiBitSet != 0) {
			op1 ^= GF2m::fieldMod_;
		}
		op2 >>= 1;
	}

	return result & (GF2m::highBitPoly_ - 1);*/
	thrust::device_vector<uint32_t> dv(bitsPerElement_);
	for (uint32_t i = 0; i < bitsPerElement_; i++) {
		dv[i] = i;
	}
	shift<uint32_t> f(op1, op2, highBitPoly_, fieldMod_);
	thrust::transform(dv.begin(), dv.end(), dv.begin(), f);
	uint32_t result = thrust::reduce(dv.begin(), dv.end(), (uint32_t)0, thrust::bit_xor<uint32_t>());

	return result;
}

std::string GF2m::element_to_bit_str(uint32_t value) {
	std::string result = "";
	uint32_t mask = highBitPoly_ >> 1;

	for (uint32_t i = 0; i < bitsPerElement_; i++) {
		result += ((mask & value) != 0) ? "1" : "0";
		mask >>= 1;
	}

	return result;
}

bool GF2m::is_element_by_value(uint32_t value) {
	for (auto it = elements_.begin(); it != elements_.end(); ++it) {
		if (it->second == value) {
			return true;
		}
	}

	return false;
}

GF2m::GF2m(uint32_t irreduciblePolynomial) {
	uint32_t primitiveElement = 2,
		fieldElement;

	highBitPoly_ = get_high_bit(irreduciblePolynomial);
	bitsPerElement_ = (uint32_t)log2(highBitPoly_);
	fieldMod_ = irreduciblePolynomial ^ highBitPoly_;

	while (elements_.size() != highBitPoly_) {
		elements_.clear();
		elementsToPrint_.clear();
		elements_.insert(std::pair<int, uint32_t>(-1, 0));
		elements_.insert(std::pair<int, uint32_t>(0, 1));
		elementsToPrint_.push_back(std::pair<int, uint32_t>(-1, 0));
		elementsToPrint_.push_back(std::pair<int, uint32_t>(0, 1));

		fieldElement = primitiveElement;
		for (uint32_t i = 0; i < highBitPoly_ - 2; i++) {
			if (is_element_by_value(fieldElement)) {
				break;
			}
			elements_.insert(std::pair<int, uint32_t>(i + 1, fieldElement));
			elementsToPrint_.push_back(std::pair<int, uint32_t>(i + 1, fieldElement));
			fieldElement = elements_mult(fieldElement, primitiveElement);
		}

		primitiveElement++;
	}
}

std::ostream& operator<<(std::ostream &os, GF2m& field) {
	if (field.elements_.size() == 0) {
		os << "empty field" << std::endl;
		return os;
	}

	std::string result = "";
	for (auto element : field.elementsToPrint_) {
		if (element.first == -1) {
			result += "0 " + field.element_to_bit_str(element.second) + "\n";
		}
		else if (element.first == -1) {
			result += "1 " + field.element_to_bit_str(element.second) + "\n";
		}
		else {
			result += "a" + std::to_string(element.first) + " " + field.element_to_bit_str(element.second) + "\n";
		}
	}
	result.erase(result.length() - 1);

	os << "GF(2^" << field.bitsPerElement_ << ") is:" << std::endl << result << std::endl;

	return os;
}