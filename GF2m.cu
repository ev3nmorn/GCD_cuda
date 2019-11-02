#include "GF2m.cuh"
#include <cmath>
#include <iostream>
#include <list>
#include <string>
#include <unordered_map>

uint32_t GF2m::fieldMod_ = 0;
uint32_t GF2m::bitsPerElement_ = 0;
uint32_t GF2m::highBitPoly_ = 0;
std::unordered_map<std::string, uint32_t> GF2m::elements_ = *new std::unordered_map<std::string, uint32_t>();
std::list<std::pair<std::string, uint32_t>> GF2m::elementsToPrint_ = *new std::list<std::pair<std::string, uint32_t>>();

uint32_t GF2m::getHighBit(uint32_t value) {
	uint32_t result = 2147483648;

	while (value < result) {
		result >>= 1;
	}

	return result;
}

uint32_t GF2m::elementsMult(uint32_t op1, uint32_t op2) {
	uint32_t result = 0,
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
	
	return result & (GF2m::highBitPoly_ - 1);
}

std::string GF2m::elementToBitStr(uint32_t value) {
	std::string result = "";
	uint32_t mask = highBitPoly_ >> 1;

	for (uint32_t i = 0; i < bitsPerElement_; i++) {
		result += ((mask & value) != 0) ? "1" : "0";
		mask >>= 1;
	}
	
	return result;
}

bool GF2m::isElementByValue(uint32_t value) {
	for (auto it = elements_.begin(); it != elements_.end(); ++it) {
		if (it->second == value) {
			return true;
		}
	}

	return false;
}

void GF2m::build(uint32_t irreduciblePolynomial) {
	uint32_t primitiveElement = 2,
		fieldElement;

	highBitPoly_ = getHighBit(irreduciblePolynomial);
	bitsPerElement_ = (uint32_t)log2(highBitPoly_);
	fieldMod_ = irreduciblePolynomial ^ highBitPoly_;

	while (elements_.size() != highBitPoly_) {
		elements_.clear();
		elementsToPrint_.clear();
		elements_.insert(std::pair<std::string, uint32_t>("0", 0));
		elements_.insert(std::pair<std::string, uint32_t>("1", 1));
		elementsToPrint_.push_back(std::pair<std::string, uint32_t>("0", 0));
		elementsToPrint_.push_back(std::pair<std::string, uint32_t>("1", 1));

		fieldElement = primitiveElement;
		for (uint32_t i = 0; i < highBitPoly_ - 2; i++) {
			if (isElementByValue(fieldElement)) {
				break;
			}
			elements_.insert(std::pair<std::string, uint32_t>("a" + std::to_string(i + 1), fieldElement));
			elementsToPrint_.push_back(std::pair<std::string, uint32_t>("a" + std::to_string(i + 1), fieldElement));
			fieldElement = elementsMult(fieldElement, primitiveElement);
		}

		primitiveElement++;
	}
}

void GF2m::print() {
	if (elements_.size() == 0) {
		std::cout << "empty field" << std::endl;
		return;
	}

	std::string result = "";
	for(auto element : elementsToPrint_) {
		result += element.first + " " + elementToBitStr(element.second) + "\n";
	}
	result.erase(result.length() - 1);

	std::cout << "GF(2^" << bitsPerElement_ << ") is:" << std::endl << result << std::endl;
}