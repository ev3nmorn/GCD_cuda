#pragma once

#include <iostream>
#include <list>
#include <string>
#include <unordered_map>


class GF2m {
private:
	static std::unordered_map<std::string, uint32_t> elements_;
	static std::list<std::pair<std::string, uint32_t>> elementsToPrint_;
	static uint32_t fieldMod_, bitsPerElement_, highBitPoly_;

	GF2m() {}
	~GF2m() {}
	static uint32_t getHighBit(uint32_t value);
	static uint32_t elementsMult(uint32_t op1, uint32_t op2);
	static std::string elementToBitStr(uint32_t value);
	static bool isElementByValue(uint32_t value);

public:
	static void build(uint32_t irreduciblePolynomial);
	static uint32_t getElementsCount() { return highBitPoly_; }
	static std::unordered_map<std::string, uint32_t> getElements() { return elements_; }
	static void print();
};