#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "GF2m.cuh"
#include "Polynomial.cuh"
#include "Term.cuh"
#include <iostream>
#include <regex>
#include <string>
#include <vector>
#include <sstream>

Polynomial GCD(Polynomial p1, Polynomial p2) {
	Polynomial q,
		max = (*p1.getTerms().begin()).getPower() >= (*p2.getTerms().begin()).getPower() ? p1 : p2,
		min = (*p1.getTerms().begin()).getPower() < (*p2.getTerms().begin()).getPower() ? p1 : p2;
	std::vector<Polynomial> row1 = { *new Polynomial(*new std::vector<std::string> { "1" }), *new Polynomial(), max },
		row2 = { *new Polynomial(), *new Polynomial(*new std::vector<std::string>{ "1" }), min },
		temp;

	while ((*row2[2].getTerms().begin()).getCoefficient() != "0") {
		temp = { *new Polynomial(row2[0].getTerms()), *new Polynomial(row2[1].getTerms()), *new Polynomial(row2[2].getTerms()) };

		q = row1[2] / row2[2];

		row2[0] = row1[0] - (row2[0] * q);
		row2[1] = row1[1] - (row2[1] * q);
		row2[2] = row1[2] - (row2[2] * q);
		row1 = { *new Polynomial(temp[0].getTerms()), *new Polynomial(temp[1].getTerms()), *new Polynomial(temp[2].getTerms()) };
	}

	return row1[2];
}

Polynomial multiGCD(std::vector<Polynomial> polynomials) {
	Polynomial result;

	for (int i = 0; i < polynomials.size() - 1; i++) {
		result = (i == 0) ? GCD(polynomials[i], polynomials[i + 1]) : GCD(result, polynomials[i + 1]);
	}

	return result;
}

uint32_t binaryStrToUInt(std::string str) {
	uint32_t result = 0,
		temp = 1;

	for (size_t i = str.length(); i > 0; i--) {
		if (str[i - 1] == '1') {
			result += temp;
		}
		temp *= 2;
	}

	return result;
}

std::vector<std::string> split(const std::string& s, char delimiter) {
	std::vector<std::string> tokens;
	std::string token;
	std::istringstream tokenStream(s);
	while (std::getline(tokenStream, token, delimiter)) {
		tokens.push_back(token);
	}
	return tokens;
}

bool checkTerms(std::vector<std::string> terms) {
	std::regex regex1 = *new std::regex(".+(?=(\\*x(\\^[0-9]+)?))"),
		regex2 = *new std::regex("^x(\\^[0-9]+)?");
	std::smatch match;
	auto elemets = GF2m::getElements();

	for (auto term : terms) {
		if (std::regex_search(term, match, regex1)) {
			if (elemets.find(match.str()) == elemets.end()) {
				return false;
			}
			continue;
		}

		if (std::regex_match(term, regex2))	{
			continue;
		}

		if (elemets.find(term) == elemets.end()) {
			return false;
		}
	}

	return true;
}

int main() {
	std::cout << "Enter irreducible polynomial: ";
	std::string irreduciblePolynomialStr;
	std::cin >> irreduciblePolynomialStr;
	if (!std::regex_match(irreduciblePolynomialStr, *new std::regex("^1[0,1]*1$")) ||
		irreduciblePolynomialStr.size() > 32) {
		std::cout << "Invalid input" << std::endl;
		return 0;
	}

	GF2m::build(binaryStrToUInt(irreduciblePolynomialStr));
	GF2m::print();

	std::cout << "Enter number of polynomials: " << std::endl;
	int n;
	std::cin >> n;
	if (n < 2) {
		std::cout << "At least 2 polynomials" << std::endl;
		return 0;
	}

	std::cout << "Enter " << n << " polynomials:" << std::endl;
	std::string polynomial;
	std::vector<std::string> terms;
	std::vector<Polynomial> polynomials;
	for (int i = 0; i < n; ++i)	{
		std::cin >> polynomial;
		
		terms = split(polynomial, '+');
		if (!checkTerms(terms)) {
			std::cout << "Invalid input" << std::endl;
			return 0;
		}
		polynomials.push_back(*new Polynomial(terms));
	}

	std::cout << "gcd is: " << multiGCD(polynomials) << std::endl;

	return 0;
}