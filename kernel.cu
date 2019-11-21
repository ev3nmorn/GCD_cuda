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

std::vector<Term> to_terms(std::vector<std::string> strTerms, GF2m field) {
	std::vector<Term> terms_;
	std::string coefficient = "", powerStr = "";

	for (auto term : strTerms) {
		if (term.find("x^") != std::string::npos) {
			powerStr = term.substr(term.find("x^") + 2, term.size() - term.find("x^") - 2);

			if (term.find('*') != std::string::npos) {
				coefficient = term.substr(0, term.find("*x^"));
			}
			else {
				coefficient = term.substr(0, term.find("x^"));
			}


			if (coefficient == "") {
				terms_.push_back(Term(0, std::stoi(powerStr), field));
			}
			else {
				terms_.push_back(Term(std::stoi(coefficient.substr(1, coefficient.size() - 1)), std::stoi(powerStr), field));
			}
		}
		else if (term.find("x") != std::string::npos) {
			if (term.find('*') != std::string::npos) {
				coefficient = term.substr(0, term.find("*x"));
			}
			else {
				coefficient = term.substr(0, term.find("x"));
			}

			if (coefficient == "") {
				terms_.push_back(Term(0, 1, field));
			}
			else {
				terms_.push_back(Term(std::stoi(coefficient.substr(1, coefficient.size() - 1)), 1, field));
			}
		}
		else {
			if (term == "1") {
				terms_.push_back(Term(0, 0, field));
			}
			else {
				terms_.push_back(Term(std::stoi(term.substr(1, term.size() - 1)), 0, field));
			}
		}
	}

	return terms_;
}

Polynomial gcd(Polynomial p1, Polynomial p2, GF2m field) {
	Polynomial q(field),
		max = (*p1.get_terms().begin()).get_power() >= (*p2.get_terms().begin()).get_power() ? p1 : p2,
		min = (*p1.get_terms().begin()).get_power() < (*p2.get_terms().begin()).get_power() ? p1 : p2;
	std::vector<Polynomial> row1 = { Polynomial(to_terms(std::vector<std::string> { "1" }, field), field), Polynomial(field), max },
		row2 = { Polynomial(field), Polynomial(to_terms(std::vector<std::string>{ "1" }, field), field), min },
		temp;

	while ((*row2[2].get_terms().begin()).get_coefficient() != -1) {
		temp = { Polynomial(row2[0].get_terms(), field), Polynomial(row2[1].get_terms(), field), Polynomial(row2[2].get_terms(), field) };

		q = row1[2] / row2[2];

		row2[0] = row1[0] - (row2[0] * q);
		row2[1] = row1[1] - (row2[1] * q);
		row2[2] = row1[2] - (row2[2] * q);
		row1 = { Polynomial(temp[0].get_terms(), field), Polynomial(temp[1].get_terms(), field), Polynomial(temp[2].get_terms(), field) };
	}

	return row1[2];
}

Polynomial multi_gcd(std::vector<Polynomial> polynomials, GF2m field) {
	Polynomial result(field);

	for (int i = 0; i < polynomials.size() - 1; i++) {
		result = (i == 0) ? gcd(polynomials[i], polynomials[i + 1], field) : gcd(result, polynomials[i + 1], field);
	}

	return result;
}

uint32_t binary_str_to_uint(std::string str) {
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

bool check_terms(std::vector<std::string> terms, GF2m field) {
	std::regex regex1(".+(?=(\\*x(\\^[0-9]+)?))"),
		regex2("^x(\\^[0-9]+)?");
	std::smatch match;
	auto elemets = field.get_elements();

	for (auto term : terms) {
		if (std::regex_search(term, match, regex1)) {
			if (elemets.find(std::stoi(match.str().substr(1, match.str().size() - 1))) == elemets.end()) {
				return false;
			}
			continue;
		}

		if (std::regex_match(term, regex2)) {
			continue;
		}

		if (term == "1") {
			continue;
		}

		if (elemets.find(std::stoi(term.substr(1, term.size() - 1))) == elemets.end()) {
			return false;
		}
	}

	return true;
}

int main() {
	std::cout << "Enter irreducible polynomial: ";
	std::string irreduciblePolynomialStr;
	std::cin >> irreduciblePolynomialStr;
	if (!std::regex_match(irreduciblePolynomialStr, std::regex("^1[0,1]*1$")) ||
		irreduciblePolynomialStr.size() > 32) {
		std::cout << "Invalid input" << std::endl;
		return 0;
	}

	GF2m field(binary_str_to_uint(irreduciblePolynomialStr));
	std::cout << field;

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
	for (int i = 0; i < n; ++i) {
		std::cin >> polynomial;

		terms = split(polynomial, '+');
		if (!check_terms(terms, field)) {
			std::cout << "Invalid input" << std::endl;
			return 0;
		}
		polynomials.push_back(Polynomial(to_terms(terms, field), field));
	}

	std::cout << "gcd is: " << multi_gcd(polynomials, field) << std::endl;

	return 0;
}