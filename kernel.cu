#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "GF2m.cuh"
#include "Polynomial.cuh"
#include <iostream>
#include <regex>
#include <string>
#include <vector>
#include <sstream>

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

int* to_terms(std::vector<std::string> strTerms, int& size, GF2m field) {
	int* terms;
	std::string coefficient = "", powerStr = "";
	std::vector<std::pair<int, int>> terms_vector;

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
				terms_vector.push_back(std::make_pair(1, std::stoi(powerStr)));
			}
			else {
				terms_vector.push_back(std::make_pair(field.get_field_elements()[std::stoi(coefficient.substr(1, coefficient.size() - 1))], std::stoi(powerStr)));
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
				terms_vector.push_back(std::make_pair(1, 1));
			}
			else {
				terms_vector.push_back(std::make_pair(field.get_field_elements()[std::stoi(coefficient.substr(1, coefficient.size() - 1))], 1));
			}
		}
		else {
			if (term == "1") {
				terms_vector.push_back(std::make_pair(1, 0));
			}
			else if (term == "0") {
				terms_vector.push_back(std::make_pair(0, 0));
			}
			else {
				terms_vector.push_back(std::make_pair(field.get_field_elements()[std::stoi(term.substr(1, term.size() - 1))], 0));
			}
		}
	}

	std::sort(terms_vector.begin(), terms_vector.end(), [](std::pair<int, int> t1, std::pair<int, int> t2) { return t1.second > t2.second;  });
	
	size = terms_vector[0].second + 1;
	terms = new int[size];
	for (int i = 0; i < size; i++) {
		terms[i] = 0;
	}
	for (int i = 0; i < terms_vector.size(); i++) {
		terms[terms_vector[i].second] = terms_vector[i].first;
	}

	return terms;
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

Polynomial gcd(Polynomial p1, Polynomial p2, GF2m field) {
	int size = 0;
	Polynomial q,
		max = p1.get_array_size() >= p2.get_array_size() ? p1 : p2,
		min = p1.get_array_size() < p2.get_array_size() ? p1 : p2;

	std::vector<Polynomial> row1 = { Polynomial(to_terms(std::vector<std::string> { "1" }, size, field), size, field), Polynomial(field), max },
		row2 = { Polynomial(field), Polynomial(to_terms(std::vector<std::string>{ "1" }, size, field), size, field), min },
		temp;

	while (row2[2].get_terms_array()[row2[2].get_array_size() - 1] != 0) {
		temp = { Polynomial(row2[0].get_terms_array(), row2[0].get_array_size(), field),
			Polynomial(row2[1].get_terms_array(), row2[1].get_array_size(), field),
			Polynomial(row2[2].get_terms_array(), row2[2].get_array_size(), field) };

		q = row1[2] / row2[2];

		row2[0] = row1[0] - (row2[0] * q);
		row2[1] = row1[1] - (row2[1] * q);
		row2[2] = row1[2] - (row2[2] * q);
		row1 = { Polynomial(temp[0].get_terms_array(), temp[0].get_array_size(), field),
			Polynomial(temp[1].get_terms_array(), temp[1].get_array_size(), field),
			Polynomial(temp[2].get_terms_array(), temp[2].get_array_size(), field) };
	}

	return row1[2];
}

Polynomial multi_gcd(std::vector<Polynomial> polynomials, GF2m field) {
	Polynomial result;

	for (int i = 0; i < polynomials.size() - 1; i++) {
		result = (i == 0) ? gcd(polynomials[i], polynomials[i + 1], field) : gcd(result, polynomials[i + 1], field);
	}

	return result;
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

	int size = 0;
	std::cout << "Enter " << n << " polynomials:" << std::endl;
	std::string polynomial;
	std::vector<std::string> terms;
	std::vector<Polynomial> polynomials;
	for (int i = 0; i < n; ++i) {
		std::cin >> polynomial;

		terms = split(polynomial, '+');
		polynomials.push_back(Polynomial(to_terms(terms, size, field), size, field));
	}

	std::cout << "gcd is: " << multi_gcd(polynomials, field) << std::endl;

	return 0;
}