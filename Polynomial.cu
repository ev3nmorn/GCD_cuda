#include "Polynomial.cuh"

void Polynomial::sortByPower() {
	terms_.sort([](Term& t1, Term& t2) { return t1.getPower() > t2.getPower();  });
}

void Polynomial::deleteZeroElements() {
	terms_.remove_if([](Term t) { return t.getCoefficient() == "0"; });

	if (terms_.size() == 0) {
		terms_.push_back(*new Term("0", 0));
	}
}

bool Polynomial::hasEqualPower(Term term) {
	for (Term t : terms_) {
		if (t.getPower() == term.getPower()) {
			return true;
		}
	}

	return false;
}

void Polynomial::addToEqualPower(Term term) {
	for (auto it = terms_.begin(); it != terms_.end(); it++) {
		if ((*it).getPower() == term.getPower()) {
			(*it) = (*it) + term;
		}
	}
}

void Polynomial::addUnique(Term term) {
	if (hasEqualPower(term)) {
		addToEqualPower(term);
	}
	else {
		terms_.push_back(term);
	}
}

Polynomial::Polynomial() {
	terms_ = *new std::list<Term>();
	terms_.push_back(*new Term("0", 0));
}

Polynomial::Polynomial(std::list<Term> terms) {
	terms_ = *new std::list<Term>();
	terms_.assign(terms.begin(), terms.end());
	deleteZeroElements();
	sortByPower();
}

Polynomial::Polynomial(std::vector<std::string> terms) {
	terms_ = *new std::list<Term>();
	std::string coefficient = "", powerStr = "";

	for (auto term : terms) {
		if (term.find("x^") != std::string::npos) {
			powerStr = term.substr(term.find("x^") + 2, term.size() - term.find("x^") - 2);

			if (term.find('*') != std::string::npos) {
				coefficient = term.substr(0, term.find("*x^"));
			}
			else {
				coefficient = term.substr(0, term.find("x^"));
			}


			if (coefficient == "") {
				terms_.push_back(*new Term("1", std::stoi(powerStr)));
			}
			else {
				terms_.push_back(*new Term(coefficient, std::stoi(powerStr)));
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
				terms_.push_back(*new Term("1", 1));
			}
			else {
				terms_.push_back(*new Term(coefficient, 1));
			}
		}
		else {
			terms_.push_back(*new Term(term, 0));
		}
	}

	deleteZeroElements();
	sortByPower();
}

/*__global__ void Polynomial::parallelColumns(std::vector<Polynomial> row1, std::vector<Polynomial> row2, Polynomial q) {
	row2[threadIdx.x] = row1[threadIdx.x] - (row2[threadIdx.x] * q);
}

Polynomial Polynomial::GCD(Polynomial p1, Polynomial p2) {
	Polynomial q,
		max = (*p1.terms_.begin()).getPower() >= (*p2.terms_.begin()).getPower() ? p1 : p2,
		min = (*p1.terms_.begin()).getPower() < (*p2.terms_.begin()).getPower() ? p1 : p2;
	std::vector<Polynomial> row1 = { *new Polynomial(*new std::vector<std::string> { "1" }), *new Polynomial(), max },
		row2 = { *new Polynomial(), *new Polynomial(*new std::vector<std::string>{ "1" }), min },
		temp;
	std::vector<Polynomial> *dev_row1, *dev_row2;
	Polynomial *dev_q;

	while ((*row2[2].terms_.begin()).getCoefficient() != "0") {
		temp = { *new Polynomial(row2[0].terms_), *new Polynomial(row2[1].terms_), *new Polynomial(row2[2].terms_) };

		q = row1[2] / row2[2];

		cudaMalloc((void **)&dev_row1, sizeof(row1));
		cudaMalloc((void **)&dev_row2, sizeof(row2));
		cudaMalloc((void **)&dev_q, sizeof(q));
		
		parallelColumns<<<3, 1>>>(dev_row1, dev_row2, dev_q);

		row2[0] = row1[0] - (row2[0] * q);
		row2[1] = row1[1] - (row2[1] * q);
		row2[2] = row1[2] - (row2[2] * q);
		row1 = { *new Polynomial(temp[0].terms_), *new Polynomial(temp[1].terms_), *new Polynomial(temp[2].terms_) };

	}

	return row1[2];
}

Polynomial Polynomial::multiGCD(std::vector<Polynomial> polynomials) {
	Polynomial result;

	for (int i = 0; i < polynomials.size() - 1; i++) {
		result = (i == 0) ? GCD(polynomials[i], polynomials[i + 1]) : GCD(result, polynomials[i + 1]);
	}

	return result;
}*/

Polynomial operator*(Polynomial p1, Polynomial p2) {
	Polynomial result;

	for (auto t1 : p1.terms_) {
		for (auto t2 : p2.terms_) {
			result.addUnique(t1 * t2);
		}
	}

	result.deleteZeroElements();
	result.sortByPower();

	return result;
}

Polynomial operator/(Polynomial p1, Polynomial p2) {
	std::list<Term> result = *new std::list<Term>(), buf;
	Term nextResultTerm;
	Polynomial temp;

	while ((*p1.terms_.begin()).getPower() >= (*p2.terms_.begin()).getPower()) {
		nextResultTerm = (*p1.terms_.begin()) / (*p2.terms_.begin());
		result.push_back(nextResultTerm);

		buf = *new std::list<Term>();
		buf.push_back(nextResultTerm);
		temp = *new Polynomial(buf) * p2;

		p1 = p1 - temp;
		if (p1.terms_.size() == 1 && (*p1.terms_.begin()).getCoefficient() == "0") {
			break;
		}
	}

	return *new Polynomial(result);
}

Polynomial operator+(Polynomial p1, Polynomial p2) {
	Polynomial result = *new Polynomial(p1.terms_);

	for (auto term : p2.terms_) {
		result.addUnique(term);
	}

	result.deleteZeroElements();
	result.sortByPower();

	return result;
}

Polynomial operator-(Polynomial p1, Polynomial p2) {
	return p1 + p2;
}

std::ostream& operator <<(std::ostream& os, const Polynomial& p) {
	std::string result = "";

	for(Term term : p.terms_) {
		result += term.to_string() + "+";
	}

	return os << result.erase(result.size() - 1, 1);
}

