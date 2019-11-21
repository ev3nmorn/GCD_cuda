#include "Polynomial.cuh"

void Polynomial::sort_by_power() {
	std::sort(terms_.begin(), terms_.end(), [](Term& t1, Term& t2) { return t1.get_power() > t2.get_power();  });
}

void Polynomial::delete_zero_elements() {
	terms_.erase(std::remove_if(terms_.begin(), terms_.end(), [](Term t) { return t.get_coefficient() == -1; }), terms_.end());

	if (terms_.size() == 0) {
		terms_.push_back(*new Term(-1, 0, field));
	}
}

bool Polynomial::has_equal_power(Term term) {
	for (Term t : terms_) {
		if (t.get_power() == term.get_power()) {
			return true;
		}
	}

	return false;
}

void Polynomial::add_to_equal_power(Term term) {
	for (auto it = terms_.begin(); it != terms_.end(); it++) {
		if ((*it).get_power() == term.get_power()) {
			(*it) = (*it) + term;
		}
	}
}

void Polynomial::add_unique(Term term) {
	if (has_equal_power(term)) {
		add_to_equal_power(term);
	}
	else {
		terms_.push_back(term);
	}
}

Polynomial::Polynomial(GF2m field) {
	this->field = field;
	terms_ = *new std::vector<Term>();
	terms_.push_back(*new Term(-1, 0, field));
}

Polynomial::Polynomial(std::vector<Term> terms, GF2m field) {
	this->field = field;
	terms_ = *new std::vector<Term>();
	terms_ = terms;
	delete_zero_elements();
	sort_by_power();
}

Polynomial operator*(Polynomial p1, Polynomial p2) {
	Polynomial result(p1.field);

	for (auto t1 : p1.terms_) {
		for (auto t2 : p2.terms_) {
			result.add_unique(t1 * t2);
		}
	}

	result.delete_zero_elements();
	result.sort_by_power();

	return result;
}

Polynomial operator/(Polynomial p1, Polynomial p2) {
	std::vector<Term> result = *new std::vector<Term>(), buf;
	Term nextResultTerm(p1.field);
	Polynomial temp(p1.field);

	while ((*p1.terms_.begin()).get_power() >= (*p2.terms_.begin()).get_power()) {
		nextResultTerm = (*p1.terms_.begin()) / (*p2.terms_.begin());
		result.push_back(nextResultTerm);

		buf = *new std::vector<Term>();
		buf.push_back(nextResultTerm);
		temp = *new Polynomial(buf, p1.field) * p2;

		p1 = p1 - temp;
		if (p1.terms_.size() == 1 && (*p1.terms_.begin()).get_coefficient() == -1) {
			break;
		}
	}

	return *new Polynomial(result, p1.field);
}

Polynomial operator+(Polynomial p1, Polynomial p2) {
	Polynomial result = *new Polynomial(p1.terms_, p1.field);

	for (auto term : p2.terms_) {
		result.add_unique(term);
	}

	result.delete_zero_elements();
	result.sort_by_power();

	return result;
}

Polynomial operator-(Polynomial p1, Polynomial p2) {
	return p1 + p2;
}

std::ostream& operator<<(std::ostream &os, Polynomial& poly) {
	std::string result = "";

	for (Term term : poly.get_terms()) {
		result += term.to_string() + "+";
	}

	return os << result.erase(result.size() - 1, 1);
}