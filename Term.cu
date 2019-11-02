#include "Term.cuh"

Term::Term(std::string coefficient, int power) {
	coefficient_ = coefficient;
	power_ = power;
}

Term operator *(Term t1, Term t2) {
	if (t1.coefficient_ == "0" || t2.coefficient_ == "0") {
		return *new Term("0", 0);
	}

	if (t1.coefficient_ == "1") {
		return *new Term(t2.coefficient_, t1.power_ + t2.power_);
	}

	if (t2.coefficient_ == "1") {
		return *new Term(t1.coefficient_, t1.power_ + t2.power_);
	}

	int indT1 = t1.coefficient_.find('a'),
		degreeT1 = std::stoi(t1.coefficient_.substr(indT1 + 1, t1.coefficient_.size() - indT1 - 1)),
		indT2 = t2.coefficient_.find('a'),
		degreeT2 = std::stoi(t2.coefficient_.substr(indT2 + 1, t2.coefficient_.size() - indT2 - 1)),
		newDegree = (degreeT1 + degreeT2) % (int)(GF2m::getElementsCount() - 1);

	if (newDegree == 0)
	{
		return *new Term("1", t1.power_ + t2.power_);
	}
	else
	{
		return *new Term("a" + std::to_string(newDegree), t1.power_ + t2.power_);
	}
}

Term operator /(Term t1, Term t2) {
	if (t1.coefficient_ == "0") {
		return *new Term("0", 0);
	}

	if (t2.coefficient_ == "0") {
		throw new std::exception("divide by zero element");
	}

	if (t2.coefficient_ == "1") {
		return *new Term(t1.coefficient_, t1.power_ - t2.power_);
	}

	if (t1.coefficient_ == "1") {
		int ind = t2.coefficient_.find('a'),
			degree = std::stoi(t2.coefficient_.substr(ind + 1, t2.coefficient_.size() - ind - 1));
		return *new Term("a" + std::to_string(GF2m::getElementsCount() - degree - 1), t1.power_ - t2.power_);
	}



	int indT1 = t1.coefficient_.find('a'),
		degreeT1 = std::stoi(t1.coefficient_.substr(indT1 + 1, t1.coefficient_.size() - indT1 - 1)),
		indT2 = t2.coefficient_.find('a'),
		degreeT2 = std::stoi(t2.coefficient_.substr(indT2 + 1, t2.coefficient_.size() - indT2 - 1)),
		newDegree = degreeT1 - degreeT2;

	if (newDegree < 0) {
		newDegree = GF2m::getElementsCount() - 1 + newDegree;
	}

	if (newDegree == 0) {
		return *new Term("1", t1.power_ - t2.power_);
	}
	else {
		return *new Term("a" + std::to_string(newDegree), t1.power_ - t2.power_);
	}
}

Term operator +(Term t1, Term t2) {
	if (t1.power_ != t2.power_) {
		throw std::exception("Not equal powers in sum operator");
	}

	uint32_t newCoefficientValue = GF2m::getElements()[t1.coefficient_] ^ GF2m::getElements()[t2.coefficient_];
	std::string newCoefficient = "";

	for(auto element : GF2m::getElements()) {
		if (newCoefficientValue == element.second) {
			newCoefficient = element.first;
			break;
		}
	}

	return *new Term(newCoefficient, t1.power_);
}

Term operator -(Term t1, Term t2) {
	return t1 + t2;
}

std::string Term::to_string() {
	std::string result = "";

	if (coefficient_ == "1") {
		if (power_ == 0) {
			result = coefficient_;
		}
		else if (power_ == 1) {
			result = "x";
		}
		else {
			result = "x^" + std::to_string(power_);
		}
	}
	else {
		if (power_ == 0) {
			result = coefficient_;
		}
		else if (power_ == 1) {
			result = coefficient_ + "*x";
		}
		else {
			result = coefficient_ + "*x^" + std::to_string(power_);
		}
	}

	return result;
}