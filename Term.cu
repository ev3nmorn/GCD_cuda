#include "Term.cuh"

Term::Term(int coefficient, int power, GF2m field) {
	this->field = field;
	coefficient_ = coefficient;
	power_ = power;
}

Term operator *(Term t1, Term t2) {
	if (t1.coefficient_ == -1 || t2.coefficient_ == -1) {
		return Term(-1, 0, t1.field);
	}

	if (t1.coefficient_ == 0) {
		return Term(t2.coefficient_, t1.power_ + t2.power_, t2.field);
	}

	if (t2.coefficient_ == 0) {
		return Term(t1.coefficient_, t1.power_ + t2.power_, t1.field);
	}

	int newDegree = (t1.coefficient_ + t2.coefficient_) % (int)(t1.field.get_elements_count() - 1);

	if (newDegree == 0) {
		return Term(0, t1.power_ + t2.power_, t1.field);
	}
	else {
		return Term(newDegree, t1.power_ + t2.power_, t1.field);
	}
}

Term operator /(Term t1, Term t2) {
	if (t1.coefficient_ == -1) {
		return Term(-1, 0, t1.field);
	}
	if (t2.coefficient_ == 0) {
		return Term(t1.coefficient_, t1.power_ - t2.power_, t1.field);
	}
	if (t1.coefficient_ == 0) {
		return Term(t1.field.get_elements_count() - t2.coefficient_ - 1, t1.power_ - t2.power_, t1.field);
	}

	int newDegree = t1.coefficient_ - t2.coefficient_;

	if (newDegree < 0) {
		newDegree = t1.field.get_elements_count() - 1 + newDegree;
	}

	return Term(newDegree, t1.power_ - t2.power_, t1.field);
}

Term operator +(Term t1, Term t2) {
	uint32_t newCoefficientValue = t1.field.get_elements()[t1.coefficient_] ^ t1.field.get_elements()[t2.coefficient_];
	int newCoefficient;

	for (auto element : t1.field.get_elements()) {
		if (newCoefficientValue == element.second) {
			newCoefficient = element.first;
			break;
		}
	}

	return Term(newCoefficient, t1.power_, t1.field);
}

Term operator -(Term t1, Term t2) {
	return t1 + t2;
}

std::string Term::to_string() {
	std::string result = "";

	if (coefficient_ == 0) {
		if (power_ == 0) {
			result = "1";
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
			result = "a" + std::to_string(coefficient_);
		}
		else if (power_ == 1) {
			result = "a" + std::to_string(coefficient_) + "*x";
		}
		else {
			result = "a" + std::to_string(coefficient_) + "*x^" + std::to_string(power_);
		}
	}

	return result;
}