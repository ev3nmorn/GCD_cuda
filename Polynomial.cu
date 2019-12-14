#include "Polynomial.cuh"

Polynomial::Polynomial(GF2m field) {
	terms_array_ = new int[1]; terms_array_[0] = 0;
	array_size_ = 1;
	this->field = field;
}

Polynomial::Polynomial(int* terms, int size, GF2m field) {
	this->field = field;

	array_size_ = size;
	terms_array_ = new int[array_size_];
	for (int i = 0; i < array_size_ + 1; i++) {
		terms_array_[i] = terms[i];
	}
}

__global__ void polynomial_mult(int* p1, int* p2, int* result, int* elements, int* elements_size) {
	if (p1[threadIdx.x] == 0 || p2[threadIdx.y] == 0) {
		return;
	}

	int p1_power, p2_power, result_power, result_element;

	for (int i = 0; i < *elements_size; i++) {
		if (elements[i] == p1[threadIdx.x]) {
			p1_power = i;
			break;
		}
	}

	for (int i = 0; i < *elements_size; i++) {
		if (elements[i] == p2[threadIdx.y]) {
			p2_power = i;
			break;
		}
	}

	result_power = (p1_power + p2_power) % *elements_size;
	result_element = elements[result_power];
	atomicXor(&result[threadIdx.x + threadIdx.y], result_element);
}

__global__ void polynomial_sum(int* result_p, int* p2) {
	result_p[threadIdx.x] ^= p2[threadIdx.x];
}

int Polynomial::dirive_elements(int el1, int el2) {
	if (el1 == 0) {
		return 0;
	}
	if (el2 == 1) {
		return el1;
	}

	int el1_power, el2_power;

	for (int i = 0; i < field.get_array_size(); i++) {
		if (field.get_field_elements()[i] == el1) {
			el1_power = i;
			break;
		}
	}

	for (int i = 0; i < field.get_array_size(); i++) {
		if (field.get_field_elements()[i] == el2) {
			el2_power = i;
			break;
		}
	}

	int new_degree = el1_power - el2_power;

	if (new_degree < 0) {
		new_degree = field.get_array_size() + new_degree;
	}

	return field.get_field_elements()[new_degree];
}

int Polynomial::is_normalized() {
	if (terms_array_[array_size_ - 1] != 0) {
		return -1;
	}

	for (int i = array_size_; i >= 0; i--) {
		if (terms_array_[i - 1] != 0) {
			return i;
		}
	}

	return 0;
}

void Polynomial::normalize(int size) {
	if (size == 0) {
		array_size_ = 1;
		delete[] terms_array_;
		terms_array_ = new int[1];
		terms_array_[0] = 0;
		return;
	}

	int* new_terms_array = new int[size];

	for (int i = 0; i < size; i++) {
		new_terms_array[i] = terms_array_[i];
	}
	delete[] terms_array_;

	terms_array_ = new_terms_array;
	array_size_ = size;
}

int* Polynomial::copy_terms() {
	int* result = new int[array_size_];

	for (int i = 0; i < array_size_; i++) {
		result[i] = terms_array_[i];
	}

	return result;
}

Polynomial operator*(Polynomial p1, Polynomial p2) {
	int result_size = p1.array_size_ + p2.array_size_ - 1, size_elements = p1.field.get_array_size();
	int *result_p_host, *p1_host, *p2_host;

	int *result_p_device, *p1_device, *p2_device,
		*field_elements_device, *size_elements_device;

	p1_host = p1.terms_array_;
	p2_host = p2.terms_array_;
	result_p_host = new int[result_size];
	for (int i = 0; i < result_size; i++) {
		result_p_host[i] = 0;
	}

	cudaMalloc((void**)&result_p_device, result_size * sizeof(int));
	cudaMalloc((void**)&p1_device, p1.array_size_ * sizeof(int));
	cudaMalloc((void**)&p2_device, p2.array_size_ * sizeof(int));
	cudaMalloc((void**)&field_elements_device, size_elements * sizeof(unsigned int));
	cudaMalloc((void**)&size_elements_device, sizeof(int));

	cudaMemcpy(result_p_device, result_p_host, result_size * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(p1_device, p1_host, p1.array_size_ * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(p2_device, p2_host, p2.array_size_ * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(field_elements_device, p1.field.get_field_elements(), size_elements * sizeof(unsigned int), cudaMemcpyHostToDevice);
	cudaMemcpy(size_elements_device, &size_elements, sizeof(int), cudaMemcpyHostToDevice);

	dim3 threads_per_block(p1.array_size_, p2.array_size_);
	polynomial_mult << < 1, threads_per_block >> > (p1_device, p2_device, result_p_device, field_elements_device, size_elements_device);

	cudaMemcpy(result_p_host, result_p_device, result_size * sizeof(int), cudaMemcpyDeviceToHost);

	cudaFree(p1_device); cudaFree(p2_device); cudaFree(result_p_device);
	cudaFree(field_elements_device);

	return *new Polynomial(result_p_host, result_size, p1.field);
}

Polynomial operator/(Polynomial p1, Polynomial p2) {
	int result_size = p1.array_size_ - p2.array_size_ + 1;
	int* result = new int[result_size];
	for (int i = 0; i < result_size; i++) {
		result[i] = 0;
	}

	int new_degree;
	Polynomial temp_poly;
	while (p1.array_size_ >= p2.array_size_) {
		new_degree = p1.array_size_ - p2.array_size_;
		result[new_degree] = p1.dirive_elements(p1.terms_array_[p1.array_size_ - 1], p2.terms_array_[p2.array_size_ - 1]);

		int* temp = new int[new_degree + 1];
		for (int i = 0; i < new_degree; i++) {
			temp[i] = 0;
		}
		temp[new_degree] = result[new_degree];
		Polynomial buf(temp, new_degree + 1, p2.field);
		temp_poly = buf * p2;

		p1 = p1 - temp_poly;
		if (p1.array_size_ == 1 && p1.terms_array_[p1.array_size_ - 1] == 0) {
			break;
		}
	}

	return Polynomial(result, result_size, p2.field);
}

Polynomial operator+(Polynomial p1, Polynomial p2) {
	int result_size, p2_size;
	int *result_p_host, *p2_host;

	int *result_p_device, *p2_device;

	if (p1.array_size_ > p2.array_size_) {
		result_p_host = p1.copy_terms();
		result_size = p1.array_size_;
		p2_host = p2.terms_array_;
		p2_size = p2.array_size_;
	}
	else {
		result_p_host = p2.copy_terms();
		result_size = p2.array_size_;
		p2_host = p1.terms_array_;
		p2_size = p1.array_size_;
	}

	cudaMalloc((void**)&result_p_device, result_size * sizeof(int));
	cudaMalloc((void**)&p2_device, p2_size * sizeof(int));

	cudaMemcpy(result_p_device, result_p_host, result_size * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(p2_device, p2_host, p2_size * sizeof(int), cudaMemcpyHostToDevice);

	polynomial_sum << < 1, p2_size >> > (result_p_device, p2_device);

	cudaMemcpy(result_p_host, result_p_device, result_size * sizeof(int), cudaMemcpyDeviceToHost);

	cudaFree(result_p_device); cudaFree(p2_device);

	Polynomial result_polynomial(result_p_host, result_size, p1.field);

	int is_normalized = result_polynomial.is_normalized();
	if (is_normalized != -1) {
		result_polynomial.normalize(is_normalized);
	}

	return result_polynomial;
}

Polynomial operator-(Polynomial p1, Polynomial p2) {
	return p1 + p2;
}

int Polynomial::get_power_by_value(int coefficient) {
	for (int i = 0; i < field.get_array_size(); i++) {
		if (coefficient == field.get_field_elements()[i]) {
			return i;
		}
	}

	return -1;
}

std::string Polynomial::term_to_string(int coefficient, int power) {
	std::string result = "";

	if (coefficient == 1) {
		if (power == 0) {
			result = "1";
		}
		else if (power == 1) {
			result = "x";
		}
		else {
			result = "x^" + std::to_string(power);
		}
	}
	else {
		if (power == 0) {
			result = "a" + std::to_string(get_power_by_value(coefficient));
		}
		else if (power == 1) {
			result = "a" + std::to_string(get_power_by_value(coefficient)) + "*x";
		}
		else {
			result = "a" + std::to_string(get_power_by_value(coefficient)) + "*x^" + std::to_string(power);
		}
	}

	return result;
}

std::ostream& operator<<(std::ostream &os, Polynomial& poly) {
	std::string result = "";

	for (int i = poly.array_size_ - 1; i >= 0; i--) {
		if (poly.terms_array_[i] != 0) {
			result += poly.term_to_string(poly.terms_array_[i], i) + "+";
		}
	}

	if (result == "") {
		result = "0";
	}
	else {
		result = result.erase(result.size() - 1, 1);
	}

	return os << result;
}

//a2*x^3+a10*x^2+a5*x+a8
//a14*x^3+a10*x^2+a5*x+a3