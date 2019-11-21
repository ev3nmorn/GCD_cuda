#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

template<class T>
struct shift {
	T operand1_, operand2_, highBitPoly_, fieldMod_;

	shift(T operand1, T operand2, T highBitPoly, T fieldMod) {
		operand1_ = operand1;
		operand2_ = operand2;
		highBitPoly_ = highBitPoly;
		fieldMod_ = fieldMod;
	}
	__host__ __device__ T operator()(T &x) const {
		if ((operand2_ & (1 << x)) == 0) {
			return 0;
		}

		unsigned int result = operand1_, hiBitSet;
		for (int i = 0; i < x; i++) {
			hiBitSet = result & (highBitPoly_ >> 1);
			result <<= 1;
			if (hiBitSet != 0) {
				result ^= fieldMod_;
			}
		}
		return result & (highBitPoly_ - 1);
	}
};