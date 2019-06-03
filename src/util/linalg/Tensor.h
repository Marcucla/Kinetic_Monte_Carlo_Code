#ifndef __TENSOR_H
#define __TENSOR_H

#include "Vector3.h"
#include "VectorN.h"
#include "Matrix3.h"

/// A first order tensor is just a one-dimensional vector.
typedef Vector3 Tensor1;

/// A second order tensor is just a two-dimensional matrix.
typedef Matrix3 Tensor2;

class SymmetricTensor2;		// defined below
class SymmetricTensor4;		// defined below

/******************************************************************************
* A symmetric 2nd order tensor.
* Stores only the lower left part of the 3x3 matrix.
******************************************************************************/
class SymmetricTensor2
{
private:
	VectorN<CAFloat, 6> a;	// 6 independent entries.

	// Private vector constructor.
	SymmetricTensor2(const VectorN<CAFloat, 6>& v) : a(v) {}

public:

	/// Default constructor.
	/// Does NOT initialize the tensor entries.
	SymmetricTensor2() {}

	/// Constant value constructor.
	explicit SymmetricTensor2(CAFloat v) : a(v) {}

	/// Element access.
	CAFloat& operator()(size_t row, size_t column) {
		SIMULATION_ASSERT(row<3 && column<3);
		if(row < column) std::swap(row, column);
		switch(row-column) {
			case 0: return a(row);
			case 1: return a(row+2);
			case 2: return a(5);
			default: SIMULATION_ASSERT(false); return a(0);
		}
	}

	/// Element access (constant version).
	const CAFloat& operator()(size_t row, size_t column) const {
		SIMULATION_ASSERT(row<3 && column<3);
		if(row < column) std::swap(row, column);
		switch(row-column) {
			case 0: return a(row);
			case 1: return a(row+2);
			case 2: return a(5);
			default: SIMULATION_ASSERT(false); return a(0);
		}
	}

	/// Element read access.
	CAFloat get(size_t row, size_t col) const { return (*this)(row, col); }
	/// Element write access.
	void set(size_t row, size_t col, CAFloat value) { (*this)(row, col) = value; }

	/// Returns a column of the tensor matrix.
	Vector3 column(size_t col) const {
		return Vector3((*this)(0, col), (*this)(1, col), (*this)(2, col));
	}

	/// Returns a row of the tensor matrix.
	Vector3 row(size_t row) const {
		return Vector3((*this)(row, 0), (*this)(row, 1), (*this)(row, 2));
	}

	/// Subtracts the identity from the tensor.
	SymmetricTensor2 operator-(const IdentityMatrix& IDENTITY) {
		SymmetricTensor2 t(*this);
		t.a(0) -= 1; t.a(1) -= 1; t.a(2) -= 1;
		return t;
	}
	/// Adds another symmetric tensors.
	SymmetricTensor2& operator+=(const SymmetricTensor2& B) {
		for(size_t i=0; i<6; i++) a[i] += B.a[i];
		return *this;
	}

	/// Assigns a symmetric tensor to a full tensor.
	operator Tensor2() const {
		Tensor2 t;
		for(size_t i=0; i<3; i++)
			for(size_t j=0; j<3; j++)
				t(i,j) = (*this)(i,j);
		return t;
	}

	// Multiplies the tensor with a scalar value.
	SymmetricTensor2 operator*(const CAFloat s) const { return SymmetricTensor2(a * s); }
	// Divides the tensor by a scalar value.
	SymmetricTensor2 operator/(const CAFloat s) const { return SymmetricTensor2(a / s); }
	// Divides the tensor by a scalar value.
	SymmetricTensor2& operator/=(const CAFloat s) { a /= s; return *this; }
	// Multiplies a 2nd order tensor with a vector.
	Tensor1 operator*(const Tensor1& v) const {
		return Tensor1(
			v.X*a(0) + v.Y*a(3) + v.Z*a(5),
			v.X*a(3) + v.Y*a(1) + v.Z*a(4),
			v.X*a(5) + v.Y*a(4) + v.Z*a(2));
	}

	/// Computes the eigenvalues and optionally eigenvectors of this symmetric tensor.
	void eigenvalues(CAFloat lambdas[3], Matrix3* eigenvectors = NULL) const;
	/// Finds the maximum absolute eigenvalue of the symmetric tensor.
	CAFloat maxEigenvalue() const;
	/// Finds the minimum absolute eigenvalue of the symmetric tensor.
	CAFloat minEigenvalue() const;

	/// Computes the second invariant of the tensor.
	CAFloat secondInvariant() const { return a(0)*a(1) + a(1)*a(2) + a(0)*a(2) - square(a(3)) - square(a(4)) - square(a(5)); }

	/// Returns a pointer to the internal data array of the tensor where the 6 independent elements are stored (Voigt notation).
	CAFloat* data() { return a.v; }
	/// Returns a pointer to the internal data array of the tensor where the 6 independent elements are stored (Voigt notation).
	const CAFloat* constData() const { return a.v; }

	/// Returns a reference to the internal representation of the tensor as a 6-component vector (Voigt notation).
	VectorN<CAFloat, 6>& voigtVector() { return a; }
	/// Returns a constant reference to the internal representation of the tensor as a 6-component vector (Voigt notation).
	const VectorN<CAFloat, 6>& voigtVector() const { return a; }

private:
	/// Perform Givens rotation of this symmetric tensor.
	void Givens(size_t p, size_t q, Matrix3* V);
};

// Multiplies the tensor with a scalar value.
inline SymmetricTensor2 operator*(const CAFloat s, const SymmetricTensor2& t) { return t * s; }

/// Multiplies a tensor with a symmetric tensor.
inline Tensor2 operator*(const Tensor2& A, const SymmetricTensor2& B)
{
	Tensor2 T2;
	for(size_t i=0; i<3; i++)
		for(size_t j=0; j<3; j++) {
			CAFloat b = 0;
			for(size_t k=0; k<3; k++) b += A(i,k) * B(k,j);
			T2(i,j) = b;
		}
	return T2;
}

/// Computes the sum of two symmetric tensors.
inline SymmetricTensor2 operator+(const SymmetricTensor2& A, const SymmetricTensor2& B)
{
	SymmetricTensor2 C;
	for(size_t i=0; i<6; i++) C.data()[i] = A.constData()[i] + B.constData()[i];
	return C;
}

// Specialized multiplications with symmetric result:

/// Computes A^t * A.
inline SymmetricTensor2 Product_AtA(const Tensor2& A)
{
	SymmetricTensor2 S;
	for(size_t i=0; i<3; i++)
		for(size_t j=0; j<=i; j++) {
			CAFloat b = 0;
			for(size_t k=0; k<3; k++) b += A(k,i) * A(k,j);
			S(i,j) = b;
		}
	return S;
}

/// Computes A * A^t.
inline SymmetricTensor2 Product_AAt(const Tensor2& A)
{
	SymmetricTensor2 S;
	for(size_t i=0; i<3; i++)
		for(size_t j=0; j<=i; j++) {
			CAFloat b = 0;
			for(size_t k=0; k<3; k++) b += A(i,k)*A(j,k);
			S(i,j) = b;
		}
	return S;
}

/// Computes A * S * A^t.
inline SymmetricTensor2 TripleProduct_ASAt(const Tensor2& A, const SymmetricTensor2& S)
{
	Tensor2 AS = A * S;
	SymmetricTensor2 R;
	for(size_t i=0; i<3; i++)
		for(size_t j=0; j<=i; j++) {
			CAFloat b = 0;
			for(size_t k=0; k<3; k++) b += AS(i,k)*A(j,k);
			R(i,j) = b;
		}
	return R;
}

/// Compute the double contraction of two tensors (A : B)
inline CAFloat DoubleContraction(const SymmetricTensor2& A, const SymmetricTensor2& B)
{
	CAFloat d = 0;
	for(size_t i=0; i<3; i++)
		d += A.constData()[i]*B.constData()[i];
	for(size_t i=3; i<6; i++)
		d += 2.0 * A.constData()[i]*B.constData()[i];
	return d;
}

/// Prints the symmetric tensor to a stream.
inline std::ostream& operator<<(std::ostream &os, const SymmetricTensor2& m)
{
	for(size_t row = 0; row < 3; row++)
		os << m(row, 0) << " " << m(row, 1) << " " << m(row, 2) << std::endl;
	return os;
}

/******************************************************************************
* A 4th order tensor with symmetry in minor and major index.
*
* A_ijkl = A_klij = A_jikl
*
* Stores only the lower left part for every symmetry (i>j,k>l,ij>kl)
******************************************************************************/
class SymmetricTensor4
{
private:
	VectorN<CAFloat, 21> a;	// 21 independent entries.

public:

	/// Default constructor.
	/// This does NOT initialize the tensor elements.
	SymmetricTensor4() {}

	/// Constant value constructor.
	explicit SymmetricTensor4(CAFloat v) : a(v) {}

	/// Element access.
	CAFloat& operator()(size_t i, size_t j, size_t k, size_t l) {
		SIMULATION_ASSERT(i<3 && j<3 && k<3 && l<3);
		if(j > i) std::swap(i, j);
		if(l > k) std::swap(k, l);
		size_t I, J;
		switch(i-j) {
			case 0: I = i; break;
			case 1: I = i+2; break;
			case 2: I = 5; break;
		}
		switch(k-l) {
			case 0: J = k; break;
			case 1: J = k+2; break;
			case 2: J = 5; break;
		}
		if(J>I) std::swap(I, J);
		switch(I-J) {
			case 0: return a(I);
			case 1: return a(I+5);
			case 2: return a(I+9);
			case 3: return a(I+12);
			case 4: return a(I+14);
			case 5: return a(20);
			default: SIMULATION_ASSERT(false); return a(0);
		}
	}

	/// Element access (constant version).
	const CAFloat& operator()(size_t i, size_t j, size_t k, size_t l) const {
		SIMULATION_ASSERT(i<3 && j<3 && k<3 && l<3);
		if(j > i) std::swap(i, j);
		if(l > k) std::swap(k, l);
		size_t I, J;
		switch(i-j) {
			case 0: I = i; break;
			case 1: I = i+2; break;
			case 2: I = 5; break;
		}
		switch(k-l) {
			case 0: J = k; break;
			case 1: J = k+2; break;
			case 2: J = 5; break;
		}
		if(J>I) std::swap(I, J);
		switch(I-J) {
			case 0: return a(I);
			case 1: return a(I+5);
			case 2: return a(I+9);
			case 3: return a(I+12);
			case 4: return a(I+14);
			case 5: return a(20);
			default: SIMULATION_ASSERT(false); return a(0);
		}
	}

	/// Element access using reduced 6x6 Voigt notation.
	CAFloat& operator()(size_t I, size_t J) {
		SIMULATION_ASSERT(I<6 && J<6);
		if(J>I) std::swap(I, J);
		switch(I-J) {
			case 0: return a(I);
			case 1: return a(I+5);
			case 2: return a(I+9);
			case 3: return a(I+12);
			case 4: return a(I+14);
			case 5: return a(20);
			default: SIMULATION_ASSERT(false); return a(0);
		}
	}

	/// Element access using reduced 6x6 Voigt notation.
	const CAFloat& operator()(size_t I, size_t J) const {
		SIMULATION_ASSERT(I<6 && J<6);
		if(J>I) std::swap(I, J);
		switch(I-J) {
			case 0: return a(I);
			case 1: return a(I+5);
			case 2: return a(I+9);
			case 3: return a(I+12);
			case 4: return a(I+14);
			case 5: return a(20);
			default: SIMULATION_ASSERT(false); return a(0);
		}
	}

	/// Element read access.
	CAFloat get(size_t I, size_t J) const { return (*this)(I, J); }
	/// Element write access.
	void set(size_t I, size_t J, CAFloat value) { (*this)(I, J) = value; }

	/// Returns a pointer to the internal data array of the tensor where the 21 independent elements are stored.
	CAFloat* data() { return a.v; }
	/// Returns a pointer to the internal data array of the tensor where the 21 independent elements are stored.
	const CAFloat* constData() const { return a.v; }

	/// Returns a reference to the internal representation of the tensor as a 21-component vector.
	VectorN<CAFloat, 21>& vector() { return a; }
	/// Returns a constant reference to the internal representation of the tensor as a 21-component vector.
	const VectorN<CAFloat, 21>& vector() const { return a; }
};

/// Computes the double contraction (A : B) of a double symmetric 4th order tensor A
/// and a symmetric 2nd order tensor B.
inline SymmetricTensor2 operator*(const SymmetricTensor4& A, const SymmetricTensor2 &B)
{
	// Off-diagonal elements occur twice.
	SymmetricTensor2 BB(B);
	BB.data()[3] *= 2.0;
	BB.data()[4] *= 2.0;
	BB.data()[5] *= 2.0;

	// The result.
	SymmetricTensor2 R;
	const CAFloat *Aj = A.constData();
	const CAFloat *bjp = BB.data(), *bjm;
	CAFloat *rjp = R.data(), *rjm;

	for(size_t i=0; i<6; i++)
		*rjp++ = (*Aj++)*(*bjp++);
	for(size_t i=1; i<6; i++) {
		rjp=&R.data()[i];
		rjm=&R.data()[0];
		bjp=&BB.data()[i];
		bjm=&BB.data()[0];
		for(size_t j=i; j<6; j++) {
			*rjm++ += *Aj*(*bjp++);
			*rjp++ += *Aj*(*bjm++);
			Aj++;
		}
	}
	return R;
}

/// Performs a push-forward operation on a 4th order tensor.
/// Calculates B_ijkl = F_ip F_jq F_kr F_ls A_pqrs
inline SymmetricTensor4 PushForward(const SymmetricTensor4& A, const Tensor2& F)
{
	SymmetricTensor4 B;

	for(size_t i=0; i<3; i++)
		for(size_t j=0; j<=i; j++) {
			size_t I;
			switch(i-j) {
				case 0: I = i; break;
				case 1: I = i+2; break;
				case 2: I = 5; break;
				default: SIMULATION_ASSERT(false);
			}

			for(size_t k=0; k<3; k++)
				for(size_t l=0; l<=k; l++) {
					size_t J;
					switch(k-l) {
						case 0: J = k; break;
						case 1: J = k+2; break;
						case 2: J = 5; break;
						default: SIMULATION_ASSERT(false);
					}
					if(J > I) continue;

					CAFloat value = 0.0;
					for(size_t p=0; p<3; p++)
						for(size_t q=0; q<3; q++)
							for(size_t r=0; r<3; r++)
								for(size_t s=0; s<3; s++)
									value += F(i,p) * F(j,q) * F(k,r) * F(l,s) * A(p,q,r,s);

					B(i,j,k,l) = value;
				}
		}
	return B;
}

/// Writes the symmetric tensor to an output stream.
inline std::ostream& operator<<(std::ostream &os, const SymmetricTensor4& m) {
	for(size_t row = 0; row < 6; row++) {
		for(size_t col = 0; col < 6; col++) {
			os << m(row, col) << " ";
		}
		os << std::endl;
	}
	return os;
}

#endif // __TENSOR_H
