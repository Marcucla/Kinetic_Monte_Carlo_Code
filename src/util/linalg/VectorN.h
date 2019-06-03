#ifndef __VECTORN_H
#define __VECTORN_H

#include "Vector3.h"

/**
 * A vector with an arbitrary but constant number of components.
 */
template<typename T, std::size_t N>
class VectorN
{
public:
	/// The components of the vector.
	T v[N];

	/////////////////////////////// Constructors /////////////////////////////////

	/// \brief Constructs a vector without initializing its component.
	/// \note All components are left uninitialized by this constructor and will therefore have a random value!
	VectorN() {}

	/// \brief Initializes all components of the vector with the given value.
	/// \param val The value that is assigned to all components of the vector.
	explicit VectorN(T val) { for(size_t k=0; k<N; k++) v[k] = val; }

	/// \brief Initializes the vector to the null vector.
	/// All vector components are set to zero.
    VectorN(NullVector NULL_VECTOR) { for(int k=0; k<N; k++) v[k] = (T)0; }

	/// \brief Initializes the components of the vector with the values in the given array.
    VectorN(T val[N]) { for(size_t k=0; k<N; k++) v[k] = val[k]; }

	/// \brief Initializes the four components of the vector with the given values.
	/// \note This constructor may only be used if the vector template class has exactly four components.
    VectorN(T x, T y, T z, T w) { 
    	SIMULATION_ASSERT_MSG(N==4, "VectorN constructor", "This constructor may only be used for the Vector4 class.");
    	v[0] = x; v[1] = y; v[2] = z; v[3] = w;
    }

    ///////////////////////////// Component access ///////////////////////////////

	/// \brief Returns a reference to the i-th component of the vector.
	/// \param i The index specifying the component to return.
	/// \return A reference to the i-th component that can be used to change the component's value. 
	T& operator[](size_t i) {
		SIMULATION_ASSERT(i<size());
		return v[i]; 
	}

	/// \brief Returns a reference to the i-th component of the vector.
	/// \param i The index specifying the component to return.
	/// \return The i-th component of the vector. 
	const T& operator[](size_t i) const {
		SIMULATION_ASSERT(i<size());
		return v[i]; 
	}

	/// \brief Returns a reference to the i-th component of the vector.
	/// \param i The index specifying the component to return.
	/// \return A reference to the i-th component that can be used to change the component's value. 
	T& operator()(size_t i) {
		SIMULATION_ASSERT(i<size());
		return v[i]; 
	}

	/// \brief Returns a reference to the i-th component of the vector.
	/// \param i The index specifying the component to return.
	/// \return The i-th component of the vector. 
	const T& operator()(size_t i) const {
		SIMULATION_ASSERT(i<size());
		return v[i]; 
	}
	
	/// \brief Returns a pointer to the first element of the vector.
	/// \sa constData() 
	T* data() { 		
		return v;
	}

	/// \brief Returns a pointer to the first element of the vector for read-only access.
	/// \sa data()
	const T* constData() const {
		return v;
	}

    /////////////////////////////// Unary operators //////////////////////////////

	/// \brief Returns the inverse of the vector.
	/// \return A vector with negated components.
	VectorN<T, N> operator-() const { 
		VectorN<T, N> n; 
        for(size_t k=0; k<N; k++) n.v[k] = -v[k];
		return n;
	}

	///////////////////////////// Assignment operators ///////////////////////////

	/// \brief Adds another vector to this vector and stores the result in this vector.
	VectorN<T, N>& operator+=(const VectorN<T, N>& w) { for(size_t k=0; k<N; k++) v[k] += w.v[k]; return *this; }

	/// \brief Subtracts another vector from this vector and stores the result in this vector.
	VectorN<T, N>& operator-=(const VectorN<T, N>& w) { for(size_t k=0; k<N; k++) v[k] -= w.v[k]; return *this; }

	/// \brief Multiplies each component of the vector with a scalar value and stores the result in this vector.
	VectorN<T, N>& operator*=(T s) { for(size_t k=0; k<N; k++) v[k] *= s; return *this; }

	/// \brief Divides each component of the vector by a scalar value and stores the result in this vector.
	VectorN<T, N>& operator/=(T s) { for(size_t k=0; k<N; k++) v[k] /= s; return *this; }

	////////////////////////////////// Comparison ////////////////////////////////

	/// \brief Compares two vector for equality.
	/// \return true if each of the components are equal; false otherwise.
	bool operator==(const VectorN<T, N>& w) const { 
		for(size_t k=0; k<N; k++) if(v[k] != w.v[k]) return false;
		return true;
	}

	/// \brief Compares two vector for inequality.
	/// \return true if any of the components are not equal; false if all are equal.
	bool operator!=(const VectorN<T, N>& w) const { 
		for(size_t k=0; k<N; k++) if(v[k] != w.v[k]) return true;
		return false;
	}

	/// \brief Checks whether the vector is the null vector, i.e. all components are zero.
	/// \return true if all of the components are zero; false otherwise
	bool operator==(const NullVector& NULL_VECTOR) const { 
		for(size_t k=0; k<N; k++) if(v[k] != 0) return false; 
		return true;
	}

	/// \brief Checks whether the vector is not a null vector, i.e. any of the components is nonzero.
	/// \return true if any of the components is nonzero; false if this is the null vector otherwise
	bool operator!=(const NullVector& NULL_VECTOR) const { 
		for(size_t k=0; k<N; k++) if(v[k] != 0) return true;
		return false;
	}

	/// \brief Checks whether two vectors are equal within a given tolerance.
	/// \param w The vector that should be compared to this vector.
	/// \param tolerance A non-negative threshold for the equality test. The two vectors are considered equal when 
	///        the component-wise differences are all smaller than this tolerance value.
	/// \return true if this vector is equal to the given vector within the given tolerance.
	bool equals(const VectorN<T, N>& w, T tolerance) const {
		for(size_t k=0; k<N; k++) if(abs(v[k] - w.v[k]) > tolerance) return false;
		return true;
	}

	/////////////////////////////// Binary operators /////////////////////////////

	/// \brief Computes the sum of two vectors.
	/// \return The sum of two vectors.
	VectorN<T, N> operator+(const VectorN<T, N>& w) const { 
        VectorN<T, N> r(*this);
        r += w;
		return r;
	}

	/// \brief Computes the difference of two vectors.
	/// \return The difference of two vectors.
	VectorN<T, N> operator-(const VectorN<T, N>& v) const { 
        VectorN<T, N> r(*this);
        r -= v;
		return r;
	}

	/// \brief Computes the product of a vector and a scalar value. All
	///        components of the vector are multiplied by the scalar.
	VectorN<T, N> operator*(T s) const { 
        VectorN<T, N> r(*this);
        r *= s;
		return r;
	}

	/// \brief Computes the division of a vector by a scalar value. All
	///        components of the vector are divided by the scalar.
	VectorN<T, N> operator/(T s) const { 
        VectorN<T, N> r(*this);
        r /= s;
		return r;
	}

	/// \brief Returns the scalar (inner or dot product) of two vectors.
	T operator*(const VectorN<T, N>& w) const { 
		T s(0);
		for(size_t k=0; k<N; k++) s += v[k] * w.v[k];
		return s;
	}

	///////////////////////////////// Information ////////////////////////////////
	
	/// \brief Returns the number of components in this vector (the constant N).
	size_t size() const { return N; }

	/// \brief Returns the dimension of this vector (the constant N).
	size_t dimension() const { return N; }
};

/// \brief Computes the scalar product of two vectors.
template<typename T, int N>
inline T DotProduct(const VectorN<T, N>& a, const VectorN<T, N>& b) { 
	return a * b;
}

/// \brief Returns the squared length of a vector.
template<typename T, int N>
inline T LengthSquared(const VectorN<T, N>& a) {
	return a * a;
}

/// \brief Returns the length of a vector.
template<typename T, int N>
inline T Length(const VectorN<T, N>& a) {
	return (T)sqrt(LengthSquared(a));
}

/// \brief Returns the vector divided by its length.
template<typename T, int N>
inline VectorN<T, N> Normalize(const VectorN<T, N>& a) {
	SIMULATION_ASSERT(a != NULL_VECTOR);
	return a / Length(a);
}

/// \brief Writes the vector to a text output stream.
template<typename T, int N>
inline std::ostream& operator<<(std::ostream& os, const VectorN<T, N>& w) {
	for(size_t k=0; k<N; k++) {
		os << w[k];
		if(k != N - 1) os << ' ';
	}
	return os;
}

/** 
 * \fn typedef Vector4
 * \brief Template class instance of the VectorN template used for floating-point vectors. 
 */
typedef VectorN<CAFloat, 4> Vector4;

/// Writes the vector to an output stream.
inline std::ostream& operator<<(std::ostream& os, const Vector4& w) {	
	return os << w[0] << ' ' << w[1] << ' ' << w[2] << ' ' << w[3];
}

#endif

