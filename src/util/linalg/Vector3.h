#ifndef __VECTOR3_H
#define __VECTOR3_H

// Empty tag class.
class NullVector {};
// This dummy variable is needed to initialize the Vector3 class to default values.
extern NullVector NULL_VECTOR;

/**
 * \brief A vector with three components X, Y, Z.
 *
 * This is one of the basic vector algebra classes. It represents a three
 * dimensional vector in space. There are two instances of this template
 * vector class: \c Vector3 is for floating-point vectors and \c Vector3I is
 * for integer vectors with three components.
 *
 * Note that there is also a class called Point_3 that is used for points
 * in a three dimensional coordinate system.
 */
template<typename TX, typename TY, typename TZ>
class Vector_3
{
public:
	/// \brief The X component of the vector.
	TX X;
	/// \brief The Y component of the vector.
	TY Y;
	/// \brief The Z component of the vector.
	TZ Z;

	/////////////////////////////// Constructors /////////////////////////////////

	/// \brief Constructs a vector without initializing its component.
	/// \note All components are left uninitialized by this constructor and will therefore have a random value!
	Vector_3() {}

	/// \brief Constructs a vector with all components set to the given value.
	/// \param val The value to be assigned to each of the vector's components.
	//explicit Vector_3(ValueType val) { X = Y = Z = val; }

	/// \brief Initializes the components of the vector with the given component values.
	/// \param x The X component of the new vector.
	/// \param y The Y component of the new vector.
	/// \param z The Z component of the new vector.
    Vector_3(TX x, TY y, TZ z) : X(x), Y(y), Z(z) {}

	/// \brief Initializes the components of the vector with the values in the given array.
	/// \param val The array that contains the vector components.
    //Vector_3(ValueType val[3]) : X(val[0]), Y(val[1]), Z(val[2]) {}

	/// \brief Initializes the vector to the null vector. All components are set to zero.
	/// \param NULL_VECTOR A dummy parameter to distinguish this overloaded constructor from the others.
	///        When using this constructor, just use the special value \c NULL_VECTOR here.
    Vector_3(NullVector NULL_VECTOR) : X(0), Y(0), Z(0) {}

	/// \brief Casts the vector to a vector with another data type.
	template<typename TX2, typename TY2, typename TZ2>
	operator Vector_3<TX2,TY2,TZ2>() const { return Vector_3<TX2,TY2,TZ2>((TX2)X, (TY2)Y, (TZ2)Z); }

    /////////////////////////////// Unary operators //////////////////////////////

	/// \brief Returns the inverse of the vector.
	/// \return A vector with negated components: (-X, -Y, -Z).
	Vector_3<TX,TY,TZ> operator-() const { return(Vector_3<TX,TY,TZ>(-X, -Y, -Z)); }

	///////////////////////////// Assignment operators ///////////////////////////

	/// \brief Adds another vector to this vector and stores the result in this vector.
	/// \param v The vector to add to this vector.
	/// \return A reference to \c this vector, which has been changed.
	Vector_3<TX,TY,TZ>& operator+=(const Vector_3<TX,TY,TZ>& v) { X += v.X; Y += v.Y; Z += v.Z; return *this; }

	/// \brief Subtracts another vector from this vector and stores the result in this vector.
	/// \param v The vector to subtract from this vector.
	/// \return A reference to \c this vector, which has been changed.
	Vector_3<TX,TY,TZ>& operator-=(const Vector_3<TX,TY,TZ>& v) { X -= v.X; Y -= v.Y; Z -= v.Z; return *this; }

	/// \brief Multiplies each component of the vector with a scalar value and stores the result in this vector.
	/// \param s The scalar value to multiply this vector with.
	/// \return A reference to \c this vector, which has been changed.
	template<typename T>
	Vector_3<TX,TY,TZ>& operator*=(T s) { X *= s; Y *= s; Z *= s; return *this; }

	/// \brief Divides each component of the vector by a scalar value and stores the result in this vector.
	/// \param s The scalar value.
	/// \return A reference to \c this vector, which has been changed.
	template<typename T>
	Vector_3<TX,TY,TZ>& operator/=(T s) { X /= s; Y /= s; Z /= s; return *this; }

	//////////////////////////// Component read access //////////////////////////

	/// \brief Returns the value of the X component of this vector.
	const TX& x() const { return X; }

	/// \brief Returns the value of the Y component of this vector.
	const TY& y() const { return Y; }

	/// \brief Returns the value of the Z component of this vector.
	const TZ& z() const { return Z; }

	//////////////////////////// Component write access //////////////////////////

	/// \brief Sets the X component of this vector to a new value.
	/// \param value The new value that is assigned to the vector component.
	void setx(const TX& value) { X = value; }

	/// \brief Sets the Y component of this vector to a new value.
	/// \param value The new value that is assigned to the vector component.
	void sety(const TY& value) { Y = value; }

	/// \brief Sets the Z component of this vector to a new value.
	/// \param value The new value that is assigned to the vector component.
	void setz(const TZ& value) { Z = value; }

	////////////////////////////////// Comparison ////////////////////////////////

	/// \brief Compares two vectors for equality.
	/// \return true if each of the components are equal; false otherwise.
	bool operator==(const Vector_3<TX,TY,TZ>& v) const { return (v.X==X) && (v.Y==Y) && (v.Z==Z); }

	/// \brief Compares two vectors for inequality.
	/// \return true if any of the components are not equal; false if all are equal.
	bool operator!=(const Vector_3<TX,TY,TZ>& v) const { return (v.X!=X) || (v.Y!=Y) || (v.Z!=Z); }

	/// \brief Checks whether the vector is the null vector, i.e. all components are zero.
	/// \return true if all of the components are zero; false otherwise
	bool operator==(const NullVector& NULL_VECTOR) const { return (X==(TX)0) && (Y==(TY)0) && (Z==(TZ)0); }

	/// \brief Checks whether the vector is not a null vector, i.e. any of the components is nonzero.
	/// \return true if any of the components is nonzero; false if this is the null vector otherwise
	bool operator!=(const NullVector& NULL_VECTOR) const { return (X!=(TX)0) || (Y!=(TY)0) || (Z!=(TZ)0); }

	/////////////////////////////// Binary operators /////////////////////////////

	/// \brief Computes the sum of two vectors.
	/// \param v The second operand.
	/// \return The sum of two vectors.
	Vector_3<TX,TY,TZ> operator+(const Vector_3<TX,TY,TZ>& v) const { return Vector_3<TX,TY,TZ>(X + v.X, Y + v.Y, Z + v.Z); }

	/// \brief Computes the difference of two vectors.
	/// \param v The second operand.
	/// \return The difference of two vectors.
	Vector_3<TX,TY,TZ> operator-(const Vector_3<TX,TY,TZ>& v) const { return Vector_3<TX,TY,TZ>(X - v.X, Y - v.Y, Z - v.Z); }

	/// \brief Computes the product of a vector and a scalar value. All
	///        components of the vector are multiplied by the scalar \a s.
	/// \param s The second operand.
	/// \return The product of the vector with a scalar \a s: (s*X, s*Y, s*Z).
	template<typename T>
	Vector_3<TX,TY,TZ> operator*(T s) const { return Vector_3<TX,TY,TZ>(X*s, Y*s, Z*s); }

	/// \brief Computes the division of a vector by a scalar value. All
	///        components of the vector are divided by the scalar \a s.
	/// \param s The second operand.
	/// \return The division of the vector by a scalar \a s: (X/s, Y/s, Z/s).
	template<typename T>
	Vector_3<TX,TY,TZ> operator/(T s) const { return Vector_3<TX,TY,TZ>(X/s, Y/s, Z/s); }
};

/// \brief Returns the product of the vector with a scalar: (v.X*s, v.Y*s, v.Z*s).
template<typename TX, typename TY, typename TZ, typename T>
inline Vector_3<TX,TY,TZ> operator*(const T& s, const Vector_3<TX,TY,TZ>& v) { return v * s; }

/// \brief Computes the scalar product of two vectors.
/// \param a The first vector.
/// \param b The second vector.
/// \return The scalar vector product: a.X*b.X + a.Y*b.Y + a.Z*b.Z
template<typename T>
inline T DotProduct(const Vector_3<T,T,T>& a, const Vector_3<T,T,T>& b) {
	return a.X*b.X + a.Y*b.Y + a.Z*b.Z;
}

template<typename T>
inline bool isDotProductPositive(const Vector_3<T,T,T>& a, const Vector_3<T,T,T>& b) {
	return a.X*b.X + a.Y*b.Y + a.Z*b.Z > (T)0;
}

/// \brief Computes the cross product of two vectors.
/// \param a The first vector.
/// \param b The second vector.
/// \return The cross product \a a x \a b.
template<typename T>
inline Vector_3<T,T,T> CrossProduct(const Vector_3<T,T,T>& a, const Vector_3<T,T,T>& b) {
	return Vector_3<T,T,T>(a.Y * b.Z - a.Z * b.Y,
				       a.Z * b.X - a.X * b.Z,
				       a.X * b.Y - a.Y * b.X);
}

/// \brief Returns the squared length of a vector.
/// \param a The input vector.
/// \return The squared length of vector \a a: X*X + Y*Y + Z*Z
/// \sa Length()
template<typename T>
inline T LengthSquared(const Vector_3<T,T,T>& a) {
	return a.X*a.X + a.Y*a.Y + a.Z*a.Z;
}

/// \brief Returns the length of a vector.
/// \param a The input vector.
/// \return The length of vector \a a: sqrt(X*X + Y*Y + Z*Z)
/// \sa LengthSquared(), Normalize()
template<typename T>
inline T Length(const Vector_3<T,T,T>& a) {
	return (T)sqrt(LengthSquared(a));
}

/// \brief Normalizes a vector to unit length.
/// \param a The input vector.
/// \return The vector \a a divided by its length.
/// \note If \a is the null vector then an assertion message is generated in debug builds. In release builds the behavior is undefined.
template<typename T>
inline Vector_3<T,T,T> Normalize(const Vector_3<T,T,T>& a) {
	SIMULATION_ASSERT_MSG(a != NullVector(), "Normalize(const Vector3&)", "Cannot normalize a null vector.");
	return a / Length(a);
}

/// \brief Normalizes a vector to unit length only if it is non-zero.
/// \param a The input vector (can be the null vector).
/// \return The vector \a a divided by its length if \a was non-zero;
///         The null vector if \a is close to zero within a small threshold.
template<typename T>
inline Vector_3<T,T,T> NormalizeSafely(const Vector_3<T,T,T>& a) {
	if(a.equals(NullVector(), CAFLOAT_EPSILON)) return NullVector();
	return a / Length(a);
}

/// \brief Finds the index of the component with the maximum value.
/// \param a The input vector.
/// \return The index (0-2) with the largest value.
template<typename T>
inline int MaxComponent(const Vector_3<T,T,T>& a) {
    return ((a.X >= a.Y) ? ((a.X >= a.Z) ? 0 : 2) : ((a.Y >= a.Z) ? 1 : 2));
}

/// \brief Finds the index of the component with the minimum value.
/// \param a The input vector.
/// \return The index (0-2) with the smallest value.
template<typename T>
inline int MinComponent(const Vector_3<T,T,T>& a) {
    return ((a.X <= a.Y) ? ((a.X <= a.Z) ? 0 : 2) : ((a.Y <= a.Z) ? 1 : 2));
}

/// \brief Finds the index of the component with the maximum absolute value.
/// \param a The input vector.
/// \return The index (0-2) with the largest absolute value.
template<typename T>
inline int MaxAbsComponent(const Vector_3<T,T,T>& a) {
    return ((fabs(a.X) >= fabs(a.Y)) ? ((fabs(a.X) >= fabs(a.Z)) ? 0 : 2) : ((fabs(a.Y) >= fabs(a.Z)) ? 1 : 2));
}

/// \brief Writes the vector to a text output stream.
/// \param os The output stream.
/// \param v The vector to write to the output stream \a os.
/// \return The output stream \a os.
template<typename TX, typename TY, typename TZ>
inline std::ostream& operator<<(std::ostream &os, const Vector_3<TX,TY,TZ> &v) {
	return os << '(' << v.X << ' ' << v.Y  << ' ' << v.Z << ')';
}

/**
 * \fn typedef Vector3
 * \brief Template class instance of the Vector_3 class used for floating-point vectors.
 */
typedef Vector_3<CAFloat,CAFloat,CAFloat>		Vector3;

/**
 * \fn typedef Vector3I
 * \brief Template class instance of the Vector_3 class used for integer vectors.
 */
typedef Vector_3<int,int,int>					Vector3I;

/// The three unit vectors.
extern Vector3 unitVectors[3];

#endif // __VECTOR3_H

