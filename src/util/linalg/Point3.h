#ifndef __POINT3_H
#define __POINT3_H

#include "Vector3.h"

// Empty tag class.
class Origin {};
extern Origin ORIGIN;

/**
 * \brief A point in three-dimensional space.
 *
 * This is one of the basic vector algebra classes. It represents a three
 * dimensional point in space defined by the coordinates X, Y and Z.
 * There are two instances of this template point class: \c Point3 is for floating-point points
 * and \c Point3I is for integer points with three components.
 *
 * Note that there is also a class called Vector_3 that is used for vector
 * in a three dimensional coordinate system.
 */
template<typename TX, typename TY, typename TZ>
class Point_3
{
public:
	/// \brief The X coordinate of the point.
	TX X;
	/// \brief The Y coordinate of the point.
	TY Y;
	/// \brief The Z coordinate of the point.
	TZ Z;

	/////////////////////////////// Constructors /////////////////////////////////

	/// \brief Constructs a point without initializing its coordinates.
	/// \note All coordinates are left uninitialized by this constructor and will therefore have a random value!
	Point_3() {}

	/// \brief Constructs a point with all coordinates set to the given value.
	/// \param val The value to be assigned to each of the points's coordinates.
	//Point_3(T val) { X = Y = Z = val; }

	/// \brief Initializes the coordinates of the point with the given values.
	/// \param x The X coordinate of the new point.
	/// \param y The Y coordinate of the new point.
	/// \param z The Z coordinate of the new point.
	Point_3(TX x, TY y, TZ z) : X(x), Y(y), Z(z) {}

	/// \brief Initializes the coordinates of the point with the values in the given array.
	/// \param val An array of three coordinates.
	//Point_3(T val[3]) : X(val[0]), Y(val[1]), Z(val[2]) {}

	/// \brief Initializes the point to the origin. All coordinates are set to zero.
	/// \param ORIGIN A dummy parameter to distinguish this overloaded constructor from the others.
	///        When using this constructor, just use the special value \c ORIGIN here.
	Point_3(Origin ORIGIN) : X((TX)0), Y((TY)0), Z((TZ)0) {}

	/// \brief Casts the point to a point with another data type.
	template<typename TX2, typename TY2, typename TZ2>
	operator Point_3<TX2,TY2,TZ2>() const { return Point_3<TX2,TY2,TZ2>((TX2)X, (TY2)Y, (TZ2)Z); }

	///////////////////////////// Assignment operators ///////////////////////////

	/// \brief Adds a vector to this point and stores the result in this point.
	/// \param v The vector to add to this point.
	/// \return A reference to \c this point, which has been changed.
	Point_3<TX,TY,TZ>& operator+=(const Vector_3<TX,TY,TZ>& v) { X += v.X; Y += v.Y; Z += v.Z; return *this; }

	/// \brief Subtracts a vector from this point and stores the result in this point.
	/// \param v The vector to subtract from this point.
	/// \return A reference to \c this point, which has been changed.
	Point_3<TX,TY,TZ>& operator-=(const Vector_3<TX,TY,TZ>& v) { X -= v.X; Y -= v.Y; Z -= v.Z; return *this; }

	/// \brief Adds a point to this point.
	/// \param v The point to add to this point.
	/// \return A reference to \c this point, which has been changed.
	Point_3<TX,TY,TZ>& operator+=(const Point_3<TX,TY,TZ>& p) { X += p.X; Y += p.Y; Z += p.Z; return *this; }

	//////////////////////////// Coordinate read access //////////////////////////

	/// Returns the X coordinate of the point.
	const TX& x() const { return X; }

	/// Returns the Y coordinate of the point.
	const TY& y() const { return Y; }

	/// Returns the Z coordinate of the point.
	const TZ& z() const { return Z; }

	//////////////////////////// Coordinate write access //////////////////////////

	/// \brief Sets the X coordinate of this point to a new value.
	/// \param value The new value that is assigned to the point coordinate.
	void setx(const TX& value) { X = value; }

	/// \brief Sets the Y coordinate of this point to a new value.
	/// \param value The new value that is assigned to the point coordinate.
	void sety(const TY& value) { Y = value; }

	/// \brief Sets the Z coordinate of this point to a new value.
	/// \param value The new value that is assigned to the point coordinate.
	void setz(const TZ& value) { Z = value; }

	////////////////////////////////// Comparison ////////////////////////////////

	/// \brief Compares two points for equality.
	/// \param p The point this point should be compared to.
	/// \return \c true if each of the coordinates are equal; \c false otherwise.
	bool operator==(const Point_3<TX,TY,TZ>& p) const { return (p.X==X) && (p.Y==Y) && (p.Z==Z); }

	/// \brief Compares two points for inequality.
	/// \param p The point this point should be compared to.
	/// \return \c true if any of the coordinates are not equal; \c false if all are equal.
	bool operator!=(const Point_3<TX,TY,TZ>& p) const { return (p.X!=X) || (p.Y!=Y) || (p.Z!=Z); }

	/// \brief Checks whether the point is the origin, i.e. all coordinates are zero.
	/// \param ORIGIN Just use the special value \c ORIGIN here.
	/// \return \c true if all of the null vector are zero; \c false otherwise
	bool operator==(const Origin& ORIGIN) const { return (X==(TX)0) && (Y==(TY)0) && (Z==(TZ)0); }

	/// \brief Checks whether the point is not the origin, i.e. any of the coordinates is nonzero.
	/// \param ORIGIN Just use the special value \c ORIGIN here.
	/// \return \c true if any of the coordinates is nonzero; \c false if this is the null vector otherwise
	bool operator!=(const Origin& ORIGIN) const { return (X!=(TX)0) || (Y!=(TY)0) || (Z!=(TZ)0); }

	/////////////////////////////// Binary operators /////////////////////////////

	/// \brief Computes the sum of two points.
	/// \param p The second operand.
	/// \return The sum of two points.
	Point_3<TX,TY,TZ> operator+(const Point_3<TX,TY,TZ>& p) const { return Point_3<TX,TY,TZ>(X + p.X, Y + p.Y, Z + p.Z); }

	/// \brief Computes the sum of a point and a vector.
	/// \param v The second operand.
	/// \return The new point.
	Point_3<TX,TY,TZ> operator+(const Vector_3<TX,TY,TZ>& v) const { return Point_3<TX,TY,TZ>(X + v.X, Y + v.Y, Z + v.Z); }

	/// \brief Computes the difference of two points.
	/// \param p The second operand.
	/// \return The vector connecting the two points.
	Vector_3<TX,TY,TZ> operator-(const Point_3<TX,TY,TZ>& p) const { return Vector_3<TX,TY,TZ>(X - p.X, Y - p.Y, Z - p.Z); }

	/// \brief Substracts a vector from a point.
	/// \param v The second operand.
	/// \return The new point.
	Point_3<TX,TY,TZ> operator-(const Vector_3<TX,TY,TZ>& v) const { return Point_3<TX,TY,TZ>(X - v.X, Y - v.Y, Z - v.Z); }

	/// \brief Converts the point to a vector.
	/// \param ORIGIN Just use the special value \c ORIGIN here.
	/// \return A vector with its components equal to the coordinates of this point.
	const Vector_3<TX,TY,TZ>& operator-(Origin ORIGIN) const {
		SIMULATION_STATIC_ASSERT(sizeof(Point_3<TX,TY,TZ>) == sizeof(Vector_3<TX,TY,TZ>));
		return *reinterpret_cast< const Vector_3<TX,TY,TZ>* >(this);
	}

	/// \brief Multiplies all coordinates of a point with a scalar.
	/// \param s The scalar value.
	/// \return A new point with the scaled coordinates.
	template<typename T>
	Point_3<TX,TY,TZ> operator*(T s) const { return Point_3<TX,TY,TZ>(X*s, Y*s, Z*s); }

	/// \brief Divides all coordinates of a point by a scalar.
	/// \param s The scalar value.
	/// \return A new point with the scaled coordinates.
	template<typename T>
	Point_3<TX,TY,TZ> operator/(T s) const { return Point_3<TX,TY,TZ>(X/s, Y/s, Z/s); }
};

/// \brief Converts a vector to a point.
/// \param ORIGIN Just use the special value \c ORIGIN here.
/// \param v The vector to convert.
/// \return A point with its coordinates equal to the components of the vector \a v.
template<typename TX, typename TY, typename TZ>
Point_3<TX,TY,TZ> operator+(const Origin& ORIGIN, const Vector_3<TX,TY,TZ>& v) {
	return Point_3<TX,TY,TZ>(v.X, v.Y, v.Z);
}

/// \brief Converts a vector to a point.
/// \param ORIGIN Just use the special value \c ORIGIN here.
/// \param v The negative vector to convert.
/// \return A point with its coordinates equal to the negative components of the vector \a v.
template<typename TX, typename TY, typename TZ>
Point_3<TX,TY,TZ> operator-(const Origin& ORIGIN, const Vector_3<TX,TY,TZ>& v) {
	return Point_3<TX,TY,TZ>(-v.X, -v.Y, -v.Z);
}

/// \brief Calculates the squared distance between two points.
/// \param a The first point.
/// \param b The second point.
/// \return The squared distance between \a a and \a b.
template<typename T>
inline T DistanceSquared(const Point_3<T,T,T>& a, const Point_3<T,T,T>& b) {
	return square(a.X - b.X) + square(a.Y - b.Y) + square(a.Z - b.Z);
}

/// \brief Calculates the distance between two points.
/// \param a The first point.
/// \param b The second point.
/// \return The distance between \a a and \a b.
template<typename T>
inline T Distance(const Point_3<T,T,T>& a, const Point_3<T,T,T>& b) {
	return (T)sqrt(DistanceSquared(a, b));
}

/// Returns the index of the coordinate with the maximum value.
/// Post-Condition: 0 <= Return Value < 3
template<typename T>
inline int MaxComponent(const Point_3<T,T,T>& a) {
    return ((a.X >= a.Y) ? ((a.X >= a.Z) ? 0 : 2) : ((a.Y >= a.Z) ? 1 : 2));
}

/// Returns the index of the coordinate with the minimum value.
/// Post-Condition: 0 <= Return Value < 3
template<typename T>
inline int MinComponent(const Point_3<T,T,T>& a) {
    return ((a.X <= a.Y) ? ((a.X <= a.Z) ? 0 : 2) : ((a.Y <= a.Z) ? 1 : 2));
}

/// Returns the index of the component with the maximum absolute value.
/// Post-Condition: 0 <= Return Value < 3
template<typename T>
inline int MaxAbsComponent(const Point_3<T,T,T>& a) {
    return ((fabs(a.X) >= fabs(a.Y)) ? ((fabs(a.X) >= fabs(a.Z)) ? 0 : 2) : ((fabs(a.Y) >= fabs(a.Z)) ? 1 : 2));
}

/// \brief Writes the point to a text output stream.
/// \param os The output stream.
/// \param p The point to write to the output stream \a os.
/// \return The output stream \a os.
template<typename TX, typename TY, typename TZ>
inline std::ostream& operator<<(std::ostream &os, const Point_3<TX,TY,TZ> &p) {
	return os << '(' << p.X << ' ' << p.Y  << ' ' << p.Z << ')';
}

/**
 * \fn typedef Point3
 * \brief Template class instance of the Point_3 class used for floating-point points.
 */
typedef Point_3<CAFloat,CAFloat,CAFloat>	Point3;

/**
 * \fn typedef Point3I
 * \brief Template class instance of the Point_3 class used for integer points.
 */
typedef Point_3<int,int,int>				Point3I;

typedef Point_3<double,double,double>		Point3D;

#endif // __POINT3_H
