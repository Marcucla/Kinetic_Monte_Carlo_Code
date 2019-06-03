#ifndef __MATRIX3_H
#define __MATRIX3_H

#include "Vector3.h"

// Empty tag class.
class IdentityMatrix {};
extern IdentityMatrix IDENTITY;

// Empty tag class.
class NullMatrix {};
extern NullMatrix NULL_MATRIX;

/**
 * \brief A 3x3 matrix class.
 *
 * This class stores an array of 3 times 3 floating-point value.
 * It is therefore a matrix with 3 rows and 3 columns.
 *
 * Such a matrix is used to describe affine transformations in 3d space
 * like rotation, shearing and scaling.
 */
class Matrix3
{
public:
	/// The 3 x 3 elements of the matrix.
	/// Elements are stored in column-major order, i.e. the first
	/// array index specifies the columns and the second index the row.
	CAFloat m[3][3];

public:

	/// \brief Constructs a matrix without initializing its elements.
	/// \note All elements are left uninitialized by this constructor and have therefore a random value!
	Matrix3() {}

	/// \brief Constructor that initializes all 9 elements of the matrix to the given values.
	/// \note Values are given in row-major order, i.e. row by row.
	Matrix3(CAFloat m11, CAFloat m12, CAFloat m13,
			CAFloat m21, CAFloat m22, CAFloat m23,
			CAFloat m31, CAFloat m32, CAFloat m33)
	{
		m[0][0] = m11; m[0][1] = m21; m[0][2] = m31;
		m[1][0] = m12; m[1][1] = m22; m[1][2] = m32;
		m[2][0] = m13; m[2][1] = m23; m[2][2] = m33;
	}

	/// \brief Constructor that initializes the three columns.
	explicit Matrix3(const Vector3& col1, const Vector3& col2, const Vector3& col3)
	{
		m[0][0] = col1.X; m[0][1] = col1.Y; m[0][2] = col1.Z;
		m[1][0] = col2.X; m[1][1] = col2.Y; m[1][2] = col2.Z;
		m[2][0] = col3.X; m[2][1] = col3.Y; m[2][2] = col3.Z;
	}

	/// \brief Constructor that initializes the three columns.
	explicit Matrix3(const Vector3 columns[3]) {
		memcpy(m, columns, sizeof(m));
	}

	/// \brief Initializes the matrix to the null matrix.
	/// All matrix elements are set to zero by this constructor.
	Matrix3(const NullMatrix& NULL_MATRIX) {
		m[0][0] = 0; m[0][1] = 0; m[0][2] = 0;
		m[1][0] = 0; m[1][1] = 0; m[1][2] = 0;
		m[2][0] = 0; m[2][1] = 0; m[2][2] = 0;
	}

	/// \brief Initializes the matrix to the identity matrix.
	/// All diagonal elements are set to one and all off-diagonal elements are set to zero.
	Matrix3(const IdentityMatrix& IDENTITY) {
		m[0][0] = 1.0; m[0][1] = 0.0; m[0][2] = 0.0;
		m[1][0] = 0.0; m[1][1] = 1.0; m[1][2] = 0.0;
		m[2][0] = 0.0; m[2][1] = 0.0; m[2][2] = 1.0;
	}

	/// \brief Returns the value of a matrix element.
	/// \param row The row of the element to return (0-2).
	/// \param col The column of the element to return (0-2).
	/// \return The value of the matrix element.
	CAFloat operator()(int row, int col) const {
		SIMULATION_ASSERT_MSG(row >= 0 && row<3, "Matrix3", "Row index out of range");
		SIMULATION_ASSERT_MSG(col >= 0 && col<3, "Matrix3", "Column index out of range");
		return m[col][row];
	}

	/// \brief Returns a reference to a matrix element.
	/// \param row The row of the element to return (0-2).
	/// \param col The column of the element to return (0-2).
	CAFloat& operator()(int row, int col) {
		SIMULATION_ASSERT_MSG(row >= 0 && row<3, "Matrix3", "Row index out of range");
		SIMULATION_ASSERT_MSG(col >= 0 && col<3, "Matrix3", "Column index out of range");
		return m[col][row];
	}

	/// \brief Returns a column from the matrix.
	/// \param i The column to return (0-2).
	/// \return The i-th column of the matrix as a vector.
	const Vector3& column(int i) const {
		SIMULATION_ASSERT_MSG(i >= 0 && i<3, "Matrix3::column()", "Column index out of range.");
		SIMULATION_STATIC_ASSERT(sizeof(m[i]) == sizeof(Vector3));
		return *reinterpret_cast<const Vector3*>(&m[i][0]);
	}

	/// \brief Returns a reference to a column in the matrix.
	/// \param i The column to return (0-2).
	/// \return The i-th column of the matrix as a vector reference. Modifying the vector modifies the matrix.
	Vector3& column(int i) {
		SIMULATION_ASSERT_MSG(i >= 0 && i<3, "Matrix3::column()", "Column index out of range.");
		SIMULATION_STATIC_ASSERT(sizeof(m[i]) == sizeof(Vector3));
		return *reinterpret_cast<Vector3*>(&m[i][0]);
	}

	/// \brief Returns a column from the matrix.
	/// \param i The column to return (0-2).
	/// \return The i-th column of the matrix as a vector.
	const Vector3& getColumn(int i) const { return column(i); }

	/// \brief Sets all elements in one column of the matrix.
	/// \param i The column to set (0-2).
	/// \param c The new element values as a vector.
	void setColumn(int i, const Vector3& c) { column(i) = c; }

	/// \brief Returns a row from the matrix.
	/// \param i The row to return (0-2).
	/// \return The i-th row of the matrix as a vector.
	Vector3 row(int i) const {
		SIMULATION_ASSERT_MSG(i >= 0 && i<3, "Matrix3::row()", "Row index out of range.");
		return Vector3(m[0][i], m[1][i], m[2][i]);
	}

	/// \brief Sets all elements in one row of the matrix.
	/// \param i The row to set (0-2).
	/// \param r The new element values as a vector.
	void setRow(int i, const Vector3& r) {
		SIMULATION_ASSERT_MSG(i >= 0 && i<3, "Matrix3::setRow()", "Row index out of range.");
		m[0][i] = r.X; m[1][i] = r.Y; m[2][i] = r.Z;
	}

	/// \brief Computes the inverse of the matrix.
	///
	/// Determinant must be non-zero.
	Matrix3 inverse() const {
		const CAFloat det = determinant();
		SIMULATION_ASSERT_MSG(fabs(det) > 1e-30, "Matrix3::inverse()", "Singular matrix cannot be inverted: determinant is zero.");
		return Matrix3( (m[1][1]*m[2][2] - m[1][2]*m[2][1])/det,
						(m[2][0]*m[1][2] - m[1][0]*m[2][2])/det,
						(m[1][0]*m[2][1] - m[1][1]*m[2][0])/det,
						(m[2][1]*m[0][2] - m[0][1]*m[2][2])/det,
						(m[0][0]*m[2][2] - m[2][0]*m[0][2])/det,
						(m[0][1]*m[2][0] - m[0][0]*m[2][1])/det,
						(m[0][1]*m[1][2] - m[1][1]*m[0][2])/det,
						(m[0][2]*m[1][0] - m[0][0]*m[1][2])/det,
						(m[0][0]*m[1][1] - m[1][0]*m[0][1])/det);
	}

	/// \brief Computes the determinant of the matrix.
	CAFloat determinant() const {
		return((m[0][0]*m[1][1] - m[0][1]*m[1][0])*(m[2][2])
			  -(m[0][0]*m[1][2] - m[0][2]*m[1][0])*(m[2][1])
			  +(m[0][1]*m[1][2] - m[0][2]*m[1][1])*(m[2][0]));
	}

	/// \brief Returns the transpose of this matrix.
	/// \return A new matrix with columns and rows swapped.
	Matrix3 transposed() const {
		return Matrix3(
			m[0][0], m[0][1], m[0][2],
			m[1][0], m[1][1], m[1][2],
			m[2][0], m[2][1], m[2][2]);
	}

	/// \brief Tests whether this matrix is a pure rotation matrix.
	/// \return \c If the matrix is a pure rotation matrix; \c false otherwise.
	///
	/// Matrix A is a pure rotation matrix if:
	///   (1) det(A) = 1  and
	///   (2) A * A^T = I
	bool isRotationMatrix(CAFloat epsilon = CAFLOAT_EPSILON) const {
		if(fabs(m[0][0]*m[1][0] + m[0][1]*m[1][1] + m[0][2]*m[1][2]) > epsilon) return false;
		if(fabs(m[0][0]*m[2][0] + m[0][1]*m[2][1] + m[0][2]*m[2][2]) > epsilon) return false;
		if(fabs(m[1][0]*m[2][0] + m[1][1]*m[2][1] + m[1][2]*m[2][2]) > epsilon) return false;
		if(fabs(m[0][0]*m[0][0] + m[0][1]*m[0][1] + m[0][2]*m[0][2] - 1.0) > epsilon) return false;
		if(fabs(m[1][0]*m[1][0] + m[1][1]*m[1][1] + m[1][2]*m[1][2] - 1.0) > epsilon) return false;
		if(fabs(m[2][0]*m[2][0] + m[2][1]*m[2][1] + m[2][2]*m[2][2] - 1.0) > epsilon) return false;
		return(fabs(determinant() - 1.0) <= epsilon);
	}

	/// \brief Tests whether this matrix is a pure rotation and/or reflection matrix.
	/// \return \c If the matrix is a pure rotation reflection matrix; \c false otherwise.
	///
	/// Matrix A is a pure rotation reflection matrix if:
	///   (1) |det(A)| = 1  and
	///   (2) A * A^T = I
	bool isRotationReflectionMatrix(CAFloat epsilon = CAFLOAT_EPSILON) const {
		if(fabs(m[0][0]*m[1][0] + m[0][1]*m[1][1] + m[0][2]*m[1][2]) > epsilon) return false;
		if(fabs(m[0][0]*m[2][0] + m[0][1]*m[2][1] + m[0][2]*m[2][2]) > epsilon) return false;
		if(fabs(m[1][0]*m[2][0] + m[1][1]*m[2][1] + m[1][2]*m[2][2]) > epsilon) return false;
		if(fabs(m[0][0]*m[0][0] + m[0][1]*m[0][1] + m[0][2]*m[0][2] - 1.0) > epsilon) return false;
		if(fabs(m[1][0]*m[1][0] + m[1][1]*m[1][1] + m[1][2]*m[1][2] - 1.0) > epsilon) return false;
		if(fabs(m[2][0]*m[2][0] + m[2][1]*m[2][1] + m[2][2]*m[2][2] - 1.0) > epsilon) return false;
		return(fabs(fabs(determinant()) - 1.0) <= epsilon);
	}

	/// \brief Checks whether two matrices are equal within a given tolerance.
	/// \param m The matrix that should be compared to this matrix.
	/// \param tolerance A non-negative threshold for the equality test. The two matrices are considered equal when
	///        the differences in the components are all smaller than this tolerance value.
	/// \return \c true if this matrix is equal to the given matrix \a m within the given tolerance.
	bool equals(const Matrix3& other, CAFloat tolerance = CAFLOAT_EPSILON) const {
		for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
				if(fabs(other.m[i][j] - m[i][j]) > tolerance)
					return false;
		return true;
	}

	/// Returns true if this is the null matrix, i.e. all its elements are exactly zero.
	bool operator==(const NullMatrix& NULL_MATRIX) const {
		for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
				if(m[i][j] != 0.0)
					return false;
		return true;
	}

	/// Returns true if this is not the null matrix, i.e. at least one of its elements is non-zero.
	bool operator!=(const NullMatrix& NULL_MATRIX) const {
		for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
				if(m[i][j] != 0.0)
					return true;
		return false;
	}

	/// Returns true if this is the identity matrix.
	bool operator==(const IdentityMatrix& IDENTITY) const {
		for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
				if(m[i][j] != (i==j))
					return false;
		return true;
	}

	/// Returns true if this is not the identity matrix.
	bool operator!=(const IdentityMatrix& IDENTITY) const {
		for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
				if(m[i][j] != (i==j))
					return true;
		return false;
	}

	/// \brief Adds a 3x3 matrix.
	inline Matrix3& operator+=(const Matrix3& b) {
		for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
				m[i][j] += b.m[i][j];
		return *this;
	}

	/// \brief Subtracts a 3x3 matrix.
	inline Matrix3& operator-=(const Matrix3& b) {
		for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
				m[i][j] -= b.m[i][j];
		return *this;
	}

	/// \brief Multiplies all elements of the matrix with the given scalar.
	inline Matrix3& operator*=(CAFloat s) {
		for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
				m[i][j] *= s;
		return *this;
	}

	/// \brief Divides all elements of the matrix by the given scalar.
	inline Matrix3& operator/=(CAFloat s) {
		for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
				m[i][j] /= s;
		return *this;
	}

	/// Returns a copy of this matrix, where small matrix element close to zero have been
	/// rounded to exactly zero.
	Matrix3 roundZeroElements(CAFloat epsilon = CAFLOAT_EPSILON) const {
		Matrix3 roundedMatrix;
		for(int col = 0; col < 3; col++) {
			for(int row = 0; row < 3; row++) {
				if(fabs(m[col][row]) <= epsilon)
					roundedMatrix.m[col][row] = 0;
				else
					roundedMatrix.m[col][row] = m[col][row];
			}
		}
		return roundedMatrix;
	}

	/// \brief Generates a scaling matrix.
	static Matrix3 scaling(const Vector3& scalingFactors) {
		return Matrix3(scalingFactors.X, 0.0, 0.0,
							0.0, scalingFactors.Y, 0.0,
							0.0, 0.0, scalingFactors.Z);
	}
};

/// \brief Multiplies a 3x3 matrix with a Vector3.
template<typename TX, typename TY, typename TZ>
inline Vector3 operator*(const Matrix3& a, const Vector_3<TX,TY,TZ>& v)
{
	return Vector3(a.m[0][0]*v.X + a.m[1][0]*v.Y + a.m[2][0]*v.Z,
				   a.m[0][1]*v.X + a.m[1][1]*v.Y + a.m[2][1]*v.Z,
				   a.m[0][2]*v.X + a.m[1][2]*v.Y + a.m[2][2]*v.Z);
}

/// \brief Multiplies a 3x3 matrix with a Point3.
template<typename TX, typename TY, typename TZ>
inline Point3 operator*(const Matrix3& a, const Point_3<TX,TY,TZ>& v)
{
	return Point3(a.m[0][0]*v.X + a.m[1][0]*v.Y + a.m[2][0]*v.Z,
				  a.m[0][1]*v.X + a.m[1][1]*v.Y + a.m[2][1]*v.Z,
				  a.m[0][2]*v.X + a.m[1][2]*v.Y + a.m[2][2]*v.Z);
}

/// \brief Multiplies a 3x3 matrix with a 3x3 Matrix.
inline Matrix3 operator*(const Matrix3& a, const Matrix3& b)
{
	Matrix3 m;
	for(int i=0; i<3; i++) {
		for(int j=0; j<3; j++) {
			CAFloat v(0);
			for(int k=0; k<3; k++)
				v += a(i, k) * b(k, j);
			m(i, j) = v;
		}
	}
	return m;
}

/// \brief Multiplies a 3x3 matrix with a scalar value.
/// Each element of the matrix is multiplied by the scalar value.
inline Matrix3 operator*(const Matrix3& a, CAFloat s)
{
	Matrix3 b;
	for(int i=0; i<3; i++)
		for(int j=0; j<3; j++)
			b.m[i][j] = a.m[i][j] * s;
	return b;
}
/// \brief Multiplies a 3x3 matrix with a scalar value.
/// Each element of the matrix is multiplied by the scalar value.
inline Matrix3 operator*(CAFloat s, const Matrix3& a) { return a * s; }

/// \brief Subtracts a 3x3 matrix.
inline Matrix3 operator-(const Matrix3& a, const Matrix3& b)
{
	Matrix3 m;
	for(int i=0; i<3; i++)
		for(int j=0; j<3; j++)
			m(i, j) = a(i,j) - b(i,j);
	return m;
}

/// \brief Adds a 3x3 matrix.
inline Matrix3 operator+(const Matrix3& a, const Matrix3& b)
{
	Matrix3 m;
	for(int i=0; i<3; i++)
		for(int j=0; j<3; j++)
			m(i, j) = a(i,j) + b(i,j);
	return m;
}

/// \brief Writes the matrix to a text output stream.
inline std::ostream& operator<<(std::ostream &os, const Matrix3& m) {
	return os << m.row(0) << std::endl << m.row(1) << std::endl << m.row(2) << std::endl;
}

#endif // __MATRIX3_H
