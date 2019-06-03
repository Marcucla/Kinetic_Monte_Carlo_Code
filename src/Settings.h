#ifndef __SETTINGS_H
#define __SETTINGS_H

/// The code version number.
#define PROGRAM_VERSION_STRING	 					"0.1.0"

// Use 64-bit double as default floating point type.
typedef double CAFloat;

// This tells the program that we're using 64-bit floating point.
#define FLOATTYPE_DOUBLE

/// A small epsilon value for the CAFloat.
#define CAFLOAT_EPSILON	((CAFloat)1e-6)

/// The maximum value for floating point variables.
#define CAFLOAT_MAX	(std::numeric_limits<CAFloat>::max())

// The Pi constant
#define CAFLOAT_PI	((CAFloat)M_PI)

/// Computes the square of a number.
template<typename T> inline constexpr T square(const T& f) { return f*f; }

/// Converts an angle from radians to degrees.
constexpr inline CAFloat rad2deg(CAFloat angle) { return angle * (CAFloat(180) / CAFloat(M_PI)); }

/// Converts an angle from degrees to radians.
constexpr inline CAFloat deg2rad(CAFloat angle) { return angle * (CAFloat(M_PI) / CAFloat(180)); }


#endif // __SETTINGS_H

