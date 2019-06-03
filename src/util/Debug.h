#ifndef __DEBUG_H
#define __DEBUG_H

/// This points to a function that is called when an assertion check fails.
extern void (*fatalAssertionErrorHandler)(const char* file, int line, const char* location, const char* msg);

// This is for debugging purposes. The SIMULATION_ASSERT() macros are used in the code to check
// if everything runs as expected.

#ifdef DEBUG_SIMULATION

	inline void simulation_assertion_noop() {}

	#define SIMULATION_ASSERT(cond) ((!(cond)) ? fatalAssertionErrorHandler(__FILE__,__LINE__,NULL,NULL) : simulation_assertion_noop())
	#define SIMULATION_ASSERT_MSG(cond,location,msg) ((!(cond)) ? fatalAssertionErrorHandler(__FILE__,__LINE__,location,msg) : simulation_assertion_noop())
	#define SIMULATION_CHECK_VALUE(v) SIMULATION_ASSERT((v) == (v))

#else

	#define SIMULATION_ASSERT(cond)
	#define SIMULATION_ASSERT_MSG(cond,location,msg)
	#define SIMULATION_CHECK_VALUE(v)

#endif

/// Used for compile-time checks.
/// Use dynamic assertion check if static version is not available.
#ifdef BOOST_STATIC_ASSERT
	#define SIMULATION_STATIC_ASSERT BOOST_STATIC_ASSERT
#else
	#define SIMULATION_STATIC_ASSERT SIMULATION_ASSERT
#endif

#endif	// __DEBUG_H
