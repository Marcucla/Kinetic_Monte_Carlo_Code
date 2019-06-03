#include "../StandardIncludes.h"
#include "Debug.h"

/********************************************************************************
 * This function aborts the program after printing an error message.
 *******************************************************************************/
void handleFatalAssertionError(const char* file, int line, const char* location, const char* msg)
{
	printf("\n====================================================================================\n");
	printf("ASSERTION FAILURE!\n");
	if(location && msg)
		printf("Source location: %s, line %i, %s\n%s\n", file, line, location, msg);
	else
		printf("Source location: %s, line %i\n", file, line);
	printf("The program stumbled over something unexpected. Please contact the author.\n");
	printf("The bug report should include the program output up to this point and, if possible,\n");
	printf("the input file(s) and command line options used to run the program.\n");
	printf("====================================================================================\n");

	// This invokes the debugger:
	//__builtin_trap();

	// Abort program.
	exit(1);
}

/// Global pointer to the handler function that is called on an assertion error.
void (*fatalAssertionErrorHandler)(const char* file, int line, const char* location, const char* msg) = handleFatalAssertionError;
