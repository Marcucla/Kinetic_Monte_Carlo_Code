///////////////////////////////////////////////////////////////////////////////
//
//  Copyright 2011, Alexander Stukowski
//  All rights reserved. See README.txt for more information.
//
///////////////////////////////////////////////////////////////////////////////

#include "Context.h"

using namespace std;

/******************************************************************************
* Context constructor.
*****************************************************************************/
Context::Context() : _verbosityLevel(VERBOSITY_NORMAL)
{
	// Send all log and error messages to stdout by default.
	_msgLogger = &cout;
}

/*********************************************************************
* Initializes the context object.
**********************************************************************/
void Context::init(int argc, char* argv[])
{
}

/*********************************************************************
* Aborts the analysis and reports an error before exiting the program.
**********************************************************************/
void Context::error(const char* errorFormatString, ...)
{
	va_list ap;
	va_start(ap, errorFormatString);
	char buffer[2048];
	vsprintf(buffer, errorFormatString, ap);
	va_end(ap);
	errorImpl(buffer);
}

/*********************************************************************
* Aborts the analysis and reports an error before exiting the program.
**********************************************************************/
void Context::errorImpl(const char* errorMessage)
{
	throw runtime_error(errorMessage);
}
