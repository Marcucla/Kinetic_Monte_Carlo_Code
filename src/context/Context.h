#ifndef __CONTEXT_H
#define __CONTEXT_H

#include "../StandardIncludes.h"
#include "../Settings.h"

#define VERBOSITY_NONE		0
#define VERBOSITY_NORMAL	1
#define VERBOSITY_HIGH		2
#define VERBOSITY_DEBUG		3

/*
 * Global object that provides basic routines for error handling and message logging.
 */
class Context
{
public:

	/******************************** Initialization **************************************/

	/// Default constructor.
	Context();

	/// Initializes the context object.
	void init(int argc, char* argv[]);

	/******************************** Error handling **************************************/

	/// Aborts the program and reports an error to the user before exiting.
	void error(const char* errorFormatString, ...);

	/************************************ Logging ****************************************/

	/// Sets the output stream to which log messages are sent.
	void setMsgLogger(std::ostream& stream) { _msgLogger = &stream; }

	/// Returns the output stream to which log messages are sent.
	std::ostream& msgLogger() const { return *_msgLogger; }

	/// Returns the output stream to which log messages are sent.
	std::ostream& msgLogger(int verbosity) const { return (verbosity <= _verbosityLevel) ? *_msgLogger : const_cast<Context&>(*this)._nullStream; }

	/// Returns the current output verbosity level.
	int verbosityLevel() const { return _verbosityLevel; }

	/// Sets the output verbosity level. Larger values mean more output.
	void setVerbosityLevel(int level) { _verbosityLevel = level; }

protected:

	/// Aborts the analysis and reports an error before exiting the program.
	void errorImpl(const char* errorMessage);

	/// An STL stream class that sends its data to nowhere.
	class NullStream : public std::ostream {
	public:
		NullStream() : std::ostream(NULL) {}
	};

protected:

	/// The output stream to which log messages are sent.
	std::ostream* _msgLogger;

	/// The stream that does nothing.
	NullStream _nullStream;

	/// Output verbosity level.
	int _verbosityLevel;
};

/*
 * Helper base class that stores a reference to the global context object.
 */
class UsingContext
{
public:

	/// Constructor.
	UsingContext(Context& context) : _context(context) {}

	/// Returns a reference to the global context object.
	Context& context() { return _context; }

private:

	/// Reference to the global context object.
	Context& _context;
};

#endif // __CONTEXT_H
