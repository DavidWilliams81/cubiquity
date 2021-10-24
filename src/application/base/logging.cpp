#include "logging.h"

#include "base.h"
#include "utility.h"

#include <stdarg.h>
#include <stdio.h>

LogLevel g_verbosity;

void setVerbosity(LogLevel verbosity) { g_verbosity = verbosity; }

void log(LogLevel severity, std::ostringstream& oss)
{
    if (severity >= g_verbosity) {
        FILE* stream = severity >= LogLevel::Warning ? stderr : stdout;

        // Set the colour for the text. Note that I'm not quite sure how widely supported this
        // is, therefore I avoid manipulating the colour codes for 'Info' (the most common case)
        // so that the output should be largely readable if colors are not supported.
        // See https://stackoverflow.com/a/2616912
        switch (severity) {
        case LogLevel::Debug:
            fprintf(stream, "\033[36m"); // Cyan
            break;
        case LogLevel::Info:
            // Nothing to do 
            break;
        case LogLevel::Warning:
            fprintf(stream, "\033[33m"); // Yellow
            break;
        case LogLevel::Error:
            fprintf(stream, "\033[1;31m"); // Red (in bold)
            break;
        }

        // Print the text
        fprintf(stream, oss.str().c_str());

        // Reset the colour manipulators if we changed them
        if (severity != LogLevel::Info) {
            fprintf(stream, "\033[0m");
        }

        // Append a newline and flush
        fprintf(stream, "\n");
        fflush(stream);
    }
}

void cubiquityLogHandler(int severity, const char* message)
{
    // This cubiquity application does not offer as many log level as the core 
    // library (for simplicity). Generally users of the library are not interested 
    // in messages from the library. Therefore the trace, debug, and info severities
    // are all mapped to 'debug' here, while warning and error are propergated.
    switch(severity) {
    case Cubiquity::ERR:
        log(Error, message);
        break;
    case Cubiquity::WARN:
        log(Warning, message);
        break;
    default:
        log(Debug, message);
    }
}
