#ifndef CUBIQUITY_APP_LOGGING_H
#define CUBIQUITY_APP_LOGGING_H

#include <sstream>

enum LogLevel { Debug, Info, Warning, Error };

void setVerbosity(LogLevel verbosity);

void log(LogLevel severity, std::ostringstream& oss);

template <typename T, typename ...Rest>
void log(LogLevel severity, std::ostringstream& oss, T&& t, Rest&&... rest)
{
	oss << std::forward<T>(t);
	log(severity, oss, std::forward<Rest>(rest)...);
}

template <typename ...Args>
void log(LogLevel severity, Args&&... args)
{
	std::ostringstream oss;
	log(severity, oss, std::forward<Args>(args)...);
}

void cubiquityLogHandler(int severity, const char* message);

#endif // CUBIQUITY_APP_LOGGING_H