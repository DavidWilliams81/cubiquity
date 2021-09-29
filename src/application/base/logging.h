#ifndef CUBIQUITY_APP_LOGGING_H
#define CUBIQUITY_APP_LOGGING_H

enum LogLevel { Debug, Info, Warning, Error };

void setVerbosity(LogLevel verbosity);
void log(LogLevel severity, const char* format, ...);

void cubiquityLogHandler(int severity, const char* message);

#endif // CUBIQUITY_APP_LOGGING_H