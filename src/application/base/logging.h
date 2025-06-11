#ifndef CUBIQUITY_APP_LOGGING_H
#define CUBIQUITY_APP_LOGGING_H

#define FMT_UNICODE 0 // Otherwise Visual C++ needs more linker flags
#define FMT_HEADER_ONLY
#include "fmt/color.h"
#include "fmt/format.h"
#include "fmt/std.h"

// Support for printing maths types
#include "types.h"

// According to the libfmt docs it should be possible to wrap a type's ostream
// operator<< but I couldn't get it to compile on GCC (it worked on VS2022).
// Providing these format_as() functions id the next easiest option.
namespace linalg {
    std::string format_as( vec3 v);
    std::string format_as(dvec3 v);
    std::string format_as(ivec3 v);
    std::string format_as(uvec3 v);
    std::string format_as(bvec3 v);
}

void set_color_enabled(bool color_enabled);

enum class log_level { debug, info, note, warning, error };
void set_verbosity(log_level threshold);

// Override print function to prevent accidental writing to stdout,
// which is reserved for real program output (e.g. volume data, not text).
// Always prints (no verbosity control) and no newline appended.
template <typename... T>
void print(fmt::format_string<T...> fmt, T&&... args) {
  fmt::vprint(stderr, fmt, fmt::make_format_args(args...));
}

// Logging functions always append a newline, and again always write to stderr.
// Based on model here: https://fmt.dev/11.0/api/#type-erasure
void vlog(                          log_level severity,
                                    fmt::string_view fmt,
                                    fmt::format_args args);
void vlog_with_newline(             log_level severity,
                                    fmt::string_view fmt,
                                    fmt::format_args args);
void vlog_with_newline_and_style(   log_level severity,
                                    const fmt::text_style& ts,
                                    fmt::string_view fmt,
                                    fmt::format_args args);
void vlog_with_newline_and_style_if(bool condition,
                                    log_level severity,
                                    const fmt::text_style& ts,
                                    fmt::string_view fmt,
                                    fmt::format_args args);

// Debug
template <typename... T>
void log_debug(fmt::format_string<T...> fmt, T&&... args) {
  vlog_with_newline_and_style(log_level::debug,
                              fg(fmt::terminal_color::bright_black),
                              fmt, fmt::make_format_args(args...));
}

template <typename... T>
void log_debug_if(bool condition,
                  fmt::format_string<T...> fmt, T&&... args) {
  vlog_with_newline_and_style_if(condition, log_level::debug,
                                 fg(fmt::terminal_color::bright_black),
                                 fmt, fmt::make_format_args(args...));
}

// Info uses default color (don't assume black/white, in case user has changed their terminal).
template <typename... T>
void log_info(fmt::format_string<T...> fmt, T&&... args) {
  vlog_with_newline(log_level::info, fmt, fmt::make_format_args(args...));
}

// Info without newline
template <typename... T>
void log_info_no_newline(fmt::format_string<T...> fmt, T&&... args) {
  vlog(log_level::info, fmt, fmt::make_format_args(args...));
}

// Note
template <typename... T>
void log_note(fmt::format_string<T...> fmt, T&&... args) {
  vlog_with_newline_and_style(log_level::note,
                              fg(fmt::terminal_color::cyan),
                              fmt, fmt::make_format_args(args...));
}

template <typename... T>
void log_note_if(bool condition,
                 fmt::format_string<T...> fmt, T&&... args) {
  vlog_with_newline_and_style_if(condition, log_level::note,
                                 fg(fmt::terminal_color::cyan),
                                 fmt, fmt::make_format_args(args...));
}

// Warning
template <typename... T>
void log_warning(fmt::format_string<T...> fmt, T&&... args) {
  vlog_with_newline_and_style(log_level::warning,
                              fg(fmt::terminal_color::bright_yellow),
                              fmt, fmt::make_format_args(args...));
}

template <typename... T>
void log_warning_if(bool condition,
                    fmt::format_string<T...> fmt, T&&... args) {
  vlog_with_newline_and_style_if(condition, log_level::warning,
                                 fg(fmt::terminal_color::bright_yellow),
                                 fmt, fmt::make_format_args(args...));
}

// Error
template <typename... T>
void log_error(fmt::format_string<T...> fmt, T&&... args) {
  vlog_with_newline_and_style(log_level::error,
                              fg(fmt::terminal_color::red),
                              fmt, fmt::make_format_args(args...));
}

template <typename... T>
void log_error_if(bool condition,
                  fmt::format_string<T...> fmt, T&&... args) {
  vlog_with_newline_and_style_if(condition, log_level::error,
                                 fg(fmt::terminal_color::red),
                                 fmt, fmt::make_format_args(args...));
}

#endif // CUBIQUITY_APP_LOGGING_H
