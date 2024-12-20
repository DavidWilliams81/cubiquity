#include "logging.h"

#include "base.h"

log_level g_threshold = log_level::info;

void set_verbosity(log_level threshold) { g_threshold = threshold; }

void vlog(log_level severity,
          fmt::string_view fmt,
          fmt::format_args args) {
    if(severity >= g_threshold) {
        fmt::print(stderr, "{}", fmt::vformat(fmt, args));
    }
}

void vlog_with_newline(log_level severity,
                       fmt::string_view fmt,
                       fmt::format_args args) {
    if(severity >= g_threshold) {
        fmt::print(stderr, "{}\n", fmt::vformat(fmt, args));
    }
}

// Log using user-provided color/style.
void vlog_with_newline_and_style(log_level severity,
                                 const fmt::text_style& ts,
                                 fmt::string_view fmt,
                                 fmt::format_args args) {
    if(severity >= g_threshold) {
        fmt::print(stderr, ts, "{}\n", fmt::vformat(fmt, args));
    }
}
