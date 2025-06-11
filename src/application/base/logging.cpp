#include "logging.h"

#include "base.h"

namespace linalg {
    std::string format_as( vec3 v) {
        return fmt::format("({},{},{})", v.x, v.y, v.z);
    }
    std::string format_as(dvec3 v) {
        return fmt::format("({},{},{})", v.x, v.y, v.z);
    }
    std::string format_as(ivec3 v) {
        return fmt::format("({},{},{})", v.x, v.y, v.z);
    }
    std::string format_as(uvec3 v) {
        return fmt::format("({},{},{})", v.x, v.y, v.z);
    }
    std::string format_as(bvec3 v) {
        return fmt::format("({},{},{})", v.x, v.y, v.z);
    }
}

bool g_color_enabled = true;
void set_color_enabled(bool color_enabled) { g_color_enabled = color_enabled; }

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
        if (g_color_enabled) {
            fmt::print(stderr, ts, "{}\n", fmt::vformat(fmt, args));
        } else {
            fmt::print(stderr,     "{}\n", fmt::vformat(fmt, args));
        }
    }
}

// Conditional log using user-provided color/style.
void vlog_with_newline_and_style_if(bool condition,
                                    log_level severity,
                                    const fmt::text_style& ts,
                                    fmt::string_view fmt,
                                    fmt::format_args args) {
    if(condition && (severity >= g_threshold)) {
        if (g_color_enabled) {
            fmt::print(stderr, ts, "{}\n", fmt::vformat(fmt, args));
        } else {
            fmt::print(stderr,     "{}\n", fmt::vformat(fmt, args));
        }
    }
}
