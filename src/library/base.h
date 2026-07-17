/***************************************************************************************************
* Cubiquity - A micro-voxel engine for games and other interactive applications                    *
*                                                                                                  *
* Written in 2019 by David Williams                                                                *
*                                                                                                  *
* To the extent possible under law, the author(s) have dedicated all copyright and related and     *
* neighboring rights to this software to the public domain worldwide. This software is distributed *
* without any warranty.                                                                            *
*                                                                                                  *
* You should have received a copy of the CC0 Public Domain Dedication along with this software.    *
* If not, see http://creativecommons.org/publicdomain/zero/1.0/.                                   *
***************************************************************************************************/
#ifndef CUBIQUITY_BASE_H
#define CUBIQUITY_BASE_H

#include <cstdint>
#include <cstring>
#include <type_traits>

// Used by logging to provide parameter checking if supported
#if defined(__GNUC__) || defined(__clang__)
#define PRINTF_LIKE(fmt, first) \
		__attribute__((format(printf, fmt, first)))
#else
#define PRINTF_LIKE(fmt, first)
#endif

namespace Cubiquity
{
	// Integer typedefs (should these be in Interals namespace?)
	typedef int8_t  i8;
	typedef int16_t i16;
	typedef int32_t i32;
	typedef int64_t i64;

	typedef uint8_t  u8;
	typedef uint16_t u16;
	typedef uint32_t u32;
	typedef uint64_t u64;

	typedef unsigned int uint;

	namespace Internals
	{
		// Utility functions
		bool is_power_of_2(u32 uInput);
		u32 log_base_2(u64 value);
		u32 next_power_of_2(u32 value);

		// Hashing and mixing
		u64 fnv1a(const void* data, i64 length, u64 seed = UINT64_C(0xcbf29ce484222325));
		u64 bit_mix(u64 bits);
		u64 hash_combine(u64 hash_1, u64 hash_2);

		// Simple alternative to std::hash which is consistent across compilers
		// and standard libraries. Note this is still byte-order dependent.
		template <typename T>
		u64 hash_value(const T& value)
		{
			static_assert(std::is_trivially_copyable_v<T>);

			T canonical = value;
			if constexpr (std::is_same_v<T, float> || std::is_same_v<T, double>) {
				// Ensure hash is the same for +0.0 vs -0.0.
				if (canonical == T(0)) canonical = T(0);
			} else {
				static_assert(std::has_unique_object_representations_v<T>);
			}

			if constexpr (sizeof(T) <= sizeof(u64)) {
				u64 bits = 0;
				std::memcpy(&bits, &canonical, sizeof(canonical));
				return bit_mix(bits);
			} else {
				return fnv1a(&canonical, sizeof(canonical));
			}
		}

		// Main functions used for logging internally.
		void log_debug(const char* fmt, ...) PRINTF_LIKE(1, 2);
		void log_warning(const char* fmt, ...) PRINTF_LIKE(1, 2);

		// Function used internally for reporting progress on long-running tasks.
		void reportProgress(const char* task_desc, int first_step, int current_step, int last_step);

		// Inherit from this to prevent copying. See https://stackoverflow.com/a/22495199
		class NonCopyable
		{
		public:
			NonCopyable() {}

			//Non copyable
			NonCopyable(NonCopyable const&) = delete;
			NonCopyable& operator=(NonCopyable const&) = delete;

			// But still movable
			NonCopyable(NonCopyable&&) = default;
			NonCopyable& operator=(NonCopyable&&) = default;
		};
	}

	typedef void (*MessageHandlerPtr)(const char* message);
	void set_debug_handler(MessageHandlerPtr debug_handler);
	void set_warning_handler(MessageHandlerPtr warning_handler);

	typedef void (*ProgressHandlerPtr)(const char* task_desc, int first_step, int current_step, int last_step);
	void setProgressHandler(ProgressHandlerPtr progress_handler);
}

#endif // CUBIQUITY_BASE_H
