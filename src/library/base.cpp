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
#include "base.h"

#include <algorithm>
#include <cassert>
#include <cstdarg>
#include <cstdio>

namespace Cubiquity
{
	namespace Internals
	{
		bool is_power_of_2(u32 uInput)
		{
			if (uInput == 0)
				return false;
			else
				return ((uInput & (uInput - 1)) == 0);
		}

		// See https://graphics.stanford.edu/~seander/bithacks.html#IntegerLogObvious
		// Faster but more complex approaches also available on that page.
		u32 log_base_2(u64 value)
		{
			assert(value != 0); // Log of zero is undefined

			u32 result = 0;
			while (value >>= static_cast<u64>(1)) {
				result++;
			}
			return result;
		}

		// See https://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2
		// Alternatively there are intrinsic instructions to count leading zeros.
		u32 next_power_of_2(u32 value)
		{
			assert(value != 0); // Would silently wrap around to 0

			value--;
			value |= value >> 1;
			value |= value >> 2;
			value |= value >> 4;
			value |= value >> 8;
			value |= value >> 16;
			value++;
			return value;
		}

		// Implementation of FNV-1a hash function (64-bit version)
		u64 fnv1a(const void* data, i64 length, u64 seed)
		{
			for (i64 i = 0; i < length; i++) {
				seed ^= static_cast<const u8*>(data)[i];
				seed *= UINT64_C(0x00000100000001B3);
			}
			return seed;
		}

		// Jon Maiga's 'xmxmx' mixer as used by Boost
		u64 bit_mix(u64 bits)
		{
			bits = ((bits >> 32) ^ bits) * UINT64_C(0x0e9846af9b1a615d);
			bits = ((bits >> 32) ^ bits) * UINT64_C(0x0e9846af9b1a615d);
			return ((bits >> 28) ^ bits);
		}

		// Matches behaviour of Boost hash_combine()
		u64 hash_combine(u64 hash_1, u64 hash_2)
		{
			return bit_mix(hash_1 + 0x9e3779b9 + hash_2);
		}
	}

	MessageHandlerPtr g_debug_handler = nullptr;
	MessageHandlerPtr g_warning_handler = nullptr;

	void set_debug_handler(MessageHandlerPtr debug_handler)
	{
		g_debug_handler = debug_handler;
	}

	void set_warning_handler(MessageHandlerPtr warning_handler)
	{
		g_warning_handler = warning_handler;
	}

	void Internals::log_debug(const char* fmt, ...)
	{
		char buffer[1024];

		va_list args;
		va_start(args, fmt);
		std::vsnprintf(buffer, sizeof(buffer), fmt, args);
		va_end(args);

		if (g_debug_handler) {
			(*g_debug_handler)(buffer);
		}
	}

	void Internals::log_warning(const char* fmt, ...)
	{
		char buffer[1024];

		va_list args;
		va_start(args, fmt);
		std::vsnprintf(buffer, sizeof(buffer), fmt, args);
		va_end(args);

		if (g_warning_handler) {
			(*g_warning_handler)(buffer);
		}
	}

	ProgressHandlerPtr g_progress_handler = NULL;
	void setProgressHandler(ProgressHandlerPtr progress_handler)
	{
		g_progress_handler = progress_handler;
	}

	// When using this function keep in mind that a task with N steps will have N + 1 states.
	// One for no steps completed, and then one for each of the steps. However, it is usually
	// more convenient to call this function only N times(in the loop body), and this is still
	// sufficient to give an impression of the progress while not duplicating the function call.
	// Take care as to whether the loop starts at 0 or 1 and goes to N or N - 1, as this may
	// depend on the type of loop and whether the function call is placed at the start or end.
	//
	// In any case, ensure that the first call for a given string has current step set to first
	// step and that the final call has current step set to last step as the user may depend
	// on these behaviours to drive their progress display.
	void Internals::reportProgress(const char* task_desc, int first_step, int current_step, int last_step)
	{
		assert(first_step < last_step && "First step is greater than last step");
		assert(current_step >= first_step && "Current step is less than first step");
		assert(current_step <= last_step && "Current step is greater than last step");

		if (g_progress_handler) {
			// Make sure that we call the handler for the first and last step
			// (as the user may do something special with these events), as
			// well as an appropriate number of interim steps (but not too many).
			const int updates = 100;
			const int interval = ::std::max((last_step - first_step) / updates, 1);
			if (current_step == first_step || current_step == last_step || current_step % interval == 0) {
				(*g_progress_handler)(task_desc, first_step, current_step, last_step);
			}
		}
	}
}
