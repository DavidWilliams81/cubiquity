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

namespace Cubiquity
{
	namespace Internals
	{
		////////////////////////////////////////////////////////////////////////////////////////////
		//   _   _ _   _ _ _ _                                                                    //
		//  | | | | | (_) (_) |                                                                   //
		//  | | | | |_ _| |_| |_ _   _                                                            //
		//  | | | | __| | | | __| | | |                                                           //
		//  | |_| | |_| | | | |_| |_| |                                                           //
		//   \___/ \__|_|_|_|\__|\__, |                                                           //
		//                        __/ |                                                           //
		//                       |___/                                                            //
		////////////////////////////////////////////////////////////////////////////////////////////
		//                                                                                        //
		// Various useful low-level/bitwise operations                                            //
		//                                                                                        //
		////////////////////////////////////////////////////////////////////////////////////////////

		// See https://stackoverflow.com/a/1898487
		bool isAligned(const void *ptr, unsigned int alignment)
		{
			return reinterpret_cast<uintptr_t>(ptr) % alignment == 0;
		}

		bool isPowerOf2(uint32_t uInput)
		{
			if (uInput == 0)
				return false;
			else
				return ((uInput & (uInput - 1)) == 0);
		}

		// Impementation of GLSL's findMSB() for *unsigned* parameter.
		// The behaviour of the signed version is more complex.
		int findMSB(uint32 value)
		{
			int result = -1;
			while (value) {
				result++;
				value >>= 1;
			}
			return result;
		}

		// Simple brute-force approach to computing the log base 2. There are
		// faster but more complex approaches with public domain implementations
		// at https://graphics.stanford.edu/~seander/bithacks.html if needed.
		uint32 logBase2(uint64 value)
		{
			assert(value != 0); // Log of zero is undefined

			uint32 result = 0;
			while (value >>= static_cast<uint64>(1))
			{
				result++;
			}
			return result;
		}

		// See https://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2
		// Alternatively there are intrinsic instructions to count leading zeros.
		uint32 roundUpToPowerOf2(uint32 value)
		{
			value--;
			value |= value >> 1;
			value |= value >> 2;
			value |= value >> 4;
			value |= value >> 8;
			value |= value >> 16;
			value++;
			return value;
		}

		////////////////////////////////////////////////////////////////////////////////////////////
		//  ___  ___                                _   _           _      _____                  //
		//  |  \/  |                               | | | |         | |    |____ |                 //
		//  | .  . |_   _ _ __ _ __ ___  _   _ _ __| |_| | __ _ ___| |__      / /                 //
		//  | |\/| | | | | '__| '_ ` _ \| | | | '__|  _  |/ _` / __| '_ \     \ \                 //
		//  | |  | | |_| | |  | | | | | | |_| | |  | | | | (_| \__ \ | | |.___/ /                 //
		//  \_|  |_/\__,_|_|  |_| |_| |_|\__,_|_|  \_| |_/\__,_|___/_| |_|\____/                  //
		//                                                                                        //
		////////////////////////////////////////////////////////////////////////////////////////////
		//                                                                                        //
		// The hashing functionality below is from MurmurHash3 by Austin Appleby. This function   //
		// was chosen because it's one of the fastest available and its (relatively simple)       //
		// reference implementation is in the public domain. There may be faster alternatives     //
		// which can be evaluated in the future if the need arises.                               //
		//                                                                                        //
		// The original MurmurHash3 reference implementation is avaiable here:                    //
		//                                                                                        //
		//     https://github.com/aappleby/smhasher                                               //
		//                                                                                        //
		////////////////////////////////////////////////////////////////////////////////////////////

		// Supress warnings in third-party (MurmurHash3) code below
#if defined(__GNUC__)
		#pragma GCC diagnostic push
		#pragma GCC diagnostic ignored "-Wimplicit-fallthrough"
#endif // __GNUC__

		// Platform-specific functions and macros

		// Microsoft Visual Studio

		#if defined(_MSC_VER)

		#define FORCE_INLINE	__forceinline

		#include <stdlib.h>

		#define ROTL32(x,y)	_rotl(x,y)
		#define ROTL64(x,y)	_rotl64(x,y)

		#define BIG_CONSTANT(x) (x)

		// Other compilers

		#else	// defined(_MSC_VER)

		#define	FORCE_INLINE inline __attribute__((always_inline))

		inline uint32_t rotl32(uint32_t x, int8_t r)
		{
			return (x << r) | (x >> (32 - r));
		}

		inline uint64_t rotl64(uint64_t x, int8_t r)
		{
			return (x << r) | (x >> (64 - r));
		}

		#define	ROTL32(x,y)	rotl32(x,y)
		#define ROTL64(x,y)	rotl64(x,y)

		#define BIG_CONSTANT(x) (x##LLU)

		#endif // !defined(_MSC_VER)

		//-----------------------------------------------------------------------------
		// Block read - if your platform needs to do endian-swapping or can only
		// handle aligned reads, do the conversion here

		FORCE_INLINE uint32_t getblock32(const uint32_t * p, int i)
		{
			return p[i];
		}

		FORCE_INLINE uint64_t getblock64(const uint64_t * p, int i)
		{
			return p[i];
		}

		//-----------------------------------------------------------------------------
		// Finalization mix - force all bits of a hash block to avalanche

		FORCE_INLINE uint32_t fmix32(uint32_t h)
		{
			h ^= h >> 16;
			h *= 0x85ebca6b;
			h ^= h >> 13;
			h *= 0xc2b2ae35;
			h ^= h >> 16;

			return h;
		}

		//----------

		FORCE_INLINE uint64_t fmix64(uint64_t k)
		{
			k ^= k >> 33;
			k *= BIG_CONSTANT(0xff51afd7ed558ccd);
			k ^= k >> 33;
			k *= BIG_CONSTANT(0xc4ceb9fe1a85ec53);
			k ^= k >> 33;

			return k;
		}

		//-----------------------------------------------------------------------------

		void MurmurHash3_x86_32(const void * key, int len,
			uint32_t seed, void * out)
		{
			const uint8_t * data = (const uint8_t*)key;
			const int nblocks = len / 4;

			uint32_t h1 = seed;

			const uint32_t c1 = 0xcc9e2d51;
			const uint32_t c2 = 0x1b873593;

			//----------
			// body

			const uint32_t * blocks = (const uint32_t *)(data + nblocks * 4);

			for (int i = -nblocks; i; i++)
			{
				uint32_t k1 = getblock32(blocks, i);

				k1 *= c1;
				k1 = ROTL32(k1, 15);
				k1 *= c2;

				h1 ^= k1;
				h1 = ROTL32(h1, 13);
				h1 = h1 * 5 + 0xe6546b64;
			}

			//----------
			// tail

			const uint8_t * tail = (const uint8_t*)(data + nblocks * 4);

			uint32_t k1 = 0;

			switch (len & 3)
			{
			case 3: k1 ^= tail[2] << 16;
			case 2: k1 ^= tail[1] << 8;
			case 1: k1 ^= tail[0];
				k1 *= c1; k1 = ROTL32(k1, 15); k1 *= c2; h1 ^= k1;
			};

			//----------
			// finalization

			h1 ^= len;

			h1 = fmix32(h1);

			*(uint32_t*)out = h1;
		}

		// Restore warnings for our own code.
#if defined(__GNUC__)
		#pragma GCC diagnostic pop
#endif // __GNUC__

		// Slightly nicer C++ wrapper functions (not part of the MurmurHash3 reference impementation)
		uint32 mixBits(uint32 value)
		{
			return fmix32(value); // Force inlined
		}

		uint32 murmurHash3(const void * key, int len, uint32 seed)
		{
			uint32 result;
			MurmurHash3_x86_32(key, len, seed, &result);
			return result;
		}
	}

	void defaultLogHandler(int severity, const char* message)
	{
		switch (severity) {
		case TRACE:
			printf("Trace: "); break;
		case DBG:
			printf("Debug: "); break;
		case INF:
			printf("Info: "); break;
		case WARN:
			printf("Warning: "); break;
		case ERR:
			printf("Error: "); break;
		}
		printf(message);
		printf("\n");
		fflush(stdout);
	}

	LogHandlerPtr gLogHandler = &defaultLogHandler;
	void setLogHandler(LogHandlerPtr logHandler)
	{
		gLogHandler = logHandler;
	}

	void Internals::log(Severity severity, std::ostringstream& oss)
	{
		if (gLogHandler)
		{
			// Forward the print to our log handler.
			(*gLogHandler)(severity, oss.str().c_str());
		}
	}

	ProgressHandlerPtr gProgressHandler = NULL;
	void setProgressHandler(ProgressHandlerPtr progressHandler)
	{
		gProgressHandler = progressHandler;
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
	void Internals::reportProgress(const char* taskDesc, int firstStep, int currentStep, int lastStep)
	{
		assert(firstStep < lastStep && "First step is greater than last step");
		assert(currentStep >= firstStep && "Current step is less than first step");
		assert(currentStep <= lastStep && "Current step is greater than last step");

		if (gProgressHandler) {
			// Make sure that we call the handler for the first and last step
			// (as the user may do something special with these events), as
			// well as an appropriate number of interim steps (but not too many).
			const int updates = 100;
			const int interval = ::std::max((lastStep - firstStep) / updates, 1);
			if (currentStep == firstStep || currentStep == lastStep || currentStep % interval == 0) {
				(*gProgressHandler)(taskDesc, firstStep, currentStep, lastStep);
			}
		}
	}
}
