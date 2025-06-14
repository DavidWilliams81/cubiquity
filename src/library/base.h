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

// Features can be activated by uncommenting the defines here, defining them before
// including any Cubiquity files, or by telling the compiler to define them.
//#define CUBIQUITY_USE_AVX

#include <cstdint>
#include <string>

namespace Cubiquity
{
	// Integer typedefs
	typedef int8_t  i8;
	typedef int16_t i16;
	typedef int32_t i32;
	typedef int64_t i64;

	typedef uint8_t  u8;
	typedef uint16_t u16;
	typedef uint32_t u32;
	typedef uint64_t u64;

	typedef unsigned int uint;

	// Having to use shortened names here for now, as otherwise they  conflict
	// with something in instance_list.cpp in the wrapper appliation...
	enum Severity { TRACE, DBG, INF, WARN, ERR };

	//! An integer type used to represent a material.
	/*!
	    An environment in Cubiquity is stored as a 3D grid (see Volume) of material identifiers, with each one specifying the material
		at a given point in space. That is, MaterialId is the type used for each voxel. The material identifiers are simple integer
		values which the user is expected to map so a more complete material description via some application-specific method. For example,
		an application might use [1='rock', 2='wood', 3='glass', ...], or might map striaght to colours [1='red', 2='yellow', 3='blue', ...].

		The MaterialId is a 16-bit value but it is *not* expected that you use the entire 16-bit range. Most applications will only have a
		small number of materials, perhaps up to a few hundred. The compression method used in Cubiquity works by identifying identical
		regions in the scene and storing them only once, and using more materials will generally result in larger file sizes and/or higher
		memory usage.

		Given the above, it is reasonable to ask why an 8-bit data type is not used instead? There are a few potential resons why a 16 bit type is still useful:

		* 256 is too few materials
		* Encoding colours or other sparse values
		* Perhaps spatially large volumes aren't needed.

		Note typedef should not be changed to differnt integer type.
	*/
	typedef u8 MaterialId;

	typedef void (*LogFuncPtr)(const char* message);
	void setLogDebugFunc(LogFuncPtr logDebugFunc);
	void setLogWarningFunc(LogFuncPtr logDebugFunc);

	typedef void (*ProgressHandlerPtr)(const char* taskDesc, int firstStep, int currentStep, int lastStep);
	void setProgressHandler(ProgressHandlerPtr progressHandler);

	namespace Internals
	{
		// Utility functions
		bool isAligned(const void *ptr, unsigned int alignment);
		bool isPowerOf2(u32 uInput);
		int findMSB(u32 value);
		u32 logBase2(u64 value);
		u32 roundUpToPowerOf2(u32 value);

		// Hashing
		u32 mixBits(u32 value);
		u32 murmurHash3(const void * key, int len, u32 seed = 0);

		// Function object to hash a given type with MurmurHash3, allowing it to be used with
		// std::map. Should only be used on simple types as it does not follow pointers, etc.
		template <typename Key>
        struct MurmurHash3
		{
			u32 operator()(const Key& key) const
			{
				return murmurHash3(&key, sizeof(key));
			}
		};

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
 
		// Main functions used for logging internally.
		void log_debug(const std::string& msg);
		void log_warning(const std::string& msg);

		// Function used internally for reporting progress on long-running tasks.
		void reportProgress(const char* taskDesc, int firstStep, int currentStep, int lastStep);
	}
}

#endif // CUBIQUITY_BASE_H
