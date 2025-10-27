#include "hash.h"

// Implementation of FNV-1a hash function
u64 fnv1a_64(const u8* data, size_t size)
{
	u64 hash = UINT64_C(0xcbf29ce484222325);
	for (size_t i = 0; i < size; i++) {
		hash ^= data[i];
		hash *= UINT64_C(0x00000100000001B3);
	}
	return hash;
}