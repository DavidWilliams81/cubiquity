#ifndef CUBIQUITY_APP_HASH_H_
#define CUBIQUITY_APP_HASH_H_

#include "base/types.h"

u64 fnv1a_64(const u8* data, size_t size);

template <typename T>
u64 fnv1a_64(const T& t)
{
	return fnv1a_64(reinterpret_cast<const u8*>(&t), sizeof(T));
}

#endif // CUBIQUITY_APP_HASH_H_