#ifndef CUBIQUITY_APP_POSITION_GENERATOR_H
#define CUBIQUITY_APP_POSITION_GENERATOR_H

#include "base/types.h"

#include <cassert>
#include <numeric>
#include <random>

// Iterates over all positions in the specified inclusive bounds. Iteration
// can be sequential or pseudorandom (if the 'shuffle' flag is set). Shuffling
// is implemented via a modular multiplicative mapping. See:
// https://stackoverflow.com/a/956103 and https://stackoverflow.com/a/79378698
class PositionGenerator
{
public:
	class iterator
	{
	public:
		// Implements a C++ input iterator which can be used with sequential
		// for_each. Supporting *parallel* for_each requires a forward iterator
		// so it can do multiple passes. This can probably be done in the
		// future, but there may be some subtleties with returning references/
		// proxys as we generate values on the fly and don't store them.
		// See https://stackoverflow.com/a/55304579
		using value_type = ivec3;
		using pointer = const ivec3*;
		using reference = const ivec3&;
		using difference_type = std::ptrdiff_t;
		using iterator_category = std::input_iterator_tag;

		// Constructor for iterator
		iterator(const PositionGenerator* pos_gen, u64 counter)
			:m_pos_gen(pos_gen), m_counter(counter) {}

		// I think the default equality comparison is sufficient.
		bool operator==(const iterator& other) const = default;

		// Dereference operator would normally return by refernce, but I believe
		// this is an exception: https://stackoverflow.com/a/15454341
		// Note that no 'structure dereference operator' (->) is provided. This
		// is tricky for generator-style iterators because it must return a
		// pointer-like object but the values are not stored. Returning a light-
		// weight proxy is a workaround but we don't need '->' in our case.
		ivec3 operator*() const {
			assert(m_counter < m_pos_gen->m_voxel_count &&
				   "Cannot dereference end iter");

			// Apply modular multiplicative mapping to shuffle the counter. Add
			// a random offset otherwise the sequence will always start at zero.
			u64 position = (m_counter *
				            m_pos_gen->m_scale +
				            m_pos_gen->m_offset) %
				            m_pos_gen->m_voxel_count;

			// Convert into local (unsigned) 3D position
			u64 x = (position % m_pos_gen->m_dims.x);
			u64 y = (position / m_pos_gen->m_dims.x) % m_pos_gen->m_dims.y;
			u64 z = (position / m_pos_gen->m_dims.x) / m_pos_gen->m_dims.y;

			// Adjust for specified lower bound
			return m_pos_gen->m_lower + ivec3(x, y, z);
		}

		// Pre-increament
		iterator& operator++() {
			assert(m_counter < m_pos_gen->m_voxel_count &&
				   "Cannot increment end iter");
			++m_counter;
			return *this;
		}

		// Post-increament
		iterator operator++(int) {
			iterator temp = *this;
			++(*this);
			return temp;
		}

	private:

		const PositionGenerator* m_pos_gen = nullptr;
		u64 m_counter = 0;

	}; // End nested iterator

	// Constructor for PositionGenerator
	PositionGenerator(const ivec3& lower, const ivec3& upper,
		              bool shuffle = false, u64 seed = 0)
		:m_lower(lower), m_dims((upper - lower) + ivec3(1, 1, 1)) {

		// Promote components before multiplication to avoid overflow.
		m_voxel_count = static_cast<u64>(m_dims.x) *
			static_cast<u64>(m_dims.y) *
			static_cast<u64>(m_dims.z);

		if (shuffle) {
			// Mersenne Twister gives required 64-bit values and is consistent
			// across all platforms. I'm not so sure about the seed (it uses
			// seed_seq internally) but it seems to match between GCC and VC++.
			std::mt19937_64 rng(seed);

			// Pick candidate scale values at random until we find one which is
			// coprime with the voxel count. We also map them to the middle part
			// of the valid range as this gives a lower-discrepancy sequence.
			do {
				m_scale = rng() % m_voxel_count;
				m_scale = (m_scale / 2) + (m_voxel_count / 4);
			} while (std::gcd(m_scale, m_voxel_count) != 1); // Not coprime

			m_offset = rng(); // Offset can be anything (just shifts sequence).
		}
	}

	// Iterator access
	iterator begin() const { return iterator(this, 0); }
	iterator end()   const { return iterator(this, m_voxel_count); }

private:
	// Limits
	ivec3 m_lower;
	ivec3 m_dims;
	u64   m_voxel_count = 0;

	// Sequence parameters
	u64 m_scale  = 1; // 1 for sequential order, large coprime for random order.
	u64 m_offset = 0; // Sets starting point of the sequence.	
};

static_assert(std::input_iterator<PositionGenerator::iterator>);

#endif // CUBIQUITY_APP_POSITION_GENERATOR_H
