// fms_grassmann.h - Grassmann algebra
#pragma once
#include <algorithm>
#include <bit>
#include <bitset>
#include <compare>
#include <concepts>
#include <ios>
#include <map>
#include <string>
#include <type_traits>

namespace fms::grassmann {

	// element of Grassmann algebra
	template<std::unsigned_integral I = unsigned, typename A = double>
		requires std::is_arithmetic_v<A>
	class element
	{
		// total number of permutations to order i,j
		static inline constexpr int perm(const I i, I j) noexcept
		{
			int s = 0;

			// check if popcount > size/2???

			for (auto k = std::bit_floor(j); 0 != std::popcount(k); k = std::bit_floor(j)) {
				s += std::popcount(i & (k - 1));
				j &= ~k;
			}

			return s;
		}
		// sign of permutation
		static inline constexpr int sign(const I i, I j) noexcept
		{
			return perm(i, j) & 1 ? -1 : 1;
		}

#ifdef _DEBUG
	public:
		static int test_bits()
		{
			{
				I a(0b1), b(0b10);
				assert(0 == perm(0b1u, a));
				assert(1 == perm(a, b));
				assert(-1 == sign(a, b));
				assert(0 == perm(b, a));
				assert(1 == sign(b, a));
			}

			return 0;
		}
	private:
#endif // _DEBUG

		struct cmp {
			constexpr bool operator()(const I i, const I j) const noexcept
			{
				auto cmp = std::popcount(i) - std::popcount(j);

				return cmp < 0 ? true : cmp == 0 ? i < j : false;
			}
		};

		using map = std::map<I, A, cmp>;
		map e;

		using T = map::value_type;
		static T join(const T& a, const T& b)
		{
			return (a.first & b.first) ? T(0u,0) : T(a.first | b.first, sign(a.first, b.first) * a.second* b.second);
		}
		static T meet(const T& a, const T& b)
		{
			return (a.first & b.first) ? T(a.first & b.first, a.second * b.second) : T(0u, 0);
		}


	public:	
		using key_type = typename map::key_type;
		using value_type = typename map::value_type;
		using mapped_type = typename map::mapped_type;
		using reference = typename map::reference;
		using iterator = typename map::iterator;
		using const_iterator = typename map::const_iterator;

		element() = default;
		element(const element&) = default;
		element& operator=(const element&) = default;
		element(element&&) = default;
		element& operator=(element&&) = default;
		~element() = default;

		constexpr element(I i, A a = 1)
			: e({ std::pair(i, a) })
		{ }
		// real number b
		template<typename B>
			requires std::is_arithmetic_v<B>
		constexpr element(const B b)
			: element(0u, b)
		{ }
		constexpr explicit element(const map& e)
			: e(e)
		{ }

		auto operator<=>(const element& f) const = default;

		auto size() const
		{
			return e.size();
		}
		auto contains(const key_type& k) const
		{
			return e.contains(k);
		}
		auto operator[](const key_type& k) const
		{
			return e[k];
		}
		auto& operator[](const key_type& k)
		{
			return e[k];
		}
		auto begin()
		{
			return e.begin();
		}
		auto begin() const
		{
			return e.begin();
		}
		auto end()
		{
			return e.end();
		}
		auto end() const
		{
			return e.end();
		}

		// remove elements with 0 coefficients
		element& trim()
		{
			std::erase_if(e, [](const auto& b) { return b.second == 0; });

			return *this;
		}

		// unary minus
		element& operator-()
		{
			std::for_each(e.begin(), e.end(), [](auto& b) { b.second = -b.second; });

			return *this;
		}

		// dagger
		element& operator~()
		{
			map _e;
			for (const auto& a : e) {
				_e[~a.first] = sign(a.first, ~a.first) * a.second;
			}
			e = _e;

			return *this;
		}

		template<class B>
			requires std::is_arithmetic_v<B>
		element& operator+=(const B b)
		{
			e[0] += b;
		}
		template<class B>
			requires std::is_arithmetic_v<B>
		element& operator-=(const B b)
		{
			e[0] -= b;
		}

		element& operator+=(const value_type& b)
		{
			e[b.first] += b.second;

			return *this;
		}
		element& operator+=(const element& f)
		{
			for (const auto& b : f) {
				operator+=(b);
			}

			return *this;
		}

		element& operator-=(const value_type& b)
		{
			e[b.first] -= b.second;

			return *this;
		}
		element& operator-=(const element& f)
		{
			for (const auto& b : f) {
				operator-=(b);
			}

			return *this;
		}

		// scalar multiplication
		template<class B>
			requires std::is_arithmetic_v<B>
		element& operator*=(const B b)
		{
			for (auto& a : e) {
				a.second *= b;
			}

			return *this;
		}
		// scalar division
		template<class B>
			requires std::is_arithmetic_v<B>
		element& operator/=(const B b)
		{
			for (auto& a : e) {
				a.second /= b;
			}

			return *this;
		}

		// quotient of compatible extensors
		A operator/(const value_type& b) const
		{
			A q = 0;

			for (const auto& a : *this) {
				if (a.first == b.first) {
					q += a.second / b.second;
				}
				else {
					q = std::numeric_limits<A>::quiet_NaN();
					break;
				}
			}

			return q;
		}
		A operator/(const element& f) const
		{
			if (1 != f.size()) {
				return std::numeric_limits<A>::quiet_NaN();
			}

			return operator/(*f.begin());
		}

		// progressive product
		element& operator|=(const element& f)
		{
			map e_;
			for (const auto& a : e) {
				for (const auto& b : f) {
					auto ab = join(a, b);
					if (ab.first and ab.second) {
						e_[ab.first] += ab.second;
					}
				}
			}
			e = e_;

			return *this;
		}

		// regressive product
		element& operator&=(const element& f)
		{
			map e_;
			for (const auto& a : e) {
				for (const auto& b : f) {
					auto ab = meet(a, b);
					if (ab.first and ab.second) {
						e_[ab.first] += ab.second;
					}
				}
			}
			e = e_;

			return *this;
		}

	};

}

template<std::unsigned_integral I, typename A>
	requires std::is_arithmetic_v<A>
inline std::ostream& operator<<(std::ostream& os, const fms::grassmann::element<I, A>& e)
{
	for (const auto& [i, a] : e) {
		os << std::showpos << a << "*" << std::hex << i << " ";
	}

	return os;
}

template<std::unsigned_integral I, typename A>
	requires std::is_arithmetic_v<A>
inline fms::grassmann::element<I, A> operator+(fms::grassmann::element<I, A> e, const fms::grassmann::element<I, A>& f)
{
	return e += f;
}
template<std::unsigned_integral I, typename A>
	requires std::is_arithmetic_v<A>
inline fms::grassmann::element<I, A> operator-(fms::grassmann::element<I, A> e, const fms::grassmann::element<I, A>& f)
{
	return e -= f;
}

template<std::unsigned_integral I, typename A, typename B>
	requires std::is_arithmetic_v<A> && std::is_arithmetic_v<B>
inline fms::grassmann::element<I, A> operator*(fms::grassmann::element<I, A> e, const B b)
{
	return e *= b;
}
template<std::unsigned_integral I, typename A, typename B>
	requires std::is_arithmetic_v<A> && std::is_arithmetic_v<B>
inline fms::grassmann::element<I, A> operator*(const B b, fms::grassmann::element<I, A> e)
{
	return e *= b;
}

template<std::unsigned_integral I, typename A, typename B>
	requires std::is_arithmetic_v<A> && std::is_arithmetic_v<B>
inline fms::grassmann::element<I, A> operator/(fms::grassmann::element<I, A> e, const B b)
{
	return e *= 1/b;
}
template<std::unsigned_integral I, typename A>
	requires std::is_arithmetic_v<A>
inline A operator/(const fms::grassmann::element<I, A>& e, const fms::grassmann::element<I, A>& f)
{
	return e / f;
}

template<std::unsigned_integral I, typename A>
	requires std::is_arithmetic_v<A>
inline fms::grassmann::element<I, A> operator|(fms::grassmann::element<I, A> e, const fms::grassmann::element<I, A>& f)
{
	return e |= f;
}

template<std::unsigned_integral I, typename A>
	requires std::is_arithmetic_v<A>
inline fms::grassmann::element<I, A> operator&(fms::grassmann::element<I, A> e, const fms::grassmann::element<I, A>& f)
{
	return e &= f;
}

// dagger
template<std::unsigned_integral I, typename A>
inline fms::grassmann::element<I, A> operator~(fms::grassmann::element<I, A> e)
{
	return ~e;
}

