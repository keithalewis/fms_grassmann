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

	// graded sort
	template<std::unsigned_integral I>
	inline constexpr auto less = [](const I i, const I j) noexcept
	{
		auto cmp = std::popcount(i) - std::popcount(j);

		return cmp < 0 ? true : cmp == 0 ? i < j : false;
	};

	// total number of permutations or order i,j
	template<std::unsigned_integral I>
	inline constexpr int perm(const I i, I j) noexcept
	{
		int s = 0;

		for (auto k = std::bit_floor(j); 0 != std::popcount(k); k = std::bit_floor(j)) {
			s += std::popcount(i & (k - 1));
			j &= ~k;
		}

		return s;
	}
	// sign of permutation
	template<std::unsigned_integral I>
	inline constexpr int sign(const I i, I j) noexcept
	{
		return perm(i, j) & 1 ? -1 : 1;
	}

#if 0
	//
	// blade
	//

	// weighted single extensor
	template<std::unsigned_integral I, typename A>
		requires std::is_arithmetic_v<A>
	using blade = std::pair<I, A>;

	// unary minus
	template<std::unsigned_integral I, typename A>
		requires std::is_arithmetic_v<A>
	inline constexpr blade<I, A> operator-(blade<I, A> a) noexcept
	{
		a.second = -a.second;

		return a;
	}

	// scalar multiplication
	template<std::unsigned_integral I, typename A, typename B>
		requires std::is_arithmetic_v<A> && std::is_arithmetic_v<B>
	inline constexpr blade<I, A>& operator*=(blade<I, A>& a, const B b) noexcept
	{
		a.second *= b;

		return a;
	}
	template<std::unsigned_integral I, typename A, typename B>
		requires std::is_arithmetic_v<A>&& std::is_arithmetic_v<B>
	inline constexpr blade<I, A> operator*(blade<I, A> a, const B b) noexcept
	{
		return a *= b;
	}
	template<std::unsigned_integral I, typename A, typename B>
		requires std::is_arithmetic_v<A> && std::is_arithmetic_v<B>
	inline constexpr blade<I, A> operator*(const B b, blade<I, A> a) noexcept
	{
		return a *= b;
	}

	// scalar right division
	template<std::unsigned_integral I, typename A, typename B>
		requires std::is_arithmetic_v<A> && std::is_arithmetic_v<B>
	inline constexpr blade<I, A>& operator/=(blade<I, A>& a, const B b) noexcept
	{
		a.second /= b;

		return a;
	}
	template<std::unsigned_integral I, typename A, typename B>
		requires std::is_arithmetic_v<A>&& std::is_arithmetic_v<B>
	inline constexpr blade<I, A> operator/(blade<I, A> a, const B b) noexcept
	{
		return a /= b;
	}

	// division of compatible extensors
	template<std::unsigned_integral I, typename A>
		requires std::is_arithmetic_v<A>
	inline constexpr A operator/(const blade<I, A>& a, const blade<I, A>& b) noexcept
	{
		return a.first == b.first ? a.second / b.second : std::numeric_limits<A>::quiet_NaN();
	}

	// progressive product
	template<std::unsigned_integral I, typename A>
		requires std::is_arithmetic_v<A>
	inline constexpr blade<I, A>& operator|(blade<I, A>& a, const blade<I, A>& b) noexcept
	{
		if (a.first & b.first) {
			a = blade<I, A>{};
		}
		else {
			a.first |= b.first;
			a.second = sign(a.first, b.first) * a.second * b.second;
		}

		return a;
	}

	// regressive product 
	template<std::unsigned_integral I, typename A>
		requires std::is_arithmetic_v<A>
	inline constexpr blade<I, A>& operator&(blade<I, A>& a, const blade<I, A>& b) noexcept
	{
		if (a.first & b.first) {
			a.first &= b.first;
			a.second *= b.second;
		}
		else {
			a = blade<I, A>{};
		}
	
		return a;
	}
#endif // 0
	//
	// element of Grassmann algebra
	// 
	template<std::unsigned_integral I = unsigned, typename A = double>
		requires std::is_arithmetic_v<A>
	class element
	{
		using map = std::map<I, A, decltype(fms::grassmann::less<I>)>;
		using key_type = typename map::key_type;
		using value_type = typename map::value_type;
		using iterator = typename map::iterator;
		using const_iterator = typename map::const_iterator;
		map e;
	public:
		// real number b
		template<typename B>
			requires std::is_arithmetic_v<B>
		element(const B b = 0)
			: e({ {0u, b} })
		{ }
		element(const value_type& b)
			: e({ b })
		{ }
		element() = default;
		element(const element&) = default;
		element& operator=(const element&) = default;
		~element() = default;

		auto operator<=>(const element& f) const
		{
			return e <=> f.e;
		}

		auto size() const
		{
			return e.size();
		}
		auto contains(const key_type& k) const
		{
			return e.contains(k);
		}
		auto operator[](const key_type& k)
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

		// progressive product
		element& operator|=(const value_type& b)
		{
			std::erase_if(e, [&b](const auto& a) { return a.first & b.first; });

			map e_;
			for (const auto& a : e) {
				e_.insert({ a.first | b.first, sign(a.first, b.first) * a.second * b.second });
			}
			//e = e_; !!! fails !!!
			//e.swap(e_);
	
			return *this;
		}
		element& operator|=(const element& f)
		{
			for (const auto& b : f) {
				operator|=(b);
			}

			return *this;
		}

		// regressive product
		element& operator&=(const value_type& b)
		{
			std::erase_if(e, [&b](const auto& a) { return !(a.first & b.first); });

			map e_;
			for (const auto& a : e) {
				e_.insert({ a.first & b.first, a.second * b.second });
			}
			e = e_;

			return *this;
		}
		element& operator&=(const element& f)
		{
			for (const auto& b : f) {
				operator&=(b);
			}

			return *this;
		}
	};

}

template<std::unsigned_integral I, typename A>
	requires std::is_arithmetic_v<A>
inline std::ostream& operator<<(std::ostream& os, const fms::grassmann::element<I, A>& e)
{
	for (const auto& [i, a] : e) {
		os << std::showpos << a << ":" << std::hex << i << " ";
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

template<std::unsigned_integral I, typename A>
	requires std::is_arithmetic_v<A>
inline fms::grassmann::element<I, A> operator*(fms::grassmann::element<I, A> e, const fms::grassmann::element<I, A>& f)
{
	return e *= f;
}
template<std::unsigned_integral I, typename A, typename B>
	requires std::is_arithmetic_v<A> && std::is_arithmetic_v<B>
inline fms::grassmann::element<I, A> operator*(fms::grassmann::element<I, A> e, const B b)
{
	return e *= b;
}
template<std::unsigned_integral I, typename A, typename B>
	requires std::is_arithmetic_v<A>&& std::is_arithmetic_v<B>
inline fms::grassmann::element<I, A> operator*(const B b, fms::grassmann::element<I, A> e)
{
	return e *= b;
}
template<std::unsigned_integral I, typename A, typename B>
	requires std::is_arithmetic_v<A>&& std::is_arithmetic_v<B>
inline fms::grassmann::element<I, A> operator/(fms::grassmann::element<I, A> e, const B b)
{
	return e *= 1/b;
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
/*
// dagger
template<std::unsigned_integral I, typename A>
inline fms::grassmann::element<N, A> operator~(const fms::grassmann::element<N, A>& e)
{
	fms::grassmann::element<N, A> _e;

	for (const auto& [i, a] : e) {
		_e[~i] = a;
	}

	return _e;
}
*/
