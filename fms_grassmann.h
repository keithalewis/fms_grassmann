// fms_grassmann.h - Grassmann algebra
#pragma once
#include <algorithm>
#include <bit>
#include <bitset>
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

	// total number of permutations
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

	// weighted single extensor
	template<std::unsigned_integral I, typename A>
		requires std::is_arithmetic_v<A>
	using blade = std::pair<I, A>;

	// unary minus
	template<std::unsigned_integral I, typename A>
		requires std::is_arithmetic_v<A>
	inline constexpr blade<I,A> operator-(blade<I,A> a) noexcept
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
		return a * b;
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

	// division of compatible polyhedra
	template<std::unsigned_integral I, typename A>
		requires std::is_arithmetic_v<A>
	inline constexpr A operator/(const blade<I, A>& a, const blade<I, A>& b) noexcept
	{
		return a.first == b.first ? a.second / b.second : std::numeric_limits<A>::quiet_NaN();
	}

	// progressive product
	template<std::unsigned_integral I, typename A>
		requires std::is_arithmetic_v<A>
	inline constexpr blade<I, A> operator|(const blade<I, A>& a, const blade<I, A>& b) noexcept
	{
		return (a.first & b.first) 
			? blade<I, A>{} 
			: blade<I, A>({ a.first | b.first, sign(a.first, b.first) * a.second * b.second });
	}

	// regressive product !!! sign is probably not right
	template<std::unsigned_integral I, typename A>
		requires std::is_arithmetic_v<A>
	inline constexpr blade<I, A> operator&(const blade<I, A>& a, const blade<I, A>& b) noexcept
	{
		return (a.first & b.first)
			? blade<I, A>({ a.first & b.first, sign(a.first, b.first) * a.second * b.second })
			: blade<I, A>{};
	}

	template<std::unsigned_integral I = unsigned, typename A = double>
		requires std::is_arithmetic_v<A>
	using element = std::map<I, A, decltype(fms::grassmann::less<I>)>;
	
	// unary minus
	template<std::unsigned_integral I, typename A = double>
		requires std::is_arithmetic_v<A>
	element<I,A>& operator-(element<I,A>& e)
	{
		std::for_each(e.begin(), e.end(), [](blade<I, A>& b) { b.second = -b.second; });

		return e;
	}
	template<std::unsigned_integral I, typename A = double>
		requires std::is_arithmetic_v<A>
	element<I, A> operator-(element<I, A> e)
	{
		return -e;
	}

	template<std::unsigned_integral I, typename A = double>
			requires std::is_arithmetic_v<A>
	inline element<I,A>& operator+(element<I, A>& a, const blade<I,A>& b)
	{
		a[b.first] += b.second;

		return a;
	}
	template<std::unsigned_integral I, typename A = double>
		requires std::is_arithmetic_v<A>
	inline element<I, A> operator+(const element<I, A>& a, const blade<I, A>& b)
	{
		return a + b;
	}
	template<std::unsigned_integral I, typename A = double>
		requires std::is_arithmetic_v<A>
	inline element<I, A>& operator+(const blade<I, A>& b, element<I, A>& a)
	{
		return a + b;
	}
	template<std::unsigned_integral I, typename A = double>
		requires std::is_arithmetic_v<A>
	inline element<I, A>& operator+(element<I, A>& e, const element<I, A>& f)
	{
		for (const auto& b : f) {
			e += b;
		}

		return e;
	}
	template<std::unsigned_integral I, typename A = double>
		requires std::is_arithmetic_v<A>
	inline element<I, A> operator+(const element<I, A>& e, const element<I, A>& f)
	{
		return e + f;
	}

	template<std::unsigned_integral I, typename A = double>
	inline element<I, A>& trim(element<I, A>& e)
	{
		std::erase_if(e, [](const auto& i) { return i.second == 0; });

		return e;
	}


}

template<std::unsigned_integral I, typename A>
inline std::ostream& operator<<(std::ostream& os, const fms::grassmann::element<I, A>& e)
{
	for (const auto& [i, a] : e) {
		os << std::showpos << a << ":" << std::hex << i << " ";
	}

	return os;
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
