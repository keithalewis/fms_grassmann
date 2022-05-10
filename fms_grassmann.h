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

	template<std::size_t N>
	inline bool operator<(const std::bitset<N>& i, const std::bitset<N>& j)
	{
		static_assert(N <= 8 * sizeof(unsigned long long));

		return i.to_ullong() < j.to_ullong();
	};
	template<std::size_t N>
	inline auto less = [](const std::bitset<N>& i, const std::bitset<N>& j) { return i < j; };

	// total number of permutations
	template<std::size_t N>
	inline int perm(const std::bitset<N>& i, std::bitset<N> j)
	{
		int s = 0;

		unsigned long long ui = i.to_ullong(), uj = j.to_ullong();

		for (auto k = std::bit_floor(uj); 0 != std::popcount(k); k = std::bit_floor(uj)) {
			s += std::popcount(ui & (k - 1));
			uj &= ~k;
		}

		return s;
	}
	template<std::size_t N>
	inline int sign(const std::bitset<N>& i, std::bitset<N> j)
	{
		return perm(i, j) & 1 ? -1 : 1;
	}

	// weighted single polyhedron
	template<std::size_t N, typename A>
		requires std::is_arithmetic_v<A>
	using blade = std::pair<std::bitset<N>, A>;

	// unary minus
	template<std::size_t N, typename A>
		requires std::is_arithmetic_v<A>
	inline blade<N,A>& operator-(blade<N,A>& b)
	{
		b.second = -b.second;

		return b;
	}

	// scalar multiplication
	template<std::size_t N, typename A, typename B>
		requires std::is_arithmetic_v<A> && std::is_arithmetic_v<B>
	inline blade<N, A> operator*(const blade<N, A>& a, const B& b)
	{
		return blade<N, A>({ a.first, a.second * b });
	}
	template<std::size_t N, typename A, typename B>
		requires std::is_arithmetic_v<A> && std::is_arithmetic_v<B>
	inline blade<N, A> operator*(const B& b, const blade<N, A>& a)
	{
		return a * b;
	}

	// scalar right division
	template<std::size_t N, typename A, typename B>
		requires std::is_arithmetic_v<A> && std::is_arithmetic_v<B>
	inline blade<N, A> operator/(const blade<N, A>& a, const B& b)
	{
		return blade<N, A>({ a.first, a.second/b });
	}
	// division of compatible polyhedra
	template<std::size_t N, typename A>
		requires std::is_arithmetic_v<A>
	inline A operator/(const blade<N, A>& a, const blade<N, A>& b)
	{
		return a.first == b.first ? a.second / b.second : std::numeric_limits<A>::quiet_NaN();
	}

	// progressive product
	template<std::size_t N, typename A>
		requires std::is_arithmetic_v<A>
	inline blade<N, A> operator|(const blade<N, A>& a, const blade<N, A>& b)
	{
		return (a.first & b.first) 
			? blade<N, A>{} 
			: blade<N, A>({ a.first | b.first, sign(a.first, b.first) * a.second * b.second });
	}

	// regressive product !!! sign is probably not right
	template<std::size_t N, typename A>
		requires std::is_arithmetic_v<A>
	inline blade<N, A> operator&(const blade<N, A>& a, const blade<N, A>& b)
	{
		return (a.first & b.first)
			? blade<N, A>({ a.first & b.first, sign(a.first, b.first) * a.second * b.second })
			: blade<N, A>{};
	}

	template<std::size_t N, typename A = double>
		requires std::is_arithmetic_v<A>
	using element = std::map<std::bitset<N>, A, decltype(less<N>)>;
	
	// unary minus
	template<std::size_t N, typename A = double>
		requires std::is_arithmetic_v<A>
	element<N,A>& operator-(element<N,A>& e)
	{
		std::for_each(e.begin(), e.end(), [](blade<N, A>& b) { b.second = -b.second; });

		return e;
	}
	template<std::size_t N, typename A = double>
		requires std::is_arithmetic_v<A>
	element<N, A> operator-(element<N, A> e)
	{
		return -e;
	}

	template<std::size_t N, typename A = double>
			requires std::is_arithmetic_v<A>
	inline element<N,A>& operator+(element<N, A>& a, const blade<N,A>& b)
	{
		a[b.first] += b.second;

		return a;
	}
	template<std::size_t N, typename A = double>
		requires std::is_arithmetic_v<A>
	inline element<N, A> operator+(const element<N, A>& a, const blade<N, A>& b)
	{
		return a + b;
	}
	template<std::size_t N, typename A = double>
		requires std::is_arithmetic_v<A>
	inline element<N, A>& operator+(const blade<N, A>& b, element<N, A>& a)
	{
		return a + b;
	}
	template<std::size_t N, typename A = double>
		requires std::is_arithmetic_v<A>
	inline element<N, A>& operator+(element<N, A>& e, const element<N, A>& f)
	{
		for (const auto& b : f) {
			e += b;
		}

		return e;
	}
	template<std::size_t N, typename A = double>
		requires std::is_arithmetic_v<A>
	inline element<N, A> operator+(const element<N, A>& e, const element<N, A>& f)
	{
		return e + f;
	}

	template<std::size_t N, typename A = double>
	inline element<N, A>& trim(element<N, A>& e)
	{
		std::erase_if(e, [](const auto& i) { return i.second == 0; });

		return e;
	}


}

template<std::size_t N, typename A>
inline std::ostream& operator<<(std::ostream& os, const fms::grassmann::element<N, A>& e)
{
	for (const auto& [i, a] : e) {
		os << std::showpos << a << ":" << i << " ";
	}

	return os;
}
/*
// dagger
template<std::size_t N, typename A>
inline fms::grassmann::element<N, A> operator~(const fms::grassmann::element<N, A>& e)
{
	fms::grassmann::element<N, A> _e;

	for (const auto& [i, a] : e) {
		_e[~i] = a;
	}

	return _e;
}
*/
