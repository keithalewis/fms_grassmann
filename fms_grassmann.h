// fms_grassmann.h - Grassmann algebra
#pragma once
#include <bit>
#include <bitset>
#include <concepts>
#include <ios>
#include <map>
#include <string>
#include <type_traits>

namespace fms::grassmann {

	template<class I>
	inline auto less = [](I i, I j)
	{
		auto pi = std::popcount(i);
		auto pj = std::popcount(j);

		return pi < pj || pi == pj && i < j;
	};

	template<typename I = uint64_t, typename A = double>
	using element = std::map<I, A, decltype(less<I>)>;

	template<typename I = uint64_t, typename A = double>
	inline element<I, A>& trim(element<I, A>& e)
	{
		std::erase_if(e, [](const auto& i) { return i.second == 0; });

		return e;
	}


}

template<typename I, typename A>
inline std::ostream& operator<<(std::ostream& os, const fms::grassmann::element<I, A>& e)
{
	for (const auto& [i, a] : e) {
		os << std::showpos << a << "`" << std::hex << i << " ";
			//std::bitset<8*sizeof(I)>(i) << ">";
	}

	return os;
}

// dagger
template<typename I, typename A>
inline fms::grassmann::element<I, A> operator~(const fms::grassmann::element<I, A>& e)
{
	fms::grassmann::element<I, A> _e;

	for (const auto& [i, a] : e) {
		_e[~i] = a;
	}

	return _e;
}

template<typename I, typename A>
inline fms::grassmann::element<I, A> operator+(fms::grassmann::element<I, A> e, const fms::grassmann::element<I, A>& f)
{
	for (const auto& [j, b] : f) {
		A a = e[j] += b;
		if (0 == a) {
			e.erase(j);
		}
		
	}

	return e;
}

template<typename I, typename A, typename E>
	requires std::is_arithmetic_v<E>
inline fms::grassmann::element<I, A> operator*(E e, fms::grassmann::element<I, A> f)
{
	for (auto& j : f) {
		j.second *= e;
	}

	return f;
}

template<typename I, typename A, typename F>
	requires std::is_arithmetic_v<F>
inline fms::grassmann::element<I, A> operator*(fms::grassmann::element<I, A> e, F f)
{
	for (auto& i : e) {
		i.second *= f;
	}

	return e;
}

template<typename I, typename A, typename F>
	requires std::is_arithmetic_v<F>
inline fms::grassmann::element<I, A> operator/(fms::grassmann::element<I, A> e, F f)
{
	for (auto& i : e) {
		i.second /= f;
	}

	return e;
}

template<typename I, typename A>
	requires std::is_arithmetic_v<F>
inline A operator/(const fms::grassmann::element<I, A>& e, const fms::grassmann::element<I, A>& f)
{
	A nan = std::numeric_limits<A>::quiet_NaN();

	if (1 != e.size() || 1 != f.size()) {
		return nan;
	}
	const auto& [i, a] = *e.begin();
	const auto& [j, b] = *f.begin();
	
	return i == j ? a / b : nan;
}

template<typename I, typename A>
inline fms::grassmann::element<I, A> operator-(fms::grassmann::element<I, A> e, const fms::grassmann::element<I, A>& f)
{
	for (const auto& [j, b] : f) {
		A a = e[j] -= b;
		if (0 == a) {
			e.erase(j);
		}
	}

	return e;
}

// progressive product
template<typename I, typename A>
inline fms::grassmann::element<I, A> operator|(const fms::grassmann::element<I, A>& e, const fms::grassmann::element<I, A>& f)
{
	fms::grassmann::element<I, A> ef;

	for (const auto& [i, a] : e) {
		for (const auto& [j, b] : f) {
			if (!(i & j)) {
				ef[i | j] = a * b;
			}
		}
	}

	return ef;
}

// regressive product
template<typename I, typename A>
inline fms::grassmann::element<I, A> operator&(const fms::grassmann::element<I, A>& e, const fms::grassmann::element<I, A>& f)
{
	fms::grassmann::element<I, A> ef;

	for (const auto& [i, a] : e) {
		for (const auto& [j, b] : f) {
			if ((i & j)) {
				ef[i & j] = a * b;
			}
		}
	}

	return ef;
}