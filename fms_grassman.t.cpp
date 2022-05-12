#include <cassert>
#include <iostream>
#include "fms_grassmann.h"

using namespace fms::grassmann;

using E = fms::grassmann::element<>;

int test_bits()
{
	{
		unsigned a(0b1), b(0b10);
		assert(0 == perm(0b1u, a));
		assert(1 == perm(a, b));
		assert(-1 == sign(a, b));
		assert(0 == perm(b, a));
		assert(1 == sign(b, a));
	}

	return 0;
}

int test_blade()
{
	{
		auto a = blade(0b11u, 1);
		a = -a;
		assert(-1 == a.second);
		a *= 2;
		assert(-2 == a.second);
		auto b = a * 3;
		assert(-6 == b.second);
		b = -4 * b;
		assert(24 == b.second);
		b /= 6;
		assert(4 == b.second);
		b = b / 2;
		assert(2 == b.second);
		auto q = a / b;
		assert(q == a.second / b.second);
	}

	return 0;
}

int main()
{
	test_bits();
	test_blade();
	{
		E A{ {0b1u,0} };
		assert(1 == A.size());
		trim(A);
		assert(0 == A.size());
	}
	/*
	{
		E A{ {0b1u, 2} };
		E B{ {0b10u, 3} };
		std::cout << (A + B)/2 << std::endl;
		std::cout << (A - B) << std::endl;
		std::cout << (A | B) << std::endl;
		std::cout << (A & A) << std::endl;
		std::cout << (~A) << std::endl;	
	}
	{
		E P1{ { 0b1, 1.} };
		E P2{ { 0b10, 1.} };
		E P3{ { 0b100, 1.} };
		auto P12 = P1 | P2;
		auto P23 = P2 | P3;
		assert(P2 == (P12 & P23));
		auto P = P1 | P2 | P3 | (P1 + P2 + P3);
		assert(0 == P.size());
	}
	*/

	return 0;
}