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
int main()
{
	test_bits();
	{
		E A;
		assert(1 == A.size());
		assert(A.contains(0));
		assert(0 == A[0]);
		A.trim();
		assert(0 == A.size());
	}
	{
		E A(2.);
		assert(1 == A.size());
		assert(A.contains(0));
		assert(2 == A[0]);
		A.trim();
		assert(1 == A.size());
	}
	{
		E A{ {0b1u,0} };
		assert(1 == A.size());
		A.trim();
		assert(0 == A.size());
	}
	{
		E A{ {0b1u, 2} };
		E B{ {0b10u, 3} };
		std::cout << (A + B) << std::endl;
		std::cout << (A - B) << std::endl;
		std::cout << (A | B) << std::endl;
		//std::cout << (A & A) << std::endl;
		//std::cout << (~A) << std::endl;	
	}
	/*
	std::map<I, A, decltype(fms::grassmann::less<I>)>	{
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