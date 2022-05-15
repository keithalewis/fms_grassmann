#include <cassert>
#include <iostream>
#include "fms_grassmann.h"

using namespace fms::grassmann;

using E = fms::grassmann::element<>;

int test_bits()
{

	return 0;
}

int main()
{
	E::test_bits();
	{
		// empty
		E A;
		assert(0 == A.size());
	}
	{
		// number 2
		E A(2.);
		assert(1 == A.size());
		assert(A.contains(0));
		assert(2 == A[0]);
		A.trim();
		assert(1 == A.size());
		A[0u] = 0;
		A.trim();
		assert(0 == A.size());
	}
	// unit extensor of i
	auto P = [](unsigned i) { return E(1u << i); };

	{
		auto P3 = P(3);
		assert(1 == P3.size());
		auto p3 = *P3.begin();
		assert(0b1000u == p3.first);
		assert(1 == p3.second);
		P3[8] = 0;
		P3.trim();
		assert(0 == P3.size());
	}
	{
		auto P01 = (P(0) | P(1));
		assert(1 == P01.size());
		auto p01 = *P01.begin();
		assert(0b11u == p01.first);
		assert(-1 == p01.second);

		P01 = P(1) | P(0);
		assert(1 == P01.size());
		auto p10 = *P01.begin();
		assert(0b11u == p10.first);
		assert(1 == p10.second);
	}
	{
		E A(0*P(0));
		assert(1 == A.size());
		E A0 = 0. * P(0);
		assert(A0 == A);
		A.trim();
		assert(0 == A.size());
	}
	{
		E A{ 2 * P(0) };
		E B{ 3 * P(1) };
		assert((A + B) == 2 * P(0) + 3 * P(1));
		assert((B + A) == 2 * P(0) + 3 * P(1));
		assert((A - B) == 2 * P(0) - 3 * P(1));
		assert((A - B) == 2 * P(0) + -3 * P(1));
		assert((A | B) == -6 * E(1u + 2u));
		assert((B | A) == 6 * E(2u + 1u));
		assert((A & A) == 2 * A);
		assert(~A == E(~1, -2));
	}
	{
		E P1 = P(0);
		E P2 = P(1);
		E P3 = P(2);
		auto P12 = P1 | P2;
		auto P23 = P2 | P3;
		assert(P2 == (P12 & P23));
		auto P = P1 | P2 | P3 | (P1 + P2 + P3);
		assert(0 == P.size());
	}
	{
		auto P1 = P(0);
		auto P2 = P(1);
		auto Q = 2 * P1 + 3 * P2;
		assert(2 == (Q | P2) / (P1 | P2));
		assert(3 == (P1 | Q) / (P1 | P2));
	}

	return 0;
}