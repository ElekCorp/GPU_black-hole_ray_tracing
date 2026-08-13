#ifndef POWN_H
#define POWN_H

// Integer powers, as the symbolic generators in tools/ emit them: sympy's
// ccode maps Pow with an integer exponent onto pown() (see USER_FUNCS there).
//
// This used to live at the bottom of cuda_ray.h.  It moved here when a second
// generated header (hamiltonian.h) started emitting pown() too, so that the
// canonical-form geodesic code does not have to pull in the whole renderer to
// get one small function.

template <class FP>
inline FP pown_gen(FP const x, int const n);

template <class FP>
inline FP pown(FP const x, int const n)
{
    switch(n)
    {
        case 0: return FP(1);
        case 1: return x;
        case 2: return x*x;
        case 3: return x*x*x;
        case 4: { FP x2=x*x; return x2*x2; }
        case 5: { FP x2=x*x; return x2*x2*x; }
        case 6: { FP x3=x*x*x; return x3*x3; }
    }

    // fallback generic
    return pown_gen(x,n);
}
template <class FP>
inline FP pown_gen(FP const x, int const n)
{
    if (n == 0) return FP(1); // Handle negative exponent safely (including INT_MIN)
    bool neg = (n < 0);
    long long nl = (long long)n;
    unsigned long long exp = neg ? (unsigned long long)(-nl) : (unsigned long long)nl;
    FP result = FP(1);
    FP base = x;
    while (exp)
    {
         if (exp & 1ull) result = result * base;
	 base = base * base;
	 exp >>= 1ull;
    }
    if (neg) return FP(1) / result;
    return result;
}

#endif // POWN_H
