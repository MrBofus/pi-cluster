// deep_sieve.cpp

#include "deep_sieve.h"

#include <cmath>
#include <cstring>
#include <stdexcept>

namespace {

// Build the table of odd primes <= bound using a segmented bit sieve.
// One byte per odd number is wasteful but trivially correct; a real
// segmented sieve trims memory dramatically. For bound=1e8 this needs
// ~50 MB transiently, which fits even on a Pi Zero 2 W (512 MB RAM)
// but you'd want to switch to segmented if you push toward 1e9.
std::vector<uint32_t> build_primes_up_to(uint64_t bound) {
    if (bound < 3) return {};
    if (bound > UINT32_MAX) {
        throw std::invalid_argument(
            "DeepSieve: bound exceeds uint32_t storage. "
            "Switch primes_ to uint64_t for larger bounds.");
    }

    // is_composite[i] corresponds to the odd number 2i + 3.
    // i.e. index 0 -> 3, 1 -> 5, 2 -> 7, ...
    const uint64_t n_odds = (bound - 1) / 2;  // count of odds in [3, bound]
    std::vector<uint8_t> is_composite(n_odds + 1, 0);

    const uint64_t limit = static_cast<uint64_t>(std::sqrt((double)bound)) + 1;

    for (uint64_t p = 3; p <= limit; p += 2) {
        const uint64_t idx = (p - 3) / 2;
        if (is_composite[idx]) continue;
        // Mark p*p, p*p+2p, ... as composite.
        // Map odd q to index (q-3)/2.
        for (uint64_t q = p * p; q <= bound; q += 2 * p) {
            is_composite[(q - 3) / 2] = 1;
        }
    }

    std::vector<uint32_t> primes;
    primes.reserve(static_cast<size_t>(bound / std::log((double)bound) * 1.05));
    for (uint64_t i = 0; i <= n_odds; ++i) {
        if (!is_composite[i]) {
            primes.push_back(static_cast<uint32_t>(2 * i + 3));
        }
    }
    return primes;
}

}  // namespace

DeepSieve::DeepSieve(uint64_t bound)
    : primes_(build_primes_up_to(bound)) {}

bool DeepSieve::fails(const mpz_t n) const {
    // mpz_fdiv_ui returns n mod d as an unsigned long. We use this rather
    // than mpz_divisible_ui_p because it's marginally faster on the hot
    // path and we don't need the boolean indirection.
    //
    // For a 50,000-digit mpz_t, each mpz_fdiv_ui is roughly one pass over
    // the limbs (~6,000 limbs at 64-bit). At ~1 ns per limb on an A53 with
    // good cache behavior, that's ~6 us per prime tested. For bound=1e8
    // (~5.7M primes), worst-case all-the-way-through is ~30s — but in
    // practice the vast majority of composites are caught in the first
    // few hundred primes, so average cost is dominated by the small primes.
    for (uint32_t p : primes_) {
        if (mpz_fdiv_ui(n, p) == 0) {
            return true;  // p divides n, fails sieve
        }
    }
    return false;  // survived; send to BPSW
}



// ************************************************************ //

// ----- form-aware sieve ----------------------------------------------
//
// For n = k·2^m + s  (s = ±1), and a small odd prime p (p != 2, p does
// not divide 2^m which is automatic since p is odd):
//
//     n ≡ 0 (mod p)
//   ⟺ k·2^m ≡ -s (mod p)
//   ⟺ k ≡ -s · (2^m)^(-1) (mod p)
//
// We precompute that residue once per (prime, sign).
//
// Modular inverse via Fermat: since p is prime, a^(-1) ≡ a^(p-2) (mod p).
// 2^m mod p computed by fast exponentiation in 64-bit ints — m fits, p
// fits, intermediate products fit in __int128.

static uint64_t mulmod(uint64_t a, uint64_t b, uint64_t m) {
    return (uint64_t)((__uint128_t)a * b % m);
}

static uint64_t powmod(uint64_t base, uint64_t exp, uint64_t m) {
    uint64_t r = 1 % m;
    base %= m;
    while (exp) {
        if (exp & 1) r = mulmod(r, base, m);
        base = mulmod(base, base, m);
        exp >>= 1;
    }
    return r;
}

// Precomputed bad-k residues for one sign.
struct SievePlan {
    vector<uint32_t> primes;
    vector<uint32_t> bad_k_minus;  // residue of k mod p that kills n-1
    vector<uint32_t> bad_k_plus;   // residue of k mod p that kills n+1
};

static SievePlan build_sieve_plan(uint64_t bound, unsigned long m) {
    SievePlan sp;
    sp.primes = build_primes_up_to(bound);
    sp.bad_k_minus.resize(sp.primes.size());
    sp.bad_k_plus.resize(sp.primes.size());

    for (size_t i = 0; i < sp.primes.size(); ++i) {
        uint64_t p = sp.primes[i];
        // 2^m mod p
        uint64_t pow2m = powmod(2, m, p);
        // (2^m)^(-1) mod p  via Fermat
        uint64_t inv = powmod(pow2m, p - 2, p);

        // n - 1 = k·2^m - 1 ≡ 0 ⟺ k ≡  +inv (mod p)
        // n + 1 = k·2^m + 1 ≡ 0 ⟺ k ≡  -inv ≡ p - inv (mod p)
        sp.bad_k_minus[i] = (uint32_t)(inv % p);
        sp.bad_k_plus[i]  = (uint32_t)((p - inv) % p);
    }
    return sp;
}

// Given a candidate k (as a uint64_t — we sweep ranges in 64-bit), test
// whether the sieve rejects k for the given sign. Returns true if k is
// killed by some small prime (composite n), false if k survives.
//
// For ranges, you'd switch this to a strike-out over a bit array — much
// faster than per-k checking. This per-k version is the simple drop-in.
static bool sieve_kills(const SievePlan& sp, uint64_t k, char sign) {
    const auto& bad = (sign == '-') ? sp.bad_k_minus : sp.bad_k_plus;
    for (size_t i = 0; i < sp.primes.size(); ++i) {
        uint32_t p = sp.primes[i];
        if (k % p == bad[i]) {
            // Edge case: if k itself equals p (so n is exactly p·2^m + s),
            // then n is divisible by p — and unless n == p, that's still
            // composite. n == p is impossible at our sizes.
            return true;
        }
    }
    return false;
}

// ----- Proth's theorem -----------------------------------------------
//
// For n = k·2^m + 1 with k odd and k < 2^m:
//   n is prime  ⟺  there exists a with (a/n) = -1 and a^((n-1)/2) ≡ -1 (mod n).
//
// Standard procedure: try small primes a = 3, 5, 7, 11, ... For an actual
// prime n, at least half of all a values are witnesses, so we expect to
// succeed within a couple of tries. For composite n, the test fails for
// every a — we declare composite as soon as one a with (a/n) = -1 fails.
//
// Returns:
//    1  if n is definitely prime
//    0  if n is definitely composite
//   -1  if no a with (a/n) = -1 was found in the tries (shouldn't happen
//       for n with k < 2^m and odd k > 1; flag as inconclusive).

static int proth_is_prime(const mpz_t n) {
    static const unsigned long trial_a[] = {3, 5, 7, 11, 13, 17, 19, 23};
    static const int n_trial = sizeof(trial_a) / sizeof(trial_a[0]);

    mpz_t exp, r, a;
    mpz_inits(exp, r, a, NULL);

    // exp = (n-1)/2
    mpz_sub_ui(exp, n, 1);
    mpz_fdiv_q_2exp(exp, exp, 1);

    int verdict = -1;
    for (int i = 0; i < n_trial; ++i) {
        mpz_set_ui(a, trial_a[i]);
        // Jacobi symbol (a/n). For n prime, (a/n) = -1 means a is a
        // quadratic non-residue, which is what Proth requires.
        int j = mpz_jacobi(a, n);
        if (j != -1) continue;          // wrong type of a, skip

        // r = a^((n-1)/2) mod n
        mpz_powm(r, a, exp, n);

        // Compare r against -1 mod n, i.e. n - 1.
        mpz_t nm1;
        mpz_init(nm1);
        mpz_sub_ui(nm1, n, 1);
        int is_neg_one = (mpz_cmp(r, nm1) == 0);
        mpz_clear(nm1);

        verdict = is_neg_one ? 1 : 0;
        break;  // one good a is enough either way
    }

    mpz_clears(exp, r, a, NULL);
    return verdict;
}