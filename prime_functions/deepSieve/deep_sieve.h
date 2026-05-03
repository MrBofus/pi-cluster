// deep_sieve.h
//
// Drop-in replacement for a 30-wheel pre-filter. Given an mpz_t candidate n,
// returns true if n has a small factor (i.e., FAILS the sieve and should be
// skipped), false if it survived the sieve and warrants a BPSW test.
//
// Usage:
//     DeepSieve sieve(/*bound=*/100'000'000ULL);   // primes up to 10^8
//     if (sieve.fails(n)) continue;                // composite, skip
//     if (mpz_probab_prime_p(n, 0)) { ... found ... }
//
// Construction is one-time and ~slow (a few seconds for bound=1e8). The
// resulting prime table is ~5M primes for bound=1e8, ~20MB at 4 bytes each.
// For larger bounds use uint64_t storage or a delta encoding.
//
// Thread-safety: fails() is const and reads only the prime table, so it's
// safe to call from multiple threads on the same DeepSieve instance.

#pragma once

#include <cstdint>
#include <vector>
#include <gmp.h>

class DeepSieve {
public:
    // Build a sieve that trial-divides by all primes p with 2 <= p <= bound.
    // Skips p=2 internally (assume caller passes only odd candidates).
    explicit DeepSieve(uint64_t bound);

    // Returns true if n is divisible by any prime in the table (composite,
    // FAIL the sieve). Returns false if no small factor was found (PASS,
    // proceed to BPSW). For odd n only — even n always returns true.
    bool fails(const mpz_t n) const;

    // Number of primes in the table (excluding 2). Useful for diagnostics.
    size_t prime_count() const { return primes_.size(); }

    // Largest prime in the table.
    uint64_t max_prime() const {
        return primes_.empty() ? 0 : primes_.back();
    }

private:
    // Odd primes only, in ascending order. We store as uint32_t to keep
    // the table small; bound up to ~4e9 fits.
    std::vector<uint32_t> primes_;
};