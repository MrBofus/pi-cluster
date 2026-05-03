/****************************************************************************************************************************/
// import libraries
//
// libraries imported:
//    prime_functions.h  --  library importing gmp and gmpxx

#include "prime_functions.h"
#include "prime_functions_art.h"


/****************************************************************************************************************************/
// function: trial composite

//  (to use in our 30-wheel filter below)
//      valid residues mod 30 for primes > 5:
//      1, 7, 11, 13, 17, 19, 23, 29
static const bool primeResidue[30] = {
    /* 0 */ false, /* 1 */ true,  /* 2 */ false, /* 3 */ false, /* 4 */ false,
    /* 5 */ false, /* 6 */ false, /* 7 */ true,  /* 8 */ false, /* 9 */ false,
    /*10 */ false, /*11 */ true,  /*12 */ false, /*13 */ true,  /*14 */ false,
    /*15 */ false, /*16 */ false, /*17 */ true,  /*18 */ false, /*19 */ true,
    /*20 */ false, /*21 */ false, /*22 */ false, /*23 */ true,  /*24 */ false,
    /*25 */ false, /*26 */ false, /*27 */ false, /*28 */ false, /*29 */ true
};

bool checkLastDigit_mpz_fast(mpz_t value){

    // first, check the last digit for obvious failed candidates
    unsigned long lastDigit = mpz_fdiv_ui(value, 10);
    if (lastDigit == 0 || lastDigit == 2 || lastDigit == 4 ||
        lastDigit == 5 || lastDigit == 6 || lastDigit == 8)

        // return true, meaning test failed
        return true;
    
    // next, check if the sum of digits is divisible by 3
	char *str = mpz_get_str(NULL, 10, value);
	int sum = 0;
	for (char *p = str; *p; ++p)
		// add up the digits of the value
		sum += (*p - '0');
	free(str);

	// if the sum is divisible by 3, return and fail the test
	if (sum % 3 == 0) return true;

    
    // finally, use the 30-wheel filter
    unsigned long mod30 = mpz_fdiv_ui(value, 30);
    if (!primeResidue[mod30])
        return true;
    
    

    // if we made it here, the value is a candidate for primality
    return false;
}



/****************************************************************************************************************************/
// function: trial composite

bool trial_composite_fast(mpz_t t1, mpz_t t2, mpz_t t3,
                          mpz_t nprime, mpz_t two,
                          mpz_t a, mpz_t d, mpz_t n, unsigned long s){

    // cout << "here" << endl;
    mpz_t iterator;
    mpz_init(iterator);
    mpz_set_ui(iterator, 0);

	if (mpz_cmp_ui(t1, 1) == 0){
		return false;
	}

    while (mpz_cmp_ui(iterator, s) < 0){

        mpz_powm(t2, two, iterator, n);
        mpz_mul(t2, t2, d);
        mpz_powm(t3, a, t2, n);

        if (mpz_cmp(t3, nprime) == 0){
            mpz_clear(iterator);
            return false;
        }

		mpz_add_ui(iterator, iterator, 1);
    }
    mpz_clear(iterator);
    return true;
}

// x:   starts as a^d mod n, and will be modified in-place
// n1:  n - 1
// n:   the candidate
// s:   exponent of 2 in n-1 (i.e., n-1 = d * 2^s)
//
// Returns true  => 'a' is a witness that n is composite.
// Returns false => this base 'a' does NOT prove compositeness.
bool trial_composite_fast_light(mpz_t x, const mpz_t nprime, const mpz_t n, unsigned long s)
{
    // If x == 1 or x == n-1 then this round passes
    if (mpz_cmp_ui(x, 1) == 0 || mpz_cmp(x, nprime) == 0) {
        return false;   // not a witness → n still maybe prime
    }

    // Repeat squaring up to s-1 times
    for (unsigned long r = 1; r < s; ++r) {
        // x = x^2 mod n
        mpz_mul(x, x, x);
        mpz_mod(x, x, n);

        if (mpz_cmp(x, nprime) == 0) {
            return false;   // found n-1 → this base passes
        }

        if (mpz_cmp_ui(x, 1) == 0) {
            // We found 1 *not* at the start: nontrivial square root of 1 mod n
            // ⇒ composite
            return true;
        }
    }

    // Never hit n-1 → 'a' is a witness for compositeness
    return true;
}


/****************************************************************************************************************************/
// function: is prime

void mpz_powm_fast(mpz_t result, mpz_t base, mpz_t exponent, mpz_t mod){
    cout << "here" << endl;
    mpz_t iterator;
    mpz_init(iterator);
    mpz_set_ui(iterator, 2);

    mpz_mul(result, base, base);
    cout << "here" << endl;

    while(mpz_cmp(iterator, exponent) < 0){

        mpz_mul(result, result, base);
        mpz_mod(result, result, mod);

        mpz_add_ui(iterator, iterator, 1);
    }
}


/****************************************************************************************************************************/

bool isPrime_mpz_fast(mpz_t value, gmp_randstate_t rstate, unsigned int ntests){

	mpz_t a, d, dprime;
	mpz_inits(a, d, dprime, NULL);

    mpz_t t1, t2, t3, nprime, two;
    mpz_inits(t1, t2, t3, nprime, two, NULL);

    mpz_set_ui(two, 2);
    mpz_sub_ui(nprime, value, 1);

    // cout << "\t\tvalidating primality...\t0% complete";
	// cout.flush();

	if (checkLastDigit_mpz_fast(value)){
		// cout << "\r\t\tvalidating primality...\t100% complete";
		// cout << endl;
		// cout.flush();

        mpz_clears(a, d, dprime, NULL);
        mpz_clears(t1, t2, t3, nprime, two, NULL);
		return false;
	}

    /*
    mpz_set_ui(s, 0);
    mpz_sub_ui(d, value, 1);

    mpz_mod_ui(dprime, d, 2);
	while (mpz_cmp_ui(dprime, 0) == 0){

		mpz_rshift(d, d, 1);
		mpz_add_ui(s, s, 1);

        mpz_mod_ui(dprime, d, 2);
	}
    */

    mpz_sub_ui(d, value, 1);         // d = n - 1

    unsigned long s = mpz_scan1(d, 0); // s = index of first 1 bit starting from bit 0
    mpz_tdiv_q_2exp(d, d, s);        // d = d / 2^s  (now d is odd)

	for (int i = 0; i < ntests; i++){
        cout << "test " << i << " of " << ntests << "..." << flush;

		// float percent = 100*i/12;
		// cout << "\r\t\tvalidating primality...\t" << percent << "% complete";
		// cout.flush();

        mpz_urandomm(a, rstate, value);
        if (mpz_cmp_ui(a, 0) == 0) { mpz_set_ui(a, 1); }
        mpz_powm(t1, a, d, value);

        if (trial_composite_fast(t1, t2, t3,
                                nprime, two,
                                a, d, value, s)){

            // cout << "\r\t\tvalidating primality...\t100% complete";
            // cout << endl;
            // cout.flush();

            cout << " failed." << endl;

            mpz_clears(a, d, dprime, NULL);
            mpz_clears(t1, t2, t3, nprime, two, NULL);
            return false;
		}
        cout << " passed." << endl;
	}

    // cout << "\r\t\tvalidating primality...\t100% complete";
    // cout << endl;
	// cout.flush();

	mpz_clears(a, d, dprime, NULL);
   	mpz_clears(t1, t2, t3, nprime, two, NULL);

	return true;
}


bool isPrime_mpz_fast_light(mpz_t value, gmp_randstate_t rstate, unsigned int ntests){
    mpz_t d, n1;
    mpz_inits(d, n1, NULL);

    mpz_sub_ui(n1, value, 1);    // n1 = n - 1
    mpz_set(d, n1);              // d  = n - 1

    if (checkLastDigit_mpz_fast(value)){
        mpz_clears(d, n1, NULL);
		return false;
	}

    unsigned long s = mpz_scan1(d, 0);  // number of trailing zeros
    mpz_tdiv_q_2exp(d, d, s);           // d = (n-1) / 2^s (now odd)

    mpz_t a, x, range;
    mpz_inits(a, x, range, NULL);

    // range = n - 3, so a ∈ [0, n-4], then we shift to [2, n-2]
    mpz_sub_ui(range, value, 3);

    for (int i = 0; i < ntests; ++i) {
        // Choose random base a in [2, n-2]
        mpz_urandomm(a, rstate, range); // a ∈ [0, n-4]
        mpz_add_ui(a, a, 2);            // a ∈ [2, n-2]

        // x = a^d mod n
        mpz_powm(x, a, d, value);

        // If 'a' is a Miller–Rabin witness, n is composite
        if (trial_composite_fast_light(x, n1, value, s)) {
            mpz_clears(a, x, range, d, n1, NULL);
            return false;
        }
    }

    mpz_clears(a, x, range, d, n1, NULL);
    return true;
}
