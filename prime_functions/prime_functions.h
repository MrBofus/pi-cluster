#include <iostream>
#include <gmp.h>
#include <gmpxx.h>
#include <ctime>

using namespace std;

/****************************************************************************************************************************/

inline void mpz_rshift(mpz_t rop, mpz_t op1, mp_bitcnt_t op2) {
	mpz_fdiv_q_2exp(rop, op1, op2);
}

/****************************************************************************************************************************/

bool checkLastDigit(mpz_t value){

	// first, check the last digit

	// initialize last digit as a mpz_t
    mpz_t last_digit;
    mpz_init(last_digit);

	// determine last digit by taking the value % 10 and cast it as an unsigned long
	mpz_mod_ui(last_digit, value, 10);
	unsigned long ld = mpz_get_ui(last_digit);

	// check if any of the last digits fail the test
    if (ld == 0 || ld == 2 || ld == 4 
		  || ld == 5 || ld == 6 || ld == 8){
		
		// clear the last digit and return true (meaning test failed)
		mpz_clear(last_digit);
        return true;
    }

	// clear the last digit before proceeding
    mpz_clear(last_digit);


	// next, check if the sum of digits is divisible by 3
	char *str = mpz_get_str(NULL, 10, value);
	int sum = 0;
	for (char *p = str; *p; ++p)
		// add up the digits of the value
		sum += (*p - '0');
	free(str);

	// if the sum is divisible by 3, return and fail the test
	if (sum % 3 == 0) return true;


	// finally, check if the value is of the form 6k ± 1
	unsigned long mod6 = mpz_fdiv_ui(value, 6);
    if (mod6 != 1 && mod6 != 5) {
		// return true, meaning the test failed
        return true;
    }

	// if we made it here, the value is a candidate for primality
    return false;
}


/****************************************************************************************************************************/


char* return_seed(mpz_t value, unsigned int pl){
	char* mpz_seed_str;

	mpz_t mpz_seed, mpz_seed_key, base;
	mpz_inits(mpz_seed, mpz_seed_key, base, NULL);
	
	mpz_set_ui(base, 10);
	mpz_pow_ui(mpz_seed_key, base, pl);
	
	mpz_mod(mpz_seed, value, mpz_seed_key);
	
	mpz_seed_str = mpz_get_str(NULL, 62, mpz_seed);
	
	mpz_clears(mpz_seed, mpz_seed_key, NULL);
	return mpz_seed_str;
}


/****************************************************************************************************************************/


void print_gmp(mpz_t val){
	// printf("val: ");
	mpz_out_str(stdout, 10, val);
	printf("\n");
}


/****************************************************************************************************************************/



bool LucasLehmer(unsigned int c){
	mpz_t candidate, base, s, zero;
	mpz_inits(candidate, base, s, zero, NULL);

	mpz_set_ui(candidate, 0);
	mpz_set_str(base, "2", 10);

	printf("  |---- base\t\t(should be 2): ");
	print_gmp(base);

	mpz_pow_ui(candidate, base, c);
	mpz_sub_ui(candidate, candidate, 1);

	mpz_set_ui(s, 4);

	printf("  |---- initial s\t(should be 4): ");
	print_gmp(s);

	size_t length = mpz_sizeinbase(candidate, 10);
	cout << "  |---- candidate is " << length << " digits\n";
	printf("\n");

	if (checkLastDigit(candidate)){
		cout << "\rvalidating primality...\t100% complete (0/0)";
		cout.flush();
		mpz_clears(candidate, base, s, zero, NULL);
		return false;
	}
	
	for(int i = 1; i < c-1; i++){
		if (i%15==0){
			float percent = 100*i/(c-1);
			// printf("\rvalidating primality...\t%.0f%% complete (%d/%d)", percent, i, (c-1));
			cout << "\rvalidating primality...\t" << percent << "% complete (" << i << "/" << c-1 << ")";
			cout.flush();
		}
        // printf("checking... %d\n", i%15);

		mpz_mul(s, s, s);
		mpz_sub_ui(s, s, 2);

		// printf("s^2 - 2 = ");
		// print_gmp(s);

		mpz_mod(s, s, candidate);

		// printf("s mod x = ");
		// print_gmp(s);
	}
	printf("\n");
	mpz_set_ui(zero, 0);

	// printf("zero: ");
	// print_gmp(zero);

	// printf("maybe zero? ");
	// print_gmp(s);

	int check = mpz_cmp(s, zero);

	if (check == 0){
		mpz_clears(candidate, base, s, zero, NULL);
		return true;
	}

	mpz_clears(candidate, base, s, zero, NULL);
	return false;
}



/****************************************************************************************************************************/




bool trial_composite(mpz_t a, mpz_t d, mpz_t n, mpz_t s){
	mpz_t t1, t2, t3, nprime, one, two, iterator;
	mpz_inits(t1, t2, t3, nprime, one, two, iterator, NULL);

	mpz_set_ui(one, 1);
	mpz_set_ui(two, 2);
	mpz_sub_ui(nprime, n, 1);

	mpz_powm(t1, a, d, n);

	if (mpz_cmp(t1, one) == 0){
		mpz_clears(t1, t2, t3, nprime, one, two, iterator, NULL);
		return false;
	}

	mpz_set_ui(iterator, 0);
	while (true){
		mpz_powm(t2, two, iterator, n);
		mpz_mul(t2, t2, d);

		mpz_powm(t3, a, t2, n);
		if (mpz_cmp(t3, nprime) == 0){
			mpz_clears(t1, t2, t3, nprime,  one, two, iterator, NULL);
			return false;
		}

		mpz_add_ui(iterator, iterator, 1);
        if (mpz_cmp(iterator, s) == 0) { break; }
	}

	mpz_clears(t1, t2, t3, nprime, one, two, iterator, NULL);
	return true;
}



/****************************************************************************************************************************/


bool isPrime(unsigned int value_){
	unsigned int l = value_%10;

	if (l == 2 || l == 4 || l == 5 || l == 6 || l == 8 || l == 0){
		return false;
	}

        mpz_t value, a, d, dprime, s, zero;
        mpz_inits(value, a, d, dprime, s, zero, NULL);

	mpz_set_ui(value, value_);
	mpz_set_ui(s, 0);
	mpz_set_ui(d, value_-1);
	mpz_set_ui(zero, 0);

	while (true){
		mpz_mod_ui(dprime, d, 2);
		if (mpz_cmp(dprime, zero) != 0){
			break;
		}

		mpz_rshift(d, d, 1);
		mpz_add_ui(s, s, 1);
	}

	gmp_randstate_t rstate;
    	gmp_randinit_mt(rstate);

	unsigned long seed;
    	gmp_randseed_ui(rstate, seed);

	for (int i = 0; i < 8; i++){
		mpz_urandomm(a, rstate, value);
		if (trial_composite(a, d, value, s)){
			gmp_randclear(rstate);
			mpz_clears(value, a, d, dprime, s, zero, NULL);
			return false;
		}
	}

	gmp_randclear(rstate);
	mpz_clears(value, a, d, dprime, s, zero, NULL);
	return true;
}



bool isPrime_mpz(mpz_t value){
	// printf("entered test\n");
	mpz_t l, a, d, dprime, s, zero;
	// printf("defined items\n");
	mpz_inits(l, a, d, dprime, s, zero, NULL);
	// printf("initialized to zero\n");

	// printf("\n");
	// printf("computing mod...");

	mpz_mod_ui(l, value, 10);

	// printf("validating primality...\t%.0f%% complete", 0.0);

	if (mpz_cmp_ui(l, 2) == 0 || mpz_cmp_ui(l, 4) == 0 || mpz_cmp_ui(l, 5) == 0 || mpz_cmp_ui(l, 6) == 0 || mpz_cmp_ui(l, 8) == 0 || mpz_cmp_ui(l, 0) == 0){
		printf("\rvalidating primality...\t%.0f%% complete", 100.0);
		printf("\n");
		return false;
	}

	mpz_set_ui(s, 0);
	mpz_sub_ui(d, value, 1);
	mpz_set_ui(zero, 0);

	while (true){
		mpz_mod_ui(dprime, d, 2);
		if (mpz_cmp(dprime, zero) != 0){
			break;
		}

		mpz_rshift(d, d, 1);
		mpz_add_ui(s, s, 1);
	}

	gmp_randstate_t rstate;
    	gmp_randinit_mt(rstate);

	unsigned long seed;
    	gmp_randseed_ui(rstate, seed);

	for (int i = 0; i < 8; i++){
		float percent = 100*i/8;
		printf("validating primality...\t%.0f%% complete", percent);

		mpz_urandomm(a, rstate, value);
		if (trial_composite(a, d, value, s)){
			gmp_randclear(rstate);
			mpz_clears(l, a, d, dprime, s, zero, NULL);
			printf("\n");
			return false;
		}
	}

	gmp_randclear(rstate);
	mpz_clears(l, a, d, dprime, s, zero, NULL);
	printf("\n");
	return true;
}
