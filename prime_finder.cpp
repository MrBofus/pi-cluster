// compile using this:
// g++ -o prime_finder prime_finder.cpp -lgmp
// or, use this:
/*
first, optimize for specific hardware:
sudo apt install m4

wget https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz
tar xf gmp-6.3.0.tar.xz
cd gmp-6.3.0

./configure \
  --prefix=$HOME/local/gmp-fast \
  CFLAGS="-O3 -march=native -mtune=native" \
  CXXFLAGS="-O3 -march=native -mtune=native"

make -j$(nproc)
make check
make install


compile with:
g++ -O3 -march=native -mtune=native -DNDEBUG \
    prime_finder.cpp \
    -I$HOME/local/gmp-fast/include \
    -L$HOME/local/gmp-fast/lib \
    -Wl,-rpath,$HOME/local/gmp-fast/lib \
    -lgmp -lgmpxx \
    -o prime_finder
*/

// import libraries
#include <iomanip>
#include <fstream>
#include <sstream>
#include <math.h>
#include <cstdlib>
#include <cmath>
#include <cstdlib>

#include "prime_functions/prime_functions_fast.h"
#include "prime_functions/deepSieve/deep_sieve.cpp"

// initialize mode
int mode = 1;

// initialize file to write primes to
FILE* primefile;

// function to write primes to file
void write_primes_to_file(ostringstream& buffer, const char* filename){
	ofstream outfile(filename, ios::app);
	outfile << buffer.str();
	outfile.close();
	buffer.str("");
	buffer.clear();
}

// main loop
int main(void) {
	// buffer to hold primes before writing them to a file
	const size_t BUFFER_SIZE = 1;
	ostringstream prime_buffer;

	// file name to store prime candidates
	char* str_lower = NULL;
	const char* file = ".out/prime_candidates.txt";

	// initialize the sieve before anything starts
	// our sieve will check for 10^8 primes before sending
	// it to the primality test
	cout << "building sieve..." << endl;
	// static const DeepSieve g_sieve(100'000'000ULL);
	static const DeepSieve g_sieve(1000000);
	cout << "sieve built" << endl;

	// get from user the desired number of digits
	// to test candidates for primality
	unsigned int val = 0;
	cout << "number of digits: ";
	cin >> val;

	// ---- compute m so 2^m has ~val digits ----
	// log2(10) ≈ 3.3219280948873623
	const unsigned long m =
		(unsigned long) floor((val - 1) * 3.3219280948873623);

	cout << "m = " << m << "  (2^m has ~" << val << " digits)\n";

	// random number for k
	gmp_randstate_t rstate;
	gmp_randinit_mt(rstate);

	srand((unsigned int)time(0));
	unsigned long seed = rand()%100000000000;

	gmp_randseed_ui(rstate, seed);

	// initialize values before starting
	mpz_t n, n_minus, n_plus, step, k_start;
	mpz_inits(n, n_minus, n_plus, step, k_start, NULL);

	// step = 2^(m+1): what we add to n each time k advances by 2
	mpz_set_ui(step, 1);
	mpz_mul_2exp(step, step, m + 1);

	// Random odd starting k. 64 random bits is plenty of starting
	// entropy for a single node; on the cluster, the coordinator hands
	// each worker a disjoint k-range instead.
	mpz_urandomb(k_start, rstate, 64);
	mpz_setbit(k_start, 0);                 // force odd

	// n = k_start · 2^m
	mpz_mul_2exp(n, k_start, m);
	// n - 1
	mpz_sub_ui(n_minus, n, 1);

	//`````````````````````````````````````````````````````//

	// initialize counter and prime number counter
	unsigned int counter = 0;
	unsigned int pcounter = 0;
	unsigned int pcounter_check = 0;

	// begin random search
	while (true) {

		// check value against sieve to see if it's worth checking
		while(g_sieve.fails(n_minus)) {
			// iterate counter
			counter++;
			mpz_add_ui(n_minus, n_minus, 1);
		}

		cout << "\r  analyzing candidate #" << counter;
		cout << " --- total of " << pcounter << " found" << flush;

		if (mpz_probab_prime_p(n_minus, 0)){

			
			// if it is prime, increment the prime counter
			pcounter++;
			pcounter_check++;

			// write prime to text file
			str_lower = mpz_get_str(NULL, 10, n_minus);
			prime_buffer << "\n" << str_lower << "\n";
			if (pcounter_check >= BUFFER_SIZE){
				write_primes_to_file(prime_buffer, file);
				pcounter_check = 0;
			}
		}
		else {}

		// increment candidate by one
		mpz_add_ui(n_minus, n_minus, 1);
	}

	// (unreachable in practice...)
	gmp_randclear(rstate);
	mpz_clears(n, n_minus, n_plus, step, k_start, NULL);
	return 0;
}
