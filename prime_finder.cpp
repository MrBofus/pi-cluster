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

// initialize buffer to communicate with
FILE* combuffer;

// function to write primes to file
void write_to_file(ostringstream& buffer, const char* filename){
	ofstream outfile_(filename, ios::app);
	outfile_ << buffer.str();
	outfile_.close();
	buffer.str("");
	buffer.clear();
}


void write_counter(ostringstream& buffer, const char* filename){
	ofstream outfile_(filename, ios::trunc);
	outfile_ << buffer.str();
	outfile_.close();
	buffer.str("");
	buffer.clear();
}

// main loop
int main(void) {
	unsigned int chunksize;
	unsigned int startingint;

	// get chunk size from tasking service
	if (!(cin >> chunksize)) { cerr << "bad chunksize" << endl; exit(1); }

	// get seed to sync all worker nodes
	string seed_str;
	if (!(cin >> seed_str)) { cerr << "bad seed" << endl; exit(1); }

	mpz_t seed_mpz;
	mpz_init(seed_mpz);
	if (mpz_set_str(seed_mpz, seed_str.c_str(), 10) != 0) {
		cerr << "could not parse seed: " << seed_str << endl;
	}

	//````````````````````````````````````````````````````````````````````````````````````````````````````````//
	// initialize the random state
	gmp_randstate_t rng_state;
	gmp_randinit_mt(rng_state);
	gmp_randseed(rng_state, seed_mpz);
	mpz_clear(seed_mpz);

	//```````````````````````````````````````````````````````````````````````````````````````````````````````//
	// initialize the number
	const unsigned long N_DIGITS = 1000;

	mpz_t base, lower, upper, range;
	mpz_inits(base, lower, upper, range);

	mpz_ui_pow_ui(lower, 10, N_DIGITS - 1); // 10^49,999
	mpz_ui_pow_ui(upper, 10, N_DIGITS); 	// 10^50,000

	mpz_sub(range, upper, lower); 		// 9 * 10^49,999

	mpz_urandomm(base, rng_state, range); 	// [0, 9 * 10^49,999)
	mpz_add(base, base, lower); 		// [10^49,999, 10^50,000)

	mpz_clears(lower, upper, range, NULL); 	// clear the intermediate variables

	//```````````````````````````````````````````````````````````````````````````````````````````````````````//
	// checkouts on the pc side
	char* base_str = mpz_get_str(NULL, 10, base);
	size_t blen = strlen(base_str);
	char head[13] = {0}, tail [13] = {0};

	strncpy(head, base_str, 12);
	strncpy(tail, base_str + blen - 12, 12);

	cout << "ready: base " << head << "..." << tail
		<< " (" << blen << " digits)" << endl;
	cout.flush();
	free(base_str);

	//````````````````````````````````````````````````````````````````````````````````````````````````````````//
	// buffer to hold primes before writing them to a file
	const size_t COUNTER_BUFFER_SIZE = 500;
	ostringstream counter_buffer;

	const size_t PRIME_BUFFER_SIZE = 1;
	ostringstream prime_buffer;

	// file name to store prime candidates
	char* str_lower = NULL;
	const char* pfile = ".out/prime_candidates.txt";
	const char* cfile = ".comms/combuffer.txt";

	// initialize the sieve before anything starts
	// our sieve will check for 10^8 primes before sending
	// it to the primality test
	// static const DeepSieve g_sieve(100'000'000ULL);
	// static const DeepSieve g_sieve(1000000);


	// initialize values before starting
	mpz_t n, offset;
	mpz_inits(n, offset, NULL);


	//````````````````````````````````````````````````````````````````````````````````````````````````````````//

	// initialize counter and prime number counter
	// unsigned int counter = 0;
	// unsigned int counter_check = 0;

	unsigned int pcounter = 0;
	unsigned int pcounter_check = 0;

	// begin random search
	while (true) {
		string offset_str;
		if (!(cin >> offset_str)) { exit(0); }
		if (offset_str == "12121212") { exit(0); }

		if (mpz_set_str(offset, offset_str.c_str(), 10) != 0) {
			cerr << "bad offset: " << offset_str << endl;
		}

		// n = base + offset
		mpz_add(n, base, offset);

		unsigned int counter = 0;
		unsigned int counter_check = 0;


		// reset the comms buffer
		write_counter(counter_buffer, cfile);

		// start clock to time the batch
		clock_t t_to_start = clock();
		while ( counter <= chunksize ) {

			// check value against sieve to see if it's worth checking
			while( checkLastDigit_mpz_fast(n) ) {
				// iterate counter
				counter++;
				counter_check++;

				if (counter >= chunksize) { break; }

				mpz_add_ui(n, n, 1);
			}

			if (mpz_probab_prime_p(n, 0)) {
				// if it exited the g_sieve, the counter didn't increment,
				// so let's increment it here
				counter++;
				counter_check++;

				// if it is prime, increment the prime counter
				pcounter++;
				pcounter_check++;

				// write prime to text file
				str_lower = mpz_get_str(NULL, 10, n);
				prime_buffer << "\n" << str_lower << "\n";
				free(str_lower);
				if (pcounter_check >= PRIME_BUFFER_SIZE) {
					write_to_file(prime_buffer, pfile);
					pcounter_check = 0;
				}


				if (counter >= chunksize) { break; }
			}
			else {
				// if it exited the g_sieve, the counter didn't increment,
				// so let's increment it here
				counter++;
				counter_check++;

				if (counter >= chunksize) { break; }
			}

			// write counter to communications buffer so the head node knows what's going on
			counter_buffer << "\n" << counter << "\n";
			if (counter_check >= COUNTER_BUFFER_SIZE) {
				write_counter(counter_buffer, cfile);
				counter_check = 0;
			}

			// increment candidate by one
			mpz_add_ui(n, n, 1);
		}

		// stop clock to time the batch
		clock_t t_to_end = clock();

		// determine time it took to validate candidate
		float delta_t = t_to_end - t_to_start;
		delta_t = delta_t/1000000;

		// cout << "took " << delta_t << "s to validate" << endl;
		cout << "finished chunk after: " << counter << " trials, found: " << pcounter << " primes, took: " << delta_t << " seconds to validate" << endl;
		cout.flush();

		// make sure to write the final counter when everything is done
		counter_buffer.str("");
		counter_buffer.clear();
		counter_buffer << "\n" << counter << "\t" << delta_t << "\n";
		write_counter(counter_buffer, cfile);
	}

	// (unreachable in practice...)
	mpz_clears(n, offset, NULL);
	mpz_clear(base);
	gmp_randclear(rng_state);
	return 0;
}
