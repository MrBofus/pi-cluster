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
	ofstream outfile(filename, ios::app);
	outfile << buffer.str();
	outfile.close();
	buffer.str("");
	buffer.clear();
}

// main loop
int main(void) {
	unsigned int chunksize;
	unsigned int startingint;

	// get chunk size from tasking service
	cin >> chunksize;

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
	mpz_t n, step;
	mpz_inits(n, step, NULL);

	mpz_set_ui(step, chunksize);


	//````````````````````````````````````````````````````````````````````````````````````````````````````````//

	// initialize counter and prime number counter
	unsigned int counter = 0;
	unsigned int counter_check = 0;

	unsigned int pcounter = 0;
	unsigned int pcounter_check = 0;

	// begin random search
	while (true) {

		cin >> startingint;
		if (startingint == 12121212) { exit(0); }
		mpz_set_ui(n, startingint);

		clock_t t_to_start = clock();
		while ( counter < startingint+chunksize ) {
			// check value against sieve to see if it's worth checking
			while( checkLastDigit_mpz_fast(n) ) {
				// iterate counter
				counter++;
				counter_check++;
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
				if (pcounter_check >= PRIME_BUFFER_SIZE) {
					write_to_file(prime_buffer, pfile);
					pcounter_check = 0;
				}
			}
			else {
				// if it exited the g_sieve, the counter didn't increment,
				// so let's increment it here
				counter++;
				counter_check++;
			}

			// write counter to communications buffer so the head node knows what's going on
			counter_buffer << "\n" << counter << "\n";
			if (counter_check >= COUNTER_BUFFER_SIZE) {
				ofstream outfile(cfile);
				outfile << counter_buffer.str();
				outfile.close();

				counter_buffer.str("");
				counter_buffer.clear();

				counter_check = 0;
			}

			// increment candidate by one
			mpz_add_ui(n, n, 1);
		}

		clock_t t_to_end = clock();

		// determine time it took to validate candidate
		float delta_t = t_to_end - t_to_start;
		delta_t = delta_t/1000000;

		cout << "took " << delta_t << "s to validate" << endl;

		// make sure to write the final counter when everything is done
		counter_buffer.str("");
		counter_buffer.clear();
		counter_buffer << "\n" << counter << "\t" << delta_t << "\n";
		ofstream outfile(cfile);
		outfile << counter_buffer.str();
		outfile.close();

		counter_buffer.str("");
		counter_buffer.clear();
	}

	// (unreachable in practice...)
	mpz_clears(n,step, NULL);
	return 0;
}
