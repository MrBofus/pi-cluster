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

	// present user with menu 
	cout << "choose a mode: \n\t(1) - for lucas lehmer";
	cout << "\n\t(2) - for random miller-rabin";
	cout << "\n\t(3) - specific value primality test\n\n";
	cout << "remember to have fun out there !! :)\n";
	cout << "mode: ";
	cin >> mode;

	// ``````````````````````````````````````````````````````````````````````````````````````` //
	// Lucas-Lehmer search mode:
	if (mode == 1) {
		
		// initialize lower and upper bounds for powers to check
		unsigned int lower = 65123456;
		unsigned int upper = 65223456;

		// initialize power variable
		unsigned int power = lower;

		// determine ntrials and initialize counter to zero
		unsigned int ntrials = upper-lower;
		unsigned int counter = 0;

		// LL search loop
		while (counter < ntrials){
			// clear the screen
			printf("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n");
			printf("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n");
			printf("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n");

			// show full art logo + some info
			std::cout << full_logo << '\n';
			printf("              **Lucas-Lehmer Variant**\n");
			printf("\n\n\n`````````````````````````````````````````````````````\n");
			printf("finding new candidate...\n");

			// find a new power candidate (power must be prime)
			while (power < upper){
				// if power is not prime, increment by one
				power += 1;
				if (isPrime(power)){
					// if power is prime, break
					break;
				}
			}

			// determine number of digits of power + display it
			float ndigits = log10(power);
			printf("\npower is %.0f digits (%d)\n", ndigits, power);

			// check if mersenne number is prime given power
			bool val = LucasLehmer(power);

			// if the mersenne number prime, enter here
			if (val){
				printf("\ncandidate is prime\n");

				// initialize mpz number to 0
				mpz_t result, base;
				mpz_inits(result, base, NULL);

				// compute prime mersenne number given power
				mpz_set_str(base, "2", 10);
				mpz_pow_ui(result, base, power);
				mpz_sub_ui(result, result, 1);

				// write prime mersenne number to file
				primefile = fopen("text_files/test.txt", "a");
				fputs("\n", primefile);
				fprintf(primefile, "%u", power);
				fputs(" -- ", primefile);
				mpz_out_str(primefile, 10, result);
				fputs("\n", primefile);
				fclose(primefile);

				// clear the mpz numbers

				mpz_clears(result, base, NULL);
			}
			else{
				printf("\ncandidate is not prime\n");
			}

			// if all numbers in range have been checked, break
			if (power > upper){
				break;
			}

			// iterate counter
			counter++;
		}
	}


	// ``````````````````````````````````````````````````````````````````````````````````````` //
	// miller-rabin random search mode:
	else if (mode == 2) {

		// initialize time counter and number of digits
		float delta_t = 0;
		unsigned int val = 0;
		char* mpz_seed_str;

		// define and initialize random number generator
		gmp_randstate_t rstate;
		gmp_randinit_mt(rstate);

		// define and initialize random seed for rng
		unsigned long seed = time(0);
		gmp_randseed_ui(rstate, seed);

		// define and initialize the mpz numbers to zero
		mpz_t lower, upper, base, random;
		mpz_inits(lower, upper, base, random, NULL);

		// intialize 'base' to 10
		mpz_set_ui(base, 10);

		// display info + ask user how many digits for prime number
		cout << "\n\n\n\n" << endl;
		cout << "beginning random number prime search !! :)" << endl;
		cout << "number of digits: ";
		cin >> val;

		// define lower and upper bound values given number of digits
		mpz_pow_ui(lower, base, val-1);
		mpz_pow_ui(upper, base, val);

		// generate random number in range
		mpz_urandomm(random, rstate, upper);

		// add random number to values being checked
		mpz_add(lower, random, lower);
		mpz_add(upper, random, upper);
		
		//`````````````````````````````````````````````````````//

		// clear screen and display info
		cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n" << endl;
		cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n" << endl;
		cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n" << endl;
		cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n" << endl;
		cout << "\n\n\n`````````````````````````````````````````````````````" << endl;
		cout << "beginning search for " << val << " digit prime..." << endl;
		cout << "\n`````````````````````````````````````````````````````" << endl;
		cout << "\n\n\n\n" << endl;

		// initialize counter and prime number counter
		unsigned int counter = 0;
		unsigned int pcounter = 0;

		// begin miller-rabin random search
		while (true){

			// iterate counter
			counter++;

			// check the last digit of number to see if it's worth checking
			while(checkLastDigit_mpz_fast(lower)){
				mpz_add_ui(lower, lower, 1);
			}

			// determine number of digits in candidate and chance of being prime
			size_t length = mpz_sizeinbase(lower, 10);
			float chance = 100.00/length;
			
			// show logo and info to user
			cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n" << endl;
			
			/*
			std::cout << full_logo << '\n';
			printf("            **Random Miller-Rabin Search**\n");
			cout << "\n````````````````````````````````";
			cout << "``````````````````````````" << endl;
			mpz_seed_str = return_seed(lower, 64);
			cout << "seed: " << mpz_seed_str << endl;
			cout << "\n````````````````````````````````";
			cout << "``````````````````````````" << endl;
			// print_gmp(lower);
			cout << "\tanalyzing candidate #" << counter << endl;
			cout << "\t | --- candidate has " << length << " digits" << endl;
			cout << "\t | --- prev. candidate took " << delta_t << "s to validate" << endl;
			cout << "\t | --- total of " << pcounter << " found" << endl;
			cout << "\t | --- (" << chance << "% chance of being prime)" << endl;
			cout.flush();
			*/

			logoizer(return_seed(lower, 64), counter, length, delta_t, pcounter, chance);

			// start timer and check if value is prime given value and 
			// rng state
			clock_t t_to_start = clock();
			if (isPrime_mpz_fast_light(lower, rstate, 3)){

				// if it is prime, increment the prime counter
				pcounter++;

				// tell user it was prime
				cout << "\tcandidate was prime (total of " << pcounter << " found)" << endl;

				// write prime to text file
				primefile = fopen("testing_new_algo.txt", "a");
				fputs("\n", primefile);
				mpz_out_str(primefile, 10, lower);
				fputs("\n", primefile);

				fclose(primefile);

			}
			else {
				cout << "\tcandidate was not prime (total of " << pcounter << " found)" << endl;
			}

			// stop clock
			clock_t t_to_end = clock();

			// determine time it took to validate candidate
			delta_t = t_to_end - t_to_start;
			delta_t = delta_t/1000000;
			cout << "\t(took " << delta_t << "s to validate)" << endl;

			// increment candidate by one
			mpz_add_ui(lower, lower, 1);

			// if all values have been checked, break
			int check = mpz_cmp(lower, upper);
			if (check > 0){ break; }
		}

		// clear rng and mpz numbers
		gmp_randclear(rstate);
		mpz_clears(lower, upper, base, random, NULL);
	}
	

	// ``````````````````````````````````````````````````````````````````````````````````````` //
	// individual value mode:
	else if (mode == 3){

		// initialize time counter
		float delta_t = 0;

		// define and initialize random number generator
		gmp_randstate_t rstate;
		gmp_randinit_mt(rstate);

		// define and initialize random seed for rng
		unsigned long seed;
		gmp_randseed_ui(rstate, seed);

		// define and initialize the mpz numbers to zero
		mpz_t value, random;
		mpz_inits(value, random, NULL);

		// print pretty logo and request a number to validate
		// primality
		cout << "\n\n\n\n" << endl;
		cout << full_logo << "\n\n\n";
		cout << "is your value prime?" << endl;
		cout << "enter a value: ";

		// define a string and read it from the line
		std::string val;
		cin >> val;

		// determine length of string + push it into 
		// a char* array for mpz
		const int strlength = val.length(); 
		char* char_array = new char[strlength + 1];
		strcpy(char_array, val.c_str()); 

		// set the mpz number to the char array, and tell
		// it base 10
		mpz_set_str(value, char_array, 10);

		// inform user progress was made
		cout << "\n\n\n`````````````````````````````````````````````````````" << endl;
		cout << "determining primality..." << endl;
		cout << "\n`````````````````````````````````````````````````````" << endl;

		// determine information about the number and show it before validating
		size_t length = mpz_sizeinbase(value, 10);
		float chance = 100.00/length;
		cout << "\t | --- candidate has " << length << " digits" << endl;
		cout << "\t | --- (" << chance << "% chance of being prime)" << endl;

		// flush cout to make it print, sometimes it waits for later to print it
		cout.flush();

		// start the clock
		clock_t t_to_start = clock();

		// determine if the value is prime given the rng
		if (isPrime_mpz_fast(value, rstate, 12)){
			// if it's prime let the user know
			cout << "\tvalue is prime" << endl;
		}
		else {
			// if it's not prime let the user know
			cout << "\tvalue is not prime" << endl;
		}

		// stop the clock
		clock_t t_to_end = clock();

		// determine the time it took to validate primality + let user know
		delta_t = t_to_end - t_to_start;
		delta_t = delta_t/1000000;
		cout << "\t(took " << delta_t << "s to validate)\n\n" << endl;

		// clear the rng and the mpz numbers
		gmp_randclear(rstate);
		mpz_clears(value, random, NULL);
	}


	// ``````````````````````````````````````````````````````````````````````````````````````` //
	// mode 4
	else if (mode == 4) {
		
		// initialize time counter and number of digits
		float delta_t = 0;
		unsigned int val = 0;
		char* mpz_seed_str;

		// define and initialize random number generator
		gmp_randstate_t rstate;
		gmp_randinit_mt(rstate);

		// define and initialize random seed for rng
		srand((unsigned int)time(0));
		unsigned long seed = rand()%100000000000;
		gmp_randseed_ui(rstate, seed);

		// define and initialize the mpz numbers to zero
		mpz_t lower, upper, base, random;
		mpz_inits(lower, upper, base, random, NULL);

		// intialize 'base' to 10
		mpz_set_ui(base, 10);

		// display info + ask user how many digits for prime number
		cout << "\n\n\n\n" << endl;
		cout << "beginning random number prime search !! :)" << endl;
		cout << "seed: " << seed << endl;
		cout << "number of digits: ";
		cin >> val;

		// define lower and upper bound values given number of digits
		mpz_pow_ui(lower, base, val-1);
		mpz_pow_ui(upper, base, val);

		// generate random number in range
		mpz_urandomm(random, rstate, upper);

		// add random number to values being checked
		mpz_add(lower, random, lower);
		mpz_add(upper, random, upper);
		


		// buffer to hold primes before writing them to a file
		const size_t BUFFER_SIZE = 1;
		ostringstream prime_buffer;

		char* str_lower = NULL;
		const char* file = "prime_candidates_1.txt";



		//`````````````````````````````````````````````````````//

		// clear screen and display info
		cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n" << endl;
		cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n" << endl;
		cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n" << endl;
		cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n" << endl;
		mpz_seed_str = return_seed(lower, 32);
		cout << "seed: " << mpz_seed_str << endl;
		cout << "\n\n\n`````````````````````````````````````````````````````" << endl;
		cout << "beginning search for " << val << " digit prime..." << endl;
		cout << "\n`````````````````````````````````````````````````````" << endl;
		cout << "\n\n\n\n" << endl;

		// initialize counter and prime number counter
		unsigned int counter = 0;
		unsigned int pcounter = 0;
		unsigned int pcounter_check = 0;

		// begin miller-rabin random search
		while (true){

			// iterate counter
			counter++;

			// check the last digit of number to see if it's worth checking
			while(checkLastDigit_mpz_fast(lower)){
				mpz_add_ui(lower, lower, 1);
			}

			// determine number of digits in candidate and chance of being prime
			// size_t length = mpz_sizeinbase(lower, 10);
			// float chance = 100.00/length;
			
			cout << "\r  analyzing candidate #" << counter;
			cout << " --- total of " << pcounter << " found";
			// cout.flush();

			// start timer and check if value is prime given value and 
			// rng state
			if (isPrime_mpz_fast(lower, rstate, 8)){

				
				// if it is prime, increment the prime counter
				pcounter++;
				pcounter_check++;

				// write prime to text file
				str_lower = mpz_get_str(NULL, 10, lower);
				prime_buffer << "\n" << str_lower << "\n";
				if (pcounter_check >= BUFFER_SIZE){
					write_primes_to_file(prime_buffer, file);
					pcounter_check = 0;
				}
			}
			else {}

			// increment candidate by one
			mpz_add_ui(lower, lower, 1);

			// if all values have been checked, break
			int check = mpz_cmp(lower, upper);
			if (check > 0){ break; }
		}

		// clear rng and mpz numbers
		gmp_randclear(rstate);
		mpz_clears(lower, upper, base, random, NULL);
	}

	// fastest test method to date, uses built-in mpz prime finding algorithm + 30-wheel sieve
	else if (mode == 5) {
		
		// initialize time counter and number of digits
		float delta_t = 0;
		unsigned int val = 0;
		char* mpz_seed_str;

		// define and initialize random number generator
		gmp_randstate_t rstate;
		gmp_randinit_mt(rstate);

		// define and initialize random seed for rng
		srand((unsigned int)time(0));
		unsigned long seed = rand()%100000000000;
		gmp_randseed_ui(rstate, seed);

		// define and initialize the mpz numbers to zero
		mpz_t lower, upper, base, random;
		mpz_inits(lower, upper, base, random, NULL);

		// intialize 'base' to 10
		mpz_set_ui(base, 10);

		// display info + ask user how many digits for prime number
		cout << "\n\n\n\n" << endl;
		cout << "beginning random number prime search !! :)" << endl;
		cout << "seed: " << seed << endl;
		cout << "number of digits: ";
		cin >> val;

		// define lower and upper bound values given number of digits
		mpz_pow_ui(lower, base, val-1);
		mpz_pow_ui(upper, base, val);

		// generate random number in range
		mpz_urandomm(random, rstate, upper);

		// add random number to values being checked
		mpz_add(lower, random, lower);
		mpz_add(upper, random, upper);
		


		// buffer to hold primes before writing them to a file
		const size_t BUFFER_SIZE = 1;
		ostringstream prime_buffer;

		char* str_lower = NULL;
		const char* file = "prime_candidates_1.txt";



		//`````````````````````````````````````````````````````//

		// clear screen and display info
		cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n" << endl;
		cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n" << endl;
		cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n" << endl;
		cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n" << endl;
		mpz_seed_str = return_seed(lower, 32);
		cout << "seed: " << mpz_seed_str << endl;
		cout << "\n\n\n`````````````````````````````````````````````````````" << endl;
		cout << "beginning search for " << val << " digit prime..." << endl;
		cout << "\n`````````````````````````````````````````````````````" << endl;
		cout << "\n\n\n\n" << endl;

		// initialize counter and prime number counter
		unsigned int counter = 0;
		unsigned int pcounter = 0;
		unsigned int pcounter_check = 0;

		// begin miller-rabin random search
		while (true){

			// check the last digit of number to see if it's worth checking
			while(checkLastDigit_mpz_fast(lower)){
				// iterate counter
				counter++;
				mpz_add_ui(lower, lower, 1);
			}

			// determine number of digits in candidate and chance of being prime
			// size_t length = mpz_sizeinbase(lower, 10);
			// float chance = 100.00/length;
			
			cout << "\r  analyzing candidate #" << counter;
			cout << " --- total of " << pcounter << " found" << flush;

			// start timer and check if value is prime given value
			// clock_t t_to_start = clock();
			if (mpz_probab_prime_p(lower, 3)){

				
				// if it is prime, increment the prime counter
				pcounter++;
				pcounter_check++;

				// write prime to text file
				str_lower = mpz_get_str(NULL, 10, lower);
				prime_buffer << "\n" << str_lower << "\n";
				if (pcounter_check >= BUFFER_SIZE){
					write_primes_to_file(prime_buffer, file);
					pcounter_check = 0;
				}
			}
			else {}

			// stop clock
			// clock_t t_to_end = clock();

			// determine time it took to validate candidate
			// delta_t = t_to_end - t_to_start;
			// delta_t = delta_t/1000000;
			// cout << "\t(took " << delta_t << "s to validate)" << endl;

			// increment candidate by one
			mpz_add_ui(lower, lower, 1);

			// if all values have been checked, break
			int check = mpz_cmp(lower, upper);
			if (check > 0){ break; }
		}

		// clear rng and mpz numbers
		gmp_randclear(rstate);
		mpz_clears(lower, upper, base, random, NULL);
	}

	else if (mode == 6) {

		// buffer to hold primes before writing them to a file
		const size_t BUFFER_SIZE = 1;
		ostringstream prime_buffer;

		// file name to store prime candidates
		char* str_lower = NULL;
		const char* file = "prime_candidates.txt";

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

		// clear screen and display info
		cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n" << endl;
		cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n" << endl;
		cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n" << endl;
		cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n" << endl;
		cout << "seed: " << seed << endl;
		cout << "\n\n\n`````````````````````````````````````````````````````" << endl;
		cout << "beginning search for " << val << " digit prime..." << endl;
		cout << "\n`````````````````````````````````````````````````````" << endl;
		cout << "\n\n\n\n" << endl;

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
	
		// (unreachable in practice; add a termination condition if you want one)
		gmp_randclear(rstate);
		mpz_clears(n, n_minus, n_plus, step, k_start, NULL);
	}


	else if (mode == 7) {

		// buffer to hold primes before writing them to a file
		const size_t BUFFER_SIZE = 1;
		ostringstream prime_buffer;

		// file name to store prime candidates
		char* str_lower = NULL;
		const char* file = "prime_candidates.txt";

		// get from user the desired number of digits
		// to test candidates for primality
		unsigned int val = 0;
		cout << "number of digits: ";
		cin >> val;

		// ---- compute m so 2^m has ~val digits ----
		// log2(10) ≈ 3.3219280948873623
		const unsigned long m =
			(unsigned long) floor((val - 1) * 3.3219280948873623);


		// initialize the sieve before anything starts
		// our sieve will check for 10^8 primes before sending
		// it to the primality test
		cout << "building sieve..." << endl;
		const uint64_t BOUND = 100'000'000ULL;
		auto sp = build_sieve_plan(BOUND, m);
		cout << "sieve built" << endl;


		// initialize values before starting
		mpz_t mpz_k, n_minus, n_plus, two_m;
    	mpz_inits(mpz_k, n_minus, n_plus, two_m, NULL);

		mpz_set_ui(two_m, 1);
   		mpz_mul_2exp(two_m, two_m, m);       // 2^m, computed once

		// 60-bit random k: well below 2^m for our m, satisfies Proth's k < 2^m
        srand((unsigned)time(NULL));
        uint64_t k = ((uint64_t)rand() << 32) ^ (uint64_t)rand();
        k |= 1;                          // force odd
        k &= ((uint64_t)1 << 60) - 1;
        if (k < 3) k = 3;


		//`````````````````````````````````````````````````````//

		// clear screen and display info
		cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n" << endl;
		cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n" << endl;
		cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n" << endl;
		cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n" << endl;
		cout << "\n\n\n`````````````````````````````````````````````````````" << endl;
		cout << "beginning search for " << val << " digit prime..." << endl;
		cout << "\n`````````````````````````````````````````````````````" << endl;
		cout << "\n\n\n\n" << endl;

		// initialize counter and prime number counter
		unsigned int counter = 0;
		unsigned int pcounter = 0;
		unsigned int pcounter_check = 0;

		// begin random search
		while (true) {
			counter += 2;
			// we can do '-' for BPSW or '+' for Proth

			bool m_killed = sieve_kills(sp, k, '-');
			bool p_killed = sieve_kills(sp, k, '+');

			if (!m_killed) {
				cout << "\r  analyzing candidate #" << counter-1;
				cout << " --- total of " << pcounter << " found (this one is minus)" << flush;

				// Build n - 1 = k·2^m - 1
				mpz_set_ui(mpz_k, 0);
				mpz_set_ui(mpz_k, (unsigned long)k);  // k fits in u64; mpz_set_ui takes ulong
				mpz_mul(n_minus, mpz_k, two_m);
				mpz_sub_ui(n_minus, n_minus, 1);

				// BPSW (no specialized fast test cheaper than this at our m)
				if (mpz_probab_prime_p(n_minus, 0)) {

					cout << "k: " << k << "m: " << m << " -1" << endl;

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
			}

			if (!p_killed) {
				cout << "\r  analyzing candidate #" << counter;
				cout << " --- total of " << pcounter << " found (this one is plus)  " << flush;

				mpz_set_ui(mpz_k, (unsigned long)k);
				mpz_mul(n_plus, mpz_k, two_m);
				mpz_add_ui(n_plus, n_plus, 1);

				// Proth — deterministic, ~half the cost of BPSW
				int verdict = proth_is_prime(n_plus);
				if (verdict == 1) {
					cout << "k: " << k << "m: " << m << " -1" << endl;

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
				} else if (verdict == -1) {
					// No suitable a found; fall back to BPSW. Rare.
					if (mpz_probab_prime_p(n_plus, 0)) {
						cout << "k: " << k << "m: " << m << " -1" << endl;

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
				}
			}

			k += 2;
		}
	}
	return 0;
}
