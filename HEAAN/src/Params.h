/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/
#ifndef HEAAN_PARAMS_H_
#define HEAAN_PARAMS_H_

//#define NDEBUG // ifdef then all asserts and maintenance of nu,B vanish
//#define SECURE // ifdef then use secure logN,logQ, otherwise use smaller values

#include <assert.h>
#ifndef NDEBUG
/*
 * Ifndef NDEBUG then assert(<cond>) evaluates <cond> and
 * if it is false then it prints a message and aborts the run.
 * Once debugging is done, you define NDEBUG and all the assert(...)
 * commands are discarded by the preprocessor.
 */
#define CIPHERTEXT_EXTENDED
/*
 * Ifdef CIPHERTEXT_EXTENDED then a bound on the plaintexts nu and
 * a bound B on the error including noise are added to
 * classes Plaintext and Ciphertext and maintained by all functions.
 */
#else
#define SECURE
#endif

#include <NTL/ZZ.h>
using namespace NTL;

#ifdef SECURE
static const long logN = 16;
static const long logQ = 1200;
#else
// used 6,200 before adding bootstrapping
//static const long logN = 6;
//static const long logQ = 200;
// For bootstrapping we need logQ be larger.
static const long logN = 7;
static const long logQ = 2400;
#endif


static const double sigma = 3.2;
static const long h = 64;
static const long pbnd = 59.0;
static const long kbar = pbnd + 1;
static const long kbar2 = 2 * kbar;
static const long logNh = ( logN - 1 ); // This is not log(N*h) = logN + log(h)
static const long logQQ = ( 2 * logQ );
static const long N = ( 1 << logN );
static const long Nh = ( 1 << logNh );
static const long M = ( N << 1 );
static const long nprimes = ( 2 + logN + 4 * logQ + pbnd - 1 ) / pbnd;
static const long Nnprimes = ( nprimes << logN );
static const long cbnd = ( logQQ + NTL_ZZ_NBITS - 1 ) / NTL_ZZ_NBITS;
static const long bignum = 0xfffffff;
static const ZZ Q = power2_ZZ ( logQ );
static const ZZ QQ = power2_ZZ ( logQQ );

static const double sqrtN = ( 1 << ( ( logN+1 ) >>1 ) );
static const double sqrth = sqrt ( h );
static const double sqrtNh = sqrtN * sqrth;
static const double Bclean = sigma* ( sqrt ( 128 ) *N+6*sqrtN+16*sqrtNh );
static const double Brs = sqrtN*sqrt ( 1.0/3 ) * ( 3+8*sqrth );
static const double Bks = sqrt ( 64.0/3 ) *sigma*N;
#define Bmult(newq) newq*Bks+Brs;

// for Plaintext and Ciphertext, the following output the right P:
//#define P(var) (double)(1UL<<var.logp) // This works as long as logp<64.
#define P(var) exp2(var.logp) // Here, integers are too small. Faster solution?
#define Q(var) exp2(var.logq) // Here, integers are too small. Faster solution? Or does the compiler care for doing that?


#endif /* PARAMS_H_ */
