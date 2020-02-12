/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/
#ifndef HEAAN_CIPHERTEXT_H_
#define HEAAN_CIPHERTEXT_H_

#include <NTL/ZZ.h>

#include <fstream>
#include "Params.h"

using namespace std;
using namespace NTL;

class Ciphertext {
public:

	ZZ* ax = new ZZ[N];
	ZZ* bx = new ZZ[N];

	long logp;
	long logq;

	long n;
    
#ifdef CIPHERTEXT_EXTENDED
    double nu; // (public) bound on plaintext, ie. |m| < nu.
    double B; //  noise bound, ie.  |ax*sx+bx-m| < B with high probability.
#endif

	Ciphertext( long logp = 0, long logq = 0, long n = 0, double nu = FP_NAN, double B = FP_NAN );

	Ciphertext(const Ciphertext& o);

	void copyParams(Ciphertext& o);

	void copy(Ciphertext& o);

	void free();

	virtual ~Ciphertext();
	
};

/* Format ciphertext and send it to some stream. ^ostringstream to put stuff in a string
 * 
 * Ciphertext<\n
 *   logq=%d, logp=%d, n=%d|N=%d,\n
 *   |plaintext|<%.2f, error<%.2f\n
 * >\n
 * a=%dd*X^%d+...+%dd,\n
 * b=%dd*X^%d+...+%dd.
 * 
 * As the output would be tremendously long, for the time being, we only print the leading and the constant coefficient.
 */
ostream& operator<<(ostream& s, const Ciphertext& o);
/*
 * Format ciphertext and print plaintext.
 * Second version compares two ciphertexts, third a ciphertext with a plaintext.
 */
class Plaintext;
class SecretKey;
class Scheme;
ostream& operator<< ( ostream& s, const std::tuple<Ciphertext*,SecretKey*,Scheme*> o );
ostream& operator<< ( ostream& s, const std::tuple<Ciphertext*,Ciphertext*,SecretKey*,Scheme*> o );
ostream& operator<< ( ostream& s, const std::tuple<Ciphertext*,Ciphertext*,SecretKey*,Scheme*> o );
//ostream& operator<< ( ostream& s, const std::tuple<Ciphertext*,Plaintext*,SecretKey*,Scheme*> o );

#endif
