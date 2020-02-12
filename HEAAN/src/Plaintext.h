/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/
#ifndef HEAAN_PLAINTEXT_H_
#define HEAAN_PLAINTEXT_H_

#include <NTL/ZZ.h>
#include "Params.h"

using namespace std;
using namespace NTL;

class Plaintext {
public:

	ZZ* mx = new ZZ[N];

	long logp;
	long logq;
	long n;

#ifdef CIPHERTEXT_EXTENDED
    double nu; // (public) bound on plaintext, ie. |m| < nu.
    double B; //  noise bound, ie.  |ax*sx+bx-m| < B.
#endif

	Plaintext( long logp = 0, long logq = 0, long n = 0, double nu = FP_NAN, double B = FP_NAN );

	virtual ~Plaintext();
};

/* Format plaintext and send it to some stream. ^ostringstream to put stuff in a string
 * 
 * Plaintext<\n
 *   logq=%d, logp=%d, n=%d|N=%d,\n
 *   |plaintext|<%.2f, error<%.2f\n
 * >\n
 * m=%dd*X^%d+...+%dd.
 * 
 * As the output would be tremendously long, for the time being, we only print the leading and the constant coefficient.
 * 
 * The second option is to use
 *   stream << std::make_pair( plaintext, scheme )
 * which additionally prints the "slots".
 */
ostream& operator<<(ostream& s, const Plaintext& o);

class Scheme; // forward declaration
ostream& operator<<(ostream& s, const std::pair<Plaintext*,Scheme*> o);
ostream& operator<<(ostream& s, const std::tuple<Plaintext*,Plaintext*,Scheme*> o);

#endif
