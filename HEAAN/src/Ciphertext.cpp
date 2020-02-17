/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/
#include "Ciphertext.h"

#include <NTL/tools.h>
#include <assert.h>

#ifdef CIPHERTEXT_EXTENDED
Ciphertext::Ciphertext ( long logp, long logq, long n, double nu, double B )
    : logp ( logp ), logq ( logq ), n ( n ), nu ( nu ), B ( B )
{
    assert ( logp>=0 && logq>=logp && logq<=logQ && n<N );
    assert ( nu>=0 && B>=0 );
}
#else
Ciphertext::Ciphertext ( long logp, long logq, long n, double nu, double B )
    : logp ( logp ), logq ( logq ), n ( n )
{
    assert ( logp>=0 && logq>=logp && logq<=logQ && n<N );
}
#endif
/*
Ciphertext::Ciphertext(long logp, long logq, long n) : logp(logp), logq(logq), n(n) {
    assert( logp>0 && logq>logp && logq<=logQ && n<N );
#ifdef CIPHERTEXT_EXTENDED
    nu = FP_INFINITE;
    B = FP_NAN;
#endif
}//*/

Ciphertext::Ciphertext ( const Ciphertext& o ) : logp ( o.logp ), logq ( o.logq ), n ( o.n )
{
#ifdef CIPHERTEXT_EXTENDED
    nu = o.nu;
    B = o.B;
    assert ( nu>0 && B>0 );
#endif
    for ( long i = 0; i < N; ++i ) {
        ax[i] = o.ax[i];
        bx[i] = o.bx[i];
    }
}

void Ciphertext::copyParams ( Ciphertext& o )
{
    logp = o.logp;
    logq = o.logq;
    n = o.n;
#ifdef CIPHERTEXT_EXTENDED
    nu = o.nu;
    B = o.B;
    assert ( nu>0 && B>0 );
#endif
}

void Ciphertext::copy ( Ciphertext& o )
{
    copyParams ( o );
    for ( long i = 0; i < N; ++i ) {
        ax[i] = o.ax[i];
        bx[i] = o.bx[i];
    }
}

void Ciphertext::free()
{
    for ( long i = 0; i < N; ++i ) {
        clear ( ax[i] );
        clear ( bx[i] );
    }
}

Ciphertext::~Ciphertext()
{
    delete[] ax;
    delete[] bx;
}


/* Format ciphertext and send it to some stream. ^ostringstream to put stuff in a string
 *
 * Ciphertext<\n
 *   logq=%d, logp=%d, n=%d|N=%d,\n
 *   |plaintext|<=%.2f, error<=%.2f\n
 * >\n
 * a=%dd*X^%d+...+%dd,\n
 * b=%dd*X^%d+...+%dd.
 *
 * As the output would be tremendously long, for the time being, we only print the leading and the constant coefficient.
 */
inline string myformat ( const char *fmt, double x )
{
    static char s[256];
    snprintf ( s, 255, fmt, x );
    return string ( s );
}
ostream& operator<< ( ostream& s, const Ciphertext& o )
{
    long i;

    long sz = sizeof ( o )+N* ( sizeof ( o.ax[0] )+sizeof ( o.bx[0] ) );
    for ( long i=0; i<N; i++ ) sz += _ntl_gsize ( o.ax[i].rep ) +  _ntl_gsize ( o.bx[i].rep );

    s << "Ciphertext<"
      <<  " // " << sz << " Bytes " // should be the true size
      << endl;
    s << "    logq=" << o.logq << ", logp=" << o.logp << ", n=" << o.n << "|N=" << N;
#ifdef CIPHERTEXT_EXTENDED
    s << "," << endl
      << "    |plain|\u2264" << myformat ( "%.3g",o.nu )
      << ", error\u2264" << myformat ( "%.3g",o.B ) << "=" << myformat ( "2^%.2f",log ( o.B ) /log ( 2 ) );
//    s << ", pl-err\u2264" << myformat("%.2e",o.B/(double)(1UL<<o.logp)) << "=" << myformat("2^%.2f",log(o.B)/log(2)-o.logp);
#endif
    s << " >" << endl;
    if (o.n>0 && o.logq>0) {
        s << "  a=";
        for ( i=N-1; i>=0; i-- ) {
            if ( o.ax[i]!=0 ) break;
        }
        s << o.ax[i] << "*X^" << i;
        s << "+...+";
        s << o.ax[0] << "," << endl;
        s << "  b=";
        for ( i=N-1; i>=0; i-- ) {
            if ( o.bx[i]!=0 ) break;
        }
        s << o.bx[i] << "*X^" << i;
        s << "+...+";
        s << o.bx[0] << "." << endl;
    } else {
        s << "  Value not initialized" << endl;
    }
    return s;
}
#include "Plaintext.h"
#include "Scheme.h"
#include "SecretKey.h"
#define CIPHERTEXT_DISPLAYSLOTS 5
ostream& operator<< ( ostream& s, const std::tuple<Ciphertext*,SecretKey*,Scheme*> o )
{
    Ciphertext* cipher = std::get<0> ( o );
    SecretKey* sk = std::get<1> ( o );
    Scheme* scheme = std::get<2> ( o );
    
    if (!cipher) {
        s << "Ciphertext<?>" << endl;
        return s;
    }
    s << *cipher;
    
    if (!sk || !scheme) {
        s << "  Missing key or scheme." << endl;        
    }
    
    Plaintext plain;
    scheme->decryptMsg ( plain, *sk, *cipher );
    if (1) { // extended version
        long i;
        s << "  m=";
        for (i=N-1; i>=0; i-- ) {
            if ( plain.mx[i]!=0 ) break;
        }
        s << plain.mx[i] << "*X^" << i;
        if (0) { // either see all
            cout << "+" << endl;
            for ( i--; i>=1; i-- ) s << "    " << plain.mx[i] << "*X^" << i << "+" << endl;
        }
        else // or nothing in the middle
            s << "+...+";
        s << plain.mx[0] << "," << endl;
    }
    complex<double> *data = scheme->decode ( plain );
    long n = plain.n;
    double norm=0;
    for ( long i=0; i<n; i++ ) {
        norm = max(norm,abs(data[i]));
#ifndef QUIET
//#ifdef VERBOSE
        if ( i<CIPHERTEXT_DISPLAYSLOTS || i>n-2 ) {
            s << "  [" << i << "]="
              << data[i] << ","
              << endl;
        } else if ( i==CIPHERTEXT_DISPLAYSLOTS ) s << "  ..." << endl;
#endif
    }
    s << "  |plain|=" << norm << "." << endl;
    // Explicit delete needed? YES!, reflecting explicit "new []" command in scheme->decode.
    delete[] data;
    return s;
}
ostream& operator<< ( ostream& s, const std::tuple<Ciphertext*,Ciphertext*,SecretKey*,Scheme*> o )
{
    Ciphertext* cipher0 = std::get<0> ( o );
    Ciphertext* cipher1 = std::get<1> ( o );
    SecretKey* sk = std::get<2> ( o );
    Scheme* scheme = std::get<3> ( o );
    if (!cipher0 || cipher0->n<=0 || cipher0->logq<=0) {
        s << "Ciphertext<?>" << endl;
    } else if (!cipher1 || cipher1->n<=0 || cipher1->logq<=0) {
        s << "No comparison ciphertext given." << endl
          << std::make_tuple( cipher0, sk, scheme );
        return s;
    } else {
        s << *cipher0;        
    }
    if (!cipher1 || cipher1->n<=0 || cipher1->logq<=0) {
        s << "Ciphertext<?>" << endl;
    } else {
        s << *cipher1;
    }        
    if (!cipher0  || cipher0->n<=0 || cipher0->logq<=0 || !cipher1 || cipher1->n<=0 || cipher1->logq<=0)
        return s;
    if (!sk || !scheme) {
        s << "  Missing key or scheme." << endl;        
    }
    Plaintext plain0;
    Plaintext plain1;
    scheme->decryptMsg ( plain0, *sk, *cipher0 );
    scheme->decryptMsg ( plain1, *sk, *cipher1 );
    if (plain1.n != plain0.n)
        s << "  WARNING: differing #slots, taking minimum..." << endl;
    double B = plain0.B + plain1.B;
    long n = min(plain0.n,plain1.n);
    complex<double> *data0 = scheme->decode ( plain0 );
    complex<double> *data1 = scheme->decode ( plain1 );
    double norm = 0;
    double maxb = 0;
    for (long i=0; i<n; i++)
    {
        double b = abs(data0[i]-data1[i]);
        norm = max(norm,abs(data0[i]));
        if (b>maxb) maxb=b;
#ifndef QUIET
//#ifdef VERBOSE
        if (i<CIPHERTEXT_DISPLAYSLOTS || i>n-2) {
            s << "  [" << i << "]="
              << data0[i] << "~"
              << data1[i] << "\u00B1"
              << b << (( b<B) ? "\u2714":"\033[1;31m\u2620\u2620\u2620\033[0m") << ","
              << endl;
        } else if (i==CIPHERTEXT_DISPLAYSLOTS) s << "  ..." << endl;
#endif
    }
    s << "  |plain|=" << norm << ".";
    s << "  max~ = \u00B1" << maxb
      << "=2^" << log2(maxb)
      << "=B*2^" << log2(maxb/B)
      << (( maxb<B) ? "\u2714":"\033[1;31m\u2620\u2620\u2620\033[0m")
      << "." << endl;
    // Explicit deletes needed? YES!, reflecting explicit "new []" commands in scheme->decode.
    delete[] data0;
    delete[] data1;
    return s;
}
