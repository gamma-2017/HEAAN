/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/
#include "Plaintext.h"

#ifdef CIPHERTEXT_EXTENDED
Plaintext::Plaintext(long logp, long logq, long n, double nu, double B)
  : logp(logp), logq(logq), n(n), nu(nu), B(B) {

}
#else
Plaintext::Plaintext(long logp, long logq, long n, double nu, double B)
  : logp(logp), logq(logq), n(n) {

}
#endif

Plaintext::~Plaintext() {
	delete[] mx;
}

inline string myformat( const char *fmt, double x )
{
    static char s[256];
    snprintf( s, 255, fmt, x );
    return string(s);
}
#include <NTL/lip.h>
ostream& operator<<(ostream& s, const Plaintext& o)
{
    long i;
    
    long sz = sizeof(o)+N*sizeof(o.mx[0]);
    for (long i=0; i<N; i++) sz += _ntl_gsize( o.mx[i].rep );
    
    s << "Plaintext<"
      <<  " // " << sz << " Bytes " // should be the true size
      << endl
      << "    logq=" << o.logq << ", logp=" << o.logp << ", n=" << o.n << "|N=" << N;
#ifdef CIPHERTEXT_EXTENDED
    s << "," << endl
      << "    |plain|\u2264" << myformat("%.3g",o.nu) 
      << ", error\u2264" << myformat("%.3g",o.B) << "=" << myformat("2^%.2f",log(o.B)/log(2));
//    s << ", pl-err\u2264" << myformat("%.2e",o.B/(double)(1UL<<o.logp)) << "=" << myformat("2^%.2f",log(o.B)/log(2)-o.logp);
#endif
    s << " >" << endl;
    s << "  m=";
    for (i=N-1; i>=0; i--)
    {
        if (o.mx[i]!=0) break;
    }
    s << o.mx[i] << "*X^" << i;
    s << "+...+";
    s << o.mx[0] << "," << endl;
    
    return s;
}

#include "Scheme.h"

ostream& operator<<(ostream& s, const std::pair<Plaintext*,Scheme*> o)
{
    s << *o.first;
    complex<double> *res = o.second->decode(*o.first);
    double norm = 0;
    long n = o.first->n;
    for (long i=0; i<n; i++)
    {
        norm = max(norm,abs(res[i]));
#ifndef QUIET
//#ifdef VERBOSE
        if (i<5 || i>n-2) {
            s << "  [" << i << "]=" << res[i] << "," << endl;
        } else if (i==5) s << "  ..." << endl;
#endif
    }
    s << "  |plain|=" << norm << "." << endl;
    // Explicit delete needed? YES!, reflecting explicit "new []" command in scheme->decode.
    delete[] res;
    return s;
}
ostream& operator<<(ostream& s, const std::tuple<Plaintext*,Plaintext*,Scheme*> o)
{
    Plaintext* plain0 = std::get<0>(o);
    Plaintext* plain1 = std::get<1>(o);
    Scheme* scheme = std::get<2>(o);
    s << *plain0;
    s << *plain1;
    if (plain1->n != plain0->n)
        s << "  WARNING: differing #slots, taking minimum..." << endl;
    double B = plain0->B + plain1->B;
    long n = min(plain0->n,plain1->n);
    complex<double> *data0 = scheme->decode(*plain0);
    complex<double> *data1 = scheme->decode(*plain1);
    double norm = 0;
    double maxb = 0;
    for (long i=0; i<n; i++)
    {
        double b = abs(data0[i]-data1[i]);
        norm = max(norm,abs(data0[i]));
        if (b>maxb) maxb=b;
#ifndef QUIET
//#ifdef VERBOSE
        if (i<5 || i>n-2) {
            s << "  [" << i << "]="
              << data0[i] << "~"
              << data1[i] << "\u00B1"
              << b << (( b<B) ? "\u2714":"\033[1;31m\u2620\u2620\u2620\033[0m") << ","
              << endl;
        } else if (i==5) s << "  ..." << endl;
#endif
    }
    s << "  |plain|=" << norm << ".";
    s << "  max~ = \u00B1" << maxb << (( maxb<B) ? "\u2714":"\033[1;31m\u2620\u2620\u2620\033[0m")
      << "=B*2^" << log2(maxb/B)
      << "." << endl;
    // Explicit deletes needed? YES!, reflecting explicit "new []" commands in scheme->decode.
    delete[] data0;
    delete[] data1;
    return s;
}
