/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/
#include "Scheme.h"

#include "NTL/BasicThreadPool.h"
#include "StringUtils.h"
#include "SerializationUtils.h"

Scheme::Scheme ( SecretKey& secretKey, Ring& ring, bool isSerialized ) : ring ( ring ), isSerialized ( isSerialized )
{
    addEncKey ( secretKey );
    addMultKey ( secretKey );
};

Scheme::~Scheme()
{
    for ( auto const& t : keyMap )
        delete t.second;
    for ( auto const& t : leftRotKeyMap )
        delete t.second;
}

void Scheme::addEncKey ( SecretKey& secretKey )
{
    ZZ* ax = new ZZ[N];
    ZZ* bx = new ZZ[N];

    long np = ceil ( ( 1 + logQQ + logN + 2 ) / ( double ) pbnd );
    ring.sampleUniform2 ( ax, logQQ );
    ring.mult ( bx, secretKey.sx, ax, np, QQ );
    ring.subFromGaussAndEqual ( bx, QQ );

    Key* key = new Key();
    ring.CRT ( key->rax, ax, nprimes );
    ring.CRT ( key->rbx, bx, nprimes );
    delete[] ax;
    delete[] bx;

    if ( isSerialized ) {
        string path = "serkey/ENCRYPTION.txt";
        SerializationUtils::writeKey ( key, path );
        serKeyMap.insert ( pair<long, string> ( ENCRYPTION, path ) );
        delete key;
    } else {
        keyMap.insert ( pair<long, Key*> ( ENCRYPTION, key ) );
    }
}

void Scheme::addMultKey ( SecretKey& secretKey )
{
    ZZ* ax = new ZZ[N];
    ZZ* bx = new ZZ[N];
    ZZ* sxsx = new ZZ[N];

    long np = ceil ( ( 1 + logQQ + logN + 2 ) / ( double ) pbnd );
    ring.sampleUniform2 ( ax, logQQ );
    ring.mult ( bx, secretKey.sx, ax, np, QQ );
    ring.subFromGaussAndEqual ( bx, QQ );

    np = ceil ( ( 2 + logN + 2 ) / ( double ) pbnd );
    ring.mult ( sxsx, secretKey.sx, secretKey.sx, np, Q );
    ring.leftShiftAndEqual ( sxsx, logQ, QQ );
    ring.addAndEqual ( bx, sxsx, QQ );
    delete[] sxsx;

    Key* key = new Key();
    ring.CRT ( key->rax, ax, nprimes );
    ring.CRT ( key->rbx, bx, nprimes );
    delete[] ax;
    delete[] bx;
    if ( isSerialized ) {
        string path = "serkey/MULTIPLICATION.txt";
        SerializationUtils::writeKey ( key, path );
        serKeyMap.insert ( pair<long, string> ( MULTIPLICATION, path ) );
        delete key;
    } else {
        keyMap.insert ( pair<long, Key*> ( MULTIPLICATION, key ) );
    }
}

void Scheme::addConjKey ( SecretKey& secretKey )
{
    ZZ* ax = new ZZ[N];
    ZZ* bx = new ZZ[N];

    long np = ceil ( ( 1 + logQQ + logN + 2 ) / ( double ) pbnd );
    ring.sampleUniform2 ( ax, logQQ );
    ring.mult ( bx, secretKey.sx, ax, np, QQ );
    ring.subFromGaussAndEqual ( bx, QQ );

    ZZ* sxconj = new ZZ[N];
    ring.conjugate ( sxconj, secretKey.sx );
    ring.leftShiftAndEqual ( sxconj, logQ, QQ );
    ring.addAndEqual ( bx, sxconj, QQ );
    delete[] sxconj;

    Key* key = new Key();
    ring.CRT ( key->rax, ax, nprimes );
    ring.CRT ( key->rbx, bx, nprimes );
    delete[] ax;
    delete[] bx;

    if ( isSerialized ) {
        string path = "serkey/CONJUGATION.txt";
        SerializationUtils::writeKey ( key, path );
        serKeyMap.insert ( pair<long, string> ( CONJUGATION, path ) );
        delete key;
    } else {
        keyMap.insert ( pair<long, Key*> ( CONJUGATION, key ) );
    }
}

void Scheme::addLeftRotKey ( SecretKey& secretKey, long r )
{
    assert( 0<=r && r < Nh );

    ZZ* ax = new ZZ[N];
    ZZ* bx = new ZZ[N];

    long np = ceil ( ( 1 + logQQ + logN + 2 ) / ( double ) pbnd );
    ring.sampleUniform2 ( ax, logQQ );
    ring.mult ( bx, secretKey.sx, ax, np, QQ );
    ring.subFromGaussAndEqual ( bx, QQ );

    ZZ* spow = new ZZ[N];
    ring.leftRotate ( spow, secretKey.sx, r );
    ring.leftShiftAndEqual ( spow, logQ, QQ );
    ring.addAndEqual ( bx, spow, QQ );
    delete[] spow;
    
    if (0) { // MN: DEBUGGING
        Ciphertext tmp( logQ, logQ, Nh );
        tmp.logq = logQQ;
        tmp.ax = ax;
        tmp.bx = bx;
        cout << "addLeftRotKey( secretKey, " << r << "):" << endl;
        cout << std::make_tuple( &tmp, &secretKey, this ) << endl;
        cout << "storing [rax,rbx] from CRT( r?x, ?x, " << nprimes << ")" << endl;
        cout << "addLeftRotKey done." << endl;
        tmp.ax = 0; /* makes sure ax is NOT freed */
        tmp.bx = 0; /* makes sure bx is NOT freed */
    }

    Key* key = new Key();
    ring.CRT ( key->rax, ax, nprimes );
    ring.CRT ( key->rbx, bx, nprimes );
    delete[] ax;
    delete[] bx;

    if ( isSerialized ) {
        string path = "serkey/ROTATION_" + to_string ( r ) + ".txt";
        SerializationUtils::writeKey ( key, path );
        serLeftRotKeyMap.insert ( pair<long, string> ( r, path ) );
        delete key;
    } else {
        leftRotKeyMap.insert ( pair<long, Key*> ( r, key ) );
    }
}

void Scheme::addRightRotKey ( SecretKey& secretKey, long r )
{
    long idx = Nh - r;
    if ( leftRotKeyMap.find ( idx ) == leftRotKeyMap.end() && serLeftRotKeyMap.find ( idx ) == serLeftRotKeyMap.end() ) {
        addLeftRotKey ( secretKey, idx );
    }
}

void Scheme::addLeftRotKeys ( SecretKey& secretKey )
{
    for ( long i = 0; i < logN - 1; ++i ) {
        long idx = 1 << i;
        if ( leftRotKeyMap.find ( idx ) == leftRotKeyMap.end() && serLeftRotKeyMap.find ( idx ) == serLeftRotKeyMap.end() ) {
            addLeftRotKey ( secretKey, idx );
        }
    }
}

void Scheme::addRightRotKeys ( SecretKey& secretKey )
{
    for ( long i = 0; i < logN - 1; ++i ) {
        long idx = Nh - ( 1 << i );
        if ( leftRotKeyMap.find ( idx ) == leftRotKeyMap.end() && serLeftRotKeyMap.find ( idx ) == serLeftRotKeyMap.end() ) {
            addLeftRotKey ( secretKey, idx );
        }
    }
}

void Scheme::addBootKey ( SecretKey& secretKey, long logl, long logp )
{
    ring.addBootContext ( logl, logp );

    addConjKey ( secretKey );
    addLeftRotKeys ( secretKey );

    long loglh = logl/2;
    long k = 1 << loglh;
    long m = 1 << ( logl - loglh );

    for ( long i = 1; i < k; ++i ) {
        if ( leftRotKeyMap.find ( i ) == leftRotKeyMap.end() && serLeftRotKeyMap.find ( i ) == serLeftRotKeyMap.end() ) {
            addLeftRotKey ( secretKey, i );
        }
    }

    for ( long i = 1; i < m; ++i ) {
        long idx = i * k;
        if ( leftRotKeyMap.find ( idx ) == leftRotKeyMap.end() && serLeftRotKeyMap.find ( idx ) == serLeftRotKeyMap.end() ) {
            addLeftRotKey ( secretKey, idx );
        }
    }
}

void Scheme::encode ( Plaintext& plain, double* vals, long n, long logp, long logq,
                      double nu )
{
    plain.logp = logp;
    plain.logq = logq;
    plain.n = n;
#ifdef CIPHERTEXT_EXTENDED
    plain.nu = nu; // Use INFINITY if attacker knows no bound.
    plain.B = ( fpclassify ( nu ) ==FP_NORMAL ) ? ( max ( exp2 ( -logp ),nu*DBL_EPSILON ) ) : max ( exp2 ( -logp ),DBL_EPSILON );
    plain.B *= 2*N;
    for ( long i=0; i<n; i++ ) assert ( abs ( vals[i] ) <= nu );
#endif
    ring.encode ( plain.mx, vals, n, logp + logQ );
}

void Scheme::encode ( Plaintext& plain, complex<double>* vals, long n, long logp, long logq,
                      double nu )
{
//    abort();
    plain.logp = logp;
    plain.logq = logq;
    plain.n = n;
#ifdef CIPHERTEXT_EXTENDED
    plain.nu = nu; // Use INFINITY if attacker knows no bound.
    plain.B = ( fpclassify ( nu ) ==FP_NORMAL ) ? ( max ( exp2 ( -logp ),nu*DBL_EPSILON ) ) : max ( exp2 ( -logp ),DBL_EPSILON );
    plain.B *= 2*N;
    for ( long i=0; i<n; i++ ) assert ( abs ( vals[i] ) <= nu );
#endif
    ring.encode ( plain.mx, vals, n, logp + logQ );
}

complex<double>* Scheme::decode ( Plaintext& plain )
{
// cout << "Scheme::decode " << plain << "yields "; // DEBUG
    complex<double>* res = new complex<double>[plain.n];
    ring.decode ( plain.mx, res, plain.n, plain.logp, plain.logq );
// cout << res[0] << "," << res[1] << "..." << endl; // DEBUG
    return res;
}

void Scheme::encodeSingle ( Plaintext& plain, double val, long logp, long logq,
                            double nu )
{
    plain.logp = logp;
    plain.logq = logq;
    plain.n = 1;
#ifdef CIPHERTEXT_EXTENDED
    plain.nu = nu; // Use INFINITY if attacker knows no bound.
    plain.B = ( fpclassify ( nu ) ==FP_NORMAL ) ? ( max ( exp2 ( -logp ),nu*DBL_EPSILON ) ) : max ( exp2 ( -logp ),DBL_EPSILON );
    plain.B *= 2*N;
    assert ( abs ( val ) <= nu );
#endif
    plain.mx[0] = EvaluatorUtils::scaleUpToZZ ( val, logp + logQ );
}

void Scheme::encodeSingle ( Plaintext& plain, complex<double> val, long logp, long logq,
                            double nu )
{
    plain.logp = logp;
    plain.logq = logq;
    plain.n = 1;
#ifdef CIPHERTEXT_EXTENDED
    plain.nu = nu;
    plain.B = ( fpclassify ( nu ) ==FP_NORMAL ) ? ( max ( exp2 ( -logp ),nu*DBL_EPSILON ) ) : max ( exp2 ( -logp ),DBL_EPSILON );
    plain.B *= 2*N;
    assert ( abs ( val ) <= nu );
#endif
    plain.mx[0] = EvaluatorUtils::scaleUpToZZ ( val.real(), logp + logQ );
    plain.mx[Nh] = EvaluatorUtils::scaleUpToZZ ( val.imag(), logp + logQ );
}

complex<double> Scheme::decodeSingle ( Plaintext& plain )
{
    ZZ q = ring.qpows[plain.logq];

    complex<double> res;
    ZZ tmp = plain.mx[0] % q;
    if ( NumBits ( tmp ) == plain.logq ) tmp -= q;
    res.real ( EvaluatorUtils::scaleDownToReal ( tmp, plain.logp ) );

    tmp = plain.mx[Nh] % q;
    if ( NumBits ( tmp ) == plain.logq ) tmp -= q;
    res.imag ( EvaluatorUtils::scaleDownToReal ( tmp, plain.logp ) );

    return res;
}

void Scheme::encryptMsg ( Ciphertext& cipher, Plaintext& plain )
{
    cipher.logp = plain.logp;
    cipher.logq = plain.logq;
    cipher.n = plain.n;
#ifdef CIPHERTEXT_EXTENDED
    cipher.nu = plain.nu;
    cipher.B = plain.B + Bclean / P ( plain );
#endif
    ZZ qQ = ring.qpows[plain.logq + logQ];

    ZZ* vx = new ZZ[N];
    ring.sampleZO ( vx );

    Key* key = isSerialized ? SerializationUtils::readKey ( serKeyMap.at ( ENCRYPTION ) ) : keyMap.at ( ENCRYPTION );

    long np = ceil ( ( 1 + logQQ + logN + 2 ) / ( double ) pbnd );
    ring.multNTT ( cipher.ax, vx, key->rax, np, qQ );
    ring.addGaussAndEqual ( cipher.ax, qQ );

    ring.multNTT ( cipher.bx, vx, key->rbx, np, qQ );
    ring.addGaussAndEqual ( cipher.bx, qQ );
    delete[] vx;

    ring.addAndEqual ( cipher.bx, plain.mx, qQ );

    ring.rightShiftAndEqual ( cipher.ax, logQ );
    ring.rightShiftAndEqual ( cipher.bx, logQ );
    
    if (isSerialized) delete key;
}

void Scheme::decryptMsg ( Plaintext& plain, SecretKey& secretKey, Ciphertext& cipher )
{
    ZZ q = ring.qpows[cipher.logq];
    plain.logp = cipher.logp;
    plain.logq = cipher.logq;
    plain.n = cipher.n;
#ifdef CIPHERTEXT_EXTENDED
    plain.nu = cipher.nu;
    plain.B = cipher.B;
#endif
    long np = ceil ( ( 1 + cipher.logq + logN + 2 ) / ( double ) pbnd );
    ring.mult ( plain.mx, cipher.ax, secretKey.sx, np, q );
    ring.addAndEqual ( plain.mx, cipher.bx, q );
}

void Scheme::encrypt ( Ciphertext& cipher, complex<double>* vals, long n, long logp, long logq, double nu )
{
    Plaintext plain;
    assert ( 0<=logq && logq <= logQ );
    assert ( logp>=0 && logp<=logq );
    assert ( 0<=n && n<=N && 0==N%n );
    assert ( nu>=0 );
    encode ( plain, vals, n, logp, logq, nu );
    encryptMsg ( cipher, plain );
}

void Scheme::encrypt ( Ciphertext& cipher, double* vals, long n, long logp, long logq, double nu )
{
    Plaintext plain;
    assert ( 0<=logq && logq <= logQ );
    assert ( logp>=0 && logp<=logq );
    assert ( 0<=n && n<=N && 0==N%n );
    assert ( nu>=0 );
    encode ( plain, vals, n, logp, logq, nu );
    encryptMsg ( cipher, plain );
}

void Scheme::encryptZeros ( Ciphertext& cipher, long n, long logp, long logq, double nu )
{
    assert ( 0<=logq && logq <= logQ );
    assert ( logp>=0 && logp<=logq );
    assert ( 0<=n && n<=N && 0==N%n );
    assert ( nu>=0 );
    encryptSingle ( cipher, 0.0, logp, logq, nu );
    cipher.n = n;
}

complex<double>* Scheme::decrypt ( SecretKey& secretKey, Ciphertext& cipher )
{
    Plaintext plain;
    decryptMsg ( plain, secretKey, cipher );
    return decode ( plain );
}

void Scheme::encryptSingle ( Ciphertext& cipher, complex<double> val, long logp, long logq, double nu )
{
    Plaintext plain;
    assert ( 0<=logq && logq <= logQ );
    assert ( logp>=0 && logp<=logq );
    assert ( nu>=0 );
    encodeSingle ( plain, val, logp, logq, nu );
    encryptMsg ( cipher, plain );
}

void Scheme::encryptSingle ( Ciphertext& cipher, double val, long logp, long logq, double nu )
{
    Plaintext plain;
    assert ( 0<=logq && logq <= logQ );
    assert ( logp>=0 && logp<=logq );
    assert ( nu>=0 );
    encodeSingle ( plain, val, logp, logq, nu );
    encryptMsg ( cipher, plain );
}

complex<double> Scheme::decryptSingle ( SecretKey& secretKey, Ciphertext& cipher )
{
    Plaintext plain;
    decryptMsg ( plain, secretKey, cipher );
    return decodeSingle ( plain );
}

//-----------------------------------------

void Scheme::negate ( Ciphertext& res, Ciphertext& cipher )
{
    res.copyParams ( cipher );
    // res.nu, res.B are ok.
    ring.negate ( res.ax, cipher.ax );
    ring.negate ( res.bx, cipher.bx );
}

void Scheme::negateAndEqual ( Ciphertext& cipher )
{
    // cipher.nu, cipher.B are ok.
    ring.negateAndEqual ( cipher.ax );
    ring.negateAndEqual ( cipher.bx );
}

void Scheme::add ( Ciphertext& res, Ciphertext& cipher1, Ciphertext& cipher2 )
{
    assert ( cipher1.logq == cipher2.logq );
    assert ( cipher1.logp == cipher2.logp );
    assert ( cipher1.n    == cipher2.n );
    ZZ q = ring.qpows[cipher1.logq];
    res.copyParams ( cipher1 );
#ifdef CIPHERTEXT_EXTENDED
    res.nu += cipher2.nu; // by triangle inequality
    res.B  += cipher2.B; // by triangle inequality
#endif
    ring.add ( res.ax, cipher1.ax, cipher2.ax, q );
    ring.add ( res.bx, cipher1.bx, cipher2.bx, q );
}

void Scheme::addAndEqual ( Ciphertext& cipher1, Ciphertext& cipher2 )
{
    assert ( cipher1.logq == cipher2.logq );
    assert ( cipher1.logp == cipher2.logp );
    assert ( cipher1.n    == cipher2.n );
    ZZ q = ring.qpows[cipher1.logq];
#ifdef CIPHERTEXT_EXTENDED
    cipher1.nu += cipher2.nu; // by triangle inequality
    cipher1.B  += cipher2.B; // by triangle inequality
#endif
    ring.addAndEqual ( cipher1.ax, cipher2.ax, q );
    ring.addAndEqual ( cipher1.bx, cipher2.bx, q );
}

//-----------------------------------------

void Scheme::addConst ( Ciphertext& res, Ciphertext& cipher, double cnst, long logp )
{
    assert ( logp<0 || logp == cipher.logp );
    ZZ q = ring.qpows[cipher.logq];
    ZZ cnstZZ = logp < 0 ? EvaluatorUtils::scaleUpToZZ ( cnst, cipher.logp ) : EvaluatorUtils::scaleUpToZZ ( cnst, logp );
    res.copy ( cipher );
#ifdef CIPHERTEXT_EXTENDED
    res.nu += abs ( cnst );
#endif
    AddMod ( res.bx[0], res.bx[0], cnstZZ, q );
}

void Scheme::addConst ( Ciphertext& res, Ciphertext& cipher, RR& cnst, long logp )
{
    assert ( logp<0 || logp == cipher.logp );
    ZZ q = ring.qpows[cipher.logq];
    ZZ cnstZZ = logp < 0 ? EvaluatorUtils::scaleUpToZZ ( cnst, cipher.logp ) : EvaluatorUtils::scaleUpToZZ ( cnst, logp );
    res.copy ( cipher );
#ifdef CIPHERTEXT_EXTENDED
    res.nu += conv<double> ( abs ( cnst ) );
#endif
    AddMod ( res.bx[0], res.bx[0], cnstZZ, q );
}

void Scheme::addConst ( Ciphertext& res, Ciphertext& cipher, complex<double> cnst, long logp )
{
    assert ( logp<0 || logp == cipher.logp );
    ZZ q = ring.qpows[cipher.logq];
    ZZ cnstZZ = logp < 0 ? EvaluatorUtils::scaleUpToZZ ( cnst.real(), cipher.logp ) : EvaluatorUtils::scaleUpToZZ ( cnst.real(), logp );
    // MN: What happens with cnst.imag() ?
    res.copy ( cipher );
#ifdef CIPHERTEXT_EXTENDED
    res.nu += abs ( cnst );
#endif
    AddMod ( res.bx[0], cipher.bx[0], cnstZZ, q );
}

void Scheme::addConstAndEqual ( Ciphertext& cipher, double cnst, long logp )
{
    assert ( logp<0 || logp == cipher.logp );
    ZZ q = ring.qpows[cipher.logq];
    ZZ cnstZZ = logp < 0 ? EvaluatorUtils::scaleUpToZZ ( cnst, cipher.logp ) : EvaluatorUtils::scaleUpToZZ ( cnst, logp );
#ifdef CIPHERTEXT_EXTENDED
    cipher.nu += abs ( cnst );
#endif
    AddMod ( cipher.bx[0], cipher.bx[0], cnstZZ, q );
}

void Scheme::addConstAndEqual ( Ciphertext& cipher, const RR& cnst, long logp )
{
    assert ( logp<0 || logp == cipher.logp );
    ZZ q = ring.qpows[cipher.logq];
    ZZ cnstZZ = logp < 0 ? EvaluatorUtils::scaleUpToZZ ( cnst, cipher.logp ) : EvaluatorUtils::scaleUpToZZ ( cnst, logp );
#ifdef CIPHERTEXT_EXTENDED
    cipher.nu += conv<double> ( abs ( cnst ) );
#endif
    AddMod ( cipher.bx[0], cipher.bx[0], cnstZZ, q );
}

void Scheme::addConstAndEqual ( Ciphertext& cipher, complex<double> cnst, long logp )
{
    assert ( logp<0 || logp == cipher.logp );
    ZZ q = ring.qpows[cipher.logq];
    ZZ cnstrZZ = logp < 0 ? EvaluatorUtils::scaleUpToZZ ( cnst.real(), cipher.logp ) : EvaluatorUtils::scaleUpToZZ ( cnst.real(), logp );
    ZZ cnstiZZ = logp < 0 ? EvaluatorUtils::scaleUpToZZ ( cnst.imag(), cipher.logp ) : EvaluatorUtils::scaleUpToZZ ( cnst.imag(), logp );
#ifdef CIPHERTEXT_EXTENDED
    cipher.nu += abs ( cnst );
#endif
    AddMod ( cipher.bx[0], cipher.bx[0], cnstrZZ, q );
    AddMod ( cipher.bx[Nh], cipher.bx[Nh], cnstiZZ, q );
}

//-----------------------------------------

void Scheme::sub ( Ciphertext& res, Ciphertext& cipher1, Ciphertext& cipher2 )
{
    assert ( cipher1.logq == cipher2.logq );
    assert ( cipher1.logp == cipher2.logp );
    assert ( cipher1.n    == cipher2.n );
    ZZ q = ring.qpows[cipher1.logq];
    res.copyParams ( cipher1 );
#ifdef CIPHERTEXT_EXTENDED
    res.nu += cipher2.nu; // by triangle inequality
    res.B  += cipher2.B; // by triangle inequality
#endif
    ring.sub ( res.ax, cipher1.ax, cipher2.ax, q );
    ring.sub ( res.bx, cipher1.bx, cipher2.bx, q );
}

void Scheme::subAndEqual ( Ciphertext& cipher1, Ciphertext& cipher2 )
{
    assert ( cipher1.logq == cipher2.logq );
    assert ( cipher1.logp == cipher2.logp );
    assert ( cipher1.n    == cipher2.n );
    ZZ q = ring.qpows[cipher1.logq];
#ifdef CIPHERTEXT_EXTENDED
    cipher1.nu += cipher2.nu; // by triangle inequality
    cipher1.B  += cipher2.B; // by triangle inequality
#endif
    ring.subAndEqual ( cipher1.ax, cipher2.ax, q );
    ring.subAndEqual ( cipher1.bx, cipher2.bx, q );
}

/* This puts the result (cipher1-cipher2) into cipher2!
 */
void Scheme::subAndEqual2 ( Ciphertext& cipher1, Ciphertext& cipher2 )
{
    assert ( cipher1.logq == cipher2.logq );
    assert ( cipher1.logp == cipher2.logp );
    assert ( cipher1.n    == cipher2.n );
    ZZ q = ring.qpows[cipher1.logq];
#ifdef CIPHERTEXT_EXTENDED
    cipher2.nu += cipher1.nu; // by triangle inequality
    cipher2.B  += cipher1.B; // by triangle inequality
#endif
    ring.subAndEqual2 ( cipher1.ax, cipher2.ax, q );
    ring.subAndEqual2 ( cipher1.bx, cipher2.bx, q );
}

void Scheme::imult ( Ciphertext& res, Ciphertext& cipher )
{
    ZZ q = ring.qpows[cipher.logq];
    res.copyParams ( cipher );
    // res.nu, res.B are ok.
    ring.multByMonomial ( res.ax, cipher.ax, Nh );
    ring.multByMonomial ( res.bx, cipher.bx, Nh );
}

void Scheme::idiv ( Ciphertext& res, Ciphertext& cipher )
{
    ZZ q = ring.qpows[cipher.logq];
    res.copyParams ( cipher );
    // res.nu, res.B are ok.
    ring.multByMonomial ( res.ax, cipher.ax, 3 * Nh );
    ring.multByMonomial ( res.bx, cipher.bx, 3 * Nh );
}

void Scheme::imultAndEqual ( Ciphertext& cipher )
{
    // cipher.nu, cipher.B are ok.
    ring.multByMonomialAndEqual ( cipher.ax, Nh );
    ring.multByMonomialAndEqual ( cipher.bx, Nh );
}

void Scheme::idivAndEqual ( Ciphertext& cipher )
{
    // cipher.nu, cipher.B are ok.
    ring.multByMonomialAndEqual ( cipher.ax, 3 * Nh );
    ring.multByMonomialAndEqual ( cipher.bx, 3 * Nh );
}

void Scheme::mult ( Ciphertext& res, Ciphertext& cipher1, Ciphertext& cipher2 )
{
    assert ( cipher1.logq == cipher2.logq );
    assert ( cipher1.logp == cipher2.logp );
    assert ( cipher1.n == cipher2.n );
    res.copyParams ( cipher1 );
    res.logp += cipher2.logp;
#ifdef CIPHERTEXT_EXTENDED
    res.nu *= cipher2.nu;
    double Bmult = exp2 ( res.logq-logQ ) * Bks + Brs; // MN: maybe this exp2(-logQ) * Bks + Brs instead??
    res.B  = ( cipher1.nu * cipher2.B + cipher2.nu * cipher1.B + cipher1.B * cipher2.B )
             + Bmult / P ( res ); // Lemma 3
#endif

    ZZ q = ring.qpows[cipher1.logq];
    ZZ qQ = ring.qpows[cipher1.logq + logQ];

    long np = ceil ( ( 2 + cipher1.logq + cipher2.logq + logN + 2 ) / ( double ) pbnd );

    uint64_t* ra1 = new uint64_t[np << logN];
    uint64_t* rb1 = new uint64_t[np << logN];
    uint64_t* ra2 = new uint64_t[np << logN];
    uint64_t* rb2 = new uint64_t[np << logN];

    ring.CRT ( ra1, cipher1.ax, np );
    ring.CRT ( rb1, cipher1.bx, np );
    ring.CRT ( ra2, cipher2.ax, np );
    ring.CRT ( rb2, cipher2.bx, np );

    ZZ* axax = new ZZ[N];
    ZZ* bxbx = new ZZ[N];
    ZZ* axbx = new ZZ[N];
    ring.multDNTT ( axax, ra1, ra2, np, q );
    ring.multDNTT ( bxbx, rb1, rb2, np, q );

    ring.addNTTAndEqual ( ra1, rb1, np );
    ring.addNTTAndEqual ( ra2, rb2, np );
    ring.multDNTT ( axbx, ra1, ra2, np, q );

    Key* key = isSerialized ? SerializationUtils::readKey ( serKeyMap.at ( MULTIPLICATION ) ) : keyMap.at ( MULTIPLICATION );

    np = ceil ( ( cipher1.logq + logQQ + logN + 2 ) / ( double ) pbnd );
    uint64_t* raa = new uint64_t[np << logN];
    ring.CRT ( raa, axax, np );
    ring.multDNTT ( res.ax, raa, key->rax, np, qQ );
    ring.multDNTT ( res.bx, raa, key->rbx, np, qQ );
    ring.rightShiftAndEqual ( res.ax, logQ );
    ring.rightShiftAndEqual ( res.bx, logQ );

    ring.addAndEqual ( res.ax, axbx, q );
    ring.subAndEqual ( res.ax, bxbx, q );
    ring.subAndEqual ( res.ax, axax, q );
    ring.addAndEqual ( res.bx, bxbx, q );

    if (isSerialized) delete key;
    delete[] axax; // Are these not done automatically?
    delete[] bxbx;
    delete[] axbx;
    delete[] ra1;
    delete[] ra2;
    delete[] rb1;
    delete[] rb2;
    delete[] raa;
}

void Scheme::multAndEqual ( Ciphertext& cipher1, Ciphertext& cipher2 )
{
    assert ( cipher1.logq == cipher2.logq );
    assert ( cipher1.logp == cipher2.logp );
    assert ( cipher1.n == cipher2.n );
#ifdef CIPHERTEXT_EXTENDED
    cipher1.nu *= cipher2.nu;
    double Bmult = exp2 ( cipher1.logq-logQ ) * Bks + Brs; // MN: maybe this exp2(-logQ) * Bks + Brs instead??
    cipher1.B  = ( cipher1.nu * cipher2.B + cipher2.nu * cipher1.B + cipher1.B * cipher2.B )
                 + Bmult / P ( cipher2 ); // Lemma 3
#endif

    ZZ q = ring.qpows[cipher1.logq];
    ZZ qQ = ring.qpows[cipher1.logq + logQ];

    long np = ceil ( ( 2 + cipher1.logq + cipher2.logq + logN + 2 ) / ( double ) pbnd );

    uint64_t* ra1 = new uint64_t[np << logN];
    uint64_t* rb1 = new uint64_t[np << logN];
    uint64_t* ra2 = new uint64_t[np << logN];
    uint64_t* rb2 = new uint64_t[np << logN];

    ring.CRT ( ra1, cipher1.ax, np );
    ring.CRT ( rb1, cipher1.bx, np );
    ring.CRT ( ra2, cipher2.ax, np );
    ring.CRT ( rb2, cipher2.bx, np );

    ZZ* axax = new ZZ[N];
    ZZ* bxbx = new ZZ[N];
    ZZ* axbx = new ZZ[N];

    ring.multDNTT ( axax, ra1, ra2, np, q );
    ring.multDNTT ( bxbx, rb1, rb2, np, q );
    ring.addNTTAndEqual ( ra1, rb1, np );
    ring.addNTTAndEqual ( ra2, rb2, np );
    ring.multDNTT ( axbx, ra1, ra2, np, q );

    Key* key = isSerialized ? SerializationUtils::readKey ( serKeyMap.at ( MULTIPLICATION ) ) : keyMap.at ( MULTIPLICATION );

    np = ceil ( ( cipher1.logq + logQQ + logN + 2 ) / ( double ) pbnd );
    uint64_t* raa = new uint64_t[np << logN];
    ring.CRT ( raa, axax, np );
    ring.multDNTT ( cipher1.ax, raa, key->rax, np, qQ );
    ring.multDNTT ( cipher1.bx, raa, key->rbx, np, qQ );

    ring.rightShiftAndEqual ( cipher1.ax, logQ );
    ring.rightShiftAndEqual ( cipher1.bx, logQ );

    ring.addAndEqual ( cipher1.ax, axbx, q );
    ring.subAndEqual ( cipher1.ax, bxbx, q );
    ring.subAndEqual ( cipher1.ax, axax, q );
    ring.addAndEqual ( cipher1.bx, bxbx, q );

    if (isSerialized) delete key;
    delete[] axax;
    delete[] bxbx;
    delete[] axbx;
    delete[] ra1;
    delete[] ra2;
    delete[] rb1;
    delete[] rb2;
    delete[] raa;

    cipher1.logp += cipher2.logp;
}

//-----------------------------------------

void Scheme::square ( Ciphertext& res, Ciphertext& cipher )
{
    res.copyParams ( cipher );
    res.logp += cipher.logp;
#ifdef CIPHERTEXT_EXTENDED
    res.nu *= cipher.nu; // obvious
    double Bmult = exp2 ( res.logq-logQ ) * Bks + Brs; // MN: maybe this exp2(-logQ) * Bks + Brs instead??
    res.B  = ( 2*cipher.nu * cipher.B + cipher.B * cipher.B )
             + Bmult / P ( res ); // Lemma 3
#endif
    ZZ q = ring.qpows[cipher.logq];
    ZZ qQ = ring.qpows[cipher.logq + logQ];

    long np = ceil ( ( 2 * cipher.logq + logN + 2 ) / ( double ) pbnd );

    uint64_t* ra = new uint64_t[np << logN];
    uint64_t* rb = new uint64_t[np << logN];

    ring.CRT ( ra, cipher.ax, np );
    ring.CRT ( rb, cipher.bx, np );

    ZZ* axax = new ZZ[N];
    ZZ* axbx = new ZZ[N];
    ZZ* bxbx = new ZZ[N];

    ring.squareNTT ( bxbx, rb, np, q );
    ring.squareNTT ( axax, ra, np, q );
    ring.multDNTT ( axbx, ra, rb, np, q );
    ring.addAndEqual ( axbx, axbx, q );

    Key* key = isSerialized ? SerializationUtils::readKey ( serKeyMap.at ( MULTIPLICATION ) ) : keyMap.at ( MULTIPLICATION );

    np = ceil ( ( cipher.logq + logQQ + logN + 2 ) / ( double ) pbnd );
    uint64_t* raa = new uint64_t[np << logN];
    ring.CRT ( raa, axax, np );
    ring.multDNTT ( res.ax, raa, key->rax, np, qQ );
    ring.multDNTT ( res.bx, raa, key->rbx, np, qQ );

    ring.rightShiftAndEqual ( res.ax, logQ );
    ring.rightShiftAndEqual ( res.bx, logQ );

    ring.addAndEqual ( res.ax, axbx, q );
    ring.addAndEqual ( res.bx, bxbx, q );

    if (isSerialized) delete key;
    delete[] axbx;
    delete[] axax;
    delete[] bxbx;

    delete[] ra;
    delete[] rb;
    delete[] raa;
}

void Scheme::squareAndEqual ( Ciphertext& cipher )
{
#ifdef CIPHERTEXT_EXTENDED
    cipher.nu *= cipher.nu; // obvious
    double Bmult = exp2 ( cipher.logq-logQ ) * Bks + Brs; // MN: maybe this exp2(-logQ) * Bks + Brs instead??
    cipher.B  = ( 2*cipher.nu * cipher.B + cipher.B * cipher.B )
                + Bmult / P ( cipher ); // Lemma 3
#endif
    ZZ q = ring.qpows[cipher.logq];
    ZZ qQ = ring.qpows[cipher.logq + logQ];

    long np = ceil ( ( 2 + 2 * cipher.logq + logN + 2 ) / ( double ) pbnd );

    uint64_t* ra = new uint64_t[np << logN];
    uint64_t* rb = new uint64_t[np << logN];

    ring.CRT ( ra, cipher.ax, np );
    ring.CRT ( rb, cipher.bx, np );

    ZZ* axax = new ZZ[N];
    ZZ* axbx = new ZZ[N];
    ZZ* bxbx = new ZZ[N];

    ring.squareNTT ( bxbx, rb, np, q );
    ring.squareNTT ( axax, ra, np, q );

    ring.multDNTT ( axbx, ra, rb, np, q );
    ring.addAndEqual ( axbx, axbx, q );

    Key* key = isSerialized ? SerializationUtils::readKey ( serKeyMap.at ( MULTIPLICATION ) ) : keyMap.at ( MULTIPLICATION );

    np = ceil ( ( cipher.logq + logQQ + logN + 2 ) / ( double ) pbnd );

    uint64_t* raa = new uint64_t[np << logN];
    ring.CRT ( raa, axax, np );
    ring.multDNTT ( cipher.ax, raa, key->rax, np, qQ );
    ring.multDNTT ( cipher.bx, raa, key->rbx, np, qQ );

    ring.rightShiftAndEqual ( cipher.ax, logQ );
    ring.rightShiftAndEqual ( cipher.bx, logQ );

    ring.addAndEqual ( cipher.ax, axbx, q );
    ring.addAndEqual ( cipher.bx, bxbx, q );
    cipher.logp *= 2;

    if (isSerialized) delete key;
    delete[] axbx;
    delete[] axax;
    delete[] bxbx;

    delete[] ra;
    delete[] rb;
    delete[] raa;
}

//-----------------------------------------

void Scheme::multByConst ( Ciphertext& res, Ciphertext& cipher, double cnst, long logp )
{
    assert ( logp >= 0 );
    ZZ q = ring.qpows[cipher.logq];
    ZZ cnstZZ = EvaluatorUtils::scaleUpToZZ ( cnst, logp );
    res.copyParams ( cipher );
    res.logp += logp;
#ifdef CIPHERTEXT_EXTENDED
    res.nu *= abs ( cnst );
    res.B *= abs ( cnst );
#endif
    ring.multByConst ( res.ax, cipher.ax, cnstZZ, q );
    ring.multByConst ( res.bx, cipher.bx, cnstZZ, q );
}

void Scheme::multByConst ( Ciphertext& res, Ciphertext& cipher, complex<double> cnst, long logp )
{
    assert ( logp >= 0 );
    res.copy ( cipher );
    multByConstAndEqual ( res, cnst, logp ); // MN:cares for parameters
}

void Scheme::multByConstVec ( Ciphertext& res, Ciphertext& cipher, complex<double>* cnstVec, long logp )
{
    assert ( logp >= 0 );
    res.copy ( cipher );
    multByConstVecAndEqual ( res, cnstVec, logp ); // MN:cares for parameters
}

void Scheme::multByConstVecAndEqual ( Ciphertext& cipher, complex<double>* cnstVec, long logp )
{
    assert ( logp >= 0 );
    long slots = cipher.n;
    ZZ* cnstPoly = new ZZ[N];
    ring.encode ( cnstPoly, cnstVec, slots, logp );
    double nu=0;
    for ( long i=0; i<slots; i++ )
        nu = max ( nu, abs ( cnstVec[i] ) ); // compute canonical embedding max norm.
    multByPolyAndEqual ( cipher, cnstPoly, logp, nu ); // MN:cares for parameters
    delete[] cnstPoly;
}

void Scheme::multByConstAndEqual ( Ciphertext& cipher, double cnst, long logp )
{
    assert ( logp >= 0 );
    ZZ q = ring.qpows[cipher.logq];
    ZZ cnstZZ = EvaluatorUtils::scaleUpToZZ ( cnst, logp );
    ring.multByConstAndEqual ( cipher.ax, cnstZZ, q );
    ring.multByConstAndEqual ( cipher.bx, cnstZZ, q );
    cipher.logp += logp;
#ifdef CIPHERTEXT_EXTENDED
    cipher.nu *= abs ( cnst );
    cipher.B *= abs ( cnst );
#endif
}

void Scheme::multByConstAndEqual ( Ciphertext& cipher, RR& cnst, long logp )
{
    assert ( logp >= 0 );
    ZZ q = ring.qpows[cipher.logq];
    ZZ cnstZZ = EvaluatorUtils::scaleUpToZZ ( cnst, logp );
    ring.multByConstAndEqual ( cipher.ax, cnstZZ, q );
    ring.multByConstAndEqual ( cipher.bx, cnstZZ, q );
    cipher.logp += logp;
#ifdef CIPHERTEXT_EXTENDED
    cipher.nu *= conv<double> ( abs ( cnst ) );
    cipher.B *= conv<double> ( abs ( cnst ) );
#endif
}

void Scheme::multByConstAndEqual ( Ciphertext& cipher, complex<double> cnst, long logp )
{
    assert ( logp >= 0 );
    ZZ q = ring.qpows[cipher.logq];
    ZZ cnstZZReal = EvaluatorUtils::scaleUpToZZ ( cnst.real(), logp );
    ZZ cnstZZImag = EvaluatorUtils::scaleUpToZZ ( cnst.imag(), logp );

    Ciphertext tmp; // compute imagnary part
    tmp.copyParams ( cipher );
    ring.multByConst ( tmp.ax, cipher.ax, cnstZZImag, q );
    ring.multByConst ( tmp.bx, cipher.bx, cnstZZImag, q );
    ring.multByMonomialAndEqual ( tmp.ax, N / 2 );
    ring.multByMonomialAndEqual ( tmp.bx, N / 2 );


    ring.multByConstAndEqual ( cipher.ax, cnstZZReal, q );
    ring.multByConstAndEqual ( cipher.bx, cnstZZReal, q );

    ring.addAndEqual ( cipher.ax, tmp.ax, QQ );
    ring.addAndEqual ( cipher.bx, tmp.bx, QQ );

    cipher.logp += logp;
#ifdef CIPHERTEXT_EXTENDED
    cipher.nu *= abs ( cnst );
    cipher.B *= abs ( cnst );
#endif
}

void Scheme::multByPoly ( Ciphertext& res, Ciphertext& cipher, ZZ* poly, long logp, double nu )
{
    assert ( logp>=0 );
    assert ( nu>=0 );
    ZZ q = ring.qpows[cipher.logq];
    res.copyParams ( cipher );
    long bnd = ring.maxBits ( poly, N );
    long np = ceil ( ( cipher.logq + bnd + logN + 2 ) / ( double ) pbnd );
    uint64_t* rpoly = new uint64_t[np << logN];
    ring.CRT ( rpoly, poly, np );
    ring.multNTT ( res.ax, cipher.ax, rpoly, np, q );
    ring.multNTT ( res.bx, cipher.bx, rpoly, np, q );
    delete[] rpoly;
    res.logp += logp;
#ifdef CIPHERTEXT_EXTENDED
    res.nu *= nu;
    res.B *= nu;
#endif
}

void Scheme::multByPolyAndEqual ( Ciphertext& cipher, ZZ* poly, long logp, double nu )
{
    assert ( logp>=0 );
    assert ( nu>=0 );
    // assert( caninfnorm(poly) <= nu );
    ZZ q = ring.qpows[cipher.logq];
    long bnd = ring.maxBits ( poly, N );
    long np = ceil ( ( cipher.logq + bnd + logN + 2 ) / ( double ) pbnd );
    uint64_t* rpoly = new uint64_t[np << logN];
    ring.CRT ( rpoly, poly, np );
    ring.multNTTAndEqual ( cipher.ax, rpoly, np, q );
    ring.multNTTAndEqual ( cipher.bx, rpoly, np, q );
    delete[] rpoly;
    cipher.logp += logp;
#ifdef CIPHERTEXT_EXTENDED
    cipher.nu *= nu;
    cipher.B *= nu;
#endif
}

void Scheme::multByPolyNTT ( Ciphertext& res, Ciphertext& cipher, uint64_t* rpoly, long bnd, long logp, double nu )
{
    assert ( logp>=0 );
    assert ( nu>=0 );
    // assert( caninfnorm(rpoly) <= nu );
    ZZ q = ring.qpows[cipher.logq];
    res.copyParams ( cipher );
    long np = ceil ( ( cipher.logq + bnd + logN + 2 ) / ( double ) pbnd );
    ring.multNTT ( res.ax, cipher.ax, rpoly, np, q );
    ring.multNTT ( res.bx, cipher.bx, rpoly, np, q );
    res.logp += logp;
#ifdef CIPHERTEXT_EXTENDED
    res.nu *= nu;
    res.B *= nu;
#endif
}

void Scheme::multByPolyNTTAndEqual ( Ciphertext& cipher, uint64_t* rpoly, long bnd, long logp, double nu )
{
    assert ( logp>=0 );
    assert ( nu>=0 );
    // assert( caninfnorm(rpoly) <= nu );
    ZZ q = ring.qpows[cipher.logq];
    long np = ceil ( ( cipher.logq + bnd + logN + 2 ) / ( double ) pbnd );
    ring.multNTTAndEqual ( cipher.ax, rpoly, np, q );
    ring.multNTTAndEqual ( cipher.bx, rpoly, np, q );
    cipher.logp += logp;
#ifdef CIPHERTEXT_EXTENDED
    cipher.nu *= nu;
    cipher.B *= nu;
#endif
}

//-----------------------------------------

void Scheme::multByMonomial ( Ciphertext& res, Ciphertext& cipher, const long degree )
{
    // assert( degree ); // MN: Which conditions should degree satisfy?
    res.copyParams ( cipher );
    ring.multByMonomial ( res.ax, cipher.ax, degree );
    ring.multByMonomial ( res.bx, cipher.bx, degree );
}

void Scheme::multByMonomialAndEqual ( Ciphertext& cipher, const long degree )
{
    // assert( degree ); // MN: Which conditions should degree satisfy?
    ring.multByMonomialAndEqual ( cipher.ax, degree );
    ring.multByMonomialAndEqual ( cipher.bx, degree );
}

//-----------------------------------------

/* Scheme::leftShift multiplies the polynomials cipher.ax,cipher.bx by 2^bits mod 2^cipher.logq.
 */
void Scheme::leftShift ( Ciphertext& res, Ciphertext& cipher, long bits )
{
    assert ( bits>=0 );
    // MN: bits<0 executes rightShift
    // However, that only makes sense when decreasing cipher.logq at the right moment.
    // See divByPo2.
    ZZ q = ring.qpows[cipher.logq];
    res.copyParams ( cipher );
#ifdef CIPHERTEXT_EXTENDED
    res.nu *= exp2 ( bits ); // MN: assuming no modular reduction applies!
    res.B *= exp2 ( bits ); // MN: assuming no modular reduction applies!
#endif
    ring.leftShift ( res.ax, cipher.ax, bits, q );
    ring.leftShift ( res.bx, cipher.bx, bits, q );
}

void Scheme::leftShiftAndEqual ( Ciphertext& cipher, long bits )
{
    assert ( bits>=0 );
    // MN: bits<0 executes rightShift
    // However, that only makes sense when decreasing cipher.logq at the right moment.
    // See divByPo2AndEqual.
    ZZ q = ring.qpows[cipher.logq];
#ifdef CIPHERTEXT_EXTENDED
    cipher.nu *= exp2 ( bits ); // MN: assuming no modular reduction applies!
    cipher.B *= exp2 ( bits ); // MN: assuming no modular reduction applies!
#endif
    ring.leftShiftAndEqual ( cipher.ax, bits, q );
    ring.leftShiftAndEqual ( cipher.bx, bits, q );
}

void Scheme::doubleAndEqual ( Ciphertext& cipher )
{
    // MN: same as leftShiftAndEqual( cipher, 1 )
    ZZ q = ring.qpows[cipher.logq];
#ifdef CIPHERTEXT_EXTENDED
    cipher.nu *= 2; // MN: assuming no modular reduction applies!
    cipher.B *= 2; // MN: assuming no modular reduction applies!
#endif
    ring.doubleAndEqual ( cipher.ax, q );
    ring.doubleAndEqual ( cipher.bx, q );
}

void Scheme::divByPo2 ( Ciphertext& res, Ciphertext& cipher, long bits )
{
    assert ( cipher.logq>=bits );
    res.copyParams ( cipher );
#ifdef CIPHERTEXT_EXTENDED
    res.nu *= exp2 ( -bits );
    res.B *= exp2 ( -bits );
#endif
    ring.rightShift ( res.ax, cipher.ax, bits );
    ring.rightShift ( res.bx, cipher.bx, bits );
    res.logq -= bits;
}

void Scheme::divByPo2AndEqual ( Ciphertext& cipher, long bits )
{
    assert ( cipher.logq>=bits );
#ifdef CIPHERTEXT_EXTENDED
    cipher.nu *= exp2 ( -bits );
    cipher.B *= exp2 ( -bits );
    cipher.B += exp2( logN+1-cipher.logp ); /* rounding error, supposing that this transfers... */
#endif
    ring.rightShiftAndEqual ( cipher.ax, bits );
    ring.rightShiftAndEqual ( cipher.bx, bits );
    cipher.logq -= bits;
}


//-----------------------------------------

void Scheme::reScaleBy ( Ciphertext& res, Ciphertext& cipher, long dlogq )
{
    assert ( dlogq>=0 ); // MN: negative dlogq would do something but...
    assert ( cipher.logq>=dlogq );
    assert ( cipher.logp>=dlogq );
    res.copyParams ( cipher );
    ring.rightShift ( res.ax, cipher.ax, dlogq );
    ring.rightShift ( res.bx, cipher.bx, dlogq );
    res.logp -= dlogq;
    res.logq -= dlogq;
#ifdef CIPHERTEXT_EXTENDED
    // res.B ok
#endif
}

void Scheme::reScaleTo ( Ciphertext& res, Ciphertext& cipher, long logq )
{
    long dlogq = cipher.logq - logq;
    assert ( dlogq>=0 ); // MN: negative dlogq would do something but...
    assert ( cipher.logq>=dlogq );
    assert ( cipher.logp>=dlogq );
    res.copyParams ( cipher );
    ring.rightShift ( res.ax, cipher.ax, dlogq );
    ring.rightShift ( res.bx, cipher.bx, dlogq );
    res.logp -= dlogq;
#ifdef CIPHERTEXT_EXTENDED
    // res.B ok
#endif
}

void Scheme::reScaleByAndEqual ( Ciphertext& cipher, long dlogq )
{
    assert ( dlogq>=0 ); // MN: negative dlogq would do something but...
    assert ( cipher.logq>=dlogq );
    assert ( cipher.logp>=dlogq );
    ring.rightShiftAndEqual ( cipher.ax, dlogq );
    ring.rightShiftAndEqual ( cipher.bx, dlogq );
    cipher.logq -= dlogq;
    cipher.logp -= dlogq;
#ifdef CIPHERTEXT_EXTENDED
    // cipher.B ok
#endif
}

void Scheme::reScaleToAndEqual ( Ciphertext& cipher, long logq )
{
    long dlogq = cipher.logq - logq;
    assert ( dlogq>=0 ); // MN: negative dlogq would do something but...
    assert ( cipher.logq>=dlogq );
    assert ( cipher.logp>=dlogq );
    ring.rightShiftAndEqual ( cipher.ax, dlogq );
    ring.rightShiftAndEqual ( cipher.bx, dlogq );
    cipher.logq = logq;
    cipher.logp -= dlogq;
#ifdef CIPHERTEXT_EXTENDED
    // cipher.B ok
#endif
}

void Scheme::modDownBy ( Ciphertext& res, Ciphertext& cipher, long dlogq )
{
    assert ( dlogq>=0 ); // MN: negative dlogq would do something but...
    assert ( cipher.logq>=dlogq ); // MN: can decrease cipher.logq to negative bitnumber.
    assert ( cipher.logq<=logQ );
    ZZ q = ring.qpows[cipher.logq - dlogq];
    res.copyParams ( cipher );
    ring.mod ( res.ax, cipher.ax, q );
    ring.mod ( res.bx, cipher.bx, q );
    res.logq -= dlogq;
    // res.nu, res.B unaffected unless mod q destroys everything.
}

void Scheme::modDownByAndEqual ( Ciphertext& cipher, long dlogq )
{
    assert ( dlogq>=0 ); // MN: negative dlogq would do something but...
    assert ( cipher.logq>=dlogq ); // MN: can decrease cipher.logq to negative bitnumber.
    assert ( cipher.logq<=logQ );
    ZZ q = ring.qpows[cipher.logq - dlogq];
    ring.modAndEqual ( cipher.ax, q );
    ring.modAndEqual ( cipher.bx, q );
    cipher.logq -= dlogq;
    // res.nu, res.B unaffected unless mod q destroys everything.
}

void Scheme::modDownTo ( Ciphertext& res, Ciphertext& cipher, long logq )
{
    assert ( 0<=logq && logq<=cipher.logq ); // MN: can increase cipher.logq.
    assert ( cipher.logq<=logQ );
    ZZ q = ring.qpows[logq];
    res.copyParams ( cipher );
    ring.mod ( res.ax, cipher.ax, q );
    ring.mod ( res.bx, cipher.bx, q );
    res.logq = logq;
    // res.nu, res.B unaffected unless mod q destroys everything.
}

void Scheme::modDownToAndEqual ( Ciphertext& cipher, long logq )
{
    assert ( 0<=logq && logq<=cipher.logq ); // MN: can increase cipher.logq.
    assert ( cipher.logq<=logQ );
    ZZ q = ring.qpows[logq];
    cipher.logq = logq;
    ring.modAndEqual ( cipher.ax, q );
    ring.modAndEqual ( cipher.bx, q );
    // res.nu, res.B unaffected unless mod q destroys everything.
}

//----------------------------------------------------------------------------------
//   ROTATIONS & CONJUGATIONS
//----------------------------------------------------------------------------------


void Scheme::leftRotateFast ( Ciphertext& res, Ciphertext& cipher, long r )
{
    ZZ q = ring.qpows[cipher.logq];
    ZZ qQ = ring.qpows[cipher.logq + logQ];

    ZZ* bxrot = new ZZ[N];
    ZZ* axrot = new ZZ[N];

    ring.leftRotate ( bxrot, cipher.bx, r );
    ring.leftRotate ( axrot, cipher.ax, r );

    Key* key = isSerialized ? SerializationUtils::readKey ( serLeftRotKeyMap.at ( r ) ) : leftRotKeyMap.at ( r );
    res.copyParams ( cipher );
#ifdef CIPHERTEXT_EXTENDED
    res.nu = cipher.nu;
    // See Lemma 5 in [chekim18]:
    double Bmult = exp2 ( res.logq-logQ ) * Bks + Brs;
    res.B  = cipher.B + Bmult / P ( res );
#endif

    long np = ceil ( ( cipher.logq + logQQ + logN + 2 ) / ( double ) pbnd );
    uint64_t* rarot = new uint64_t[np << logN];
    ring.CRT ( rarot, axrot, np );
    ring.multDNTT ( res.ax, rarot, key->rax, np, qQ );
    ring.multDNTT ( res.bx, rarot, key->rbx, np, qQ );

    ring.rightShiftAndEqual ( res.ax, logQ );
    ring.rightShiftAndEqual ( res.bx, logQ );
    ring.addAndEqual ( res.bx, bxrot, q );

    if (isSerialized) delete key;
    delete[] bxrot;
    delete[] axrot;
    delete[] rarot;
}

void Scheme::leftRotateFastAndEqual ( Ciphertext& cipher, long r )
{
    ZZ q = ring.qpows[cipher.logq];
    ZZ qQ = ring.qpows[cipher.logq + logQ];

    ZZ* bxrot = new ZZ[N];
    ZZ* axrot = new ZZ[N];

    ring.leftRotate ( bxrot, cipher.bx, r );
    ring.leftRotate ( axrot, cipher.ax, r );
    Key* key = isSerialized ? SerializationUtils::readKey ( serLeftRotKeyMap.at ( r ) ) : leftRotKeyMap.at ( r );
#ifdef CIPHERTEXT_EXTENDED
    // See Lemma 5 in [chekim18]:
    double Bmult = exp2 ( cipher.logq-logQ ) * Bks + Brs;
    cipher.B  += Bmult / P ( cipher );
#endif
    long np = ceil ( ( cipher.logq + logQQ + logN + 2 ) / ( double ) pbnd );
    uint64_t* rarot = new uint64_t[np << logN];
    ring.CRT ( rarot, axrot, np );
    ring.multDNTT ( cipher.ax, rarot, key->rax, np, qQ );
    ring.multDNTT ( cipher.bx, rarot, key->rbx, np, qQ );

    ring.rightShiftAndEqual ( cipher.ax, logQ );
    ring.rightShiftAndEqual ( cipher.bx, logQ );

    ring.addAndEqual ( cipher.bx, bxrot, q );

    if (isSerialized) delete key;
    delete[] bxrot;
    delete[] axrot;
    delete[] rarot;
}

void Scheme::rightRotateFast ( Ciphertext& res, Ciphertext& cipher, long r )
{
    long rr = Nh - r;
    leftRotateFast ( res, cipher, rr );
}

void Scheme::rightRotateFastAndEqual ( Ciphertext& cipher, long r )
{
    long rr = Nh - r;
    leftRotateFastAndEqual ( cipher, rr );
}

void Scheme::conjugate ( Ciphertext& res, Ciphertext& cipher )
{
    ZZ q = ring.qpows[cipher.logq];
    ZZ qQ = ring.qpows[cipher.logq + logQ];

    ZZ* bxconj = new ZZ[N];
    ZZ* axconj = new ZZ[N];

    ring.conjugate ( bxconj, cipher.bx );
    ring.conjugate ( axconj, cipher.ax );

    Key* key = isSerialized ? SerializationUtils::readKey ( serKeyMap.at ( CONJUGATION ) ) : keyMap.at ( CONJUGATION );
    res.copyParams ( cipher );
#ifdef CIPHERTEXT_EXTENDED
    res.nu = cipher.nu;
    // See Lemma 5 in [chekim18]:
    double Bmult = exp2 ( res.logq-logQ ) * Bks + Brs;
    res.B  = cipher.B + Bmult / P ( res );
#endif
    long np = ceil ( ( cipher.logq + logQQ + logN + 2 ) / ( double ) pbnd );
    uint64_t* raconj = new uint64_t[np << logN];
    ring.CRT ( raconj, axconj, np );
    ring.multDNTT ( res.ax, raconj, key->rax, np, qQ );
    ring.multDNTT ( res.bx, raconj, key->rbx, np, qQ );

    ring.rightShiftAndEqual ( res.ax, logQ );
    ring.rightShiftAndEqual ( res.bx, logQ );
    ring.addAndEqual ( res.bx, bxconj, q );

    if (isSerialized) delete key;
    delete[] bxconj;
    delete[] axconj;
    delete[] raconj;
}

void Scheme::conjugateAndEqual ( Ciphertext& cipher )
{
    ZZ q = ring.qpows[cipher.logq];
    ZZ qQ = ring.qpows[cipher.logq + logQ];

    ZZ* bxconj = new ZZ[N];
    ZZ* axconj = new ZZ[N];

    ring.conjugate ( bxconj, cipher.bx );
    ring.conjugate ( axconj, cipher.ax );

    Key* key = isSerialized ? SerializationUtils::readKey ( serKeyMap.at ( CONJUGATION ) ) : keyMap.at ( CONJUGATION );
#ifdef CIPHERTEXT_EXTENDED
    // See Lemma 5 in [chekim18]:
    double Bmult = exp2 ( cipher.logq-logQ ) * Bks + Brs;
    cipher.B  += Bmult / P ( cipher );
#endif

    long np = ceil ( ( cipher.logq + logQQ + logN + 2 ) / ( double ) pbnd );
    uint64_t* raconj = new uint64_t[np << logN];
    ring.CRT ( raconj, axconj, np );
    ring.multDNTT ( cipher.ax, raconj, key->rax, np, qQ );
    ring.multDNTT ( cipher.bx, raconj, key->rbx, np, qQ );

    ring.rightShiftAndEqual ( cipher.ax, logQ );
    ring.rightShiftAndEqual ( cipher.bx, logQ );

    ring.addAndEqual ( cipher.bx, bxconj, q );

    if (isSerialized) delete key;
    delete[] bxconj;
    delete[] axconj;
    delete[] raconj;
}


///// MN: modified til here. Tests pending...

//----------------------------------------------------------------------------------
//   BOOTSTRAPPING
//----------------------------------------------------------------------------------


void Scheme::normalizeAndEqual ( Ciphertext& cipher )
{
    ZZ q = ring.qpows[cipher.logq];

    for ( long i = 0; i < N; ++i ) {
        if ( NumBits ( cipher.ax[i] ) == cipher.logq ) cipher.ax[i] -= q;
        if ( NumBits ( cipher.bx[i] ) == cipher.logq ) cipher.bx[i] -= q;
    }
}

void Scheme::coeffToSlotAndEqual ( Ciphertext& cipher )
{
    long slots = cipher.n;
    long logSlots = log2 ( slots );
    long logk = logSlots / 2;
    long k = 1 << logk;

    Ciphertext* rotvec = new Ciphertext[k];
    rotvec[0].copy ( cipher );

    NTL_EXEC_RANGE ( k - 1, first, last );
    for ( long j = first; j < last; ++j ) {
        leftRotateFast ( rotvec[j+1], rotvec[0], j + 1 );
    }
    NTL_EXEC_RANGE_END;

    BootContext* bootContext = ring.bootContextMap.at ( logSlots );

    Ciphertext* tmpvec = new Ciphertext[k];

    //for (long j=0; j<k; j++) {
    NTL_EXEC_RANGE ( k, first, last );
    for ( long j = first; j < last; ++j ) {
        multByPolyNTT ( tmpvec[j],
                        rotvec[j],
                        bootContext->rpvec[j], bootContext->bndvec[j], bootContext->logp,
                        1.0   //MN2020/01/27: TODO, replace adhoc value by reasonable one,
                        // should be a bound to the specified polynomial
                      );
    }
    NTL_EXEC_RANGE_END;

    for ( long j = 1; j < k; ++j ) {
        addAndEqual ( tmpvec[0], tmpvec[j] );
    }

    cipher.copy ( tmpvec[0] );
    for ( long ki = k; ki < slots; ki += k ) {
        NTL_EXEC_RANGE ( k, first, last );
        for ( long j = first; j < last; ++j ) {
            multByPolyNTT ( tmpvec[j],
                            rotvec[j],
                            bootContext->rpvec[j + ki], bootContext->bndvec[j + ki], bootContext->logp,
                            1.0   //MN2020/01/27: TODO, replace adhoc value by reasonable one,
                            // should be a bound to the specified polynomial
                          );
        }
        NTL_EXEC_RANGE_END;
        for ( long j = 1; j < k; ++j ) {
            addAndEqual ( tmpvec[0], tmpvec[j] );
        }
        leftRotateFastAndEqual ( tmpvec[0], ki );
        addAndEqual ( cipher, tmpvec[0] );
    }
    reScaleByAndEqual ( cipher, bootContext->logp );
    delete[] rotvec;
    delete[] tmpvec;
}

void Scheme::slotToCoeffAndEqual ( Ciphertext& cipher )
{
    long slots = cipher.n;
    long logSlots = log2 ( slots );
    long logk = logSlots / 2;
    long k = 1 << logk;

    Ciphertext* rotvec = new Ciphertext[k];
    rotvec[0].copy ( cipher );

    NTL_EXEC_RANGE ( k-1, first, last );
    for ( long j = first; j < last; ++j ) {
        leftRotateFast ( rotvec[j + 1], rotvec[0], j + 1 );
    }
    NTL_EXEC_RANGE_END;

    BootContext* bootContext = ring.bootContextMap.at ( logSlots );

    Ciphertext* tmpvec = new Ciphertext[k];

    NTL_EXEC_RANGE ( k, first, last );
    for ( long j = first; j < last; ++j ) {
        multByPolyNTT ( tmpvec[j],
                        rotvec[j],
                        bootContext->rpvecInv[j], bootContext->bndvecInv[j], bootContext->logp,
                        1.0   //MN2020/01/28: TODO, replace adhoc value by reasonable one,
                        // should be a bound to the specified polynomial
                      );
    }
    NTL_EXEC_RANGE_END;

    for ( long j = 1; j < k; ++j ) {
        addAndEqual ( tmpvec[0], tmpvec[j] );
    }
    cipher.copy ( tmpvec[0] );

    for ( long ki = k; ki < slots; ki+=k ) {
        NTL_EXEC_RANGE ( k, first, last );
        for ( long j = first; j < last; ++j ) {
            multByPolyNTT ( tmpvec[j],
                            rotvec[j],
                            bootContext->rpvecInv[j + ki], bootContext->bndvecInv[j + ki], bootContext->logp,
                            1.0   //MN2020/01/28: TODO, replace adhoc value by reasonable one,
                            // should be a bound to the specified polynomial
                          );
        }
        NTL_EXEC_RANGE_END;

        for ( long j = 1; j < k; ++j ) {
            addAndEqual ( tmpvec[0], tmpvec[j] );
        }

        leftRotateFastAndEqual ( tmpvec[0], ki );
        addAndEqual ( cipher, tmpvec[0] );
    }
    reScaleByAndEqual ( cipher, bootContext->logp );
    delete[] rotvec;
    delete[] tmpvec;
}

/*
 * Compute cipher(x) -> cipher( exp2pi(x) )
 * with exp2pi(x) = sum( (2 Pi x)^n / n!, n=0..7 )
 * approximating exp( 2 Pi x ) to degree <8.
 * 
 * cipher=x, cipher2=x², cipher4=x⁴,
 * c = 1/ ( 2*Pi );
 * cipher01=x+c;
 * c = 2*Pi;
 * cipher01*=c;
 * c = 3/ ( 2*Pi );
 * cipher23=cipher+c;
 * c = 4*Pi*Pi*Pi/3;
 * cipher23*=c;
 * cipher23*=x²;
 * cipher23+=cipher01;
 * c = 5/ ( 2*Pi );
 * cipher45=x+c;
 * c = 4*Pi*Pi*Pi*Pi*Pi/15;
 * cipher45*=c;
 * c = 7/ ( 2*Pi );
 * cipher+=c;
 * c = 8*Pi*Pi*Pi*Pi*Pi*Pi*Pi/315;
 * cipher*=c;
 * cipher*=x²;
 * cipher+=cipher45;
 * cipher*=x⁴;
 * cipher+=cipher23;
 * 
 */
void Scheme::exp2piAndEqual ( Ciphertext& cipher, long dummy )
{
    long logp = cipher.logp;
    
    Ciphertext cipher2;
    square ( cipher2, cipher );
    reScaleByAndEqual ( cipher2, logp ); // cipher2.logq : logq - logp

    Ciphertext cipher4;
    square ( cipher4, cipher2 );
    reScaleByAndEqual ( cipher4, logp ); // cipher4.logq : logq -2logp
    RR c = 1/ ( 2*Pi );
    Ciphertext cipher01;
    addConst ( cipher01, cipher, c, -1 ); // cipher01.logq : logq // MN: -1=use same prec as cipher.

    c = 2*Pi;
    multByConstAndEqual ( cipher01, c, logp );
    reScaleByAndEqual ( cipher01, logp ); // cipher01.logq : logq - logp

    c = 3/ ( 2*Pi );
    Ciphertext cipher23;
    addConst ( cipher23, cipher, c, -1 ); // cipher23.logq : logq // MN: -1=use same prec as cipher.

    c = 4*Pi*Pi*Pi/3;
    multByConstAndEqual ( cipher23, c, logp );
    reScaleByAndEqual ( cipher23, logp ); // cipher23.logq : logq - logp

    multAndEqual ( cipher23, cipher2 );
    reScaleByAndEqual ( cipher23, logp ); // cipher23.logq : logq - 2logp

    modDownToAndEqual ( cipher01, cipher23.logq );
    addAndEqual ( cipher23, cipher01 ); // cipher23.logq : logq - 2logp

    c = 5/ ( 2*Pi );
    Ciphertext cipher45;
    addConst ( cipher45, cipher, c, logp ); // cipher45.logq : logq

    c = 4*Pi*Pi*Pi*Pi*Pi/15;
    multByConstAndEqual ( cipher45, c, logp );
    reScaleByAndEqual ( cipher45, logp ); // cipher45.logq : logq - logp

    c = 7/ ( 2*Pi );
    addConstAndEqual ( cipher, c, logp ); // cipher.logq : logq

    c = 8*Pi*Pi*Pi*Pi*Pi*Pi*Pi/315;
    multByConstAndEqual ( cipher, c, logp );
    reScaleByAndEqual ( cipher, logp ); // cipher.logq : logq - logp

    multAndEqual ( cipher, cipher2 );
    reScaleByAndEqual ( cipher, logp ); // cipher.logq : logq - 2logp

    modDownByAndEqual ( cipher45, logp ); // cipher45.logq : logq - 2logp
    addAndEqual ( cipher, cipher45 ); // cipher.logq : logq - 2logp

    multAndEqual ( cipher, cipher4 );
    reScaleByAndEqual ( cipher, logp ); // cipher.logq : logq - 3logp

    modDownByAndEqual ( cipher23, logp );
    addAndEqual ( cipher, cipher23 ); // cipher.logq : logq - 3logp
}

void Scheme::evalExpAndEqual (
    Ciphertext& cipher,
    long logT,
    long logI
)
{
    long slots = cipher.n;
    long logSlots = log2 ( slots );
    BootContext* bootContext = ring.bootContextMap.at ( logSlots );
    if ( logSlots < logNh ) {
        /* MN:
         * In this branch all coefficients fit into one vector and
         * have been put into suitable positions.
         */
        Ciphertext tmp;
        /* MN:
         * cipher = 2i Im(cipher)
         */
        conjugate ( tmp, cipher );
        subAndEqual ( cipher, tmp );
        /* MN:
         * cipher /= 2^logT 2       = i Im(cipher) / 2^logT
         */
        divByPo2AndEqual ( cipher, logT + 1 ); // bitDown: logT + 1
        /* MN:
         * Approximate exp( 2 Pi cipher ) ~= exp( 2 Pi i Im(cipher) / 2^logT )
         */
        exp2piAndEqual ( cipher, bootContext->logp ); // bitDown: logT + 1 + 3(logq + logI)
#ifdef CIPHERTEXT_EXTENDED
        /* MN:
         * Disregarding errors the result should be of absolute value 1.
         * Optimistically, we thus use cipher.nu=1 as a bound to the absolute value.
         */
        cipher.nu = 1; // MN2020/02/06: slightly optimistic gues: abs( decrypted exp2pi result ) ~= 1. Putting it here, makes the subsequent error propagation estimation more realistic.
#endif
        /* MN:
         * Square the result logT times to obtain ~= exp( 2 Pi i Im(cipher) ).
         * Square logI times more to TODO Why?
         */
        long logp1 = cipher.logp;
        for ( long i = 0; i < logI + logT; ++i ) {
            squareAndEqual ( cipher ); /* doubles post period bits */
            reScaleByAndEqual ( cipher, logp1 ); /* keep logp1 post period bits */
        }
        /* MN:
         * Compute 2i times imaginary part of result.
         */
        conjugate ( tmp, cipher );
        subAndEqual ( cipher, tmp );
        /* MN:
         * Extract first half of the results and multiply by 1/(4 Pi).
         * 1/(2 Pi) is the desired scaling for the sine function, the other 2 is for the imaginary part.
         */
        multByPolyNTT ( tmp, cipher,
                        bootContext->rp1, bootContext->bnd1, bootContext->logp,
                        0.25/M_PI   //MN2020/01/28: TODO, replace adhoc value by reasonable one,
                        // should be a bound to the specified polynomial
                      );
        Ciphertext tmprot;
        leftRotateFast ( tmprot, tmp, slots );
        addAndEqual ( tmp, tmprot );
        /* MN:
         * Extract second half and...
         */
        multByPolyNTTAndEqual ( cipher,
                                bootContext->rp2, bootContext->bnd2, bootContext->logp,
                                0.25/M_PI   //MN2020/01/28: TODO, replace adhoc value by reasonable one,
                                // should be a bound to the specified polynomial
                              );
        leftRotateFast ( tmprot, cipher, slots );
        addAndEqual ( cipher, tmprot );
        /* MN:
         * Combine both halfs to full coefficient vector.
         */
        addAndEqual ( cipher, tmp );
    }
    else {
        /* MN:
         * In this branch there are twice as many coefficients as slots,
         * they are delivered as real and imaginary parts of the slots.
         * So we have to perform the entire computation twice.
         */
        Ciphertext tmp;
        conjugate ( tmp, cipher );
        Ciphertext c2;
        sub ( c2, cipher, tmp );
        addAndEqual ( cipher, tmp );
        imultAndEqual ( cipher );
        /* MN:
         * Now c2 = 2i Im(cipher) and cipher = 2i Re(cipher).
         */
        /* MN:
         * cipher~c2 /= 2^logT 2       = i (Re~Im)(cipher) / 2^logT
         */
        divByPo2AndEqual ( cipher, logT + 1 ); // cipher bitDown: logT + 1
        divByPo2AndEqual ( c2, logT + 1 ); // c2 bitDown: logT + 1 // MN: Modified
        /* MN:
         * Approximate exp( 2 Pi cipher~c2 ) ~= exp( 2 Pi i (Re~Im)(cipher) / 2^logT )
         */
        exp2piAndEqual ( cipher, bootContext->logp ); // cipher bitDown: logT + 1 + 3(logq + logI)
        exp2piAndEqual ( c2, bootContext->logp ); // c2 bitDown: logT + 1 + 3(logq + logI)
        /* MN:
         * Disregarding errors the result should be of absolute value 1.
         * Optimistically, we thus use cipher.nu=1 as a bound to the absolute value.
         */
#ifdef CIPHERTEXT_EXTENDED
        cipher.nu = 1; // MN2020/02/06: slightly optimistic guess: abs( decrypted exp2pi result ) ~= 1. Putting it here, makes the subsequent error propagation estimation more realistic.
        c2.nu = 1;
#endif
        /* MN:
         * Square the result logT times to obtain ~= exp( 2 Pi i Im(cipher) ).
         * Square logI times more to TODO Why?
         */
        long logp1 = cipher.logp;
        assert( cipher.logp == c2.logp );
        for ( long i = 0; i < logI + logT; ++i ) {
            squareAndEqual ( c2 );
            squareAndEqual ( cipher );
            reScaleByAndEqual ( c2, logp1 );
            reScaleByAndEqual ( cipher, logp1 );
        }
        /* MN:
         * Compute 2i times imaginary part of result.
         */
        conjugate ( tmp, c2 );
        subAndEqual ( c2, tmp );
        conjugate ( tmp, cipher );
        subAndEqual ( cipher, tmp );
        imultAndEqual ( cipher );
        /* MN:
         * Now cipher = -2 Im(exp(2 Pi i Re())), c2 = 2i Im(exp( 2 Pi i Im() ))
         */
        subAndEqual2 ( c2, cipher ); // cipher = c2-cipher = 2 [ Im(exp(2 Pi i Re())) + i Im(exp( 2 Pi i Im() )) ]
        RR c = 0.25/Pi;
        multByConstAndEqual ( cipher, c, bootContext->logp );
        /* MN:
         * Now approximately
         *     Re(cipher) = 1/(2 Pi) sin( 2 Pi Re(cipher) ),
         *     Im(cipher) = 1/(2 Pi) sin( 2 Pi Im(cipher) ).
         */
    }
    /* MN:
     * Move period back in place TODO: ie?
     */
    reScaleByAndEqual ( cipher, bootContext->logp + logI );
}

void Scheme::bootstrapAndEqual (
    Ciphertext& cipher, /* ciphertext to be bootstrapped */
    long logq, /* number of plaintext bits to consider, any higher bits will be destroyed */
    long logQ, /* nunber of toplevel bits */
    long logT, /* number of halvings/squarings in evalExpAndEqual/exp2piAndEqual computation */
    long logI  /* ?  in evalExpAndEqual/exp2piAndEqual computation */
)
{
    long logSlots = log2 ( cipher.n );
    long logp = cipher.logp;
#ifdef CIPHERTEXT_EXTENDED
    double nu = cipher.nu; /* incoming value bound */
    double B = cipher.B; /* incoming error bound */
#endif

    /* MN:
     * First, erase potentially available info above logq.
     */
    modDownToAndEqual ( cipher, logq ); /* reduce mod 2^logq */
    normalizeAndEqual ( cipher ); /* make sure coeffs are in -q/2..q/2-1 instead of 0..q-1 */

    /* MN:
     * ModRaise = Read ciphertext in a larger domain.
     * This requires the operation "mod q" later on,
     * which is the essential part of the bootstrapping process.
     *
     * With the later "mod q" in mind this has no "real" effect on cipher.nu, cipher.B!
     * Well, actually, it completely destroys it if referring to the plaintext m,
     * but it keeps it referring to q*I+m for some unknown integer I.
     */
    cipher.logq = logQ;
    /* MN:
     * The effect of modifying cipher.logp is that
     * the corresponding plaintext is divided by 2^(logp+4-cipher.logp).
     */
    cipher.nu *= exp2 ( cipher.logp- ( logq+4 ) );
    cipher.B *= exp2 ( cipher.logp- ( logq+4 ) );
    cipher.logp = logq + 4;  // This divides by 2^(logq+4-cipher.logp) wo losing precision.
    /* MN:
     * If not all slots are used repeat values... TODO: Is that correct?
     */
    Ciphertext rot;
    for ( long i = logSlots; i < logNh; ++i ) {
        leftRotateFast ( rot, cipher, ( 1 << i ) );
        addAndEqual ( cipher, rot );
    }
    /* MN:
     * Divide by 2^logNh.
     */
    divByPo2AndEqual ( cipher, logNh ); // bitDown: context.logNh - logSlots
    /*
     * MN:
     * CoeffToSlot: linear transform, puts coefficients into slots
     */
    coeffToSlotAndEqual ( cipher );
    /* MN:
     * ExpEval: evaluate exp( cipher ) homomorphically
     * by approximating exp to degree 8 and
     * using logT halvings/squarings to improve the domain.
     */
    evalExpAndEqual ( cipher, logT, logI ); // bitDown: context.logNh + (logI + logT + 5) * logq + (logI + logT + 6) * logI + logT + 1
    /* MN:
     * SlotToCoeff: linear transform, puts slots into coefficients
     */
    slotToCoeffAndEqual ( cipher );
    /* MN:
     * Multiply by 2^(cipher.logp-logp) to get binary point back in place.
     */
    cipher.nu *= exp2 ( cipher.logp-logp );
    cipher.B *= exp2 ( cipher.logp-logp );
    cipher.logp = logp;

#ifdef CIPHERTEXT_EXTENDED
    /* MN:
     * Estimate the error and the bootstrapped value.
     */
    // MN2020/02/05: TODO [chehan18], §5.3, describes an error bound. Use that.
    // MN2020/02/06: guessing (not yet ok in all situations, in particular if initial logp>=45 newB is too small):
    double newB = B // incoming error
                  + Brs * exp2 ( ( logq+4-logp )-cipher.logp ) // evalExp error
                  + Brs * exp2 ( logN / 2-cipher.logp ) // slotToCoeff error
                  ;
#ifndef NDEBUG
    // MN2020/02/07: Output the error/bound values following the individual operations and
    // compare to the theoretical value in newB.
    cout << "After bootstrapping iteratedly found: "
         << "B=" << cipher.B << "=2^" << log2 ( cipher.B ) << ", "
         << "nu=" << cipher.nu << ";" << endl;
    cout << "        while guessed error bound: "
         << "newB=" << newB << "=2^" << log2 ( newB ) << "=B*2^" << log2 ( newB/cipher.B ) << "." << endl;
#endif
    if ( 0 ) // Under which condition is newB correct?  How to correct?
        if ( newB<cipher.B ) cipher.B = newB; // use the iterated error or the theoretical bound whichever is better
    cipher.nu = nu + cipher.B; /* MN2020/02/11: output should be bounded same as input + err */
#endif
}
