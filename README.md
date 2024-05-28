# Efficient and constant-time Rust implementation of pairings on BLS12, BLS24 and BLS48 curves.


Pairings-based cryptography with rust: a constant-time and efficient implementation with rust for BLS12, BLS24 and BLS48 curves. 

The current library implements pairing computations for various security levels: 128-bit, 192-bit, and 256-bit, respectively, on BLS12, BLS24, and BLS48 curves. According to the standardization draft at [[1]](https://datatracker.ietf.org/doc/draft-irtf-cfrg-pairing-friendly-curves/), BLS curves are the most suitable for pairings-based cryptography due to several optimization and security benefits. We have implemented a more optimal set of parameters, particularly for the BLS48 curve, where constant-time hashing to G2 is feasible thanks to the existence of prime-order isogenies that enable the SWU mapping for all implemented curves on both G1 and G2 sub-groups. Detailed information is provided below.

## Implemented Curves

### 1. BLS12 Configuration

BLS12 curves are implemented based on the following Fp12 tower construction:

- **Fp** = GF(p)
- **GF(p²)** = Fp2<Fp, u² + 1>
- **GF(p⁶)** = Fp6<Fp2, v³ - (u + 1)>
- **GF(p¹²)** = Fp12<Fp2, w² - v>

This construction is chosen due to several optimization techniques that leverage -1 as a quadratic non-residue (QNR) for Fp and (1 + u) as a QNR for Fp2, enabling efficient butterfly multiplication.

The field GF(p) is defined by the prime p = 1/3 (x - 1)²(x⁴ - x² + 1) + x. The torsion subgroup sizes on both G1 and G2 is r = x⁴ - x² + 1, with x being a parameter of the curve. The curve on G1 is defined by y² = x³ + b, and the twist curve on G2 = GF(p²)[r] is defined by y² = x³ + b * w.

Two different sets of parameters for BLS12 are implemented:

- **BLS12-381**: Widely considered the standard for the 128-bit security level and extensively used by most existing libraries and implementations of pairing-based cryptography. However, recent attacks [[3]](https://eprint.iacr.org/2017/334.pdf) suggest that it does not truly provide 128-bit security.
- **BLS12-461**: Recommended for a guaranteed 128-bit security level.

For further details and discussions, refer to [[2]](https://github.com/zcash/zcash/issues/4065).

#### Parameters of the BLS12-381 (M-Type):

- **x** = `- 0xd201000000010000`
- **p** = `0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab`
- **r** = `0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001`
- **b** = 4
- **btw** = 4*(1+u)
- **G1 cofactor h1** = `0xd201000000010001` (simplified cofactor trick according to [https://eprint.iacr.org/2019/403.pdf](https://eprint.iacr.org/2019/403.pdf))
- **G1/G2 Isogeny order** (for SWU mapping) = 11 / 3

#### Parameters of the BLS12-461 (M-Type):

- **x** = `- 0x1ffffffbfffe00000000`
- **p** = `0x15555545554d5a555a55d69414935fbd6f1e32d8bacca47b14848b42a8dffa5c1cc00f26aa91557f00400020000555554aaaaaac0000aaaaaaab`
- **r** = `0xffffff7fffc0180017fe05fd000e801fc017ffc80001100007fefffeffffc0000000000000001`
- **b** = 4
- **btw** = 4*(1+u)
- **G1 cofactor h1** = `0x1ffffffbfffe00000001` (simplified cofactor trick according to [https://eprint.iacr.org/2019/403.pdf](https://eprint.iacr.org/2019/403.pdf))
- **G1/G2 Isogeny order** (for SWU mapping) = 7 / 3

For a full description of all related parameters, see the file (`pairings/parameters/paramlist.rs`).

### 2. BLS24 Configuration

BLS24 curves are implemented based on the following Fp24 tower construction:

- **Fp** = GF(p)
- **GF(p²)** = Fp2<Fp, u² + 1>
- **GF(p⁴)** = Fp4<Fp2, v² - (u + 1)>
- **GF(p⁸)** = Fp8<Fp4, w² - v>
- **GF(p²⁴)** = Fp24<Fp8, z³ - w>

The prime field GF(p) is defined by the prime `p = 1/3 (x-1)²(x⁸-x⁴+1)+x`, while the torsion subgroup size on both G1 and G2 is `r = x⁸ - x⁴ + 1`, with `x` being a parameter of the curve.

The curve on GF(p) is defined by `y² = x³ + b`, and the twist curve on GF(p⁴) is defined by `y² = x³ + btw`. Two different sets of parameters for the BLS24 are implemented, namely the BLS24-479 and the BLS24-559 (we recommend the BLS24-559 for a guaranteed 192-bit security level, depending on whether we consider Razvan’s attack).

#### Parameters of the BLS24-479 (D-Type):

- **x** = `- 0xfffff8f80001`
- **p** = `0x55553de5583ee8de07a5bd8d1f2b1d181350b13b2aaac13ea6db72c4cb7c8eaf4d23c01812428a6b82580149b426a18f5c662a623ec8036a82d2aaab`
- **r** = `0xffffc7c0057046b26bdc4ce1022439d3586c9e35253fdc6bb4c486b2c06e624fa8dcae61a5547eae9c43a57fe3e00001`
- **b** = 1
- **btw** = `0x2aaa9ef2ac1f746f03d2dec68f958e8c09a8589d9555609f536db96265be4757a691e00c09214535c12c00a4da1350c7ae3315311f6401b541695556*(1+u)*v`
- **G1 cofactor h1** = `0xfffff8f80002`
- **G1/G2 Isogeny order** (for SWU mapping) = 2 / 7

#### Parameters of the BLS24-559 (D-Type):

- **x** = `- 0x10007fffffffe40`
- **p** = `0x557003c04ffe89db51ef86c56921091e1c72d1e92397d41e37c255cd44eff051c38fa813a695c5d574923beee353ada1100a0b5fb2111a277d389362d77fcb581fe8b50105eb`
- **r** = `0x1004007006ff65d0bb6ec2d37308240156364161952ef8fb30d7e88dc80a18eb216a645887544b2fd0d47cff70915710376c0fff69f000001`
- **b** = `-2`
- **btw** = `(-1+u)*v`
- **G1 cofactor h1** = `0x10007fffffffe41`
- **G1/G2 Isogeny order** (for SWU mapping) = 3 / 2

### 1.3. BLS48 Configuration

BLS48 curves are implemented based on the following Fp48 tower construction:

- **Fp** = GF(p)
- **GF(p²)** = Fp2<Fp, u² + 1>
- **GF(p⁴)** = Fp4<Fp2, v² - (u + 1)>
- **GF(p⁸)** = Fp4<Fp4, w² - v>
- **GF(p²⁴)** = Fp24<Fp8, z³ - w>
- **GF(p⁴⁸)** = Fp48<Fp24, t² - z>

The prime field GF(p) is defined by the prime `p = 1/3 (x-1)²(x¹⁶-x⁸+1)+x`, while the torsion subgroup size on both G1 and G2 is `r = x¹⁶ - x⁸ + 1`, for a given x defined as a parameter of the curve.

The curve on GF(p) is defined by `y² = x³ + b`, and the twist curve on GF(p⁸) is defined by `y² = x³ + btw`. Only one set of parameters for the BLS48 is implemented, namely the BLS12-575. This is a curve we discovered using exhaustive parameters search, with respect to several optimization criteria, in particular; the existence of isogenies of prime order to hash in constant-time to G1 and G2 using SWU mapping. Unfortunately and surprisingly, the BLS48 introduced in the standarization draft [[1]](https://datatracker.ietf.org/doc/draft-irtf-cfrg-pairing-friendly-curves/) does not have an isogeny for G2 at all, and its tower construction is not adequate for optimal implementation (due to the tower construction, and much more …..). 

#### Parameters of the BLS48-575 (M-Type):

- **x** = `-0xff800801`
- **p** = `0x526222098be0d6809ec22a43ab79d1f2120ecbcfbab8934e6c19b1b30a442eea74c92f3ed3940030dbe05f531cd41c69464c793821806366ec094551c2828e22d811dcc31beb02ab`
- **r** = `0xf81e36cde28fad72337ee7497988c7a518b481d35f9cca5a5cbb9c453c85b3ef09f470219adad65475809e468f2faa5490a999c9525ef3f601cd211813004001`
- **b** = `1`
- **btw** = `w`
- **G1 cofactor h1** = `0xff800802`
- **G1/G2 Isogeny order** (for SWU mapping) = 2 / 5

For a full description of all related parameters, see the file (`pairings/parameters/paramlist.rs`).

## Constant-time Implementation

In order to avoid several side-channel attacks, a constant-time implementation of almost all operations is provided, especially for:

1. **Hashing to Elliptic Curves Sub-groups (for G1 and G2) using Simplified Shallue-van de Woestijne-Ulas Method (Simplified SWU for AB = 0)**[[4]](https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-hash-to-curve-05#section-6.6.3): This method is commonly used for hashing to G1 in several existing implementations. However, we have developed a special code to find corresponding isogenies for implemented BLS12, BLS24, and BLS48 curves (results are listed in `pairings/parameters/paramlist.rs`). To our knowledge, this is the first implementation of SWU on G2 for pairings-friendly curves.

    To explore alternative addition methods, we implemented a new mapping scheme proposed by Dmitry Khovratovich in [12](https://link.springer.com/article/10.1007/s10623-022-01012-8) for performing hashing to G1 on BLS curves. We optimized the proposed algorithm and managed to produce a constant-time version. However, the implementation is only valid for the BLS12-361 curve since it is conditioned by (p mod 9 = 4) and (p mod 27 = 1) (a condition not met by all other implemented curves). Additionally, the obtained performance is not very interesting compared to SWU, particularly due to the required long exponentiation.

3. **GLV-multiplication for Torsion Points on G1**: Here, we introduced a personalized optimization that halves the size of the precomputed look-up table (implemented using sliding window size w=3). The approach provides both constant-time and high-speed hashing to G1.

4. **GLS-multiplication on G2**: Implemented using constant-time scalars recoding and optimized cofactor-cleaning based on works in [[5]](https://eprint.iacr.org/2013/458.pdf). We implemented GLS-4, GLS-8, and GLS-16 respectively for scalar multiplication of torsion points from G2 on BLS12, BLS24, and BLS48.

5. **Fast Cofactor-cleaning using Decompositions of h2**: The cofactor of the torsion sub-group G2 is decomposed using the method from [[6]](https://ia.cr/2017/419) (Budroni-Pintore). We introduced a minor contribution that further speeds up this operator using an "inverted-endomorphism computation" (maybe will be published later…).

6. **Constant-time Multiplication for Points from E(Fp2), E(Fp4), and E(Fp8)**: Using w-sized sliding window, a fast and stable approach is used as a replacement for the "Montgomery-Ladder" one.

7. **Scalar's Recoding**: Implemented separately in a constant-time way for all involved scalar multiplications, on all curve's sub-groups. (see `pairings/tools/recoders.rs` for implementation details).

## Points and Field Elements Representation

Standard point compression and serialization are utilized to enable point encoding according to the standards defined in the ZCash serialization format [[7]](https://www.ietf.org/archive/id/draft-irtf-cfrg-pairing-friendly-curves-11.html#name-zcash-serialization-format-). The scheme is generalized to all implemented curves for elements from both G1 and G2.

Finite field element representation follows RFC standards [[8]](https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-hash-to-curve-06#name-hashing-to-a-finite-field), while the derivation also adheres to the following scheme ([RFC](https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-hash-to-curve-05#section-4.1)). Both curve points and field elements can be exported/imported in decimal representation, hexadecimal representation, byte array format, and interestingly, base64 encoding format.

## Arithmetic for Fields and Curves

Arithmetic on finite fields (on both GF(p) and GF(r) defined for handling arithmetic on curves and scalars respectively) is implemented using the Montgomery representation to provide fast and efficient computation. Special optimizations of the CIOS reduction approach [[9]](https://www.microsoft.com/en-us/research/wp-content/uploads/1998/06/97Acar.pdf) are implemented from [[10]](https://hackmd.io/@gnark/modular_multiplication#fn1) to ensure optimal arithmetic when conditions are satisfied (for both multiplication and squaring).

Arithmetic on elliptic curves is implemented using Jacobian coordinate systems to avoid inversion on affine coordinates. Optimal point addition and doubling were implemented from [[11]](https://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html#addition-madd-2007-bl).

## Implemented pairings
For implemented BLS curves, Optimal Ate pairings is a refinement of the Ate pairing on elliptic curves, optimized for efficiency. I is considered as best choice for both runtime and security consideration, compared to Weil, Tate and R-Ate variants.

##  Rust Implementation Considerations 

Existing Rust implementations of pairings-based cryptography are generally tailored to specific curves, meaning that the code works only for a given combination of parameters and hence is not extensible. One of the challenges of the proposed implementation is to permit full parametrization of the code for any new curve definition while maintaining the same level of code optimality and runtime performance. This can be achieved by simply injecting new parameters in a textual hexadecimal format (a clear illustration of this can be found in the file `pairings/parameters/paramlist.rs`.

Pairings engines are implemented as constant Rust "structurs" using the OnceCell crate. This enables fast and efficient use of the implemented curves for new implementations of pairings-based protocols and transparent manipulation of related functionalities. As an illustration, the following code demonstrates how fast and simple the implemented library can be used:

```rust
use pairings::{Bls12Curves, PairingsEngine};

fn main() {
    let engine = pairings::BLS12::_381();
    let p = engine.g1.hash_to_field("identity 1", 0);
    let q = engine.g2.hash_to_field("identity 2", 0);
    let a = engine.fr.random_element();
    let b = engine.fr.random_element();
    let e1 = engine.paire(&(a*p), &(b*q));
    let e2 = engine.paire(&(b*p), &(a*q));
    println!("Pairings Verification : {}", e1 == e2);
    println!("Non-degeneracy verification : (e1  != 1): {}, (e2 != 1): {}", e2 != engine.gt.one() ,e1 != engine.gt.one());
}
```
This code snippet demonstrates the simplicity and efficiency of the library's usage for performing pairings-based cryptography operations.

## Performance and Runtime Results

The following table illustrates the obtained performances and runtime results for various implemented operations across all curves. The tests were performed on an Intel Core i7-10700F CPU running at 2.9 GHz, with 16 GB of RAM. Some operations have timing measurements in microseconds (µs) as they take less than one millisecond. Note that for the 256-bit security level (provided by the BLS48 curve), the obtained results are significantly better than the existing ones presented [here](https://github.com/mk-math-kyushu/bls48/blob/master/README.md) using a C++ implementation of BLS48, when using a similar performance architecture.

<table>
  <thead>
    <tr>
      <th></th>
      <th>Pairings</th>
      <th>Miller Loop</th>
      <th>Final Expo.</th>
      <th>Hashing to G1</th> 
      <th>Hashing to G2</th>
      <th>Mul. on G1 (GLV)</th>
      <th>Mul. on G2 (GLS)</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>BLS12-381</th>
      <td align="center">1.65 ms</td>
      <td align="center">680.85 µs</td>
      <td align="center">958.08 µs</td>
      <td align="center">206.83 µs</td>
      <td align="center">481.27 µs</td>
      <td align="center">145.94 µs</td>
      <td align="center">366.56 µs</td>
    </tr>
    <tr>
      <th>BLS12-461</th>
      <td align="center">2.86 ms</td>
      <td align="center">1.22 ms</td>
      <td align="center">1.56 ms</td>
      <td align="center">319.69 µs</td>
      <td align="center">1.35 ms</td>
      <td align="center">262.08 µs</td>
      <td align="center">652.37 µs</td>
    </tr>
    <tr>
      <th>BLS24-479</th>
      <td align="center">8.49 ms</td>
      <td align="center">2.51 ms</td>
      <td align="center">5.99 ms</td>
      <td align="center">334.58 µs</td>
      <td align="center">4.50 ms</td>
      <td align="center">318.45 µs</td>
      <td align="center">2.20 ms</td>
    </tr>
    <tr>
    <th>BLS24-559</th>
      <td align="center">11.46 ms</td>
      <td align="center">3.43 ms</td>
      <td align="center">8.17 ms</td>
      <td align="center">280.86 µs</td>
      <td align="center">6.98 ms</td>
      <td align="center">466.58 µs</td>
      <td align="center">3.16 ms</td>
    </tr>
    <tr>
    <th>BLS48-575</th>
      <td align="center">33.82 ms</td>
      <td align="center">6.03 ms</td>
      <td align="center">27.78 ms</td>
      <td align="center">204.79 µs</td>
      <td align="center">18.12 ms</td>
      <td align="center">519.22 µs</td>
      <td align="center">10.35 ms</td>
    </tr>
  </tbody>
</table>

## Illustrative demonstrations 

The implemented curves provide conversion routines to several representation formats (Decimal, Hexadecimal, Base64, and byte arrays), making I/O integration trivial during protocol implementation. We provide two illustrative implementation code examples for BLS signatures and Identity-Based Encryption, respectively. A full "Demos" directory containing several other implementations of pairings-based cryptography protocols will be uploaded to the repository as soon as possible (feel free to try it and provide comments).


### 1. BLS Signature scheme [13](https://link.springer.com/chapter/10.1007/3-540-45682-1_30)

The BLS (Boneh-Lynn-Shacham) signature scheme is a cryptographic algorithm using pairing-based cryptography on elliptic curves. It utilizes a bilinear pairing `e: G1 × G2 → GT`, where `G1`, `G2`, and `GT` are groups of the same prime order `r`, and follows these steps:

1. **Setup**:
    - Select an elliptic curve `E` and define the groups `G1`, `G2`, and `GT` of order `p`, admitting a bilinear pairing `e: G1 × G2 → GT`.
    - Choose a generator `g ∈ G2`.

2. **Key Generation**:
    - **Private Key (sk)**: Randomly select a private key `x ∈ Fr`.
    - **Public Key (pk)**: Compute the public key `pk = x * g ∈ G2`.

3. **Signing**:
    - Compute `H(m)`, the hash of the message `m` to a point in `G1`.
    - Compute the signature `σ = x * H(m) ∈ G1`.

4. **Verification**:
    - Compute `H(m) ∈ G1` in the same way as above.
    - Verify `e(σ, g) = e(H(m), pk)`.

The following code snippets illustrate how such a scheme can easily be implemented using the present library.

```rust
use pairings::{Bls12Curves, PairingsEngine};

fn main() {    
    // Hashing to elliptic curves support two modes :  
    //      mode = 0 stands for Non-uniform-Encoding, while mode = 1 stands for Random Oracle Model encoding 

    let engine = pairings::BLS12::_461();

    // Key-paire genration :
    let sk = engine.fr.random_element();
    let pk = sk * engine.g2.default_generator();
    println!(" Secrete key (base64) = {}", sk.to_base64());
    println!(" Public Key  (base64) = {}", pk.encode_to_base64());       

    // BLS Signing : 
    let message = "This is a simple message to be signed. A message can be any arbitrary length string ....";
    let hashed_message = engine.g1.hash_to_field(&message, 0);
    let signature = sk * hashed_message;
    println!(" Signatue is (base64): {}", signature.encode_to_base64());

    // BLS Verification :
    let hashed_message = engine.g1.hash_to_field(&message, 0);
    let verification_result = engine.paire(&signature, &engine.g2.default_generator()) == engine.paire(&hashed_message, &pk);
    println!("Verification result : {}",if verification_result {"correct"} else {"incorrect"});

    // Faster way to verify using multi-pairings
    let verification_result = engine.multi_paire(&[signature,hashed_message], &[-engine.g2.default_generator(),pk]) == engine.gt.one();
    println!("Verification result : {}",if verification_result {"correct"} else {"incorrect"});       
}
```
```bash
PS C:\pairings-rust> cargo run
   Compiling pairings-rust v0.1.0 (C:\Users\Kamel\Desktop\new-pairings-rust\pairings-rust)
    Finished `dev` profile [unoptimized + debuginfo] target(s) in 1.78s
     Running `target\debug\pairings-rust.exe`
 Secrete key (base64) = D2byDX9BZ9oGoPrmLrZL6rQQnOR28f7ud8j3FkDzn5O25I1wx8qi
 Public Key  (base64) = p2jFDR75AwFg70jol1BT7v8r+PeEPy5ueGM53fdrFWyM34ipzZffmbf6x/jzALdGmewAPFVinG2URQuHxFuzZcGJVzxkuWB8kfm/Axs8kkHZUGEUNAfGiUp12NK8TUANYN2C9wiQl3gv/jhhuoZSvdWBI2g=
 Signatue is (base64): rOgxEY8EdlN1Kw/9N3lAH77SqckVVQ1CBG32ur09I2TAeFLh5CB7lpXcJgBzw/cuT8tea5adggID7w==
Verification result : correct
Verification result (using multi-pairings) : correct
PS C:\pairings-rust>
```
### 1. Identity-based encryption scheme [14](https://crypto.stanford.edu/~dabo/papers/bfibe.pdf)

The Boneh-Franklin IBE scheme is historicaly the first proposed IBE schemes, and is still widely used. We consider a bilinear pairing `e: G1 × G2 → GT`, where `G1`, `G2`, and `GT` are groups of the same prime order `r`. One of the vriants of this scheme the following mathematical steps based on elliptic curve pairings :

1. **Setup: System Initialization**:
   - The PKG selects a random master secret `s` in `Fr` and computes `MPk = s.g1`.
   - The master public key is `MPk`.
   - The master secret key is `s`.

2. **Key Generation :Private Key Generation for User**:
   - User's identity `ID` (e.g., email address) is hashed to a point on the elliptic curve using a cryptographic hash function : `QID = Hash_to_G1(ID) ∈ G1`.
   - The PKG computes the private key for `ID` as `dID = s.QID ∈ G1`.

3. **Encrypting a Message**:
    When a sender wishes to send a message `M` to a user with identity `ID`, he performs the following :
   - The sender computes `QID = Hash_to_G1(ID)`.
   - The sender selects a random `a ∈Fr' and computes:
     - Ciphertext component `C1 = r.g2 ∈ G2`.
     - Ciphertext component  `C2 = M ⊕ HDK(e(QID, MPk)^r)`, when HDK is a secure key derivation function.
   - The ciphertext is `C = (C1, C2)`.

4. **Decrypting the Ciphertext**:
   - The recipient uses their private key `dID` to compute:
     - `e(dID, C1) = e(s.QID, a.g2) = e(QID, g2)^(s.ar) = e(QID, MPk)^r`.
   - The recipient then computes the message `M` as:
     - `M = C2 ⊕ HDK(e(dID, C1))`.

This process ensures secure communication based on the hardness of the bilinear Diffie-Hellman problem and eliminates the need for a traditional public key infrastructure. 

The following is a code snippet of the Boneh-Franklin IBE's implementation using the implemented library. However, it's essential to note that the following code is just a proof of concept and requires further adjustments if a concrete real application is targeted (including especially IND-CCA updates to the scheme...).
```rust
use base64::{engine::general_purpose, Engine};
use pairings::{Bls12Curves, PairingsEngine};

fn main() {
    let engine = pairings::BLS12::_461();
    // Generation of Master Keys (Setup):
    let msk =  engine.fr.random_element();
    let mpk = msk * engine.g2.default_generator();
    println!("The Master secrete key : {} ",msk.to_base64());
    println!("The Master public key : {} \n",mpk.encode_to_base64());

    // Key extraction : generation of the user's secrete key for corresponding Identity :
    let user_identity ="ID-1";    
    let id_sk = msk * engine.g1.hash_to_field(&user_identity, 0);
    println!("User's secrete key for identity '{}' : {} \n",user_identity,id_sk.encode_to_base64()); 
    
    // Key confirmation : user can confirm the authenticity and corectness of the secrete key like follows: 
    let valide_secrete_key = engine.paire(&id_sk, &engine.g2.default_generator()) 
                                   == engine.paire(&engine.g1.hash_to_field(&user_identity, 0), &mpk);
    println!("User's secrete key confirmation : {} ",if valide_secrete_key {"Valid key\n"} else {"Invalid key\n"}); 

    //  Encryption of a message to the user using its Identity :
    let message ="This is a simple message to be signed. A message can be any arbitrary length string ....";
    println!("Plaintext message : {}\n",message);
    let message_as_bytes: Vec<u8> = message.as_bytes().to_vec();
    let a = engine.fr.random_element();
    let u = a * engine.g2.default_generator();
    let key_stream = engine.paire(&engine.g1.hash_to_field(&user_identity, 0),&mpk)
                              .pow(&a).derive_hkdf(8*message_as_bytes.len(), None);
    let encrypted_data: Vec<u8> = key_stream.iter().zip(message_as_bytes.iter()).map(|(&x1, &x2)| x1 ^ x2).collect();    
    let encrypted_message =[u.encode_to_base64(),general_purpose::STANDARD.encode(encrypted_data)];
    println!("Encrypted message : {:?}\n",encrypted_message);

    // Decryption of the message using the user's secrete key 
    let decoded_encryption = general_purpose::STANDARD.decode(&encrypted_message[1]).unwrap();
    let u = engine.g2.from_base64(&encrypted_message[0]);
    let key_stream = engine.paire(&id_sk, &u).derive_hkdf(8*decoded_encryption.len(), None); 
    let decrypted_message : Vec<u8> = key_stream.iter().zip(decoded_encryption.iter()).map(|(&x1, &x2)| x1 ^ x2).collect();    
    println!("Decrypted message : {}",std::str::from_utf8(&decrypted_message).unwrap());
}
```
```bash
PS C:\pairings-rust> cargo run
   Compiling pairings-rust v0.1.0 (/Users/neptune/Desktop/Rust-Pairings/pairings-rust)
    Finished dev [unoptimized + debuginfo] target(s) in 0.61s
     Running `target/debug/pairings-rust`
The Master secrete key : AmuxOKifYyMz691pN0TWEg9YJ3oRmlYFWeqDC87OAN//DOd7gMSK 
The Master public key : o3cYmzmzTGY+oO7nfCzo0aVbzAQA6fpKvLRxo3C1v7RFbMlA72cxwE3mQHg0fuFKKSfK/uyrmz2lWwUPD+IUqlHwSnPPK7DeHThCFF7QXJcLfLQSgaqEeP0jiUzj4J7qeCw4VfTaBTEROIZwnnJW6xJcHc4= 

User's secrete key for identity 'ID-1' : hnqW0zS1XdisCJOHrJL5v5Kr/F/w06WpPDpDoR2m8L4vZnDlH5Xji21RD0GMHrdDbB1JTjKUe0ogPg== 

User's secrete key confirmation : Valid key
 
Plaintext message : This is a simple message to be signed. A message can be any arbitrary length string ....

Encrypted message : ["rSJQ5OensV32nOuyjE59C6iAJ/8inxOSVV7Dfsw7IEUDYJggBasZq+UAkmkbJhqNw3+6BmtNopSM0xR++csiDqhAw4FeViyBNe3mkydnMShzDQIyXxTmimvcMx2zlLDab66IiHX7VHS8Ud2E1YznxhuvxgI=", "D+KS5i5y+I377rTgTkEpBKp03ITt6VPhFFTDZY/9LrcSYF62JQHC7E3Z3i0yy/ZNBLTesBrimhpNiY9PWvOpIm7aWEL1g+OUz/MA4O4/TzlVIfiekjOOmQ=="]

Decrypted message : This is a simple message to be signed. A message can be any arbitrary length string ....
PS C:\pairings-rust>
```
We plan to include a demonstration folder in the near future, which will contain implementations of several other pairing-based encryption schemes, including: Hierarchical Identity-Based Encryption (HIBE), Proxy Re-encryption, Identity-Based Signature, Signcryption, Attribute-Based Encryption, Searchable Encryption, Zero-Knowledge Proofs, and much more.


Feel free to comment, correct, or propose any additional updates to this code. Your comments are welcome at kamel_mh@yahoo.fr.






FARAOUN Kamel Mohamed.

## References 

[1] . [Pairing-Friendly Curves : draft-irtf-cfrg-pairing-friendly-curves-11.](https://datatracker.ietf.org/doc/draft-irtf-cfrg-pairing-friendly-curves/)

[2] . [Why not BLS12-461 with equivalent performance?](https://github.com/zcash/zcash/issues/4065) 

[3] . [Updating key size estimations for pairings](https://eprint.iacr.org/2017/334.pdf)

[4] . [Simplified Shallue-van de Woestijne-Ulas Method](https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-hash-to-curve-05#section-6.6.3)

[5] . [Efficient hash maps to G2 on BLS curves](https://eprint.iacr.org/2013/458.pdf) 

[6] . [Exponentiating in Pairing Groups](https://ia.cr/2017/419)

[7] . [ZCash serialization format for BLS12_381](https://www.ietf.org/archive/id/draft-irtf-cfrg-pairing-friendly-curves-11.html#name-zcash-serialization-format-)

[8] . [Hashing to a Finite Field](https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-hash-to-curve-06#name-hashing-to-a-finite-field)

[9] . [High-Speed Algorithms & Architectures For Number-Theoretic Cryptosystems](https://www.microsoft.com/en-us/research/wp-content/uploads/1998/06/97Acar.pdf)

[10] . [Faster big-integer modular multiplication for most moduli](https://hackmd.io/@gnark/modular_multiplication#fn1)
  
[11] . [Jacobian coordinates for short Weierstrass curves](https://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html#addition-madd-2007-bl).

[12] . [Koshelev, D. (2022). Indifferentiable hashing to ordinary elliptic F q-curves of j= 0 with the cost of one exponentiation in Fq.](https://link.springer.com/article/10.1007/s10623-022-01012-8) 

[13] . [Short Signatures from the Weil Pairing](https://link.springer.com/chapter/10.1007/3-540-45682-1_30)

[14] . [Identity-Based Encryption from the Weil Pairing](https://crypto.stanford.edu/~dabo/papers/bfibe.pdf)

