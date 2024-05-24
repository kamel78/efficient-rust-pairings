# rust-pairings
Pairings-based cryptography with rust: a constant-time and efficient implementation with rust for BLS12, BLS24 and BLS48 curves. 

The current library implements pairing computations for various security levels: 128-bit, 192-bit, and 256-bit, respectively, on BLS12, BLS24, and BLS48 curves. According to the standardization draft at [1], BLS curves are the most suitable for pairings-based cryptography due to several optimization and security benefits. We have implemented a more optimal set of parameters, particularly for the BLS48 curve, where constant-time hashing to G2 is feasible thanks to the existence of prime-order isogenies that enable the SWU mapping for all implemented curves on both G1 and G2 sub-groups. Detailed information is provided below.

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

- **BLS12-381**: Widely considered the standard for the 128-bit security level and extensively used by most existing libraries and implementations of pairing-based cryptography. However, recent attacks [reference] suggest that it does not truly provide 128-bit security.
- **BLS12-461**: Recommended for a guaranteed 128-bit security level.

For further details and discussions, refer to [2].

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

The curve on GF(p) is defined by `y² = x³ + b`, and the twist curve on GF(p⁸) is defined by `y² = x³ + btw`. Only one set of parameters for the BLS48 is implemented, namely the BLS12-575. This is a curve we discovered using exhaustive parameters set with respect to several optimization criteria, in particular; the existence of isogenies of prime order to hash in constant-time to G1 and G2 using SWU mapping. Unfortunately and surprisingly, the BLS48 introduced in the normalization draft [xx] does not have an isogeny for G2 at all, and its tower construction is not adequate for optimal implementation (due to the tower construction, and much more …..). 

#### Parameters of the BLS24-575 (M-Type):

- **x** = `- 0xff800801`
- **p** = `0x526222098be0d6809ec22a43ab79d1f2120ecbcfbab8934e6c19b1b30a442eea74c92f3ed3940030dbe05f531cd41c69464c793821806366ec094551c2828e22d811dcc31beb02ab`
- **r** = `0xf81e36cde28fad72337ee7497988c7a518b481d35f9cca5a5cbb9c453c85b3ef09f470219adad65475809e468f2faa5490a999c9525ef3f601cd211813004001`
- **b** = `1`
- **btw** = `w`
- **G1 cofactor h1** = `0xff800802`
- **G1/G2 Isogeny order** (for SWU mapping) = 2 / 5

For a full description of all related parameters, see the file (`pairings/parameters/paramlist.rs`).

## Constant-time Implementation

In order to avoid several side-channel attacks, a constant-time implementation of almost all operations is provided, especially for:

1. **Hashing to Elliptic Curves Sub-groups (for G1 and G2) using Simplified Shallue-van de Woestijne-Ulas Method (Simplified SWU for AB = 0)**: This method is commonly used for hashing to G1 in several existing implementations. However, we have developed a special code to find corresponding isogenies for implemented BLS12, BLS24, and BLS48 curves (results are listed in `pairings/parameters/paramlist.rs`).

2. **GLV-multiplication for Torsion Points on G1**: Here, we introduced a personalized optimization that halves the size of the precomputed look-up table (implemented using sliding window size w=3). The approach provides both constant-time and high-speed hashing to G1.

3. **GLS-multiplication on G2**: Implemented using constant-time scalars recoding and optimized cofactor-cleaning based on works in [https://eprint.iacr.org/2013/458.pdf](https://eprint.iacr.org/2013/458.pdf). We implemented GLS-4, GLS-8, and GLS-16 respectively for scalar multiplication of torsion points from G2 on BLS12, BLS24, and BLS48.

4. **Fast Cofactor-cleaning using Decompositions of h2**: The cofactor of the torsion sub-group G2 is decomposed using the method from [https://ia.cr/2017/419](https://ia.cr/2017/419) (Budroni-Pintore). We introduced a minor contribution that further speeds up this operator using an "inverted-endomorphism computation" (maybe will be published later…).

5. **Constant-time Multiplication for Points from E(Fp2), E(Fp4), and E(Fp8)**: Using w-sized sliding window, a fast and stable approach is used as a replacement for the "Montgomery-Ladder" one.

6. **Scalar's Recoding**: Implemented separately in a constant-time way for all implemented multiplications on all curves. (see `pairings/tools/recoders.rs` for implementation details).



