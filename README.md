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

## Parameters of the BLS12-381 (M-Type):

- **x** = `- 0xd201000000010000`
- **p** = `0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab`
- **r** = `0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001`
- **b** = 4
- **btw** = 4*(1+u)
- **G1 cofactor h1** = `0xd201000000010001` (simplified cofactor trick according to [https://eprint.iacr.org/2019/403.pdf](https://eprint.iacr.org/2019/403.pdf))
- **G1/G2 Isogeny order** (for SWU mapping) = 11 / 3

## Parameters of the BLS12-461 (M-Type):

- **x** = `- 0x1ffffffbfffe00000000`
- **p** = `0x15555545554d5a555a55d69414935fbd6f1e32d8bacca47b14848b42a8dffa5c1cc00f26aa91557f00400020000555554aaaaaac0000aaaaaaab`
- **r** = `0xffffff7fffc0180017fe05fd000e801fc017ffc80001100007fefffeffffc0000000000000001`
- **b** = 4
- **btw** = 4*(1+u)
- **G1 cofactor h1** = `0x1ffffffbfffe00000001` (simplified cofactor trick according to [https://eprint.iacr.org/2019/403.pdf](https://eprint.iacr.org/2019/403.pdf))
- **G1/G2 Isogeny order** (for SWU mapping) = 7 / 3

For a full description of all related parameters, see the file (`pairings/parameters/paramlist.rs`).

