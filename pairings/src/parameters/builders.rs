// Code developed by FARAOUN Kamel Mohamed.
// EEDIS-Laboratory. UDL-University. Algeria
// During May 2024.

use num_bigint::{BigInt, BigUint, ToBigInt};
use num_traits::{Num, One, ToPrimitive, Zero};
use crate::curves::curve_arithmetics::EcPoint;
use crate::curves::g1::{G1SwuIsogeniesConsts,G1Consts};
use crate::curves::g2::{G2Consts, G2SwuIsogeniesConsts};
use crate::extensions::g2_extfields::{ExtFieldG2Element, ExtG2Field};
use crate::fields::prime_fields::{FieldElement,FieldParams};
use crate::tools::arithmetic_interface::ArithmeticOperations;
use super::paramlist::CurveParams;
use crate::{extensions::ext_fields::ExFieldConsts, fields::prime_fields::PrimeField};

pub fn build_field_params<'a ,const NUMLIMBS:usize>(inputfield : &str) -> FieldParams<NUMLIMBS>
            {   
                fn mod_inverse(a: &BigInt, m: &BigInt) -> Option<BigInt> {  // Extended Euclidean Algorithm for finding modular inverse
                    let (mut old_r, mut r) = (a.clone(), m.clone());
                    let (mut old_s, mut s) = (BigInt::one(), BigInt::zero());
                    let (mut old_t, mut t) = (BigInt::zero(), BigInt::one());
                
                    while !r.is_zero() {
                        let quotient = &old_r / &r;
                        let temp_r = r.clone();
                        r = &old_r - &quotient * &r;
                        old_r = temp_r;
                
                        let temp_s = s.clone();
                        s = &old_s - &quotient * &s;
                        old_s = temp_s;
                
                        let temp_t = t.clone();
                        t = &old_t - &quotient * &t;
                        old_t = temp_t;
                    }
                
                    if old_r == BigInt::one() {
                        let mut inverse = old_s % m.clone();
                        if inverse < BigInt::zero() { inverse = inverse + m.clone();
                                                     }
                        Some(inverse)
                    } else {
                        None // No modular inverse exists
                    }}
                
                let mut mod_on_limbs = [0;NUMLIMBS];
                let bigint_modulo:BigInt = BigInt::from_str_radix(&inputfield[2..], 16).unwrap();                
                let mask           = (BigInt::one() << 64) - BigInt::one();                
                let mut num_bits      = bigint_modulo.bits();
                num_bits = ((num_bits / 64) +  ((num_bits %64 !=0) as u64)) * 64;                    
                let mut one = [0;NUMLIMBS];
                let mut rsquare = [0;NUMLIMBS];
                let mut modplus1div4 =[0;NUMLIMBS];
                let mut inv2 = [0;NUMLIMBS];
                let mut threshold = [0;NUMLIMBS];
                let mut tonelli_s = [0;NUMLIMBS];
                let mut tonelli_g = [0;NUMLIMBS];
                let r  = (BigInt::one() << num_bits) % &bigint_modulo;
                let r2 = (BigInt::one() << num_bits).pow(2) % &bigint_modulo;
                let mp1d4 = (&bigint_modulo + BigInt::one()) >> 2;
                let in2= (mod_inverse(&(BigInt::one()+BigInt::one()),&bigint_modulo).unwrap()* &r) % &bigint_modulo;
                let qp = mod_inverse(&(&(BigInt::one() << num_bits)-&bigint_modulo),&(BigInt::one() << num_bits)).unwrap();
                let th = (&bigint_modulo - BigUint::one().to_bigint().unwrap()) >> 1;
                for i in 0..NUMLIMBS {   mod_on_limbs[i] = (((&bigint_modulo & (&mask << (i * 64))) >> (i * 64)) as BigInt).to_u64().unwrap();
                                                one[i] = (((&r & (&mask << (i * 64))) >> (i * 64)) as BigInt).to_u64().unwrap();
                                                rsquare[i] = (((&r2 & (&mask << (i * 64))) >> (i * 64)) as BigInt).to_u64().unwrap();
                                                modplus1div4[i] = (((&mp1d4 & (&mask << (i * 64))) >> (i * 64)) as BigInt).to_u64().unwrap();
                                                inv2[i] = (((&in2 & (&mask << (i * 64))) >> (i * 64)) as BigInt).to_u64().unwrap();
                                                threshold[i] = (((&th & (&mask << (i * 64))) >> (i * 64)) as BigInt).to_u64().unwrap();
                                                };          
                let qprime = (&qp  & &mask).to_u128().unwrap();
                let mut sqrtid: i16 = if (&bigint_modulo % 4u8).to_u8().unwrap() == 3 {-1} else {1};
                if sqrtid ==1 { let mut s = (bigint_modulo.clone()-BigInt::one()).to_biguint().unwrap();
                                sqrtid = 0;
                                while s.bit(0) == false {s=s >>1; sqrtid = sqrtid +1}                                
                                let mut n = BigUint::one()+BigUint::one();
                                let big_modulo_uint =bigint_modulo.to_biguint().unwrap();
                                let big_m1_div2 = (&big_modulo_uint-BigUint::one()) >> 1;
                                while n.modpow(&big_m1_div2, &big_modulo_uint) != &big_modulo_uint-BigUint::one() {n = n + BigUint::one()};
                                let g = n.modpow(&s, &big_modulo_uint);
                                for i in 0..NUMLIMBS 
                                    {tonelli_s[i] = (((&((&s-BigUint::one())>>1) & (&(mask.to_biguint().unwrap()) << (i * 64))) >> (i * 64)) ).to_u64().unwrap();
                                     tonelli_g[i] = (((&g & (&(mask.to_biguint().unwrap()) << (i * 64))) >> (i * 64)) ).to_u64().unwrap();
                                     }
                                }
                let optimize_sqr = mod_on_limbs[NUMLIMBS-1]< ((1<<62)-1);
                let optimize_mul =  (mod_on_limbs[NUMLIMBS-1] & (1 <<63)) == 0;                
                FieldParams{    numlimbs: NUMLIMBS,
                                modulo:   mod_on_limbs,
                                modulo_as_strhex:inputfield,
                                one,
                                rsquare,
                                modplus1div4,
                                inv2 ,
                                qprime ,
                                sqrtid ,
                                zero:[0;NUMLIMBS],
                                optimize_sqr,
                                optimize_mul,
                                sig_theshold : threshold,
                                num_of_bits : bigint_modulo.bits() as usize,
                                tonelli_params : [tonelli_g ,tonelli_s]
                            }
            }

pub fn build_frob_params<'a ,const NUMLIMBS:usize,const NUMFROBCONSTS :usize, const MAX_ISO_COEFS_G1:usize ,const G2_ORDER:usize,const MAX_ISO_COEFS_G2:usize>
        (input : CurveParams<'a ,NUMFROBCONSTS, MAX_ISO_COEFS_G1 ,G2_ORDER, MAX_ISO_COEFS_G2>,field:&PrimeField<NUMLIMBS>) -> ExFieldConsts<NUMFROBCONSTS,NUMLIMBS>
{    
    let mut result =[field.from_str("0");NUMFROBCONSTS];
    for i in 0..NUMFROBCONSTS  { result[i] =field.from_hex_str(&input.frob_consts[i]); }
    ExFieldConsts {frobinus_consts : result, u: input.u, base_qnr :field.from_str(&input.qnr.to_string()) , fb_id:input.frob_id}
}

pub fn build_g1_params<'a ,const NUMLIMBS:usize,const NUMFROBCONSTS :usize, const MAX_ISO_COEFS_G1:usize ,const G2_ORDER:usize,const MAX_ISO_COEFS_G2:usize,const R:usize>
                        (input : CurveParams<'a ,NUMFROBCONSTS, MAX_ISO_COEFS_G1 ,G2_ORDER, MAX_ISO_COEFS_G2>,base_field:&PrimeField<NUMLIMBS>,
                                                                                           r_field:&PrimeField<R>) -> G1Consts<R,NUMLIMBS,MAX_ISO_COEFS_G1>
    {   let mut xnum : [FieldElement<NUMLIMBS>;MAX_ISO_COEFS_G1] = [base_field.one() ;MAX_ISO_COEFS_G1]; 
        let mut xden : [FieldElement<NUMLIMBS>;MAX_ISO_COEFS_G1] = [base_field.one();MAX_ISO_COEFS_G1]; 
        let mut ynum : [FieldElement<NUMLIMBS>;MAX_ISO_COEFS_G1] = [base_field.one();MAX_ISO_COEFS_G1]; 
        let mut yden : [FieldElement<NUMLIMBS>;MAX_ISO_COEFS_G1] = [base_field.one();MAX_ISO_COEFS_G1]; 
        for i in 0..MAX_ISO_COEFS_G1 { xnum[i] = base_field.from_hex_str(&input.xnum_g1[i]); 
                                           xden[i] = base_field.from_hex_str(&input.xden_g1[i]); 
                                           ynum[i] = base_field.from_hex_str(&input.ynum_g1[i]);
                                           yden[i] = base_field.from_hex_str(&input.yden_g1[i]);
                                         }
        G1Consts {  b: base_field.from_hex_str(&input.b),
                    a: base_field.from_hex_str(&input.a),
                    h1: input.h1,
                    w: base_field.from_hex_str(&input.w),
                    lambda: r_field.from_hex_str(&input.lambda),
                    security_level : input.security_level,
                    base_field_numbits : base_field.modulo_as_bigint.bits() as usize,
                    lambda_big: BigUint::from_str_radix(&input.lambda[2..], 16).unwrap(),                
                    swu_consts : G1SwuIsogeniesConsts {     z :  base_field.from_hex_str(&input.swu_z_g1),
                                                            swu_a : base_field.from_hex_str(&input.swu_a_g1),
                                                            swu_b : base_field.from_hex_str(&input.swu_b_g1),
                                                            xnum  ,
                                                            xden  ,
                                                            ynum  ,
                                                            yden  ,
                                                            inv_z : base_field.from_hex_str(&input.swu_z_g1).invert().negate(),
                                                            b_div_a : base_field.from_hex_str(&input.swu_b_g1).multiply(&base_field.from_hex_str(&input.swu_a_g1).invert()).negate()
                                                        },                              
                    default_generator :EcPoint {    x: base_field.from_hex_str(&input.g1_generator[0]),
                                                    y: base_field.from_hex_str(&input.g1_generator[1]),
                                                    z: base_field.from_hex_str(&"0x1")
                                                }
              }
}

pub fn build_g2_params<'a ,const PRAMASIZE:usize, const NUMLIMBS:usize,const NUMFROBCONSTS :usize, const MAX_ISO_COEFS_G1:usize ,const G2_ORDER:usize,const MAX_ISO_COEFS_G2:usize,const R:usize>
                        (   input : CurveParams<'a ,NUMFROBCONSTS, MAX_ISO_COEFS_G1,G2_ORDER, MAX_ISO_COEFS_G2>,
                            base_field:&ExtG2Field<NUMLIMBS,PRAMASIZE>,
                            r_field:&PrimeField<R>,
                            extfieldparams : &ExFieldConsts<PRAMASIZE, NUMLIMBS> ) -> G2Consts<PRAMASIZE,R,NUMLIMBS,MAX_ISO_COEFS_G2>
{       
    let mut xnum : [ExtFieldG2Element<NUMLIMBS,PRAMASIZE>;MAX_ISO_COEFS_G2] = [base_field.one();MAX_ISO_COEFS_G2]; 
    let mut xden : [ExtFieldG2Element<NUMLIMBS,PRAMASIZE>;MAX_ISO_COEFS_G2] = [base_field.one();MAX_ISO_COEFS_G2]; 
    let mut ynum : [ExtFieldG2Element<NUMLIMBS,PRAMASIZE>;MAX_ISO_COEFS_G2] = [base_field.one();MAX_ISO_COEFS_G2]; 
    let mut yden : [ExtFieldG2Element<NUMLIMBS,PRAMASIZE>;MAX_ISO_COEFS_G2] = [base_field.one();MAX_ISO_COEFS_G2]; 
    for i in 0..MAX_ISO_COEFS_G2 { xnum[i] = base_field.from_hex_strings(&input.xnum_g2[i]); 
                                          xden[i] = base_field.from_hex_strings(&input.xden_g2[i]); 
                                          ynum[i] = base_field.from_hex_strings(&input.ynum_g2[i]);
                                          yden[i] = base_field.from_hex_strings(&input.yden_g2[i]);
                                        }       
    G2Consts {  b: base_field.from_hex_strings(&input.btw),
                a: base_field.from_hex_strings(&input.atw),                
                lambda: r_field.from_hex_str(&input.lambda),
                security_level : input.security_level,
                extfieldparams : extfieldparams.clone(),
                base_field_numbits : base_field.basefield().modulo_as_bigint.bits() as usize,
                lambda_big: BigUint::from_str_radix(&input.lambda[2..], 16).unwrap(),    
                u : input.u,                
                swu_consts : G2SwuIsogeniesConsts {     z :  base_field.from_hex_strings(&input.swu_z_g2),
                                                        swu_a : base_field.from_hex_strings(&input.swu_a_g2),
                                                        swu_b : base_field.from_hex_strings(&input.swu_b_g2),
                                                        xnum  ,
                                                        xden  ,
                                                        ynum  ,
                                                        yden  ,
                                                        inv_z : base_field.from_hex_strings(&input.swu_z_g2).invert().negate(),
                                                        j_inv_z : base_field.from_hex_strings(&input.swu_z_g2).invert().negate().substract(&base_field.one()).invert(),
                                                        b_div_a : base_field.from_hex_strings(&input.swu_b_g2).multiply(&base_field.from_hex_strings(&input.swu_a_g2).invert()).negate()
                                                    },
                twist_type: input.twist_type,
                order: G2_ORDER,
                default_generator :EcPoint {
                    x: base_field.from_hex_strings(&input.g2_generator[0]),
                    y: base_field.from_hex_strings(&input.g2_generator[1]),
                    z: base_field.one(),
                }
            }
            
}                                                                                           
