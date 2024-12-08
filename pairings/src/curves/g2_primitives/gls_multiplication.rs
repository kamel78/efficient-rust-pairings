// Code developed by FARAOUN Kamel Mohamed.
// EEDIS-Laboratory. UDL-University. Algeria
// During May 2024.

use crate::{curves::{curve_arithmetics::EcPoint, g2::G2Element, g2_primitives::phi::phi_bls48}, 
            fields::prime_fields::FieldElement, tools::{arithmetic_interface::ArithmeticOperations, 
            recoders::{recod_scalar_gls16, recod_scalar_gls4, recod_scalar_gls8}}};
use num_bigint::BigUint;
use num_traits::{Num, One,ToPrimitive};
use std::{ops::BitAnd, str::FromStr};

use super::phi::phi_bls24;

pub fn gls4_multiply_bls12<const PRAMASIZE:usize, const R: usize, const N: usize, const MAX_COEFS_COUNT: usize>
            (input :&G2Element<PRAMASIZE, R,N,MAX_COEFS_COUNT>, scalar :&FieldElement<R>) -> G2Element<PRAMASIZE, R,N,MAX_COEFS_COUNT>
{   
    // Constant-Time multiplication for elements in G2 (4-GLS): GLS Implementation of points multiplication on G2, for the BLS12 Curve
    // Joppe W. Bos, Craig Costello, and Michael Naehrig https://eprint.iacr.org/2013/458.pdf
    
    let ff: BigUint = BigUint::from_str("255").unwrap();
    let mut code = recod_scalar_gls4(&scalar.to_big_uint(),&scalar.negate().to_big_uint(),input.consts.u.abs() as u128);
    let mut p1 = input.phi();
    let p2 = p1.phi();
    let mut p3 = p2.phi();
    let infinit = G2Element{ point : EcPoint {  x:input.point.x.one(), 
                                                                                            y :input.point.x.one(), 
                                                                                            z: input.point.x.zero() }, 
                                                                         consts: input.consts };                                                                                
    if input.consts.u < 0 {   p1 = p1.negate();
                              p3 = p3.negate()};                                                                             
    let mut lookup = [infinit;8];    
    lookup[0] = input.clone();
    lookup[1] = p1.addto(&lookup[0]);
    lookup[2] = p2.addto(&lookup[0]);
    lookup[3] = p2.addto(&lookup[1]);
    lookup[4] = p3.addto(&lookup[0]);
    lookup[5] = p3.addto(&lookup[1]);
    lookup[6] = p3.addto(&lookup[2]);
    lookup[7] = p3.addto(&lookup[3]);   
    let mut limb : u8 = (&code).bitand((ff).clone()).to_u8().unwrap();
    let out_sig : i8 =  1 - ((limb & 1) << 1) as i8;
    code = code >> 1;
    limb = (&code).bitand((ff).clone()).to_u8().unwrap();
    let mut sign : i8 = ((limb & 1) << 1) as i8 - 1 ;
    let mut idx =  ((limb & 15) >> 1) as usize;
    let mut result = if sign == 1 {lookup[idx]} else {lookup[idx].negate()};    
    code = code >> 4;
    while code != BigUint::one() {  result = result.double();
                                    limb = (&code).bitand((ff).clone()).to_u8().unwrap();
                                    sign = ((limb & 1) << 1) as i8- 1 ;
                                    idx =  ((limb & 15) >> 1) as usize;
                                    result = result.addto(& if sign == 1 {lookup[idx]} else {lookup[idx].negate()}); 
                                    code = code >> 4;
                                    }
    if out_sig == 1 {result} else {result.negate()}
}

pub fn gls8_multiply_bls24<const PRAMASIZE:usize, const R: usize, const N: usize, const MAX_COEFS_COUNT: usize>
            (input :&G2Element<PRAMASIZE, R,N,MAX_COEFS_COUNT>, scalar :&FieldElement<R>) -> G2Element<PRAMASIZE, R,N,MAX_COEFS_COUNT>
{   
    //  Constant-Time multiplication for elements in G2 (8-GLS): GLS Implementation of points multiplication on G2, for the BLS24 Curve
    //  Joppe W. Bos, Craig Costello, and Michael Naehrig https://eprint.iacr.org/2013/458.pdf    
    let ff: BigUint = BigUint::from_str("255").unwrap();
    let mut code = recod_scalar_gls8(&scalar.to_big_uint(),
                                               input.consts.u.abs() as u128,
                                               &BigUint::from_str_radix(&input.consts.lambda.fieldparams.modulo_as_strhex[2..],16).unwrap());
    let mut p1 = input.phi();
    let p2 = p1.phi();
    let mut p3 = p2.phi();
    let infinit = G2Element{ point : EcPoint {  x :input.point.x.one(), 
                                                                                            y :input.point.x.one(), 
                                                                                            z :input.point.x.zero() }, 
                                                                         consts: input.consts };                                                                                
    if input.consts.u < 0 { p1 = p1.negate();
                            p3 = p3.negate()};                                                                             
    let mut lookup1: [G2Element<PRAMASIZE, R, N, MAX_COEFS_COUNT>; 8] = [infinit;8];    
    lookup1[0] = input.clone();
    lookup1[1] = p1.addto(&lookup1[0]);
    lookup1[2] = p2.addto(&lookup1[0]);
    lookup1[3] = p2.addto(&lookup1[1]);
    lookup1[4] = p3.addto(&lookup1[0]);
    lookup1[5] = p3.addto(&lookup1[1]);
    lookup1[6] = p3.addto(&lookup1[2]);
    lookup1[7] = p3.addto(&lookup1[3]);
    let lookup2  = {
        let mut arg = [infinit; 8];
        arg.iter_mut()
            .zip(lookup1.iter().map(|x| phi_bls24(x,4)))
            .for_each(|(dest, src)| *dest = src);
        arg
    };    
    let mut limb : u8 = (&code).bitand((ff).clone()).to_u8().unwrap();
    let out_sig : i8 =  1 - ((limb & 1) << 1) as i8;    
    code = code >> 1;
    limb = (&code).bitand((ff).clone()).to_u8().unwrap();    
    let (mut c1, mut c2) = ((limb & 15) as usize , ((limb >> 4) & 15) as usize);
    let (mut sign1 , mut sign2)    = (((c1 & 1) << 1) as i8 - 1, ((c2 & 1) << 1) as i8 - 1);
    let mut result = (if sign1 == 1 {lookup1[c1 >> 1]} else {lookup1[c1 >> 1].negate()})
                                                                 .addto(&if sign2 == 1 {lookup2[c2 >> 1]} else {lookup2[c2 >> 1].negate()});    
    code = code >> 8;
    while code != BigUint::one() {  result = result.double();
                                    limb = (&code).bitand((ff).clone()).to_u8().unwrap();
                                    (c1, c2) =  ((limb & 15) as usize , ((limb >> 4) & 15) as usize);
                                    (sign1, sign2) =  (((c1 & 1) << 1) as i8 - 1, ((c2 & 1) << 1) as i8 - 1);
                                    result = result.addto(&if sign1 == 1 {lookup1[c1 >> 1]} else {lookup1[c1 >> 1].negate()})
                                                           .addto(&if sign2 == 1 {lookup2[c2 >> 1]} else {lookup2[c2 >> 1].negate()}); 
                                    code = code >> 8;
                                    }
    if out_sig == 1 {result} else {result.negate()}
}

pub fn gls16_multiply_bls48<const PRAMASIZE:usize, const R: usize, const N: usize, const MAX_COEFS_COUNT: usize>
            (input :&G2Element<PRAMASIZE, R,N,MAX_COEFS_COUNT>, scalar :&FieldElement<R>) -> G2Element<PRAMASIZE, R,N,MAX_COEFS_COUNT>
{   
    //  Constant-Time multiplication for elements in G2 (16-GLS): GLS Implementation of points multiplication on G2 for the BLS48 Curve    //  Inspired from Joppe W. Bos, Craig Costello, and Michael Naehrig https://eprint.iacr.org/2013/458.pdf    
    let ffff: BigUint = BigUint::from_str("65535").unwrap();
    let mut code = recod_scalar_gls16(&scalar.to_big_uint(),
                                               input.consts.u.abs() as u64,
                                               &BigUint::from_str_radix(&input.consts.lambda.fieldparams.modulo_as_strhex[2..],16).unwrap());                      
    let mut p1 = input.phi();
    let p2 = p1.phi();
    let mut p3 = p2.phi();
    let infinit = G2Element{ point : EcPoint {  x :input.point.x.one(), 
                                                                                            y :input.point.x.one(), 
                                                                                            z :input.point.x.zero() }, 
                                                                         consts: input.consts };                                                                                
    if input.consts.u < 0 { p1 = p1.negate();
                            p3 = p3.negate()};                                                                             
    let mut lookup1: [G2Element<PRAMASIZE, R, N, MAX_COEFS_COUNT>; 8] = [infinit;8];    
    lookup1[0] = input.clone();
    lookup1[1] = p1.addto(&lookup1[0]);
    lookup1[2] = p2.addto(&lookup1[0]);
    lookup1[3] = p2.addto(&lookup1[1]);
    lookup1[4] = p3.addto(&lookup1[0]);
    lookup1[5] = p3.addto(&lookup1[1]);
    lookup1[6] = p3.addto(&lookup1[2]);
    lookup1[7] = p3.addto(&lookup1[3]);    
    let lookup2  = {
        let mut arg = [infinit; 8];
        arg.iter_mut()
            .zip(lookup1.iter().map(|x| phi_bls48(x,4)))
            .for_each(|(dest, src)| *dest = src);
        arg
    };    
    let lookup3  = {
        let mut arg = [infinit; 8];
        arg.iter_mut()
            .zip(lookup2.iter().map(|x| phi_bls48(x,4)))
            .for_each(|(dest, src)| *dest = src);
        arg
    };    
    let lookup4  = {
        let mut arg = [infinit; 8];
        arg.iter_mut()
            .zip(lookup3.iter().map(|x| phi_bls48(x,4)))
            .for_each(|(dest, src)| *dest = src);
        arg
    };        
    let mut limb : u16 = (&code).bitand((ffff).clone()).to_u16().unwrap();
    let out_sig : i8 =  1 - ((limb & 1) << 1) as i8;    
    code = code >> 1;
    limb = (&code).bitand((ffff).clone()).to_u16().unwrap();    
    let (mut c1, mut c2,mut c3,mut c4) = ((limb & 15) as usize , ((limb >> 4) & 15) as usize,
                                                                      ((limb >> 8) & 15) as usize,((limb >> 12) & 15) as usize);
    let (mut sign1 , mut sign2, mut sign3, mut sign4)  = (((c1 & 1) << 1) as i8 - 1, ((c2 & 1) << 1) as i8 - 1,
                                                                          ((c3 & 1) << 1) as i8 - 1,((c4 & 1) << 1) as i8 - 1 );
    let mut result1 = (if sign1 == 1 {lookup1[c1 >> 1]} else {lookup1[c1 >> 1].negate()})
                                                                  .addto(&if sign2 == 1 {lookup2[c2 >> 1]} else {lookup2[c2 >> 1].negate()});    
    let mut result2 = (if sign3 == 1 {lookup3[c3 >> 1]} else {lookup3[c3 >> 1].negate()})
                                                                  .addto(&if sign4 == 1 {lookup4[c4 >> 1]} else {lookup4[c4 >> 1].negate()});    
    let mut result = result1.addto(&result2);
    code = code >> 16;
    while code != BigUint::one() {  result = result.double();
                                    limb = (&code).bitand((ffff).clone()).to_u16().unwrap();
                                    (c1,c2,c3,c4) = ((limb & 15) as usize , ((limb >> 4) & 15) as usize,((limb >> 8) & 15) as usize,((limb >> 12) & 15) as usize);
                                    (sign1,sign2,sign3,sign4) = (((c1 & 1) << 1) as i8 - 1, ((c2 & 1) << 1) as i8 - 1,((c3 & 1) << 1) as i8 - 1,((c4 & 1) << 1) as i8 - 1 );
                                    result1 = (if sign1 == 1 {lookup1[c1 >> 1]} else {lookup1[c1 >> 1].negate()})
                                               .addto(&if sign2 == 1 {lookup2[c2 >> 1]} else {lookup2[c2 >> 1].negate()});    
                                    result2 = (if sign3 == 1 {lookup3[c3 >> 1]} else {lookup3[c3 >> 1].negate()})
                                               .addto(&if sign4 == 1 {lookup4[c4 >> 1]} else {lookup4[c4 >> 1].negate()});    
                                    result = result.addto(&result1.addto(&result2)); 
                                    code = code >> 16;
                                    }
    if out_sig == 1 {result} else {result.negate()}
}