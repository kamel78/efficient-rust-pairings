// Code developed by FARAOUN Kamel Mohamed.
// EEDIS-Laboratory. UDL-University. Algeria
// During May 2024.

use std::ops::BitAnd;
use crate::fields::prime_fields::FieldElement;
use crate::fields::arithmetic;
use num_bigint::{BigUint, ToBigUint};
use num_traits::{FromPrimitive, Zero};

fn get_naf(x: BigUint) -> Vec<i8> {
    let t = &x;
    let xh: &BigUint = &(t >> 1);
    let a = t + xh;
    let b = xh.clone() ^ &a;
    let mut np: BigUint = a.bitand(&b);
    let mut nm = xh.bitand(&b);
    let mut result = Vec::<i8>::new();
    while np > BigUint::zero() {   let _npbit = np.bit(0) as i8;
                                   let _nmbit = nm.bit(0) as i8;
                                   let bit = _npbit -_nmbit;
                                   result.push(bit);
                                   np = np >> 1;
                                   nm = nm >> 1;
                                }
    result.reverse();
    result
}
pub trait Exponent<const N:usize>{
    fn to_u64_array(&self) -> Option<[u64; N]>;
    fn to_naf(&self) -> Vec<i8>; 
    fn get_len(&self)->usize;
}

impl <const N:usize,const R:usize > Exponent<N> for FieldElement<R> {
    
    fn to_u64_array(&self) -> Option<[u64; N]> {
        let mut one :[u64;R] = [0;R];
        one[0] = 1;
        let mut result :[u64;N] = [0;N];
        let from_mont  : [u64;R] = arithmetic::mul(&self.mont_limbs,&one,&self.fieldparams);
        result[0..R].copy_from_slice(&from_mont);
        Some(result)
    }
    
    fn get_len(&self)->usize { R }

    fn to_naf(&self) -> Vec<i8>{
        let mut one :[u64;R] = [0;R];
        one[0] = 1;
        let mut as_big : BigUint = BigUint::zero();    
        let from_mont  : [u64;R] = arithmetic::mul(&self.mont_limbs,&one,&self.fieldparams);
        println!("{:?}",from_mont);
        for i in 0..R {as_big += from_mont[i].to_biguint().unwrap() << (i * 64)};        
        get_naf(as_big)
        }   
}

impl <const N:usize > Exponent<N> for u8 {
    fn to_u64_array(&self) -> Option<[u64; N]> {
        let mut _a=[0;N];
        _a[0] = *self as u64;
        Some(_a)
    }
    
    fn to_naf(&self) -> Vec<i8>{
        get_naf(BigUint::from_u8(*self).unwrap())
        }

    fn get_len(&self)->usize { 1 }
}

impl <const N:usize > Exponent<N> for u32 {
    fn to_u64_array(&self) -> Option<[u64; N]> {
        let mut _a=[0;N];
        _a[0] = *self as u64;
        Some(_a)
    }

    fn to_naf(&self) -> Vec<i8>{
        get_naf(BigUint::from_u32(*self).unwrap())
        }

    fn get_len(&self)->usize { 1 }
}

impl <const N:usize >Exponent<N> for u64 {
    fn to_u64_array(&self) -> Option<[u64; N]> {
        let mut _a=[0;N];
        _a[0] = *self as u64;
        Some(_a)
    }

    fn to_naf(&self) -> Vec<i8>{
        get_naf(BigUint::from_u64(*self).unwrap())
        }

    fn get_len(&self)->usize { 1 }
}

impl  <const N:usize >Exponent<N> for u128 {
    fn to_u64_array(&self) -> Option<[u64; N]> {
        let mut _a=[0;N];
        _a[0] = *self as u64;
        _a[1] =  (*self >> 64)  as u64;
        Some(_a)
    }

    fn to_naf(&self) -> Vec<i8>{
        get_naf(BigUint::from_u128(*self).unwrap())
        }

    fn get_len(&self)->usize { 2 }
}

impl  <const N:usize >Exponent<N> for i128 {
    fn to_u64_array(&self) -> Option<[u64; N]> {
        let mut _a=[0;N];
        let x :u128 = self.abs().try_into().unwrap();
        _a[0] = x as u64;
        _a[1] =  (x >> 64)  as u64;
        Some(_a)
    }

    fn to_naf(&self) -> Vec<i8>{
        let x :u128 = self.abs().try_into().unwrap();
        get_naf(BigUint::from_u128(x).unwrap())
        }

    fn get_len(&self)->usize { 2 }
}

