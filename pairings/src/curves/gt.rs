// Code developed by FARAOUN Kamel Mohamed.
// EEDIS-Laboratory. UDL-University. Algeria
// During May 2024.

use std::{fmt, ops::{Add, Div, Mul, Neg, Sub}};
use hmac::{Hmac, Mac};
use sha2::Sha256;

use crate::{ extensions::{ext_fields::{ExtElement, ExtField}, 
             fp12::{Fp12Element, Fp12Field}, fp24::{Fp24Element, Fp24Field}, fp48::{Fp48Element, Fp48Field}}, fields::prime_fields::FieldElement, tools::{arithmetic_interface::ArithmeticOperations, exponent::Exponent}};


#[derive(Clone, Copy)]
pub enum GTElement <const N:usize, const PARAMSIZE:usize>
            {   Fp12(Fp12Element<PARAMSIZE,N>),
                Fp24(Fp24Element<PARAMSIZE,N>),
                Fp48(Fp48Element<PARAMSIZE,N>)
            }

#[derive(Debug)]
pub enum GTField <const N:usize, const PARAMSIZE:usize>
            {   Fp12(&'static Fp12Field<PARAMSIZE,N>),
                Fp24(&'static Fp24Field<PARAMSIZE,N>),
                Fp48(&'static Fp48Field<PARAMSIZE,N>)
            }

impl <const N:usize, const PARAMSIZE:usize> ArithmeticOperations for GTElement<N,PARAMSIZE>{
    fn addto(&self, other: &Self) -> Self {
        match (self, other) {   (GTElement::Fp12(x), GTElement::Fp12(y)) => x.addto(y).into(),
                                (GTElement::Fp24(x), GTElement::Fp24(y)) => x.addto(y).into(),
                                (GTElement::Fp48(x), GTElement::Fp48(y)) => x.addto(y).into(),
                                _ => unimplemented!("Addition not implemented for different types"),
                            }
    }
    fn double(&self) -> Self {
            match self {    GTElement::Fp12(x) => x.double().into(),
                            GTElement::Fp24(x)=> x.double().into(),
                            GTElement::Fp48(x)=> x.double().into(),
                        }       
    }
    fn substract(&self, other: &Self) -> Self {
        match (self, other) {   (GTElement::Fp12(x), GTElement::Fp12(y)) => x.addto(&y.negate()).into(),
                                (GTElement::Fp24(x), GTElement::Fp24(y)) => x.addto(&y.negate()).into(),
                                (GTElement::Fp48(x), GTElement::Fp48(y)) => x.addto(&y.negate()).into(),
                                _ => unimplemented!("Substractionnot implemented for different types"),
                            }
    }
    fn multiply(&self, other: &Self) -> Self {
        match (self, other) {   (GTElement::Fp12(x), GTElement::Fp12(y)) => x.multiply(y).into(),
                                (GTElement::Fp24(x), GTElement::Fp24(y)) => x.multiply(y).into(),
                                (GTElement::Fp48(x), GTElement::Fp48(y)) => x.multiply(y).into(),
                                _ => unimplemented!("Addition not implemented for different types"),
                            }
    }
    fn sqr(&self) -> Self 
        {
            match self {    GTElement::Fp12(x) => x.sqr().into(),
                            GTElement::Fp24(x)=> x.sqr().into(),
                            GTElement::Fp48(x)=> x.sqr().into(),
                        }  
        }
    fn invert(&self) -> Self 
        {
            match self {    GTElement::Fp12(x) => x.invert().into(),
                            GTElement::Fp24(x)=> x.invert().into(),
                            GTElement::Fp48(x)=> x.invert().into(),
                    }  
        }
    fn negate(&self) -> Self 
        {
            match self {    GTElement::Fp12(x) => x.negate().into(),
                            GTElement::Fp24(x)=> x.negate().into(),
                            GTElement::Fp48(x)=> x.negate().into(),
                        }  
        }
    fn equal(&self, rhs :&Self) -> bool {
        match (self, rhs) {   (GTElement::Fp12(x), GTElement::Fp12(y)) => x.equal(y),
                                (GTElement::Fp24(x), GTElement::Fp24(y)) => x.equal(y),
                                (GTElement::Fp48(x), GTElement::Fp48(y)) => x.equal(y),
                                _ => unimplemented!("Addition not implemented for different types"),
                            }
    }
    fn is_zero(&self) -> bool {
        match self {    GTElement::Fp12(x) => x.is_zero(),
                        GTElement::Fp24(x)=> x.is_zero(),
                        GTElement::Fp48(x)=> x.is_zero(),
                    }       
    }
    fn is_one(&self) -> bool {
        match self {    GTElement::Fp12(x) => x.is_one(),
                        GTElement::Fp24(x)=> x.is_one(),
                        GTElement::Fp48(x)=> x.is_one(),
                    }       
    }

    fn to_dec_string(&self) -> String {
        match self {    GTElement::Fp12(x) => x.to_a_string(),
                        GTElement::Fp24(x)=> x.to_a_string(),
                        GTElement::Fp48(x)=> x.to_a_string(),
                    }       
    }
    
    fn to_hex_string(&self) -> String {
        match self {    GTElement::Fp12(x) => x.to_hex_string(),
                        GTElement::Fp24(x)=> x.to_hex_string(),
                        GTElement::Fp48(x)=> x.to_hex_string(),
                    }       
    }
    fn one(&self) -> Self {
        match self {    GTElement::Fp12(x) => x.one().into(),
                        GTElement::Fp24(x)=> x.one().into(),
                        GTElement::Fp48(x)=> x.one().into(),
                    }       
    }
    fn zero(&self) -> Self {
        match self {    GTElement::Fp12(x) => x.zero().into(),
                        GTElement::Fp24(x)=> x.zero().into(),
                        GTElement::Fp48(x)=> x.zero().into(),
                    }       
    }
}

impl <const N:usize, const PARAMSIZE:usize>  GTElement<N,PARAMSIZE>{
    fn mulbyu8(&self, other: u8) -> Self {
        match  self  {   GTElement::Fp12(x) => x.mulbyu8(other).into(),
                         GTElement::Fp24(x) => x.mulbyu8(other).into(),
                         GTElement::Fp48(x) => x.mulbyu8(other).into(),
                     }
        }

    pub fn encode_to_base64(&self) -> String {
            match  self  {   GTElement::Fp12(x) => x.encode_to_base64(),
                             GTElement::Fp24(x) => x.encode_to_base64(),
                             GTElement::Fp48(x) => x.encode_to_base64(),
                         }
            }
    pub fn pow(&self,e : & dyn Exponent<N>) -> Self {
                match  self  {   GTElement::Fp12(x) => {if !x.conjugate().multiply(&x).is_one() {x.pow(e).into()}
                                                                                    else {x.cyclotomic_power(e, false, &None).into()}},
                                 GTElement::Fp24(x) => {if !x.conjugate().multiply(&x).is_one() {x.pow(e).into()}
                                                                                    else {x.cyclotomic_power(e, false, &None).into()}},
                                 GTElement::Fp48(x) => {if !x.conjugate().multiply(&x).is_one() {x.pow(e).into()}
                                                                                    else {x.cyclotomic_power(e, false, &None).into()}},
                             }
                }
    pub fn sparse_multiply(&self,rhs : &[&[FieldElement<N>];3]) -> Self {
                match  self  {   GTElement::Fp12(x) => x.sparse_multiply(rhs).into(),
                                 GTElement::Fp24(x) => x.sparse_multiply(rhs).into(),
                                 GTElement::Fp48(x) => x.sparse_multiply(rhs).into(),
                                }
                } 
    pub fn one(&self) -> Self {
                match  self  {   GTElement::Fp12(x) => x.one().into(),
                                 GTElement::Fp24(x) => x.one().into(),
                                 GTElement::Fp48(x) => x.one().into(),
                             }
                }       
    pub fn content(&self) -> &[FieldElement<N>] {
                    match  self  {   GTElement::Fp12(x) => &x.content,
                                     GTElement::Fp24(x) => &x.content,
                                     GTElement::Fp48(x) => &x.content,
                                 }
                    }       
    pub fn final_exponentiation(&self) -> Self {
                        match  self  {   GTElement::Fp12(x) => x.final_exponentiation(true).into(),
                                         GTElement::Fp24(x) => x.final_exponentiation(true).into(),
                                         GTElement::Fp48(x) => x.final_exponentiation(true).into(),
                                     }
                        }    
    pub fn to_byte_array(&self) -> Vec<u8> {
                            match self {    GTElement::Fp12(x) => x.to_i2osp_bytearray(),
                                            GTElement::Fp24(x)=> x.to_i2osp_bytearray(),
                                            GTElement::Fp48(x)=> x.to_i2osp_bytearray(),
                                        }       
                        }                        
    pub fn derive_hkdf(&self,sizeinbits:usize,salt :Option<&[u8]>) -> Vec<u8>
    {
        const DSIZE :usize = 16; // length of sha256 output in bytes
        let size_in_bytes = (sizeinbits / 8) + ((sizeinbits % 8)!=0) as usize;
        if size_in_bytes < DSIZE {panic!("length of the output have to be at least equal to th hash size ...")}
        if size_in_bytes > DSIZE * 255 {panic!("length of the output cannot be longer than 255 * Hashlength ...")}
        let key = self.to_byte_array();
        let size = (sizeinbits / 8) + ((sizeinbits % 8)!=0) as usize;
        let mut _salt = Vec::<u8>::with_capacity(size);
        if !salt.is_none() {_salt.resize(salt.unwrap().len(),0);
                            _salt.extend(salt.unwrap());}
        else { _salt.resize(size, 0);}
        let mut mac = Hmac::<Sha256>::new_from_slice(&_salt).expect("HMAC can take key of any size");
        mac.update(&key);
        let extracted_result = mac.finalize().into_bytes();
        let mut okm = Vec::<u8>::new();
        let mut ti = Vec::<u8>::new();
        let n = (size_in_bytes + DSIZE - 1) / DSIZE;
        for i in 0..n {  ti.push((i+1) as u8);
                                let mut mac = Hmac::<Sha256>::new_from_slice(&ti).expect("HMAC can take key of any size");
                                mac.update(&extracted_result);
                                let tmp = mac.finalize().into_bytes();
                                okm.extend(tmp);
                                ti.extend(tmp);
                            }
        okm.truncate(size_in_bytes);
        okm
    }                

}

impl<const N: usize, const PARAMSIZE: usize> From<Fp12Element<PARAMSIZE,N>> for GTElement<N, PARAMSIZE> {
    fn from(fp12: Fp12Element<PARAMSIZE,N>) -> Self {
        GTElement::Fp12(fp12)
    }
}

impl<const N: usize, const PARAMSIZE: usize> From<Fp24Element<PARAMSIZE, N>> for GTElement<N, PARAMSIZE> {
    fn from(fp24: Fp24Element<PARAMSIZE, N>) -> Self {
        GTElement::Fp24(fp24)
    }
}

impl<const N: usize, const PARAMSIZE: usize> From<Fp48Element<PARAMSIZE, N>> for GTElement<N, PARAMSIZE> {
    fn from(fp48: Fp48Element<PARAMSIZE, N>) -> Self {
        GTElement::Fp48(fp48)
    }
}

impl<'a, const N: usize, const PARAMSIZE: usize> fmt::Display for GTElement<N,PARAMSIZE> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:}", &self.to_dec_string())
    }
    }

impl  <const PARAMSIZE:usize,const N:usize> Add for GTElement<PARAMSIZE,N> {
    type Output =  GTElement<PARAMSIZE,N>;
        fn add(self, rhs: Self) -> Self::Output {   self.addto(&rhs) }
    }

impl  <const PARAMSIZE:usize,const N:usize> Sub for GTElement<PARAMSIZE,N> {
    type Output =  GTElement<PARAMSIZE,N>;
        fn sub(self, rhs: Self) -> Self::Output {   self.substract(&rhs) }
    }

impl  <const PARAMSIZE:usize,const N:usize> Neg for GTElement<PARAMSIZE,N> {
    type Output =  GTElement<PARAMSIZE,N>;
        fn neg(self) -> Self::Output { self.negate() }
    }   

impl<const PARAMSIZE:usize,const N: usize> PartialEq for GTElement<PARAMSIZE,N> {
        fn eq(&self, other: &Self) -> bool {    self.equal(other) }
    }

impl <const PARAMSIZE:usize,const N:usize> Mul<u8> for GTElement<PARAMSIZE,N> {
        type Output = GTElement<PARAMSIZE,N>;    
        fn mul(self, rhs: u8) -> Self::Output {     self.mulbyu8(rhs)   }
    }

impl<const PARAMSIZE:usize,const N: usize> Mul<GTElement<PARAMSIZE,N>> for u8 {
    type Output = GTElement<PARAMSIZE,N>;
    fn mul(self, rhs: GTElement<PARAMSIZE,N>) -> Self::Output {  rhs.mulbyu8(self)  }
    }
    

impl<const PARAMSIZE:usize,const N:usize> Mul for GTElement<PARAMSIZE,N> {
        type Output =  GTElement<PARAMSIZE,N>;
        fn mul(self, rhs: Self) -> Self::Output {  self.multiply(&rhs)  }
    }

impl<const PARAMSIZE:usize,const N:usize> Div for GTElement<PARAMSIZE,N> {
    type Output =  GTElement<PARAMSIZE,N>;
    fn div(self, rhs: Self) -> Self::Output {   self * rhs.invert()    }                                        
    }

impl<const PARAMSIZE:usize,const N: usize> Div <GTElement<PARAMSIZE,N>> for u8 {
    type Output = GTElement<PARAMSIZE,N>;
    fn div(self, rhs: GTElement<PARAMSIZE,N>) -> Self::Output {
        match  self { 0 => rhs.zero(),
                        1 =>{ rhs.invert()}
                        _ =>{ self /rhs}
                    }        
                }
    }

impl <const N:usize, const PARAMSIZE:usize> GTField<N,PARAMSIZE> {
    pub fn one(&self) -> GTElement<N,PARAMSIZE> {
        match  self  {   GTField::Fp12(x) => x.one().into(),
                         GTField::Fp24(x) => x.one().into(),
                         GTField::Fp48(x) => x.one().into(),
                     }
        }
}


