// Code developed by FARAOUN Kamel Mohamed.
// EEDIS-Laboratory. UDL-University. Algeria
// During May 2024.

use std::{str::FromStr, usize};
use std::fmt;
use std::ops::{Add,Sub,Mul,Div,Neg};
use base64::engine::general_purpose;
use base64::Engine;
use num_bigint::{BigInt, BigUint, ToBigInt, ToBigUint};
use num_traits::{ Euclid, Num, One, ToPrimitive, Zero};
use crate::tools::arithmetic_interface::ArithmeticOperations;
use crate::tools::hashs::*; 
use super::super::tools::exponent::Exponent;
use super::arithmetic;

#[derive(Debug, PartialEq, Eq)]
pub enum Endianness {
    Little,
    Big,
}

#[derive(Clone,Debug)]
pub struct FieldParams<'a, const N:usize> {
pub numlimbs : usize,
pub modulo : [u64;N],
pub modulo_as_strhex:&'a str,
pub one:[u64;N],
pub qprime:u128,
pub rsquare:[u64;N],
pub sqrtid:u8,  
pub modplus1div4:[u64;N],
pub zero:[u64;N],
pub inv2:[u64;N],
pub optimize_sqr:bool,
pub optimize_mul:bool,
pub num_of_bits:usize,
pub sig_theshold : [u64;N],
}

#[derive(Clone, Copy,Debug)]
pub struct FieldElement <const N:usize> {
pub fieldparams:&'static FieldParams<'static,N>,
pub mont_limbs:[u64;N]

}

#[derive(Clone,Debug)]
pub struct PrimeField<const N:usize>{
pub parametres :&'static FieldParams<'static , N>,
pub modulo_as_bigint: BigUint
}


impl <'a, const N:usize> PrimeField<N> {             
    pub fn new(input_params :&'static FieldParams<N>)-> PrimeField<N>
        {   let bigint_modulo : BigUint =  BigUint::from_str_radix(&input_params.modulo_as_strhex[2..], 16).unwrap(); 
            PrimeField {parametres : input_params, modulo_as_bigint: bigint_modulo}
        }   

    pub fn random_element(&self) -> FieldElement<N> 
        {   use rand::Rng;
            let mut rng = rand::thread_rng();
            let random_bytes: Vec<u8> = (0..N*8).map(|_| rng.gen()).collect();
            let randomelement  = BigUint::from_bytes_be(&random_bytes); 
            Self::from_bigint(&self,&((randomelement % self.modulo_as_bigint.clone()).to_bigint()).unwrap())
        }

    pub fn from_bigint(&self, input : &BigInt)-> FieldElement<N> 
            // Load the Fp element from a BigInt signed value, and convert to montgomery form     
        {   let mut limbs     : [u64;N] = [0;N];              
            let reduced_big   : BigUint = input.rem_euclid(&self.modulo_as_bigint.to_bigint().unwrap()).to_biguint().unwrap();
            let mask           = (BigUint::one() << 64) - BigUint::one();
            for i in 0..N { limbs[i] = (((&reduced_big & (&mask << (i * 64))) >> (i * 64)) as BigUint).to_u64().unwrap()}
            FieldElement{   fieldparams : &self.parametres,  
                            mont_limbs  : super::arithmetic::mul(&limbs, &self.parametres.rsquare, &self.parametres)
                            }
        }
    
    pub fn from_u64(&self, input : u64)-> FieldElement<N> 
        // Load the Fp element from a u8 signed value, and convert to montgomery form     
    {   let mut limbs     : [u64;N] = [0;N];              
        limbs[0] = input;
        FieldElement{   fieldparams : &self.parametres,  
                        mont_limbs  : super::arithmetic::mul(&limbs, &self.parametres.rsquare, &self.parametres)
                        }
    }

    pub fn from_str(&self,input : &str)-> FieldElement<N> 
            // Load the Fp element from a string, and convert to montgomery form
        {   let mut strin = input;
            let mut negative = false; 
            if &strin[0..1]=="-" { negative =true;
                                   strin = &strin[1..];
                                 }
            if negative {Self::from_bigint(&self,&BigInt::from_str(strin).unwrap()).negate()}
            else {Self::from_bigint(&self,&BigInt::from_str(strin).unwrap())}        
        }

    pub fn from_hex_str(&self, input : &str) -> FieldElement<N> 
            // Load the Fp element from a hex-string, and convert to montgomery form     
        {   let mut strin = input;
            let mut negative = false; 
            if &strin[0..1]=="-" { negative =true;
                                   strin = &strin[1..];
                                 }
            if &strin[0..2]=="0x" { strin = &strin[2..];}
            if negative { Self::from_bigint(&self,&BigInt::from_str_radix(&strin,16).unwrap()).negate()}
            else {Self::from_bigint(&self,&BigInt::from_str_radix(&strin,16).unwrap())}
        }

    pub fn from_byte_array(&self, source :&[u8], repre:Endianness) -> FieldElement<N>{
            let mut _source = source.to_vec();
            if repre == Endianness::Big {(_source).reverse();}
            let numbits = self.parametres.num_of_bits;                
            let sizeinbytes = (numbits >> 3) + if (numbits % 8) ==0 {0} else {1}; 
            if sizeinbytes != _source.len() {panic!("Size of input does not correspond to the field's extension ...");}                   
            Self::from_bigint(&self,&os2ip(&_source).to_bigint().unwrap())
        }                   
    
    pub fn from_base64(&self, source :&str) -> FieldElement<N>
        {
            let decoded_bytes = match general_purpose::STANDARD.decode(source) {
                Ok(bytes) => bytes,
                Err(_) => {
                    panic!("Failed to decode base64 string");
                }
            };
            self.from_byte_array(&decoded_bytes, Endianness::Little)
        }
    
    pub fn zero(&self) ->FieldElement<N>  {
            FieldElement{   fieldparams: &self.parametres, 
                            mont_limbs: [0u64; N],
                        }
                    }
    pub fn one(&self) ->FieldElement<N>  {
                        FieldElement{   fieldparams: &self.parametres, 
                                        mont_limbs: self.parametres.one,
                                    }
                                }
    pub fn hash_to_field(&self,id : &str, security_level:usize,count :usize) -> Vec<FieldElement<N>>
        {   
            hash_string_to_field(id, self, count, security_level, 1)
        }
    } 


impl  <'a, const N:usize> FieldElement<N>  {
    pub fn pow(&self, e :& dyn Exponent<N>, useladder:bool) -> FieldElement<N>
    {   
        let array = e.to_u64_array().unwrap();
        FieldElement{   mont_limbs : arithmetic::pow(&self.mont_limbs, &array, useladder, self.fieldparams),
                        fieldparams : &self.fieldparams} 
    }
    pub fn sqrt(&self) -> Option<FieldElement<N>> 
    {   let result =super::arithmetic::sqrt(&self.mont_limbs, &self.fieldparams);
        if result.is_none() { None }
        else {  Some(FieldElement{   fieldparams : &self.fieldparams, 
                                     mont_limbs  : result.unwrap()}) 
            }
    }
    
    pub fn to_big_uint(&self) -> BigUint
    {   let mut as_big : BigUint = BigUint::zero();
        let mut one :[u64;N] = [0;N];
        one[0] = 1;
        let from_mont  : [u64;N] = super::arithmetic::mul(&self.mont_limbs,&one,&self.fieldparams);
        for i in 0..N {as_big += from_mont[i].to_biguint().unwrap() << (i * 64)};
        as_big
    } 

    pub fn sign(&self) -> i8
    {   let mut one :[u64;N] = [0;N];
        one[0] = 1;
        let val = arithmetic::mul(&self.mont_limbs, &one, self.fieldparams);        
        let mut i =N-1;
        while (val[i] != self.fieldparams.sig_theshold[i]) & (i>0) { i=i-1}
        let sig = if val[i] > self.fieldparams.sig_theshold[i] {-1} else {1};        
        sig
    }

    pub fn to_i2osp_bytearray(&self) -> Vec<u8>
    {  let mut out= Vec::<u8>::new();
       let numbits = self.fieldparams.num_of_bits;                
       let sizeinbytes = (numbits >> 3) + if (numbits % 8) ==0 {0} else {1};     
       out.extend(i2osp_pf(&self, sizeinbytes));      
       out 
    }

    pub fn to_base64(&self) -> String
    {
        general_purpose::STANDARD.encode(self.to_i2osp_bytearray())
    }
}

impl <const N:usize> ArithmeticOperations for FieldElement<N> {    
    fn addto(&self, other: &Self) -> Self {
        Self {  fieldparams:self.fieldparams, 
            mont_limbs: super::arithmetic::add(&self.mont_limbs, &other.mont_limbs, &self.fieldparams)}
    }

    fn double(&self) -> Self {
        self.addto(&self)
    }

    fn substract(&self, other: &Self) -> Self {
        Self {  fieldparams:self.fieldparams, 
                mont_limbs: super::arithmetic::sub(&self.mont_limbs, &other.mont_limbs, &self.fieldparams)}
    }

    fn multiply(&self, other: &Self) -> Self {
        Self {  fieldparams:self.fieldparams, 
            mont_limbs: super::arithmetic::mul(&self.mont_limbs, &other.mont_limbs, &self.fieldparams)}
    }

    fn equal(&self, rhs :&Self) -> bool {
        super::arithmetic::equal(&self.mont_limbs, &rhs.mont_limbs)
    }

    fn to_dec_string(&self) -> String {
        self.to_big_uint().to_string()
    }
    
    fn to_hex_string(&self) -> String {
        String::from("0x") + &Self::to_big_uint(self).to_str_radix(16)
    }    
    
    fn sqr(&self) -> Self {
        FieldElement {  fieldparams : &self.fieldparams,  
            mont_limbs : super::arithmetic::sqr(&self.mont_limbs, &self.fieldparams)} 
    }

    fn invert(&self) -> Self {
        FieldElement{   fieldparams : &self.fieldparams,  
            mont_limbs : super::arithmetic::invert(&self.mont_limbs,&self.fieldparams)}        
    }

    fn negate(&self) -> Self {
        Self {  fieldparams:self.fieldparams, 
            mont_limbs: super::arithmetic::neg(&self.mont_limbs, &self.fieldparams )}
    }
    fn one(&self) -> Self {
        FieldElement{   fieldparams: &self.fieldparams, 
                        mont_limbs: self.fieldparams.one.clone(),
                    }
    }   
    fn zero(&self) -> Self {
        FieldElement{   fieldparams: &self.fieldparams, 
                        mont_limbs: self.fieldparams.zero.clone(),
                    }
    }   
    fn is_zero(&self) -> bool {
        let zero = [0u64;N];
        arithmetic::equal(&self.mont_limbs,&zero)
    }    
    fn is_one(&self) -> bool {
        arithmetic::equal(&self.mont_limbs,&self.fieldparams.one)
    }
    
}

impl  <'a, const N:usize> Add for FieldElement<N> {
    type Output =  FieldElement<N>;
    fn add(self, rhs: Self) -> Self::Output {
        self.addto(&rhs)
    }
    }
    
    impl <'a, const N:usize> Add<u64> for FieldElement<N> {
        type Output = FieldElement<N>;
        fn add(self, rhs: u64) -> Self::Output {
            let mut _rhs =[0;N];
            _rhs[0]=rhs;
            self.addto(& Self::Output{fieldparams : self.fieldparams,  
                               mont_limbs  : super::arithmetic::mul(&_rhs, &self.fieldparams.rsquare, &self.fieldparams)})
        }
    }
    
    impl  <'a, const N:usize> Sub for FieldElement<N> {
    type Output =  FieldElement<N>;
    fn sub(self, rhs: Self) -> Self::Output {
        self.substract(&rhs)
    }
    }
    
    impl  <'a, const N:usize> Mul for FieldElement<N> {
    type Output =  FieldElement<N>;
    fn mul(self, rhs: Self) -> Self::Output {
        self.multiply(&rhs)
    }
    }
    
    
    impl <'a, const N:usize> Mul<u8> for FieldElement<N> {
    type Output = FieldElement<N>;
    fn mul(self, rhs: u8) -> Self::Output {
        match  rhs { 0 => Self::Output {fieldparams:self.fieldparams, mont_limbs:self.fieldparams.zero},
                        1 => self,
                        2 => self + self,
                        3 => self + self + self,
                        4 =>{let double =self+self;
                            double + double},      
                        5 =>{let double =self+self;
                            let fourth =double+double;
                            fourth + self },
                        _ =>{let mut limbs : [u64;N] = [0;N];              
                            limbs[0] = rhs as u64;
                            limbs = super::arithmetic::mul(&limbs, &self.fieldparams.rsquare, &self.fieldparams);
                            Self::Output {fieldparams : self.fieldparams,  
                                          mont_limbs  : super::arithmetic::mul(&limbs, &self.mont_limbs, &self.fieldparams)
                                         }
                            }            
                    }        
                }
    }
    
    impl<'a, const N: usize> Mul<FieldElement<N>> for u8 {
    type Output = FieldElement<N>;
    fn mul(self, rhs: FieldElement<N>) -> Self::Output {
        match  self { 0 => Self::Output {fieldparams:rhs.fieldparams, mont_limbs:rhs.fieldparams.zero},
                        1 => rhs,
                        2 => rhs + rhs,
                        3 => rhs + rhs + rhs,
                        4 =>{let double =rhs + rhs;
                            double + double},      
                        5 =>{let double =rhs + rhs;
                            let fourth =double+double;
                            fourth + rhs },
                        _ =>{let mut limbs : [u64;N] = [0;N];              
                             limbs[0] = self as u64;
                             Self::Output { fieldparams : rhs.fieldparams,  
                                            mont_limbs  : super::arithmetic::mul(&limbs, &rhs.fieldparams.rsquare, &rhs.fieldparams)
                                         }
                            }
                    }                   
                }
    }
    
    impl  <'a, const N:usize> Div for FieldElement<N> {
    type Output =  FieldElement<N>;
    fn div(self, rhs: Self) -> Self::Output {
        self.multiply(&rhs.invert())
    }
    }
    
    impl<'a, const N: usize> PartialEq for FieldElement<N> {
    fn eq(&self, rhs: &Self) -> bool {
        self.equal(&rhs)
    }
    }
    
    impl<'a, const N: usize> Eq for FieldElement<N> {}
    
    impl  <'a, const N:usize> Neg for FieldElement<N> {
    type Output =  FieldElement<N>;
    fn neg(self) -> Self::Output {
        self.negate()
    }
    }   
    
    impl<'a, const N: usize> fmt::Display for FieldElement<N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:}", &self.to_dec_string())
    }
    }
   
