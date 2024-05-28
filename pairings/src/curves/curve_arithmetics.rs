// Code developed by FARAOUN Kamel Mohamed.
// EEDIS-Laboratory. UDL-University. Algeria
// During May 2024.

use crate::{fields::prime_fields::FieldElement, tools::arithmetic_interface::ArithmeticOperations};
use crate::tools::recoders::*;
use num_bigint::BigUint;
use std::fmt::Display;
use std::ops::BitAnd;
use std::str::FromStr;
use num_traits::{One,ToPrimitive, Zero};


#[derive(Clone, Copy, Debug)]
pub struct EcPoint<T> {
    pub x: T,
    pub y: T,
    pub z: T,
}

impl <T> EcPoint <T>
    where 
    T :ArithmeticOperations + Clone + Copy + Display,    
{      
    pub fn negate(&self) -> EcPoint<T>
    {
        Self  { x: self.x.clone(),
                y: self.y.clone().negate(),
                z: self.z.clone() }
    }

    pub fn double_jacobian(&self)-> EcPoint<T>{
                if self.z.is_zero() { self.clone() }
                else {  // https://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html#doubling-dbl-2007-bl
                        let xx   = self.x.sqr();
                        let yy   = self.y.sqr();
                        let yy_2 = yy.sqr();
                        let zz   = self.z.sqr();
                        let s    = (self.x.addto(&yy).sqr().substract(&xx).substract(&yy_2)).double();
                        let m    = xx.double().addto(&xx);
                        let t    = m.sqr().substract(&s.double());
                        Self  { x: t.clone(),
                                y: m.multiply(&s.substract(&t)).substract(&yy_2.double().double().double()),
                                z: self.y.addto(&self.z).sqr().substract(&yy).substract(&zz) 
                              }
                    }
    }

    pub fn add_jacobian(&self, rhs :&EcPoint<T>) -> EcPoint<T>{
                if self.z.is_zero() { if rhs.z.is_zero() { self.clone() }
                                      else { rhs.clone() }  
                                    }
                else { if rhs.z.is_zero() { self.clone() }
                       else {   // https://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html#addition-add-2007-bl
                                let z1_2 = self.z.sqr();
                                let z2_2 = rhs.z.sqr();
                                let u1 = self.x.multiply(&z2_2);
                                let u2 = rhs.x.multiply(&z1_2);
                                let s1 = self.y.multiply(&rhs.z).multiply(&z2_2);
                                let s2 = rhs.y.multiply(&self.z).multiply(&z1_2);
                                let h  = u2.substract(&u1);
                                let i  = h.double().sqr();
                                let j  = h.multiply(&i);
                                let r  = s2.substract(&s1).double();
                                let v  = u1.multiply(&i);
                                let x3 = r.sqr().substract(&j).substract(&v.double());
                                Self {  x: x3.clone(), 
                                        y: r.multiply(&v.substract(&x3)).substract(&s1.multiply(&j).double()),
                                        z: (self.z.addto(&rhs.z).sqr().substract(&z1_2).substract(&z2_2)).multiply(&h)}
                            }
                    }
                }

    pub fn mixed_jacobian(&self, rhs :&EcPoint<T>) -> EcPoint<T>{
                // https://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html#addition-madd-2007-bl
                let z1_2 = self.z.sqr();
                let u2 = rhs.x.multiply(&z1_2);
                let s2 = rhs.y.multiply(&self.z).multiply(&z1_2);
                let h = u2.substract(&self.x);
                let hh = h.sqr();
                let i = hh.double().double();
                let j = h.multiply(&i);
                let r = s2.substract(&self.y).double();
                let v = self.x.multiply(&i);
                let x3 = r.sqr().substract(&j).substract(&v.double());
                Self { x : x3.clone(),
                       y : r.multiply(&v.substract(&x3)).substract(&self.y.multiply(&j).double()),
                       z : self.z.addto(&h).sqr().substract(&z1_2).substract(&hh)  
                     } 
                }
                
    pub fn multiply_with_const(&self , scalar :i128) -> EcPoint<T>{
        //  not Constant-time multiplication, used when multiplying with small constant
        //  no need for resistance to side-channel attacks !, so can do faster
        let e : u128 = scalar.abs() as u128;     
        match e { 0 => {    EcPoint{ x: self.x.one(), y: self.x.one(), z: self.x.zero() }},
                  2 => {    self.double_jacobian() },
                  _ => {    let mut result = EcPoint{ x: self.x.one(), y: self.x.one(), z: self.x.zero() };
                            let bitlen = 128 - e.leading_zeros();
                            for i in (0..bitlen).rev()
                                    {   result = result.double_jacobian();
                                        if (e >> i) & 1 == 1 {result =result.add_jacobian(&self)}                
                                    }
                            if scalar>0 {result} else {result.negate()}
                        }  
                }
        }

    pub fn multiply<const N:usize>(&self , scalar :&FieldElement<N>) -> EcPoint<T>{
                // Constant-time multiplication using w-sliding window (w=3)
                let infinit = Self {x:self.x.one(), y : self.x.one(), z: self.x.zero() };
                let ff: BigUint = BigUint::from_str("255").unwrap();
                if self.z.is_zero()  {infinit}
                else {  let s = scalar.to_big_uint();
                        if s == BigUint::zero() {infinit.clone()}
                        else {  if s == BigUint::one() {self.clone()}
                                else {  if s == (BigUint::from(2u8))  { self.double_jacobian()}
                                        else {  let mut lookup = [infinit;5];
                                                lookup[0] = self.double_jacobian();
                                                lookup[1] = self.clone();
                                                lookup[2] = self.add_jacobian(&lookup[0]);  
                                                lookup[3] = lookup[0].add_jacobian(&lookup[2]);  
                                                lookup[4] = lookup[0].add_jacobian(&lookup[3]);  
                                                let mut code = recod_one_scalar(&(&s + ((&s).bitand(BigUint::one()) + BigUint::one())));                                                
                                                let mut result = lookup[(((&code).bitand((ff).clone()).to_u8().unwrap() & (WMASK >> 1)) + 1) as usize];
                                                code = code >> (WSIZE -1);
                                                while code != BigUint::one() {  let limb = (&code).bitand((ff).clone()).to_u8().unwrap();
                                                                                let sig : i8 = (2 * (limb & 1) as i8) - 1; //2 * ((limb as i8) & 1) - 1;                                                                                
                                                                                let idx: usize = (((limb & WMASK) >> 1) + 1) as usize;                                                                                
                                                                                result = result.double_jacobian(); 
                                                                                result = result.double_jacobian(); 
                                                                                result = result.double_jacobian(); 
                                                                                if sig == 1 {result = result.add_jacobian(&lookup[idx])}
                                                                                else {result = result.add_jacobian(&lookup[idx].negate())}
                                                                                code = code >> WSIZE;
                                                                             }                
                                                result.to_affine();
                                                result = result.add_jacobian(&lookup[1 - (&s).bitand(BigUint::one()).to_usize().unwrap()].negate());                                  
                                                result 
                                                } 
                                      }    
                              }
                     }
                }

    pub fn to_affine (& mut self) {
                 if ! self.z.is_zero(){ let d = self.z.invert();
                                        let d2 = d.sqr();
                                        self.x = self.x.multiply(&d2);
                                        self.y = self.y.multiply(&d.multiply(&d2));
                                        self.z = self.x.one().clone();
                                       }
                }

    pub fn to_string(&self) -> String
                {  let mut out= String::new();
                   if self.z.is_zero() {out.push_str("Infinit");}
                   else {   out.push_str("(x : ");               
                            out.push_str(&self.x.to_dec_string());
                            out.push_str(",\ny : ");       
                            out.push_str(&self.y.to_dec_string());     
                            out.push_str(",\nz: ");  
                            out.push_str(&self.z.to_dec_string());
                            out.push_str(")");
                        }
                    String::from(&out)
                }
        
    pub fn to_hex_string(&self) -> String
                {  let mut out= String::new();
                   if self.z.is_zero() {out.push_str("Infinit");}
                   else {   out.push_str("(x: ");               
                            out.push_str(&self.x.to_hex_string());
                            out.push_str(",\ny :");       
                            out.push_str(&self.y.to_hex_string());     
                            out.push_str(",\nz :"); 
                            out.push_str(&self.z.to_hex_string());
                            out.push_str(")");
                        }
                   String::from(&out)
                }
    
    pub fn equal(&self, other :&EcPoint<T>) -> bool
                {
                    let z1_2 = self.z.sqr();
                    let z1_3 = z1_2.multiply(&self.z);
                    let z2_2 = other.z.sqr();
                    let z2_3 = z2_2.multiply(&other.z);
                    self.x.multiply(&z2_2).equal(&other.x.multiply(&z1_2)) & self.y.multiply(&z2_3).equal(&other.y.multiply(&z1_3))                   
                }
    
    pub fn is_infinit(&self) -> bool
                {
                      self.z.is_zero()   
                }
    
}