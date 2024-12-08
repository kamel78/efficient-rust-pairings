// Code developed by FARAOUN Kamel Mohamed.
// EEDIS-Laboratory. UDL-University. Algeria
// During May 2024.

use std::{fmt, usize};
use super::super::super::fields::prime_fields::{FieldElement, PrimeField};
use crate::tools::arithmetic_interface::ArithmeticOperations;
use super::super::super::extensions::ext_fields::{ExtElement, ExtField,ExFieldConsts};
use super::fp2::Fp2Element;

#[derive(Clone, Copy,Debug)]
pub struct Fp4Element<const PARAMSIZE:usize,const N:usize>{
                                    pub content :[FieldElement<N>;4],
                                    pub constants :&'static ExFieldConsts<PARAMSIZE,N>
                                    }
#[derive(Clone,Debug)]
pub struct Fp4Field<const PARAMSIZE:usize,const N:usize> {
                                    base_field:PrimeField<N>,
                                    constants :&'static ExFieldConsts<PARAMSIZE,N>
                                    }

impl<const PARAMSIZE:usize,const N: usize> ExtField<PARAMSIZE,4,N> for Fp4Field<PARAMSIZE,N> {    
    type  ElementType   = Fp4Element<PARAMSIZE,N>;    
    type  BaseFieldType = PrimeField<N>;

    fn field_interface(&self)->PrimeField<N> 
        { self.base_field.clone() }

    fn new(base_field :&Self::BaseFieldType, consts :Option<&'static ExFieldConsts<PARAMSIZE,N>>) ->  Fp4Field<PARAMSIZE,N>
        { Fp4Field {base_field : (*base_field).clone(), constants: consts.unwrap()}}   
    
    fn extconsts_interface(&self) -> Option<&'static ExFieldConsts<PARAMSIZE,N>> {
        Some(self.constants)
    }
    }

impl <const PARAMSIZE:usize,const N:usize> ExtElement<PARAMSIZE,4,N> for Fp4Element<PARAMSIZE,N>{
    fn new(content :&[FieldElement<N>; 4], consts :Option<&'static ExFieldConsts<PARAMSIZE,N>>) -> Fp4Element<PARAMSIZE,N>
        {   Fp4Element{content :content.clone(), constants :consts.unwrap()}    }

    fn content_interface(&self) -> &[super::super::super::fields::prime_fields::FieldElement<N>;4]
        {   &self.content   }
    
    fn constants_interface(&self) -> Option<&'static ExFieldConsts<PARAMSIZE,N>> {
        Some(self.constants)
    }

    fn sqr(&self) -> Self {
        let a = Fp2Element{content :[self.content[0],self.content[1]], constants:self.constants };
        let b = Fp2Element{content :[self.content[2],self.content[3]], constants:self.constants };
        let v0= a.sqr();
        let v1= b.sqr();
        let b = a.addto(&b).sqr().substract(&v0).substract(&v1);
        let a = v0.addto(&v1.mul_by_u());
        Self { content : [a.content[0],a.content[1],b.content[0],b.content[1]],
               constants :self.constants  }
    }

    fn multiply(&self, rhs:&Self) -> Self {        
        let a0 = Fp2Element{content :[self.content[0],self.content[1]], constants:self.constants };
        let b0 = Fp2Element{content :[self.content[2],self.content[3]], constants:self.constants };
        let a1 = Fp2Element{content :[rhs.content[0],rhs.content[1]], constants:self.constants };
        let b1 = Fp2Element{content :[rhs.content[2],rhs.content[3]], constants:self.constants };
        let t0 = a0.multiply(&a1);
        let t1 = b0.multiply(&b1);
        let a  = t0.addto(&t1.mul_by_u());
        let b  = a0.addto(&b0).multiply(&a1.addto(&b1)).substract(&t0).substract(&t1);
        Self { content : [a.content[0],a.content[1],b.content[0],b.content[1]],
               constants :self.constants  }
    }

}

impl <const PARAMSIZE:usize,const N:usize> Fp4Element <PARAMSIZE,N>{

    pub fn conjugate(&self) -> Self {
        Self {content : [self.content[0] , self.content[1], self.content[2].negate(), self.content[3].negate()],
              constants : self.constants}
    }

    pub fn frobinus(&self) -> Self {
        Self {  content : [self.content[0] , 
                           self.content[1].negate(), 
                           self.content[2].multiply(&self.constants.frobinus_consts[0]), 
                           self.content[3].negate().multiply(&self.constants.frobinus_consts[0])],
                constants : self.constants}
    }   
 
    pub fn invert(&self) -> Self{        
        let a = Fp2Element{content :[self.content[0],self.content[1]], constants:self.constants };
        let b = Fp2Element{content :[self.content[2],self.content[3]], constants:self.constants };        
        let t = a.sqr().substract(&b.sqr().mul_by_u()).invert();               
        let a = a.multiply(&t);
        let b = b.multiply(&t).negate();
        Self { content : [a.content[0],a.content[1],b.content[0],b.content[1]],
            constants :self.constants  }      
    }

    pub fn is_qr(&self) -> bool{        
        let a = Fp2Element{content :[self.content[0],self.content[1]], constants:self.constants };
        let b = Fp2Element{content :[self.content[2],self.content[3]], constants:self.constants };
        a.sqr().substract(&b.sqr().mul_by_u()).is_qr()
        }

    pub fn sqrt(&self) -> Option<Self>{
        let outparams = self.content[0].fieldparams;
        let inv2 = FieldElement{mont_limbs: outparams.inv2, fieldparams:outparams};   
        let zero =Self{  content:[FieldElement{mont_limbs:outparams.zero,fieldparams:outparams};4], 
                                              constants : self.constants};
        let a = Fp2Element{content :[self.content[0],self.content[1]], constants:self.constants};
        let b = Fp2Element{content :[self.content[2],self.content[3]], constants:self.constants };
        let rootdelta: Option<Fp2Element<PARAMSIZE,N>> = a.sqr().substract(&b.sqr().mul_by_u()).sqrt();
        if rootdelta.is_some() {
            let rootdelta = rootdelta.unwrap();
            let mut t = Fp2Element{content :[ rootdelta.content[0].addto(&self.content[0]).multiply(&inv2),
                                                             rootdelta.content[1].addto(&self.content[1]).multiply(&inv2)],
                                                             constants:self.constants};    
            let mut r = t.sqrt();
            if r.is_none() { t = t.negate();                                     
                             r = t.sqrt();   
                             if r.is_none(){ t = t.negate();
                                             r = t.sqrt();   
                                             if r.is_none(){ t = t.substract(&rootdelta);
                                                             r = t.sqrt(); 
                                                            }
                                            }
                           }
            if r.is_none() { None }
            else { if r.unwrap().equal(&Fp2Element{ content :[FieldElement{mont_limbs:outparams.zero,fieldparams:outparams};2], 
                                                    constants:self.constants}) 
                            {  Some(zero) }
                   else {   let r= r.unwrap();
                            let t =b.multiply(&r.double().invert());
                            Some(Self{content:[ r.content[0],r.content[1],
                                                t.content[0],t.content[1]
                                               ],
                                      constants : self.constants
                                    }
                             )
                        } 
                 }   
            }
            else { None }  
        }
    
    pub fn mulby_v(&self) -> Self {
            Fp4Element {content : [self.content[3].multiply(&self.constants.base_qnr), self.content[2],self.content[0], self.content[1]],                                                                                                                                                             
                        constants : self.constants}

        }
    }

impl<const PARAMSIZE:usize,const N: usize> fmt::Display for Fp4Element<PARAMSIZE,N> {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {  write!(f, "{:}", &self.to_a_string()) }
    }

impl<const PARAMSIZE:usize,const N: usize> PartialEq for Fp4Element<PARAMSIZE,N> {
        fn eq(&self, other: &Self) -> bool {    self.equal(other) }
    }
    
