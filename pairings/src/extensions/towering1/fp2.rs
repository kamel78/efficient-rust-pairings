// Code developed by FARAOUN Kamel Mohamed.
// EEDIS-Laboratory. UDL-University. Algeria
// During May 2024.

use std::{fmt, usize};
use super::super::super::fields::prime_fields::{FieldElement, PrimeField};
use super::super::ext_fields::{ExtElement, ExtField,ExFieldConsts};
use crate::tools::arithmetic_interface::ArithmeticOperations;

#[derive(Clone, Copy,Debug)]
pub struct Fp2Element<const N:usize>{
                                    pub content :[FieldElement<N>;2]
                                    }
#[derive(Clone,Debug)]
pub struct Fp2Field<const N:usize> {
                                    base_field:PrimeField<N>,
                                    }

impl<const N: usize> ExtField<0,2,N> for Fp2Field<N> {     
    type  ElementType   = Fp2Element<N>;    
    type  BaseFieldType = PrimeField<N>;

    fn field_interface(&self)->PrimeField<N> 
        { self.base_field.clone() }
    
    fn extconsts_interface(&self) ->Option<&'static ExFieldConsts<0,N>> {
        None
    }

    fn new(base_field :&Self::BaseFieldType, _consts :Option<&ExFieldConsts<0,N>>)->  Fp2Field<N>
        { Fp2Field {base_field : (*base_field).clone()}}   
    }


impl <const N:usize> ExtElement<0,2,N> for Fp2Element<N>{
    fn new(content :&[FieldElement<N>; 2],_consts :Option<&ExFieldConsts<0,N>>) -> Fp2Element<N>
        {   Fp2Element{content :content.clone()}    }

    fn content_interface(&self) -> &[super::super::super::fields::prime_fields::FieldElement<N>;2]
        {   &self.content   }
    
    fn constants_interface(&self) -> Option<&'static ExFieldConsts<0,N>>{
        None
    }

    fn sqr(&self) -> Self {
        let v0 =self.content[0].sqr();
        let v1 =self.content[1].sqr();  
        let v2 =self.content[0].multiply(&self.content[1]);
        Self {content : [v0.substract(&v1), v2.addto(&v2)]}
    }

    fn multiply(&self, rhs:&Self) -> Self {        
        let v0 =self.content[0].multiply(&rhs.content[0]);
        let v1 =self.content[1].multiply(&rhs.content[1]);            
        Self {content : [v0.substract(&v1),
                            (self.content[0].addto(&self.content[1]).multiply(
                                            &rhs.content[0].addto(&rhs.content[1])).substract(&v0).substract(&v1))
                                    ]}
    }

}

impl <const N:usize> Fp2Element<N>{
    
    pub fn mul_by_u_p_1(&self) -> Self {
        Self {content : [self.content[0].substract(&self.content[1]),self.content[0].addto(&self.content[1])]}

    } 
    pub fn conjugate(&self) -> Self {
        Self {content : [self.content[0] , self.content[1].negate()]}
    }

    pub fn frobinus(&self) -> Self {
        // In Fp2, the frobinus is equal to the conjugate
        self.conjugate()
    }   
    pub fn invert(&self) -> Self{
        let  t = (self.content[0].sqr()).addto(&self.content[1].sqr()).invert();
        Self {content :[self.content[0].multiply(&t),(self.content[1].multiply(&t)).negate()]}
    }
    
    pub fn is_qr(&self) -> bool{
        let rootdelta = (self.content[0].sqr()).addto(&self.content[1].sqr()).sqrt(); 
        rootdelta.is_some() 
    }

    pub fn sqrt(&self) -> Option<Self>{
        let outparams = self.content[0].fieldparams;
        let inv2 = FieldElement{mont_limbs: outparams.inv2, fieldparams:outparams };   
        let zero =Self {content:[ FieldElement{mont_limbs:outparams.zero,fieldparams:outparams},
                                                    FieldElement{mont_limbs:outparams.zero,fieldparams:outparams}]};
        let rootdelta = (self.content[0].sqr()).addto(&self.content[1].sqr()).sqrt(); 
        if rootdelta.is_some() {
            let mut t1 = (self.content[0].addto(&rootdelta.unwrap())).multiply(&inv2);            
            let mut a = t1.sqrt();
            if a.is_none() { t1 = t1.substract(&rootdelta.unwrap());                             
                             a  = t1.sqrt();   
                             if a.is_none(){ t1 = t1.negate();
                                             a  = t1.sqrt();   
                                             if a.is_none(){ t1 = t1.substract(&rootdelta.unwrap());
                                                             a = t1.sqrt(); 
                                                            }
                                            }
                           }
            if a.is_none() { None }
            else { if a.unwrap().equal(&FieldElement{mont_limbs:outparams.zero,fieldparams:outparams}) 
                            {  Some(zero) }
                   else {Some(Self{content:[a.unwrap(),
                                            self.content[1].multiply(&(a.unwrap().addto(&a.unwrap())).invert())]
                                  }
                             )
                        } 
                 }   
            }
            else { None }  
        }
    }
    
impl<const N: usize> fmt::Display for Fp2Element<N> {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result { write!(f, "{:}", &self.to_a_string()) }
    }

impl<const N: usize> PartialEq for Fp2Element<N> {
        fn eq(&self, other: &Self) -> bool {    self.equal(other) }
    }



                        
        
