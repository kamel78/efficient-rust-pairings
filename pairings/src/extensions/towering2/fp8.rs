// Code developed by FARAOUN Kamel Mohamed.
// EEDIS-Laboratory. UDL-University. Algeria
// During May 2024.

use std::{fmt, usize};
use super::super::super::fields::prime_fields::{FieldElement, PrimeField};
use super::super::super::extensions::ext_fields::{ExtElement, ExtField,ExFieldConsts};
use super::fp2::Fp2Element;
use crate::tools::arithmetic_interface::ArithmeticOperations;
use super::fp4::Fp4Element;

#[derive(Clone, Copy,Debug)]
pub struct Fp8Element<const PARAMSIZE:usize,const N:usize>{
                                    pub content :[FieldElement<N>;8],
                                    pub constants :&'static ExFieldConsts<PARAMSIZE,N>
                                    }
#[derive(Clone,Debug)]
pub struct Fp8Field<const PARAMSIZE:usize,const N:usize> {
                                    base_field:PrimeField<N>,
                                    constants :&'static ExFieldConsts<PARAMSIZE,N>
                                    }

impl<const PARAMSIZE:usize,const N: usize> ExtField<PARAMSIZE,8,N> for Fp8Field<PARAMSIZE,N> {    
    type  ElementType   = Fp8Element<PARAMSIZE,N>;    
    type  BaseFieldType = PrimeField<N>;

    fn field_interface(&self)->PrimeField<N> 
        { self.base_field.clone() }

    fn new(base_field :&Self::BaseFieldType, consts :Option<&'static ExFieldConsts<PARAMSIZE,N>>)->  Fp8Field<PARAMSIZE,N>
        { Fp8Field {base_field : (*base_field).clone(), constants: consts.unwrap()}}   
    
    fn extconsts_interface(&self) ->Option<&'static ExFieldConsts<PARAMSIZE,N>> {
        Some(self.constants)
    }
    }

impl <const PARAMSIZE:usize,const N:usize> ExtElement<PARAMSIZE,8,N> for Fp8Element<PARAMSIZE,N>{
    fn new(content :&[FieldElement<N>; 8], consts :Option<&'static ExFieldConsts<PARAMSIZE,N>>) -> Fp8Element<PARAMSIZE,N>
        {   Fp8Element{content :content.clone(), constants :consts.unwrap()}    }

    fn content_interface(&self) -> &[super::super::super::fields::prime_fields::FieldElement<N>;8]
        {   &self.content   }
    
    fn constants_interface(&self) -> Option<&'static ExFieldConsts<PARAMSIZE,N>> {
        Some(self.constants)
    }

    fn sqr(&self) -> Self {
        let a = Fp4Element{ content : [ self.content[0],self.content[1],self.content[2], 
                                                          self.content[3]], constants:self.constants} ;
        let b = Fp4Element{ content : [ self.content[4],self.content[5],self.content[6],
                                                          self.content[7]], constants:self.constants};                                                       
        let t0= a.sqr();
        let t1= b.sqr();      
        let t2= a.addto(&b).sqr();  
        let a = t1.mulby_v().addto(&t0);
        let b = t2.substract(&t0).substract(&t1);
        Self {  content : [ a.content[0], a.content[1], a.content[2], a.content[3],
                            b.content[0], b.content[1], b.content[2], b.content[3]],
                constants : self.constants }                                             
    }

    fn multiply(&self, rhs:&Self) -> Self {       
        let a0 = Fp4Element{ content : [ self.content[0],self.content[1],self.content[2], 
                                                           self.content[3]], constants:self.constants} ;
        let b0 = Fp4Element{ content : [ self.content[4],self.content[5],self.content[6],
                                                           self.content[7]], constants:self.constants};                                                       
        let a1 = Fp4Element{ content : [ rhs.content[0],rhs.content[1],rhs.content[2], 
                                                           rhs.content[3]], constants:rhs.constants} ;
        let b1 = Fp4Element{ content : [ rhs.content[4],rhs.content[5],rhs.content[6],
                                                           rhs.content[7]], constants:rhs.constants};                                                       
        let t0 = a0.multiply(&a1);       
        let t1 = b0.multiply(&b1);     
        let t2 = a0.addto(&b0).multiply(&a1.addto(&b1));
        let a  = t1.mulby_v().addto(&t0);
        let b  = t2.substract(&t0).substract(&t1);
        Self {  content : [ a.content[0], a.content[1], a.content[2], a.content[3],
                            b.content[0], b.content[1], b.content[2], b.content[3]],
                constants : self.constants }  
    }
}

impl <const PARAMSIZE:usize,const N:usize> Fp8Element <PARAMSIZE,N>{
    
    pub fn conjugate(&self) -> Self {
        Self {content : [self.content[0] , self.content[1], self.content[2], self.content[3],
                         self.content[4].negate() , self.content[5].negate(), self.content[6].negate(), self.content[7].negate()
                        ],
              constants : self.constants}
    }

    fn frobinus_mode_2(&self) -> Self {
        Self {content : [self.content[0] , 
                         self.content[1].negate(), 
                         self.content[2].multiply(&self.constants.frobinus_consts[0]), 
                         self.content[3].negate().multiply(&self.constants.frobinus_consts[0]),
                         self.content[4].multiply(&self.constants.frobinus_consts[1]),
                         self.content[5].negate().multiply(&self.constants.frobinus_consts[1]),
                         self.content[6].multiply(&self.constants.frobinus_consts[2]),
                         self.content[7].negate().multiply(&self.constants.frobinus_consts[2])],
            constants : self.constants}
    }

    fn frobinus_mode_1(&self) -> Self {
        Self {content : [self.content[0] , 
                         self.content[1].negate(), 
                         self.content[2].multiply(&self.constants.frobinus_consts[0]), 
                         self.content[3].negate().multiply(&self.constants.frobinus_consts[0]),
                         self.content[5].negate().multiply(&self.constants.frobinus_consts[2]).multiply(&self.constants.base_qnr),
                         self.content[4].multiply(&self.constants.frobinus_consts[2]),
                         self.content[7].negate().multiply(&self.constants.frobinus_consts[4]).multiply(&self.constants.base_qnr),
                         self.content[6].multiply(&self.constants.frobinus_consts[4])],
            constants : self.constants}
    }

    pub fn frobinus(&self) -> Self {
        match self.constants.fb_id{
            1=> { self.frobinus_mode_1()}
            2=> { self.frobinus_mode_2()},
            _=>{panic!("Invalid frobinus mode for this construction ...")}
        }
    }    

    pub fn invert(&self) -> Self {
        let a= Fp4Element{   content : [self.content[0],self.content[1],self.content[2],self.content[3]], 
                                               constants:self.constants};
        let b= Fp4Element{   content : [self.content[4],self.content[5],self.content[6],self.content[7]], 
                                               constants:self.constants};
        let t0 = a.sqr();
        let t1 = b.sqr(); 
        let t = t0.substract(&t1.mulby_v()).invert();                                               
        let t0 = a.multiply(&t);
        let t1 = b.multiply(&t).negate();                                       
        Self {  content : [ t0.content[0], t0.content[1],t0.content[2], t0.content[3],
                            t1.content[0], t1.content[1],t1.content[2], t1.content[3] 
                          ],
                constants : self.constants}
    }

    pub fn is_qr(&self) -> bool{        
            let a = Fp4Element { content : [self.content[0],self.content[1],self.content[2],self.content[3]], 
                                                constants: self.constants};
            let b = Fp4Element { content : [self.content[4],self.content[5],self.content[6],self.content[7]], 
                                                constants: self.constants};
            a.sqr().substract(&b.sqr().mulby_v()).is_qr()            
            }

    pub fn sqrt(&self) -> Option<Self>{
        let outparams = self.content[0].fieldparams;
        let inv2 = FieldElement{mont_limbs: outparams.inv2, fieldparams:outparams};   
        let zero = Self { content:[FieldElement{mont_limbs:outparams.zero,fieldparams:outparams};8], 
                                               constants : self.constants};
        let a = Fp4Element { content : [self.content[0],self.content[1],self.content[2],self.content[3]], 
                                               constants: self.constants};
        let b = Fp4Element { content : [self.content[4],self.content[5],self.content[6],self.content[7]], 
                                               constants: self.constants};
        let t0= a.sqr();
        let t1= b.sqr();       
        let rootdelta = t0.substract(&t1.mulby_v()).sqrt();                                 
        if rootdelta.is_some() {
            let rootdelta = rootdelta.unwrap();
            let mut t = Fp4Element{content :[ rootdelta.content[0].addto(&self.content[0]).multiply(&inv2),
                                                                        rootdelta.content[1].addto(&self.content[1]).multiply(&inv2),
                                                                        rootdelta.content[2].addto(&self.content[2]).multiply(&inv2),
                                                                        rootdelta.content[3].addto(&self.content[3]).multiply(&inv2)],
                                                     constants : self.constants };    
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
            if r.is_none() { None } //Never occure ...
            else { if r.unwrap().equal(&Fp4Element{content :[FieldElement{mont_limbs:outparams.zero,fieldparams:outparams};4],
                                                   constants : self.constants }) 
                            {  Some(zero) }
                   else {   let r= r.unwrap();
                            let t = b.multiply(&(r.addto(&r)).invert());
                            Some(Self{content:[ r.content[0], r.content[1], r.content[2], r.content[3],
                                                t.content[0], t.content[1], t.content[2], t.content[3],
                                               ],
                                      constants : self.constants
                                    }
                             )
                        } 
                 }   
            }
            else { None }  //Never occure ...
        }
        
    pub fn mulby_w(&self) -> Self {
            Fp8Element {content : [self.content[7].multiply(&self.constants.base_qnr), self.content[6],
                                   self.content[4], self.content[5], self.content[0], self.content[1], self.content[2], self.content[3]],
                        constants : self.constants}

        }       
    
    pub fn sparse_multiply(&self, rhs:&[FieldElement<N>]) -> Self {                
                let a0 = Fp2Element{ content : [ self.content[0],self.content[1]], constants:self.constants} ;
                let a1 = Fp2Element{ content : [ self.content[2],self.content[3]], constants:self.constants} ;
                let a2 = Fp2Element{ content : [ self.content[4],self.content[5]], constants:self.constants} ;
                let a3 = Fp2Element{ content : [ self.content[6],self.content[7]], constants:self.constants} ;
                let b0 = Fp2Element{ content : [ rhs[0],rhs[1]], constants:self.constants} ;
                let b1 = Fp2Element{ content : [ rhs[2],rhs[3]], constants:self.constants} ;
                let t0 = a0.multiply(&b0);
                let t1 = a1.multiply(&b1);
                let r0 = t0.addto(&t1.mul_by_u());
                let r1 =a0.addto(&a1).multiply(&b0.addto(&b1)).substract(&t0).substract(&t1) ;                                      
                let t0 = a2.multiply(&b0);
                let t1 = a3.multiply(&b1);   
                let r2 = t0.addto(&t1.mul_by_u());            
                let r3 =a2.addto(&a3).multiply(&b0.addto(&b1)).substract(&t0).substract(&t1) ;                                                                                           
                Self {  content : [ r0.content[0], r0.content[1], r1.content[0], r1.content[1],
                                    r2.content[0], r2.content[1], r3.content[0], r3.content[1]],
                        constants : self.constants }  
            }
    }
    


impl<const PARAMSIZE:usize,const N: usize> fmt::Display for Fp8Element<PARAMSIZE,N> {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {  write!(f, "{:}", &self.to_a_string()) }
    }

impl<const PARAMSIZE:usize,const N: usize> PartialEq for Fp8Element<PARAMSIZE,N> {
        fn eq(&self, other: &Self) -> bool {    self.equal(other) }
    }
    