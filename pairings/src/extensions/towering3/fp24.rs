// Code developed by FARAOUN Kamel Mohamed.
// EEDIS-Laboratory. UDL-University. Algeria
// During May 2024.

use std::{fmt, usize};
use super::super::super::fields::prime_fields::{FieldElement, PrimeField};
use super::super::ext_fields::{ExtElement, ExtField,ExFieldConsts};
use crate::tools::arithmetic_interface::ArithmeticOperations;
use super::fp8::Fp8Element;

#[derive(Clone, Copy)]
pub struct Fp24Element<const PARAMSIZE:usize,const N:usize>{
                                    pub content :[FieldElement<N>;24],
                                    pub constants :&'static ExFieldConsts<PARAMSIZE,N>
                                    }
#[derive(Clone,Debug)]
pub struct Fp24Field<const PARAMSIZE:usize,const N:usize> {
                                    base_field:PrimeField<N>,
                                    constants :&'static ExFieldConsts<PARAMSIZE,N>
                                    }

fn get_slice_fp8<const PARAMSIZE:usize,const N:usize>(element : &Fp24Element<PARAMSIZE,N>, i : usize) -> Fp8Element<PARAMSIZE,N>
    {       
        let _t: [FieldElement<N>; 8] = element.content[i*8..(i+1)*8].iter()
        .map(|&element1| FieldElement {fieldparams: element1.fieldparams,
                                                       mont_limbs: element1.mont_limbs,
                                                      })
        .collect::<Vec<_>>()
        .try_into().unwrap(); 
        Fp8Element{content :_t, constants :element.constants}        
    }

impl<const PARAMSIZE:usize,const N: usize> ExtField<PARAMSIZE,24,N> for Fp24Field<PARAMSIZE,N> {    
    type  ElementType   = Fp24Element<PARAMSIZE,N>;    
    type  BaseFieldType = PrimeField<N>;

    fn field_interface(&self)->PrimeField<N> 
        { self.base_field.clone() }

    fn new(base_field :&Self::BaseFieldType, consts :Option<&'static ExFieldConsts<PARAMSIZE,N>>)->  Fp24Field<PARAMSIZE,N>
        { Fp24Field {base_field : (*base_field).clone(), constants: consts.unwrap()}}   
    
    fn extconsts_interface(&self) ->Option<&'static ExFieldConsts<PARAMSIZE,N>> {
        Some(self.constants)
    }
    }

impl <const PARAMSIZE:usize,const N:usize> ExtElement<PARAMSIZE,24,N> for Fp24Element<PARAMSIZE,N>{
    fn new(content :&[FieldElement<N>; 24], consts :Option<&'static ExFieldConsts<PARAMSIZE,N>>) -> Fp24Element<PARAMSIZE,N>
        {   Fp24Element{content :content.clone(), constants :consts.unwrap()}    }

    fn content_interface(&self) -> &[super::super::super::fields::prime_fields::FieldElement<N>;24]
        {   &self.content   }
    
    fn constants_interface(&self) -> Option<&'static ExFieldConsts<PARAMSIZE,N>> {
        Some(self.constants)
    }
    
    fn sqr(&self) -> Self {
        let a = get_slice_fp8(&self, 0);
        let b = get_slice_fp8(&self, 1);
        let c = get_slice_fp8(&self, 2);
        let t0 = a.sqr();
        let t1 = b.sqr();
        let t2 = c.sqr();
        let h0 = b.addto(&c).sqr().substract(&t1).substract(&t2).mulby_w().negate().addto(&t0);
        let h1 = t2.mulby_w().negate().addto(&a.addto(&b).sqr().substract(&t0).substract(&t1));
        let h2 = a.addto(&c).sqr().substract(&t0).substract(&t2).addto(&t1); 
        let mut result =[FieldElement{  mont_limbs:self.content[0].fieldparams.zero,
                                                               fieldparams:self.content[0].fieldparams};24];
        result[..8].copy_from_slice(&h0.content);
        result[8..16].copy_from_slice(&h1.content);
        result[16..].copy_from_slice(&h2.content);
            Self {  content : result, constants :self.constants}
        }
        
    fn multiply(&self, rhs:&Self) -> Self {          
        let a0 = get_slice_fp8(&self, 0);
        let a1 = get_slice_fp8(&rhs, 0);
        let b0 = get_slice_fp8(&self, 1);
        let b1 = get_slice_fp8(&rhs, 1);
        let c0 = get_slice_fp8(&self, 2);
        let c1 = get_slice_fp8(&rhs, 2);
        let t0 = a0.multiply(&a1);
        let t1 = b0.multiply(&b1);
        let t2 = c0.multiply(&c1);
        let h0 = b0.addto(&c0).multiply(&b1.addto(&c1)).substract(&t1).substract(&t2).mulby_w().negate().addto(&t0);        
        let h1 = t2.mulby_w().negate().addto(&a0.addto(&b0).multiply(&a1.addto(&b1)).substract(&t0).substract(&t1));
        let h2 = a0.addto(&c0).multiply(&a1.addto(&c1)).substract(&t0).substract(&t2).addto(&t1);                
        let mut result =[FieldElement{  mont_limbs:self.content[0].fieldparams.zero,
                                                               fieldparams:self.content[0].fieldparams};24];
        result[..8].copy_from_slice(&h0.content);
        result[8..16].copy_from_slice(&h1.content);
        result[16..].copy_from_slice(&h2.content);
        Self {  content : result, constants :self.constants}                                                                                                                       
    }
}

impl <const PARAMSIZE:usize,const N:usize> Fp24Element <PARAMSIZE,N>{  
    pub fn invert(&self) -> Self {
        let a = get_slice_fp8(&self, 0);
        let b = get_slice_fp8(&self, 1);
        let c = get_slice_fp8(&self, 2);
        let t0 = a.sqr().substract(&b.multiply(&c).mulby_w().negate());
        let t1 = c.sqr().mulby_w().negate().substract(&a.multiply(&b));
        let t2 = b.sqr().substract(&a.multiply(&c));                             
        let t3 = c.multiply(&t1).addto(&b.multiply(&t2)).mulby_w().negate().addto(&a.multiply(&t0)).invert();                
        let mut result =[FieldElement{  mont_limbs:self.content[0].fieldparams.zero,
                                                               fieldparams:self.content[0].fieldparams};24];
        result[..8].copy_from_slice(&t0.multiply(&t3).content);
        result[8..16].copy_from_slice(&t1.multiply(&t3).content);
        result[16..].copy_from_slice(&t2.multiply(&t3).content);
        Self {  content : result, constants :self.constants} 
    }

    pub fn mulby_z(&self) -> Self{
        let mut result =self.content.clone();
        result[0] = self.content[23].substract(&self.content[22]);
        result[1] = self.content[22].addto(&self.content[23]).negate();
        result[2] = self.content[20];
        result[3] = self.content[21];       
        result[4]=self.content[16].negate();
        result[5]=self.content[17].negate();
        result[6]=self.content[18].negate();
        result[7]=self.content[19].negate();
        result[8..24].copy_from_slice(&self.content[0..16]);
        Fp24Element {content : result,
                     constants : self.constants}
    }
    pub fn sparse_multiply_for48(&self, rhs:&[&[FieldElement<N>];3], mode :u8) -> Self {         
        match  mode { 3|0 |1 => { //  Sparse multiplication of an Fp24 element with a sparse element on two consecutive Fp8 only
                             //  Used during pairings on BLS48 (with tow possible configurations of positions M/D-Types)
                            let a0 = get_slice_fp8(&self, 0);        
                            let b0 = get_slice_fp8(&self, 1);
                            let c0 = get_slice_fp8(&self, 2);
                            let mut a1 = Fp8Element{ content: rhs[if mode!=3{2} else {0}].try_into().unwrap(), constants: self.constants };
                            let mut b1 = Fp8Element{ content: rhs[if mode!=3{1} else {1}].try_into().unwrap(), constants: self.constants };
                            if mode ==0 {a1 =a1.mulby_u();
                                         b1 = b1.mulby_u();   
                                        }
                            if mode ==3 {b1 = b1.negate();     }                                        
                            let t0 = a0.multiply(&a1);
                            let t1 = c0.multiply(&b1);
                            let t2 = b0.multiply(&b1);
                            let t3 = a0.addto(&b0).multiply(&a1.addto(&b1));
                            let res1 = t1.mulby_w().negate().addto(&t0);
                            let res2 = t3.substract(&t0).substract(&t2);
                            let res3 = t2.addto(&c0.multiply(&a1)); 
                            let mut result =[FieldElement{  mont_limbs:self.content[0].fieldparams.zero,
                                                                                   fieldparams:self.content[0].fieldparams};24];
                            result[..8].copy_from_slice(&res1.content);
                            result[8..16].copy_from_slice(&res2.content);
                            result[16..].copy_from_slice(&res3.content);
                            Self {  content : result, constants :self.constants}  
                            }      
                      4|2 => { //  Sparse multiplication of an Fp24 element with a sparse element on one Fp8 only
                             //  Used during pairings on BLS48, rhs on 4 Fp elements (with tow possible configurations of positions M/D-Types)
                            let a0 = get_slice_fp8(&self, 0);        
                            let b0 = get_slice_fp8(&self, 1);
                            let c0 = get_slice_fp8(&self, 2);
                            let y =Fp8Element {content : rhs[if mode==2 {0} else {2}].try_into().unwrap(),constants :self.constants};
                            let res1 = a0.multiply(&y);
                            let res2 = b0.multiply(&y);
                            let res3 = c0.multiply(&y);
                            let mut result =[FieldElement{  mont_limbs:self.content[0].fieldparams.zero,
                                                                                   fieldparams:self.content[0].fieldparams};24];
                            result[..8].copy_from_slice(&res1.content);
                            result[8..16].copy_from_slice(&res2.content);
                            result[16..].copy_from_slice(&res3.content);
                            Self {  content : result, constants :self.constants}  
                            }
                      _ =>{ panic!("Invalid sparse multiplication mode for FP24 .")}
                    }                                                                                                              
    }
               
}




impl<const PARAMSIZE:usize,const N: usize> fmt::Display for Fp24Element<PARAMSIZE,N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {  write!(f, "{:}", &self.to_a_string()) }
}

impl<const PARAMSIZE:usize,const N: usize> PartialEq for Fp24Element<PARAMSIZE,N> {
    fn eq(&self, other: &Self) -> bool {    self.equal(other) }
}
