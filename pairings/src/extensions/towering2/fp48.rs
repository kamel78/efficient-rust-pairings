// Code developed by FARAOUN Kamel Mohamed.
// EEDIS-Laboratory. UDL-University. Algeria
// During May 2024.

use std::{fmt, usize};
use crate::tools::exponent::Exponent;

use super::super::super::fields::prime_fields::{FieldElement, PrimeField};
use crate::tools::arithmetic_interface::ArithmeticOperations;
use super::super::super::extensions::ext_fields::{ExtElement, ExtField,ExFieldConsts};
use super::fp24::Fp24Element;
use super::fp8::Fp8Element;


#[derive(Clone, Copy)]
pub struct Fp48Element<const PARAMSIZE:usize,const N:usize>{
                                    pub content :[FieldElement<N>;48],
                                    pub constants :&'static ExFieldConsts<PARAMSIZE,N>
                                    }
#[derive(Clone,Debug)]
pub struct Fp48Field<const PARAMSIZE:usize,const N:usize> {
                                    base_field:PrimeField<N>,
                                    constants :&'static ExFieldConsts<PARAMSIZE,N>
                                    }

fn get_slice_fp24<const PARAMSIZE:usize,const N:usize>(element : &Fp48Element<PARAMSIZE,N>, i : usize) -> Fp24Element<PARAMSIZE,N>
    {       
        let _t: [FieldElement<N>; 24] = element.content[i*24..(i+1)*24].iter()
        .map(|&element1| FieldElement {fieldparams: element1.fieldparams,
                                                       mont_limbs: element1.mont_limbs,
                                                      })
        .collect::<Vec<_>>()
        .try_into().unwrap(); 
             Fp24Element{content :_t, constants :element.constants}        
    }

fn get_slice_fp8<const PARAMSIZE:usize,const N:usize>(element : &Fp48Element<PARAMSIZE,N>, i : usize) -> Fp8Element<PARAMSIZE,N>
    {       
        let _t: [FieldElement<N>; 8] = element.content[i*8..(i+1)*8].iter()
        .map(|&element1| FieldElement {fieldparams: element1.fieldparams,
                                                       mont_limbs: element1.mont_limbs,
                                                      })
        .collect::<Vec<_>>()
        .try_into().unwrap(); 
             Fp8Element{content :_t, constants :element.constants}        
    }

impl<const PARAMSIZE:usize,const N: usize> ExtField<PARAMSIZE,48,N> for Fp48Field<PARAMSIZE,N> {    
    type  ElementType   = Fp48Element<PARAMSIZE,N>;    
    type  BaseFieldType = PrimeField<N>;

    fn field_interface(&self)->PrimeField<N> 
        { self.base_field.clone() }

    fn new(base_field :&Self::BaseFieldType, consts :Option<&'static ExFieldConsts<PARAMSIZE,N>>)->  Fp48Field<PARAMSIZE,N>
        { Fp48Field {base_field : (*base_field).clone(), constants: consts.unwrap()}}   
    
    fn extconsts_interface(&self) ->Option<&'static ExFieldConsts<PARAMSIZE,N>> {
        Some(self.constants)
    }
    }

impl <const PARAMSIZE:usize,const N:usize> ExtElement<PARAMSIZE,48,N> for Fp48Element<PARAMSIZE,N>{
    fn new(content :&[FieldElement<N>; 48], consts :Option<&'static ExFieldConsts<PARAMSIZE,N>>) -> Fp48Element<PARAMSIZE,N>
        {   Fp48Element{content :content.clone(), constants :consts.unwrap()}    }

    fn content_interface(&self) -> &[super::super::super::fields::prime_fields::FieldElement<N>;48]
        {   &self.content   }
    
    fn constants_interface(&self) -> Option<&'static ExFieldConsts<PARAMSIZE,N>> {
        Some(self.constants)
    }
    
    fn sqr(&self) -> Self {
        let a = get_slice_fp24(&self, 0);
        let b = get_slice_fp24(&self, 1);
        let t0 = a.sqr();            
        let t1 = b.sqr();
        let t2 = a.addto(&b).sqr();              
        let x0 = t1.mulby_z().addto(&t0);
        let x1 = t2.substract(&t0).substract(&t1);
        let mut result:[FieldElement<N>;48] = [self.content[0];48];
        result[..24].copy_from_slice(&x0.content);
        result[24..].copy_from_slice(&x1.content);
        Self {  content : result,
                constants : self.constants}                                                                                                                                                    
    }

    fn multiply(&self, rhs:&Self) -> Self {          
        let a0 = get_slice_fp24(&self, 0);
        let b0 = get_slice_fp24(&self, 1);
        let a1 = get_slice_fp24(&rhs, 0);
        let b1 = get_slice_fp24(&rhs, 1);
        let t0 = a0.multiply(&a1);
        let t1 = b0.multiply(&b1);
        let t3 = a0.addto(&b0).multiply(&a1.addto(&b1));
        let x0 = t1.mulby_z().addto(&t0);
        let x1 = t3.substract(&t0).substract(&t1);
        let mut result:[FieldElement<N>;48] = [self.content[0];48] ;
        result[..24].copy_from_slice(&x0.content);
        result[24..].copy_from_slice(&x1.content);
        Self {  content : result,
                constants : self.constants}                                                                                                                                                    
    }
}

impl <const PARAMSIZE:usize,const N:usize> Fp48Element <PARAMSIZE,N>{

    pub fn conjugate(&self) -> Self {
        let mut result:[FieldElement<N>;48] = [self.content[0];48];
        result[..24].copy_from_slice(&self.content[0..24]);
        let b = get_slice_fp24(&self, 1).negate();
        result[24..].copy_from_slice(&b.content);
        Self {  content : result,
                constants : self.constants}                                                                                                                                                            
    }

    fn frobinus_mode_1(&self, order :u8) -> Self {      
        match  order { 1 =>{ Self { content : [ self.content[0], 
                                                self.content[1].negate(),
                                                self.content[2].multiply(&self.constants.frobinus_consts[0]),
                                                self.content[3].negate().multiply(&self.constants.frobinus_consts[0]),
                                                self.content[5].negate().multiply(&self.constants.frobinus_consts[1]),
                                                self.content[4].multiply(&self.constants.frobinus_consts[2]),
                                                self.content[7].negate().multiply(&self.constants.frobinus_consts[3]),
                                                self.content[6].multiply(&self.constants.frobinus_consts[4]),

                                                self.content[9].negate().multiply(&self.constants.frobinus_consts[6]),
                                                self.content[8].multiply(&self.constants.frobinus_consts[5]),
                                                self.content[11].negate().multiply(&self.constants.frobinus_consts[8]),
                                                self.content[10].multiply(&self.constants.frobinus_consts[7]),
                                                self.content[12].multiply(&self.constants.frobinus_consts[10]),
                                                self.content[13].negate().multiply(&self.constants.frobinus_consts[10]),
                                                self.content[14].multiply(&self.constants.frobinus_consts[11]),
                                                self.content[15].negate().multiply(&self.constants.frobinus_consts[11]),

                                                self.content[16].multiply(&self.constants.frobinus_consts[12]),
                                                self.content[17].negate().multiply(&self.constants.frobinus_consts[12]),
                                                self.content[18].multiply(&self.constants.frobinus_consts[13]),
                                                self.content[19].negate().multiply(&self.constants.frobinus_consts[13]),
                                                self.content[21].negate().multiply(&self.constants.frobinus_consts[15]),
                                                self.content[20].multiply(&self.constants.frobinus_consts[14]),
                                                self.content[23].negate().multiply(&self.constants.frobinus_consts[17]),
                                                self.content[22].multiply(&self.constants.frobinus_consts[16]),

                                                self.content[27].multiply(&self.constants.frobinus_consts[18]),
                                                self.content[26].multiply(&self.constants.frobinus_consts[19]),
                                                self.content[24].multiply(&self.constants.frobinus_consts[20]),
                                                self.content[25].negate().multiply(&self.constants.frobinus_consts[20]),
                                                self.content[30].multiply(&self.constants.frobinus_consts[21]),
                                                self.content[31].negate().multiply(&self.constants.frobinus_consts[21]),
                                                self.content[29].multiply(&self.constants.frobinus_consts[22]),
                                                self.content[28].multiply(&self.constants.frobinus_consts[23]),

                                                self.content[34].multiply(&self.constants.frobinus_consts[24]),
                                                self.content[35].negate().multiply(&self.constants.frobinus_consts[24]),
                                                self.content[33].multiply(&self.constants.frobinus_consts[25]),
                                                self.content[32].multiply(&self.constants.frobinus_consts[26]),
                                                self.content[39].multiply(&self.constants.frobinus_consts[27]),
                                                self.content[38].multiply(&self.constants.frobinus_consts[28]),                                                
                                                self.content[36].multiply(&self.constants.frobinus_consts[29]),
                                                self.content[37].negate().multiply(&self.constants.frobinus_consts[29]),

                                                self.content[43].multiply(&self.constants.frobinus_consts[30]),
                                                self.content[42].multiply(&self.constants.frobinus_consts[31]),                                                                                                
                                                self.content[40].multiply(&self.constants.frobinus_consts[32]),
                                                self.content[41].negate().multiply(&self.constants.frobinus_consts[32]),                                                                                                
                                                self.content[46].multiply(&self.constants.frobinus_consts[33]),
                                                self.content[47].negate().multiply(&self.constants.frobinus_consts[33]),
                                                self.content[45].multiply(&self.constants.frobinus_consts[34]),
                                                self.content[44].multiply(&self.constants.frobinus_consts[35])
                                                ],
                                    constants :self.constants }
                            },                    
                       2 =>{ Self { content : [ self.content[0], 
                                                self.content[1],
                                                self.content[2].negate(),
                                                self.content[3].negate(),
                                                self.content[4].negate().multiply(&self.constants.frobinus_consts[0]),
                                                self.content[5].negate().multiply(&self.constants.frobinus_consts[0]),
                                                self.content[6].multiply(&self.constants.frobinus_consts[0]),
                                                self.content[7].multiply(&self.constants.frobinus_consts[0]),

                                                self.content[8].negate().multiply(&self.constants.frobinus_consts[12]),
                                                self.content[9].negate().multiply(&self.constants.frobinus_consts[12]),
                                                self.content[10].multiply(&self.constants.frobinus_consts[12]),
                                                self.content[11].multiply(&self.constants.frobinus_consts[12]),
                                                self.content[12].multiply(&self.constants.frobinus_consts[13]),
                                                self.content[13].multiply(&self.constants.frobinus_consts[13]),
                                                self.content[14].negate().multiply(&self.constants.frobinus_consts[13]),
                                                self.content[15].negate().multiply(&self.constants.frobinus_consts[13]),

                                                self.content[16].multiply(&self.constants.frobinus_consts[10]),
                                                self.content[17].multiply(&self.constants.frobinus_consts[10]),
                                                self.content[18].negate().multiply(&self.constants.frobinus_consts[10]),
                                                self.content[19].negate().multiply(&self.constants.frobinus_consts[10]),
                                                self.content[20].negate().multiply(&self.constants.frobinus_consts[11]),
                                                self.content[21].negate().multiply(&self.constants.frobinus_consts[11]),
                                                self.content[22].multiply(&self.constants.frobinus_consts[11]),
                                                self.content[23].multiply(&self.constants.frobinus_consts[11]),

                                                self.content[25].multiply(&self.constants.frobinus_consts[8]),
                                                self.content[24].multiply(&self.constants.frobinus_consts[7]),
                                                self.content[27].negate().multiply(&self.constants.frobinus_consts[8]),
                                                self.content[26].negate().multiply(&self.constants.frobinus_consts[7]),
                                                self.content[29].multiply(&self.constants.frobinus_consts[6]),
                                                self.content[28].multiply(&self.constants.frobinus_consts[5]),
                                                self.content[31].negate().multiply(&self.constants.frobinus_consts[6]),
                                                self.content[30].negate().multiply(&self.constants.frobinus_consts[5]),

                                                self.content[33].negate().multiply(&self.constants.frobinus_consts[3]),
                                                self.content[32].negate().multiply(&self.constants.frobinus_consts[4]),
                                                self.content[35].multiply(&self.constants.frobinus_consts[3]),
                                                self.content[34].multiply(&self.constants.frobinus_consts[4]),                                                
                                                self.content[37].negate().multiply(&self.constants.frobinus_consts[1]),                                                
                                                self.content[36].negate().multiply(&self.constants.frobinus_consts[2]),                                                
                                                self.content[39].multiply(&self.constants.frobinus_consts[1]),                                                
                                                self.content[38].multiply(&self.constants.frobinus_consts[2]),                                                

                                                self.content[41].multiply(&self.constants.frobinus_consts[17]),                                                
                                                self.content[40].multiply(&self.constants.frobinus_consts[16]),                                                
                                                self.content[43].negate().multiply(&self.constants.frobinus_consts[17]),                                                
                                                self.content[42].negate().multiply(&self.constants.frobinus_consts[16]),                                                
                                                self.content[45].multiply(&self.constants.frobinus_consts[15]),                                                
                                                self.content[44].multiply(&self.constants.frobinus_consts[14]),                                                
                                                self.content[47].negate().multiply(&self.constants.frobinus_consts[15]),                                                
                                                self.content[46].negate().multiply(&self.constants.frobinus_consts[14]),                                                
                                                ],
                    constants :self.constants }
                            },                 
                       4 =>{Self { content : [self.content[0], 
                                              self.content[1],
                                              self.content[2],
                                              self.content[3],
                                              self.content[4].negate(),
                                              self.content[5].negate(),
                                              self.content[6].negate(),
                                              self.content[7].negate(),

                                              self.content[8].multiply(&self.constants.frobinus_consts[10]),
                                              self.content[9].multiply(&self.constants.frobinus_consts[10]),
                                              self.content[10].multiply(&self.constants.frobinus_consts[10]),
                                              self.content[11].multiply(&self.constants.frobinus_consts[10]),
                                              self.content[12].negate().multiply(&self.constants.frobinus_consts[10]),
                                              self.content[13].negate().multiply(&self.constants.frobinus_consts[10]),
                                              self.content[14].negate().multiply(&self.constants.frobinus_consts[10]),
                                              self.content[15].negate().multiply(&self.constants.frobinus_consts[10]),

                                              self.content[16].multiply(&self.constants.frobinus_consts[13]),
                                              self.content[17].multiply(&self.constants.frobinus_consts[13]),
                                              self.content[18].multiply(&self.constants.frobinus_consts[13]),
                                              self.content[19].multiply(&self.constants.frobinus_consts[13]),
                                              self.content[20].negate().multiply(&self.constants.frobinus_consts[13]),
                                              self.content[21].negate().multiply(&self.constants.frobinus_consts[13]),
                                              self.content[22].negate().multiply(&self.constants.frobinus_consts[13]),
                                              self.content[23].negate().multiply(&self.constants.frobinus_consts[13]),

                                              self.content[24].negate().multiply(&self.constants.frobinus_consts[12]),
                                              self.content[25].negate().multiply(&self.constants.frobinus_consts[12]),
                                              self.content[26].negate().multiply(&self.constants.frobinus_consts[12]),
                                              self.content[27].negate().multiply(&self.constants.frobinus_consts[12]),
                                              self.content[28].multiply(&self.constants.frobinus_consts[12]),
                                              self.content[29].multiply(&self.constants.frobinus_consts[12]),
                                              self.content[30].multiply(&self.constants.frobinus_consts[12]),
                                              self.content[31].multiply(&self.constants.frobinus_consts[12]),

                                              self.content[32].negate().multiply(&self.constants.frobinus_consts[0]),
                                              self.content[33].negate().multiply(&self.constants.frobinus_consts[0]),
                                              self.content[34].negate().multiply(&self.constants.frobinus_consts[0]),
                                              self.content[35].negate().multiply(&self.constants.frobinus_consts[0]),
                                              self.content[36].multiply(&self.constants.frobinus_consts[0]),
                                              self.content[37].multiply(&self.constants.frobinus_consts[0]),
                                              self.content[38].multiply(&self.constants.frobinus_consts[0]),
                                              self.content[39].multiply(&self.constants.frobinus_consts[0]),

                                              self.content[40].negate().multiply(&self.constants.frobinus_consts[11]),
                                              self.content[41].negate().multiply(&self.constants.frobinus_consts[11]),
                                              self.content[42].negate().multiply(&self.constants.frobinus_consts[11]),
                                              self.content[43].negate().multiply(&self.constants.frobinus_consts[11]),
                                              self.content[44].multiply(&self.constants.frobinus_consts[11]),
                                              self.content[45].multiply(&self.constants.frobinus_consts[11]),
                                              self.content[46].multiply(&self.constants.frobinus_consts[11]),
                                              self.content[47].multiply(&self.constants.frobinus_consts[11]),
                                              ],
                                   constants :self.constants }}                     
                       8 =>{Self { content : [self.content[0], 
                                              self.content[1],
                                              self.content[2],
                                              self.content[3],
                                              self.content[4],
                                              self.content[5],
                                              self.content[6],
                                              self.content[7],

                                              self.content[8].multiply(&self.constants.frobinus_consts[13]),
                                              self.content[9].multiply(&self.constants.frobinus_consts[13]),
                                              self.content[10].multiply(&self.constants.frobinus_consts[13]),
                                              self.content[11].multiply(&self.constants.frobinus_consts[13]),
                                              self.content[12].multiply(&self.constants.frobinus_consts[13]),
                                              self.content[13].multiply(&self.constants.frobinus_consts[13]),
                                              self.content[14].multiply(&self.constants.frobinus_consts[13]),
                                              self.content[15].multiply(&self.constants.frobinus_consts[13]),

                                              self.content[16].negate().multiply(&self.constants.frobinus_consts[10]),
                                              self.content[17].negate().multiply(&self.constants.frobinus_consts[10]),
                                              self.content[18].negate().multiply(&self.constants.frobinus_consts[10]),
                                              self.content[19].negate().multiply(&self.constants.frobinus_consts[10]),
                                              self.content[20].negate().multiply(&self.constants.frobinus_consts[10]),
                                              self.content[21].negate().multiply(&self.constants.frobinus_consts[10]),
                                              self.content[22].negate().multiply(&self.constants.frobinus_consts[10]),
                                              self.content[23].negate().multiply(&self.constants.frobinus_consts[10]),

                                              self.content[24].multiply(&self.constants.frobinus_consts[10]),
                                              self.content[25].multiply(&self.constants.frobinus_consts[10]),
                                              self.content[26].multiply(&self.constants.frobinus_consts[10]),
                                              self.content[27].multiply(&self.constants.frobinus_consts[10]),
                                              self.content[28].multiply(&self.constants.frobinus_consts[10]),
                                              self.content[29].multiply(&self.constants.frobinus_consts[10]),
                                              self.content[30].multiply(&self.constants.frobinus_consts[10]),
                                              self.content[31].multiply(&self.constants.frobinus_consts[10]),

                                              self.content[32].negate(),
                                              self.content[33].negate(),
                                              self.content[34].negate(),
                                              self.content[35].negate(),
                                              self.content[36].negate(),
                                              self.content[37].negate(),
                                              self.content[38].negate(),
                                              self.content[39].negate(),
                                              
                                              self.content[40].negate().multiply(&self.constants.frobinus_consts[13]),
                                              self.content[41].negate().multiply(&self.constants.frobinus_consts[13]),
                                              self.content[42].negate().multiply(&self.constants.frobinus_consts[13]),
                                              self.content[43].negate().multiply(&self.constants.frobinus_consts[13]),
                                              self.content[44].negate().multiply(&self.constants.frobinus_consts[13]),
                                              self.content[45].negate().multiply(&self.constants.frobinus_consts[13]),
                                              self.content[46].negate().multiply(&self.constants.frobinus_consts[13]),
                                              self.content[47].negate().multiply(&self.constants.frobinus_consts[13]),
                                            ],
                                  constants :self.constants }}
                       _=>{ panic!("Frobinus not implemented for order {}.",order  )}            
        }       
    }

    fn frobinus_mode_2(&self, order :u8) -> Self {      
        match  order { 1 =>{ Self { content : [ self.content[0], 
                                                self.content[1].negate(),
                                                self.content[2].multiply(&self.constants.frobinus_consts[0]),
                                                self.content[3].negate().multiply(&self.constants.frobinus_consts[0]),
                                                self.content[4].multiply(&self.constants.frobinus_consts[1]),
                                                self.content[5].negate().multiply(&self.constants.frobinus_consts[1]),
                                                self.content[6].multiply(&self.constants.frobinus_consts[2]),
                                                self.content[7].negate().multiply(&self.constants.frobinus_consts[2]),

                                                self.content[8].multiply(&self.constants.frobinus_consts[3]),
                                                self.content[9].negate().multiply(&self.constants.frobinus_consts[3]),
                                                self.content[10].multiply(&self.constants.frobinus_consts[4]),
                                                self.content[11].negate().multiply(&self.constants.frobinus_consts[4]),
                                                self.content[12].multiply(&self.constants.frobinus_consts[5]),
                                                self.content[13].negate().multiply(&self.constants.frobinus_consts[5]),
                                                self.content[14].multiply(&self.constants.frobinus_consts[6]),
                                                self.content[15].negate().multiply(&self.constants.frobinus_consts[6]),

                                                self.content[16].multiply(&self.constants.frobinus_consts[7]),
                                                self.content[17].negate().multiply(&self.constants.frobinus_consts[7]),
                                                self.content[18].multiply(&self.constants.frobinus_consts[8]),
                                                self.content[19].negate().multiply(&self.constants.frobinus_consts[8]),
                                                self.content[20].multiply(&self.constants.frobinus_consts[9]),
                                                self.content[21].negate().multiply(&self.constants.frobinus_consts[9]),
                                                self.content[22].multiply(&self.constants.frobinus_consts[10]),
                                                self.content[23].negate().multiply(&self.constants.frobinus_consts[10]),

                                                self.content[25].multiply(&self.constants.frobinus_consts[12]),
                                                self.content[24].multiply(&self.constants.frobinus_consts[11]),
                                                self.content[27].multiply(&self.constants.frobinus_consts[14]),
                                                self.content[26].multiply(&self.constants.frobinus_consts[13]),
                                                self.content[29].multiply(&self.constants.frobinus_consts[16]),
                                                self.content[28].multiply(&self.constants.frobinus_consts[15]),
                                                self.content[31].multiply(&self.constants.frobinus_consts[18]),
                                                self.content[30].multiply(&self.constants.frobinus_consts[17]),

                                                self.content[33].multiply(&self.constants.frobinus_consts[20]),
                                                self.content[32].multiply(&self.constants.frobinus_consts[19]),
                                                self.content[35].multiply(&self.constants.frobinus_consts[22]),
                                                self.content[34].multiply(&self.constants.frobinus_consts[21]),
                                                self.content[37].multiply(&self.constants.frobinus_consts[24]),
                                                self.content[36].multiply(&self.constants.frobinus_consts[23]),                                                
                                                self.content[39].multiply(&self.constants.frobinus_consts[26]),
                                                self.content[38].multiply(&self.constants.frobinus_consts[25]),

                                                self.content[41].multiply(&self.constants.frobinus_consts[28]),
                                                self.content[40].multiply(&self.constants.frobinus_consts[27]),                                                                                                
                                                self.content[43].multiply(&self.constants.frobinus_consts[30]),
                                                self.content[42].multiply(&self.constants.frobinus_consts[29]),                                                                                                
                                                self.content[45].multiply(&self.constants.frobinus_consts[32]),
                                                self.content[44].multiply(&self.constants.frobinus_consts[31]),
                                                self.content[47].multiply(&self.constants.frobinus_consts[34]),
                                                self.content[46].multiply(&self.constants.frobinus_consts[33])
                                                ],
                                    constants :self.constants }
                            },                      
                       2 =>{ Self { content : [ self.content[0], 
                                                self.content[1],
                                                self.content[2].negate(),
                                                self.content[3].negate(),
                                                self.content[4].multiply(&self.constants.frobinus_consts[0]),
                                                self.content[5].multiply(&self.constants.frobinus_consts[0]),
                                                self.content[6].negate().multiply(&self.constants.frobinus_consts[0]),
                                                self.content[7].negate().multiply(&self.constants.frobinus_consts[0]),

                                                self.content[8].multiply(&self.constants.frobinus_consts[7]),                                                
                                                self.content[9].multiply(&self.constants.frobinus_consts[7]),
                                                self.content[10].negate().multiply(&self.constants.frobinus_consts[7]),
                                                self.content[11].negate().multiply(&self.constants.frobinus_consts[7]),
                                                self.content[12].multiply(&self.constants.frobinus_consts[8]),
                                                self.content[13].multiply(&self.constants.frobinus_consts[8]),
                                                self.content[14].negate().multiply(&self.constants.frobinus_consts[8]),
                                                self.content[15].negate().multiply(&self.constants.frobinus_consts[8]),

                                                self.content[16].multiply(&self.constants.frobinus_consts[5]),
                                                self.content[17].multiply(&self.constants.frobinus_consts[5]),
                                                self.content[18].negate().multiply(&self.constants.frobinus_consts[5]),
                                                self.content[19].negate().multiply(&self.constants.frobinus_consts[5]),
                                                self.content[20].multiply(&self.constants.frobinus_consts[6]),
                                                self.content[21].multiply(&self.constants.frobinus_consts[6]),
                                                self.content[22].negate().multiply(&self.constants.frobinus_consts[6]),
                                                self.content[23].negate().multiply(&self.constants.frobinus_consts[6]),

                                                self.content[24].negate().multiply(&self.constants.frobinus_consts[3]),
                                                self.content[25].negate().multiply(&self.constants.frobinus_consts[3]),
                                                self.content[26].multiply(&self.constants.frobinus_consts[3]),
                                                self.content[27].multiply(&self.constants.frobinus_consts[3]),
                                                self.content[28].negate().multiply(&self.constants.frobinus_consts[4]),
                                                self.content[29].negate().multiply(&self.constants.frobinus_consts[4]),
                                                self.content[30].multiply(&self.constants.frobinus_consts[4]),
                                                self.content[31].multiply(&self.constants.frobinus_consts[4]),

                                                self.content[32].negate().multiply(&self.constants.frobinus_consts[1]),
                                                self.content[33].negate().multiply(&self.constants.frobinus_consts[1]),
                                                self.content[34].multiply(&self.constants.frobinus_consts[1]),
                                                self.content[35].multiply(&self.constants.frobinus_consts[1]),                                                
                                                self.content[36].negate().multiply(&self.constants.frobinus_consts[2]),                                                
                                                self.content[37].negate().multiply(&self.constants.frobinus_consts[2]),                                                
                                                self.content[38].multiply(&self.constants.frobinus_consts[2]),                                                
                                                self.content[39].multiply(&self.constants.frobinus_consts[2]),                                                

                                                self.content[40].negate().multiply(&self.constants.frobinus_consts[9]),                                                
                                                self.content[41].negate().multiply(&self.constants.frobinus_consts[9]),                                                
                                                self.content[42].multiply(&self.constants.frobinus_consts[9]),                                                
                                                self.content[43].multiply(&self.constants.frobinus_consts[9]),                                                
                                                self.content[44].negate().multiply(&self.constants.frobinus_consts[10]),                                                
                                                self.content[45].negate().multiply(&self.constants.frobinus_consts[10]),                                                
                                                self.content[46].multiply(&self.constants.frobinus_consts[10]),                                                
                                                self.content[47].multiply(&self.constants.frobinus_consts[10]),                                                
                                                ],
                    constants :self.constants }
                            },                                      
                       4 =>{Self { content : [self.content[0], 
                                              self.content[1],
                                              self.content[2],
                                              self.content[3],
                                              self.content[4].negate(),
                                              self.content[5].negate(),
                                              self.content[6].negate(),
                                              self.content[7].negate(),

                                              self.content[8].multiply(&self.constants.frobinus_consts[5]),
                                              self.content[9].multiply(&self.constants.frobinus_consts[5]),
                                              self.content[10].multiply(&self.constants.frobinus_consts[5]),
                                              self.content[11].multiply(&self.constants.frobinus_consts[5]),
                                              self.content[12].negate().multiply(&self.constants.frobinus_consts[5]),
                                              self.content[13].negate().multiply(&self.constants.frobinus_consts[5]),
                                              self.content[14].negate().multiply(&self.constants.frobinus_consts[5]),
                                              self.content[15].negate().multiply(&self.constants.frobinus_consts[5]),

                                              self.content[16].multiply(&self.constants.frobinus_consts[8]),
                                              self.content[17].multiply(&self.constants.frobinus_consts[8]),
                                              self.content[18].multiply(&self.constants.frobinus_consts[8]),
                                              self.content[19].multiply(&self.constants.frobinus_consts[8]),
                                              self.content[20].negate().multiply(&self.constants.frobinus_consts[8]),
                                              self.content[21].negate().multiply(&self.constants.frobinus_consts[8]),
                                              self.content[22].negate().multiply(&self.constants.frobinus_consts[8]),
                                              self.content[23].negate().multiply(&self.constants.frobinus_consts[8]),

                                              self.content[24].multiply(&self.constants.frobinus_consts[7]),
                                              self.content[25].multiply(&self.constants.frobinus_consts[7]),
                                              self.content[26].multiply(&self.constants.frobinus_consts[7]),
                                              self.content[27].multiply(&self.constants.frobinus_consts[7]),
                                              self.content[28].negate().multiply(&self.constants.frobinus_consts[7]),
                                              self.content[29].negate().multiply(&self.constants.frobinus_consts[7]),
                                              self.content[30].negate().multiply(&self.constants.frobinus_consts[7]),
                                              self.content[31].negate().multiply(&self.constants.frobinus_consts[7]),

                                              self.content[32].multiply(&self.constants.frobinus_consts[0]),
                                              self.content[33].multiply(&self.constants.frobinus_consts[0]),
                                              self.content[34].multiply(&self.constants.frobinus_consts[0]),
                                              self.content[35].multiply(&self.constants.frobinus_consts[0]),
                                              self.content[36].negate().multiply(&self.constants.frobinus_consts[0]),
                                              self.content[37].negate().multiply(&self.constants.frobinus_consts[0]),
                                              self.content[38].negate().multiply(&self.constants.frobinus_consts[0]),
                                              self.content[39].negate().multiply(&self.constants.frobinus_consts[0]),

                                              self.content[40].multiply(&self.constants.frobinus_consts[6]),
                                              self.content[41].multiply(&self.constants.frobinus_consts[6]),
                                              self.content[42].multiply(&self.constants.frobinus_consts[6]),
                                              self.content[43].multiply(&self.constants.frobinus_consts[6]),
                                              self.content[44].negate().multiply(&self.constants.frobinus_consts[6]),
                                              self.content[45].negate().multiply(&self.constants.frobinus_consts[6]),
                                              self.content[46].negate().multiply(&self.constants.frobinus_consts[6]),
                                              self.content[47].negate().multiply(&self.constants.frobinus_consts[6]),
                                              ],
                                   constants :self.constants }}                                                    
                       8 =>{Self { content : [self.content[0], 
                                              self.content[1],
                                              self.content[2],
                                              self.content[3],
                                              self.content[4],
                                              self.content[5],
                                              self.content[6],
                                              self.content[7],

                                              self.content[8].multiply(&self.constants.frobinus_consts[8]),
                                              self.content[9].multiply(&self.constants.frobinus_consts[8]),
                                              self.content[10].multiply(&self.constants.frobinus_consts[8]),
                                              self.content[11].multiply(&self.constants.frobinus_consts[8]),
                                              self.content[12].multiply(&self.constants.frobinus_consts[8]),
                                              self.content[13].multiply(&self.constants.frobinus_consts[8]),
                                              self.content[14].multiply(&self.constants.frobinus_consts[8]),
                                              self.content[15].multiply(&self.constants.frobinus_consts[8]),

                                              self.content[16].negate().multiply(&self.constants.frobinus_consts[5]),
                                              self.content[17].negate().multiply(&self.constants.frobinus_consts[5]),
                                              self.content[18].negate().multiply(&self.constants.frobinus_consts[5]),
                                              self.content[19].negate().multiply(&self.constants.frobinus_consts[5]),
                                              self.content[20].negate().multiply(&self.constants.frobinus_consts[5]),
                                              self.content[21].negate().multiply(&self.constants.frobinus_consts[5]),
                                              self.content[22].negate().multiply(&self.constants.frobinus_consts[5]),
                                              self.content[23].negate().multiply(&self.constants.frobinus_consts[5]),

                                              self.content[24].multiply(&self.constants.frobinus_consts[5]),
                                              self.content[25].multiply(&self.constants.frobinus_consts[5]),
                                              self.content[26].multiply(&self.constants.frobinus_consts[5]),
                                              self.content[27].multiply(&self.constants.frobinus_consts[5]),
                                              self.content[28].multiply(&self.constants.frobinus_consts[5]),
                                              self.content[29].multiply(&self.constants.frobinus_consts[5]),
                                              self.content[30].multiply(&self.constants.frobinus_consts[5]),
                                              self.content[31].multiply(&self.constants.frobinus_consts[5]),

                                              self.content[32].negate(),
                                              self.content[33].negate(),
                                              self.content[34].negate(),
                                              self.content[35].negate(),
                                              self.content[36].negate(),
                                              self.content[37].negate(),
                                              self.content[38].negate(),
                                              self.content[39].negate(),                         

                                              self.content[40].negate().multiply(&self.constants.frobinus_consts[8]),
                                              self.content[41].negate().multiply(&self.constants.frobinus_consts[8]),
                                              self.content[42].negate().multiply(&self.constants.frobinus_consts[8]),
                                              self.content[43].negate().multiply(&self.constants.frobinus_consts[8]),
                                              self.content[44].negate().multiply(&self.constants.frobinus_consts[8]),
                                              self.content[45].negate().multiply(&self.constants.frobinus_consts[8]),
                                              self.content[46].negate().multiply(&self.constants.frobinus_consts[8]),
                                              self.content[47].negate().multiply(&self.constants.frobinus_consts[8]),
                                            ],
                                  constants :self.constants }}
                       _=>{ panic!("Frobinus not implemented for order {}.",order  )}            
        }       
    }
    pub fn frobinus(&self, order :u8) -> Self {      
        match self.constants.fb_id {
            1 => {self.frobinus_mode_1(order)}
            2 => {self.frobinus_mode_2(order)}
            _ => {panic!("Invalid frobinus mode for this construction ....")}
        }
    }

    pub fn invert(&self) -> Self {        
        let a = get_slice_fp24(&self, 0); 
        let b = get_slice_fp24(&self, 1); 
        let t = a.sqr().substract(&b.sqr().mulby_z()).invert();
        let x0 = a.multiply(&t);
        let x1 = b.negate().multiply(&t);
        let mut result:[FieldElement<N>;48]=[self.content[0];48];
        result[..24].copy_from_slice(&x0.content);
        result[24..].copy_from_slice(&x1.content);
        Self {  content : result,
                constants : self.constants}                                                                                                                                                    
    }

    pub fn unisqr(&self) -> Self { 
        // Adapted Scott's square for 2 - 3 - 2 construction from "Choosing and generating parameters for low level pairing implementation on BN curves" section 7.2.2
        // https://eprint.iacr.org/2015/1212 
        let x0 = get_slice_fp8(&self, 0);
        let x1 = get_slice_fp8(&self, 1);
        let x2 = get_slice_fp8(&self, 2);
        let x3 = get_slice_fp8(&self, 3);
        let x4 = get_slice_fp8(&self, 4);
        let x5 = get_slice_fp8(&self, 5);
        let b02= x0.sqr();
        let b22= x1.sqr();
        let b42= x2.sqr();
        let b12= x3.sqr();
        let b32= x4.sqr();
        let b52 =x5.sqr();
        let b2pb5_2  = x1.addto(&x5).sqr();
        let b0pb3_2  = x0.addto(&x4).sqr();
        let b4pb1_2  = x2.addto(&x3).sqr();
        let a0 = b32.mulby_w().addto(&b02).mulbyu8(3u8).substract(&x0.addto(&x0));
        let a2 = b42.mulby_w().addto(&b12).mulbyu8(3u8).substract(&x1.addto(&x1));
        let a4 = b52.mulby_w().addto(&b22).mulbyu8(3u8).substract(&x2.addto(&x2));
        let a1 = b2pb5_2.substract(&b22).substract(&b52).mulbyu8(3u8).mulby_w().addto(&x3.addto(&x3));
        let a3 = b0pb3_2.substract(&b02).substract(&b32).mulbyu8(3u8).addto(&x4.addto(&x4));
        let a5 = b4pb1_2.substract(&b42).substract(&b12).mulbyu8(3u8).addto(&x5.addto(&x5)); 
        let mut result:[FieldElement<N>;48] = [self.content[0];48];
        result[..8].copy_from_slice(&a0.content);
        result[8..16].copy_from_slice(&a2.content);
        result[16..24].copy_from_slice(&a4.content);
        result[24..32].copy_from_slice(&a1.content);
        result[32..40].copy_from_slice(&a3.content);
        result[40..].copy_from_slice(&a5.content);        
        Self {  content : result,
                constants : self.constants}                                                                                                                                                    
    }
    
    pub fn cyclotomic_power(&self, e:& dyn Exponent<N>, negative :bool, naf_repre:&Option<Vec<i8>>) -> Self {
        let one = FieldElement{ mont_limbs:self.content[0].fieldparams.one,
                                                 fieldparams:self.content[0].fieldparams};
        let zero = FieldElement{ mont_limbs:[0;N],
                                                  fieldparams:self.content[0].fieldparams};        
        let mut result: [FieldElement<N>; 48] =[zero;48] ;
        result [0] = one; 
        let mut result = Self::new(&result,self.constants_interface());
        if naf_repre.is_none() { if let Some(array) = e.to_u64_array() 
                                        { let limbnum= e.get_len();
                                          for i in array[0..limbnum].as_ref().iter().rev()  {
                                                    for j in (0..64).rev(){ result = result.unisqr();                
                                                                                 if (i >> j) & 1 == 1 { result = result.multiply(&self);}                
                                                                                 }
                                                                            }                                                        
                                        }                                                                      
                            }
        else {  let selfconj=self.conjugate();
                for i in (*naf_repre).clone().unwrap() { result = result.unisqr();
                                                             if i == 1 {result = result.multiply(&self)};
                                                             if i ==-1 {result = result.multiply(&selfconj);}
                                                           }

             }
        if !negative { result} else { result.conjugate()}
    }

    pub fn final_exponentiation(&self, use_naf:bool) ->Self {      
        // It realy compute 3*(f^(p^48-1)/r) (the power of a pairings is a pairings)
        // Daiki Hayashida and Kenichiro Hayasaka and Tadanori Teruya https://eprint.iacr.org/2020/875.pdf
        let u  = self.constants.u.abs().try_into().unwrap();        
        let um1= (self.constants.u - 1).abs().try_into().unwrap();    
        let u_sign  = self.constants.u < 0; 
        let um1_sign= (self.constants.u - 1) < 0;
        let mut naf_u_repr  :Option<Vec<i8>> = None;
        let mut naf_um1_repr:Option<Vec<i8>> = None;
        if use_naf {naf_u_repr   = Some(<i128 as Exponent<N>>::to_naf(&u));
                    naf_um1_repr = Some(<i128 as Exponent<N>>::to_naf(&um1));
                   };

             //  Soft Part of the exponentiation :f^((p^24-1)*(p^8+1))                  
        let mut s = self.conjugate().multiply(&self.invert());               
        let mut t = s.frobinus(8);
        s = t.multiply(&s);
        
            // Hard part of the Exponentiation  :f^(3*(p^16-p^8+1)/r) = f^((u-1)^2*(u+p)*(u^2+p^2)*(u^4+p^4)*(u^8+p^8-1)+3) 
        let mut a = s.cyclotomic_power(&um1, um1_sign,&naf_um1_repr);               
        a = a.cyclotomic_power(&um1, um1_sign,&naf_um1_repr);
        t = a.frobinus(1);                
        a = t.multiply(&a.cyclotomic_power(&u, u_sign, &naf_u_repr)) ;               
        t = a.frobinus(2);         
        a = a.cyclotomic_power(&u, u_sign,&naf_u_repr);
        a = a.cyclotomic_power(&u, u_sign,&naf_u_repr);        
        a = t.multiply(&a);        
        t = a.frobinus(4);        
        a = a.cyclotomic_power(&u, u_sign,&naf_u_repr);
        a = a.cyclotomic_power(&u, u_sign,&naf_u_repr);
        a = a.cyclotomic_power(&u, u_sign,&naf_u_repr);
        a = a.cyclotomic_power(&u, u_sign,&naf_u_repr);
        a = t.multiply(&a);        
        t = a.frobinus(8).multiply(&a.conjugate());                
        a = a.cyclotomic_power(&u, u_sign,&naf_u_repr);
        a = a.cyclotomic_power(&u, u_sign,&naf_u_repr);
        a = a.cyclotomic_power(&u, u_sign,&naf_u_repr);
        a = a.cyclotomic_power(&u, u_sign, &naf_u_repr);
        a = a.cyclotomic_power(&u, u_sign,&naf_u_repr);
        a = a.cyclotomic_power(&u, u_sign,&naf_u_repr);
        a = a.cyclotomic_power(&u, u_sign,&naf_u_repr);
        a = a.cyclotomic_power(&u, u_sign,&naf_u_repr);                
        a = t.multiply(&a);
        a = a.multiply(&s).multiply(&s.unisqr()); 
        a
    }

    pub fn sparse_multiply(&self, rhs:&[&[FieldElement<N>];3]) -> Self {    
        let a0 = get_slice_fp24(&self, 0);
        let b0 = get_slice_fp24(&self, 1);

        let t0 = a0.sparse_multiply_for48(&rhs,2);        
        let t1 = b0.sparse_multiply_for48(&rhs,1);
        let t2 = t1.mulby_z();
        let y = [&[self.content[0].zero();8],
                                          rhs[1] ,
                                          &[rhs[0][0].addto(&rhs[2][0]),rhs[0][1].addto(&rhs[2][1]),rhs[0][2].addto(&rhs[2][2]),rhs[0][3].addto(&rhs[2][3]),
                                           rhs[0][4].addto(&rhs[2][4]),rhs[0][5].addto(&rhs[2][5]),rhs[0][6].addto(&rhs[2][6]),rhs[0][7].addto(&rhs[2][7]) ]
                                          ];
        let t3 = a0.addto(&b0).sparse_multiply_for48(&y, 1);
        let mut result:[FieldElement<N>;48] = [self.content[0];48] ;
        for i in 0..24 { result[i] = t2.content[i].addto(&t0.content[i])}
        for i in 24..48 { result[i] = t3.content[i-24].substract(&t0.content[i-24]).substract(&t1.content[i-24])}        
        Self {  content : result,
                constants : self.constants}                                                                                                                                                    
    }
    
}

impl<const PARAMSIZE:usize,const N: usize> fmt::Display for Fp48Element<PARAMSIZE,N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {  write!(f, "{:}", &self.to_a_string()) }
}

impl<const PARAMSIZE:usize,const N: usize> PartialEq for Fp48Element<PARAMSIZE,N> {
    fn eq(&self, other: &Self) -> bool {    self.equal(other) }
}

