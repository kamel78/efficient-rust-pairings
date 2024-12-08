// Code developed by FARAOUN Kamel Mohamed.
// EEDIS-Laboratory. UDL-University. Algeria
// During May 2024.

use std::{fmt, usize};
use std::ops::Add;
use crate::tools::exponent::Exponent;

use super::super::super::fields::prime_fields::{FieldElement, PrimeField};
use super::super::ext_fields::{ExtElement, ExtField,ExFieldConsts};
use crate::tools::arithmetic_interface::ArithmeticOperations;
use super::fp8::Fp8Element;
use super::fp4::Fp4Element;


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

fn get_slice_fp4<const PARAMSIZE:usize,const N:usize>(element : &Fp24Element<PARAMSIZE,N>, i : usize) -> Fp4Element<PARAMSIZE,N>
    {       
        let _t: [FieldElement<N>; 4] = element.content[i*4..(i+1)*4].iter()
        .map(|&element1| FieldElement {fieldparams: element1.fieldparams,
                                                       mont_limbs: element1.mont_limbs,
                                                      })
        .collect::<Vec<_>>()
        .try_into().unwrap(); 
        Fp4Element{content :_t, constants :element.constants}        
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
        let h0 = b.addto(&c).sqr().substract(&t1).substract(&t2).mulby_w().addto(&t0);
        let h1 = t2.mulby_w().addto(&a.addto(&b).sqr().substract(&t0).substract(&t1));
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
        let h0 = b0.addto(&c0).multiply(&b1.addto(&c1)).substract(&t1).substract(&t2).mulby_w().addto(&t0);        
        let h1 = t2.mulby_w().addto(&a0.addto(&b0).multiply(&a1.addto(&b1)).substract(&t0).substract(&t1));
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

    pub fn conjugate(&self) -> Self {
        Self {content : [self.content[0], self.content[1], self.content[2], self.content[3],
                         self.content[4].negate(), self.content[5].negate(), self.content[6].negate(), self.content[7].negate(),
                         self.content[8].negate(), self.content[9].negate(), self.content[10].negate(), self.content[11].negate(),
                         self.content[12], self.content[13], self.content[14], self.content[15], 
                         self.content[16], self.content[17],self.content[18], self.content[19], 
                         self.content[20].negate(), self.content[21].negate(), self.content[22].negate(), self.content[23].negate()],
              constants :self.constants  }
    }

    pub fn frobinus(&self, order :u8) -> Self {        
        match  order {   
          1 =>{ Self { content : [self.content[0], 
                                  self.content[1].negate(),
                                  self.content[2].addto(&self.content[3]).multiply(&self.constants.frobinus_consts[0]),
                                  self.content[2].substract(&self.content[3]).multiply(&self.constants.frobinus_consts[0]),
                                  self.content[6].negate().multiply(&self.constants.frobinus_consts[1]),
                                  self.content[7].multiply(&self.constants.frobinus_consts[1]),
                                  self.content[5].multiply(&self.constants.frobinus_consts[2]),
                                  self.content[4].multiply(&self.constants.frobinus_consts[2]),                                                
                                  self.content[11].addto(&self.content[10]).multiply(&self.constants.frobinus_consts[3]),
                                  self.content[10].substract(&self.content[11]).multiply(&self.constants.frobinus_consts[3]),
                                  self.content[8].substract(&self.content[9]).multiply(&self.constants.frobinus_consts[4]),
                                  self.content[8].addto(&self.content[9]).negate().multiply(&self.constants.frobinus_consts[4]),
                                  self.content[13].multiply(&self.constants.frobinus_consts[5]),
                                  self.content[12].multiply(&self.constants.frobinus_consts[5]),
                                  self.content[15].substract(&self.content[14]).multiply(&self.constants.frobinus_consts[6]),
                                  self.content[15].addto(&self.content[14]).multiply(&self.constants.frobinus_consts[6]),
                                  self.content[16].substract(&self.content[17]).multiply(&self.constants.frobinus_consts[7]),
                                  self.content[16].addto(&self.content[17]).negate().multiply(&self.constants.frobinus_consts[7]),
                                  self.content[18].multiply(&self.constants.frobinus_consts[5].add(1u64)),
                                  self.content[19].negate().multiply(&self.constants.frobinus_consts[5].add(1u64)),
                                  self.content[23].substract(&self.content[22]).multiply(&self.constants.frobinus_consts[8]),
                                  self.content[23].addto(&self.content[22]).multiply(&self.constants.frobinus_consts[8]),
                                  self.content[20].addto(&self.content[21]).multiply(&self.constants.frobinus_consts[9]),
                                  self.content[20].substract(&self.content[21]).multiply(&self.constants.frobinus_consts[9]),
                                  ],
                      constants :self.constants }
              },
         2 =>{Self { content : [self.content[0], 
                                self.content[1],
                                self.content[2].negate(),
                                self.content[3].negate(),
                                self.content[5].negate(),
                                self.content[4],
                                self.content[7],
                                self.content[6].negate(),                                                
                                self.content[9].negate().multiply(&self.constants.frobinus_consts[5].add(1u64)),
                                self.content[8].multiply(&self.constants.frobinus_consts[5].add(1u64)),
                                self.content[11].multiply(&self.constants.frobinus_consts[5].add(1u64)),
                                self.content[10].negate().multiply(&self.constants.frobinus_consts[5].add(1u64)),
                                self.content[12].negate().multiply(&self.constants.frobinus_consts[5].add(1u64)),
                                self.content[13].negate().multiply(&self.constants.frobinus_consts[5].add(1u64)),
                                self.content[14].multiply(&self.constants.frobinus_consts[5].add(1u64)),
                                self.content[15].multiply(&self.constants.frobinus_consts[5].add(1u64)),                                
                                self.content[16].negate().multiply(&self.constants.frobinus_consts[5]),
                                self.content[17].negate().multiply(&self.constants.frobinus_consts[5]),
                                self.content[18].multiply(&self.constants.frobinus_consts[5]),
                                self.content[19].multiply(&self.constants.frobinus_consts[5]),
                                self.content[21].multiply(&self.constants.frobinus_consts[5]),
                                self.content[20].negate().multiply(&self.constants.frobinus_consts[5]),
                                self.content[23].negate().multiply(&self.constants.frobinus_consts[5]),
                                self.content[22].multiply(&self.constants.frobinus_consts[5]),
                                  ],
                     constants :self.constants }}
         4 =>{Self { content : [self.content[0], 
                                self.content[1],
                                self.content[2],
                                self.content[3],
                                self.content[4].negate(),
                                self.content[5].negate(),
                                self.content[6].negate(),
                                self.content[7].negate(),                                                
                                self.content[8].negate().multiply(&self.constants.frobinus_consts[5]),
                                self.content[9].negate().multiply(&self.constants.frobinus_consts[5]),
                                self.content[10].negate().multiply(&self.constants.frobinus_consts[5]),
                                self.content[11].negate().multiply(&self.constants.frobinus_consts[5]),
                                self.content[12].multiply(&self.constants.frobinus_consts[5]),
                                self.content[13].multiply(&self.constants.frobinus_consts[5]),
                                self.content[14].multiply(&self.constants.frobinus_consts[5]),
                                self.content[15].multiply(&self.constants.frobinus_consts[5]),
                                self.content[16].negate().multiply(&self.constants.frobinus_consts[5].add(1u64)),
                                self.content[17].negate().multiply(&self.constants.frobinus_consts[5].add(1u64)),
                                self.content[18].negate().multiply(&self.constants.frobinus_consts[5].add(1u64)),
                                self.content[19].negate().multiply(&self.constants.frobinus_consts[5].add(1u64)),
                                self.content[20].multiply(&self.constants.frobinus_consts[5].add(1u64)),
                                self.content[21].multiply(&self.constants.frobinus_consts[5].add(1u64)),
                                self.content[22].multiply(&self.constants.frobinus_consts[5].add(1u64)),
                                self.content[23].multiply(&self.constants.frobinus_consts[5].add(1u64)),
                                  ],
                    constants :self.constants }}
         _=>{ panic!("Frobinus not implemented for order {}.",order  )}   
    }
}
    pub fn invert(&self) -> Self {
        let a = get_slice_fp8(&self, 0);
        let b = get_slice_fp8(&self, 1);
        let c = get_slice_fp8(&self, 2);
        let t0 = a.sqr().substract(&b.multiply(&c).mulby_w());
        let t1 = c.sqr().mulby_w().substract(&a.multiply(&b));
        let t2 =b.sqr().substract(&a.multiply(&c));                                    
        let t3 = c.multiply(&t1).addto(&b.multiply(&t2)).mulby_w().addto(&a.multiply(&t0)).invert();
        let mut result =[FieldElement{  mont_limbs:self.content[0].fieldparams.zero,
                                                               fieldparams:self.content[0].fieldparams};24];
        result[..8].copy_from_slice(&t0.multiply(&t3).content);
        result[8..16].copy_from_slice(&t1.multiply(&t3).content);
        result[16..].copy_from_slice(&t2.multiply(&t3).content);
        Self {  content : result, constants :self.constants} 
    }

    pub fn unisqr(&self) -> Self {
       let a0 =  get_slice_fp4(&self, 0);
       let a1 =  get_slice_fp4(&self, 1);
       let a2 =  get_slice_fp4(&self, 2);
       let a3 =  get_slice_fp4(&self, 3);
       let a4 =  get_slice_fp4(&self, 4);
       let a5 =  get_slice_fp4(&self, 5);
       let t0 = a0.sqr();
       let t1 = a1.sqr();
       let t2 = a0.addto(&a1).sqr();
       let h00= t1.mulby_v().addto(&t0).mulbyu8(3u8).substract(&a0.addto(&a0));
       let h01= t2.substract(&t0).substract(&t1).mulbyu8(3u8).addto(&a1.addto(&a1));
       let t0 = a4.sqr();
       let t1 = a5.sqr();
       let t2 = a4.addto(&a5).sqr();
       let h10= t2.substract(&t0).substract(&t1).mulby_v().mulbyu8(3u8).addto(&a2.addto(&a2));
       let h11= t1.mulby_v().addto(&t0).mulbyu8(3u8).substract(&a3.addto(&a3));
       let t0 = a2.sqr();
       let t1 = a3.sqr();
       let t2 = a2.addto(&a3).sqr();
       let h20= t1.mulby_v().addto(&t0).mulbyu8(3u8).substract(&a4.addto(&a4));
       let h21= t2.substract(&t0).substract(&t1).mulbyu8(3u8).addto(&a5.addto(&a5));
       let mut result =[FieldElement{  mont_limbs:self.content[0].fieldparams.zero,
                                                              fieldparams:self.content[0].fieldparams};24];
       result[..4].copy_from_slice(&h00.content);
       result[4..8].copy_from_slice(&h01.content);
       result[8..12].copy_from_slice(&h10.content);
       result[12..16].copy_from_slice(&h11.content);
       result[16..20].copy_from_slice(&h20.content);
       result[20..].copy_from_slice(&h21.content);
       Self {  content : result, constants :self.constants} 
    }
    
    pub fn cyclotomic_power(&self, e:& dyn Exponent<N>, negative :bool,naf_repre:&Option<Vec<i8>>) -> Self {
        let one = FieldElement{mont_limbs:self.content[0].fieldparams.one,
                                                fieldparams:self.content[0].fieldparams};
        let zero = FieldElement{mont_limbs:[0;N],
                                                 fieldparams:self.content[0].fieldparams};        
        let mut result: [FieldElement<N>; 24] =[zero;24] ;
        result [0] = one; 
        let mut result = Self::new(&result,self.constants_interface());
        if naf_repre.is_none() { if let Some(array) = e.to_u64_array() 
                                        { let limbnum= e.get_len();
                                          for i in array[0..limbnum].as_ref().iter().rev()  {
                                                    for j in (0..64).rev(){ result = result.unisqr();                
                                                                                 if (i >> j) & 1 == 1 {result = result.multiply(&self);}                
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
        if !negative { result}
        else { result.conjugate()}

    }

    pub fn final_exponentiation(&self, use_naf:bool) ->Self {
        // It actually computes 3 * (f^(p^24-1)/r) (the power of a pairings is a pairings)        
        // Daiki Hayashida, Kenichiro Hayasaka, Tadanori Teruya : Link: https://eprint.iacr.org/2020/875.pdf
        let u = self.constants.u.abs().try_into().unwrap();
        let u_sign = self.constants.u <0; 
        let um1  = (self.constants.u - 1).abs().try_into().unwrap();    
        let um1_sign = (self.constants.u - 1) < 0;
        let mut naf_u_repr:Option<Vec<i8>> = None;
        let mut naf_um1_repr:Option<Vec<i8>> = None;
        if use_naf {naf_u_repr = Some(<i128 as Exponent<N>>::to_naf(&u));
                    naf_um1_repr = Some(<i128 as Exponent<N>>::to_naf(&um1));};

                // Soft part of the exponentiation: f^((p^12-1)*(p^4+1))
        let mut s = self.conjugate().multiply(&self.invert());        
        let mut t = s.frobinus(4);
        s = t.multiply(&s);        

                // Hard part of the Exponentiation  : f^(3*(p^8-p^4+1)/r) = f^((u-1)^2*(u+p)*(u^2+p^2)*(u^4+p^4-1)+3) 
        let mut a = s.cyclotomic_power(&um1, um1_sign, &naf_um1_repr);
        a = a.cyclotomic_power(&um1, um1_sign,&naf_um1_repr);                
        t = a.frobinus(1);
        a = t.multiply(&a.cyclotomic_power(&u, u_sign,&naf_u_repr)) ;       
        t = a.frobinus(2);         
        a = a.cyclotomic_power(&u, u_sign,&naf_u_repr);
        a = a.cyclotomic_power(&u, u_sign,&naf_u_repr);
        a = t.multiply(&a);
        t = a.frobinus(4).multiply(&a.conjugate());
        a = a.cyclotomic_power(&u, u_sign, &naf_u_repr);
        a = a.cyclotomic_power(&u, u_sign,&naf_u_repr);
        a = a.cyclotomic_power(&u, u_sign,&naf_u_repr);
        a = a.cyclotomic_power(&u, u_sign, &naf_u_repr);
        a = t.multiply(&a);
        a = a.multiply(&s).multiply(&s.unisqr());
        a
    }

    pub fn mulby_z(&self) -> Self{
        let mut result =self.content.clone();
        result[0] = self.content[22].substract(&self.content[23]);
        result[1] = self.content[22].addto(&self.content[23]);
        result[2] = self.content[20];
        result[3] = self.content[21];
        result[4..8].copy_from_slice(&self.content[16..20]);
        result[8..24].copy_from_slice(&self.content[0..16]);
        Fp24Element {content : result,
                     constants : self.constants}
    }

    pub fn sparse_multiply(&self, rhs:&[&[FieldElement<N>];3],twist_type:char) -> Self {
            let mut result =[FieldElement{ mont_limbs:self.content[0].fieldparams.zero,
                                                                  fieldparams:self.content[0].fieldparams};24];        
            let res1 :  Fp8Element<PARAMSIZE, N>;     
            let res2 :  Fp8Element<PARAMSIZE, N>;     
            let res3 :  Fp8Element<PARAMSIZE, N>;     
            if twist_type =='D'{
                //  D-type sparse multiplication of an Fp24 element with a sparse element 
                //  used during pairings on BLS24 () , rhs on 12 Fp element 
            let a0 = get_slice_fp8(&self, 0);        
            let b0 = get_slice_fp8(&self, 1);
            let c0 = get_slice_fp8(&self, 2);
            let y0 = Fp8Element{ content :[rhs[0][0],rhs[0][1],rhs[0][2],rhs[0][3],rhs[1][0],rhs[1][1],rhs[1][2],rhs[1][3]], constants:self.constants};
            let y1 = Fp8Element{ content :[rhs[2][0],rhs[2][1],rhs[2][2],rhs[2][3],rhs[2][0].zero(),rhs[2][0].zero(),rhs[2][0].zero(),rhs[2][0].zero()], 
                                                            constants:self.constants};
            let t0 = a0.multiply(&y0);
            let t1 = b0.sparse_multiply(&rhs[2]);
            let t2 = c0.sparse_multiply(&rhs[2]);
            res1 = t2.mulby_w().addto(&t0);
            res2 = a0.addto(&b0).multiply(&y0.addto(&y1)).substract(&t0).substract(&t1);
            res3 = c0.multiply(&y0).addto(&t1);                       
            }

            else {  //  M-type sparse multiplication of an Fp24 element with a sparse element 
                    //  used during pairings on BLS24 () , rhs on 12 Fp element 
                    let a0 = get_slice_fp8(&self, 0);
                    let a1 = Fp8Element{ content :[rhs[0][0],rhs[0][1],rhs[0][2],rhs[0][3],rhs[2][0],rhs[2][1],rhs[2][2],rhs[2][3]], constants:self.constants};
                    let b0 = get_slice_fp8(&self, 1);
                    let c0 = get_slice_fp8(&self, 2);
                    let c1 = Fp8Element{ content :[rhs[1][0],rhs[1][1],rhs[1][2],rhs[1][3],rhs[1][0].zero(),
                                                                             rhs[1][0].zero(),rhs[1][0].zero(),rhs[1][0].zero()], constants:self.constants};
                    let t0 = a0.multiply(&a1);
                    let t2 = c0.sparse_multiply(&[rhs[1][0],rhs[1][1],rhs[1][2],rhs[1][3]]);
                    res1 = b0.sparse_multiply(&[rhs[1][0],rhs[1][1],rhs[1][2],rhs[1][3]]).mulby_w().addto(&t0);        
                    res2 = t2.mulby_w().addto(&b0.multiply(&a1));
                    res3 = a0.addto(&c0).multiply(&a1.addto(&c1)).substract(&t0).substract(&t2);                
            }
            result[..8].copy_from_slice(&res1.content);
            result[8..16].copy_from_slice(&res2.content);
            result[16..].copy_from_slice(&res3.content);
            Self {  content : result, constants :self.constants}     
    }
                
    pub fn sparse_multiply_for48(&self, rhs:&[&[FieldElement<N>];3], mode :u8) -> Self {         
        match  mode { 1 => { //  Sparse multiplication of an Fp24 element with a sparse element on two consecutive Fp8 only
                             //  Used during pairings on BLS48
                            let a0 = get_slice_fp8(&self, 0);        
                            let b0 = get_slice_fp8(&self, 1);
                            let c0 = get_slice_fp8(&self, 2);
                            let a1 = Fp8Element{ content: rhs[0].try_into().unwrap(), constants: self.constants };
                            let b1 = Fp8Element{ content: rhs[1].try_into().unwrap(), constants: self.constants };
                            let t0 = a0.multiply(&a1);
                            let t1 = c0.multiply(&b1);
                            let t2 = b0.multiply(&b1);
                            let t3 = a0.addto(&b0).multiply(&a1.addto(&b1));
                            let res1 = t1.mulby_w().addto(&t0);
                            let res2 = t3.substract(&t0).substract(&t2);
                            let res3 = t2.addto(&c0.multiply(&a1)); 
                            let mut result =[FieldElement{  mont_limbs:self.content[0].fieldparams.zero,
                                                                                   fieldparams:self.content[0].fieldparams};24];
                            result[..8].copy_from_slice(&res1.content);
                            result[8..16].copy_from_slice(&res2.content);
                            result[16..].copy_from_slice(&res3.content);
                            Self {  content : result, constants :self.constants}  
                            }      
                      2 => { //  Sparse multiplication of an Fp24 element with a sparse element on two consecutive Fp8 only
                             //  Used during pairings on BLS48, rhs on 4 Fp elements
                            let a0 = get_slice_fp8(&self, 0);        
                            let b0 = get_slice_fp8(&self, 1);
                            let c0 = get_slice_fp8(&self, 2);
                            let y =Fp8Element {content : rhs[2].try_into().unwrap(),constants :self.constants};
                            let res1 = c0.multiply(&y).mulby_w();
                            let res2 = a0.multiply(&y);
                            let res3 = b0.multiply(&y);
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
