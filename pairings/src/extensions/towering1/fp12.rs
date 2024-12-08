// Code developed by FARAOUN Kamel Mohamed.
// EEDIS-Laboratory. UDL-University. Algeria
// During May 2024.

use std::{fmt, usize};
use std::ops::Add;
use crate::tools::exponent::Exponent;
use super::super::super::fields::prime_fields::{FieldElement, PrimeField};
use super::super::ext_fields::{ExtElement, ExtField,ExFieldConsts};
use crate::tools::arithmetic_interface::ArithmeticOperations;
use super::fp2::Fp2Element;

pub mod fp6 {
    use super::super::fp2::Fp2Element;
    use super::super::super::super::fields::prime_fields::{FieldElement, PrimeField};
    use crate::tools::arithmetic_interface::ArithmeticOperations;
    use super::super::super::ext_fields::{ExtElement, ExtField,ExFieldConsts};

    #[derive(Clone, Copy)]
    pub struct Fp6Element<const PARAMSIZE:usize,const N:usize>{
                                        pub content :[FieldElement<N>;6],
                                        pub constants :Option<&'static ExFieldConsts<PARAMSIZE,N>>
                                        }
    #[derive(Clone)]
    pub struct Fp6Field<const PARAMSIZE:usize,const N:usize> {
                                        base_field:PrimeField<N>,
                                        constants :Option<&'static ExFieldConsts<PARAMSIZE,N>>
                                        }
    
    impl<const PARAMSIZE:usize,const N: usize> ExtField<PARAMSIZE,6,N> for Fp6Field<PARAMSIZE,N> {    
        type  ElementType   = Fp6Element<PARAMSIZE,N>;    
        type  BaseFieldType = PrimeField<N>;
    
        fn field_interface(&self)->PrimeField<N> 
            { self.base_field.clone() }
    
        fn new(base_field :&Self::BaseFieldType, consts :Option<&'static ExFieldConsts<PARAMSIZE,N>>)->  Fp6Field<PARAMSIZE,N>
            { Fp6Field {base_field : (*base_field).clone(), constants: consts}}   
        
        fn extconsts_interface(&self) ->Option<&'static ExFieldConsts<PARAMSIZE,N>> {
            self.constants
        }
        }       
    
    impl <const PARAMSIZE:usize,const N:usize> ExtElement<PARAMSIZE,6,N> for Fp6Element<PARAMSIZE,N>{
        fn new(content :&[FieldElement<N>; 6], consts :Option<&'static ExFieldConsts<PARAMSIZE,N>>) -> Fp6Element<PARAMSIZE,N>
            {   Fp6Element{content :content.clone(), constants :consts}    }
    
        fn content_interface(&self) -> &[FieldElement<N>;6]
            {   &self.content   }
        
        fn constants_interface(&self) -> Option<&'static ExFieldConsts<PARAMSIZE,N>> {
            self.constants
        }    
    
        fn sqr(&self) -> Self {
            let a = Fp2Element{content :[self.content[0],self.content[1]]};
            let b = Fp2Element{content :[self.content[2],self.content[3]]};
            let c = Fp2Element{content :[self.content[4],self.content[5]]};
            let t0 = a.sqr();
            let t1 = b.sqr();
            let t2 = c.sqr();
            let x0 = c.addto(&b).sqr().substract(&t1).substract(&t2).mul_by_u_p_1().addto(&t0);
            let x1 = a.addto(&b).sqr().substract(&t0).substract(&t1).addto(&t2.mul_by_u_p_1());
            let x2 = a.addto(&c).sqr().substract(&t0).substract(&t2).addto(&t1);                
            Self {content :[x0.content[0], x0.content[1], x1.content[0], x1.content[1], x2.content[0], x2.content[1]], 
                  constants : self.constants  }                                 
        }
        
    
        fn multiply(&self, rhs : &Self) -> Self {  
            let a0 = Fp2Element{content :[self.content[0],self.content[1]]};
            let b0 = Fp2Element{content :[self.content[2],self.content[3]]};
            let c0 = Fp2Element{content :[self.content[4],self.content[5]]};   
            let a1 = Fp2Element{content :[rhs.content[0],rhs.content[1]]};
            let b1 = Fp2Element{content :[rhs.content[2],rhs.content[3]]};
            let c1 = Fp2Element{content :[rhs.content[4],rhs.content[5]]};   
            let t0 = a0.multiply(&a1);
            let t1 = b0.multiply(&b1);
            let t2 = c0.multiply(&c1);
            let x0 = b0.addto(&c0).multiply(&b1.addto(&c1)).substract(&t1).substract(&t2).mul_by_u_p_1().addto(&t0);
            let x1 = a0.addto(&b0).multiply(&a1.addto(&b1)).substract(&t0).substract(&t1).addto(&t2.mul_by_u_p_1());
            let x2 = a0.addto(&c0).multiply(&a1.addto(&c1)).substract(&t0).substract(&t2).addto(&t1);
            Self {content :[x0.content[0], x0.content[1], x1.content[0], x1.content[1], x2.content[0], x2.content[1]], 
                constants : self.constants  }                                                                                                    
            }

     }
    
    impl <const PARAMSIZE:usize,const N:usize> Fp6Element <PARAMSIZE,N>{
        pub fn invert(&self) -> Self {
            let a = Fp2Element{content :[self.content[0],self.content[1]]};
            let b = Fp2Element{content :[self.content[2],self.content[3]]};
            let c = Fp2Element{content :[self.content[4],self.content[5]]}; 
            let t0 = a.sqr().substract(&b.multiply(&c).mul_by_u_p_1());
            let t1 = c.sqr().mul_by_u_p_1().substract(&a.multiply(&b));
            let t2 = b.sqr().substract(&a.multiply(&c));                
            let t = c.multiply(&t1).addto(&b.multiply(&t2)).mul_by_u_p_1().addto(&a.multiply(&t0)).invert();
            let x0 = t.multiply(&t0);
            let x1 = t.multiply(&t1); 
            let x2 = t.multiply(&t2);
            Self {content :[x0.content[0], x0.content[1], x1.content[0], x1.content[1], x2.content[0], x2.content[1]], 
                constants : self.constants  }   
            }
    
        pub fn mul_by_u(&self) -> Self{
                Self {content :[self.content[4].substract(&self.content[5]), self.content[4].addto(&self.content[5]),
                                self.content[0],self.content[1],self.content[2],self.content[3]], constants : None}
            }

        pub fn sparse_multiply(&self, rhs :&[&[FieldElement<N>];3], mode :u8) -> Self {  
            match mode { 0 => {
                                // rhs is sparse in Fp6 :(y0+y1*u)+(y2+y3*u)*v , (y4=y5=0)
                                let a0 = Fp2Element{content :[self.content[0],self.content[1]]};
                                let b0 = Fp2Element{content :[rhs[0][0],rhs[0][1]]};
                                let t0 =a0.multiply(&b0);
                                let a0 = Fp2Element{content :[self.content[2],self.content[3]]};
                                let b0 = Fp2Element{content :[rhs[1][0],rhs[1][1]]};
                                let t1 =a0.multiply(&b0);
                                let a0 = Fp2Element{content :[self.content[2].addto(&self.content[4]),self.content[3].addto(&self.content[5])]};
                                let c0 = a0.multiply(&b0).substract(&t1);
                                let res0 = c0.mul_by_u_p_1().addto(&t0);
                                let a0 = Fp2Element{content :[self.content[0].addto(&self.content[2]),self.content[1].addto(&self.content[3])]};
                                let b0 = Fp2Element{content :[rhs[0][0].addto(&rhs[1][0]),rhs[0][1].addto(&rhs[1][1])]};
                                let res1 =a0.multiply(&b0).substract(&t0).substract(&t1);
                                let a0 = Fp2Element{content :[self.content[0].addto(&self.content[4]),self.content[1].addto(&self.content[5])]};
                                let res2 = a0.multiply(&Fp2Element{content :[rhs[0][0],rhs[0][1]]}).substract(&t0).addto(&t1);
                                Self {  content :[res0.content[0], res0.content[1], res1.content[0], res1.content[1], res2.content[0], res2.content[1]], 
                                            constants : self.constants  }                                                                                                    
                                    }
                        1=> {
                              // rhs is sparse in Fp6 : (y2+y3*u)*v , (y0=y1=y4=y5=0)
                            let y2py3 = rhs[2][0].addto(&rhs[2][1]);
                            let tmp = self.content[4].addto(&self.content[5]).multiply(&y2py3);
                            let res0 = [self.content[4].multiply(&rhs[2][0]).double().substract(&tmp),
                                                              tmp.substract(&self.content[5].multiply(&rhs[2][1]).double())];
                             let t0 = self.content[0].multiply(&rhs[2][0]);
                             let t1 = self.content[1].multiply(&rhs[2][1]);
                             let res1 =[t0.substract(&t1), (self.content[0].addto(&self.content[1])).multiply(&y2py3).substract(&t0).substract(&t1)];
                             let t0 = self.content[2].multiply(&rhs[2][0]);
                             let t1 = self.content[3].multiply(&rhs[2][1]);
                             let res2 =[t0.substract(&t1), (self.content[2].addto(&self.content[3])).multiply(&y2py3).substract(&t0).substract(&t1)];
                             Self {  content :[res0[0], res0[1], res1[0], res1[1], res2[0], res2[1]], 
                                constants : self.constants  }  
                        }
                        _=> {panic!("Not implemented sparse multiplication mode for Fp12")}
            }
                
            }
            
        }

        
    }
#[derive(Clone, Copy)]
pub struct Fp12Element<const PARAMSIZE:usize,const N:usize>{
                                    pub content :[FieldElement<N>;12],
                                    pub constants :&'static ExFieldConsts<PARAMSIZE,N>
                                    }
#[derive(Clone,Debug)]
pub struct Fp12Field<const PARAMSIZE:usize,const N:usize> {
                                    base_field:PrimeField<N>,
                                    constants :&'static ExFieldConsts<PARAMSIZE,N>
                                    }

fn get_slice_fp6<const PARAMSIZE:usize,const N:usize>(element : &Fp12Element<PARAMSIZE,N>, i : usize) -> fp6::Fp6Element<PARAMSIZE,N>
    {       
        let _t: [FieldElement<N>; 6] = element.content[i*6..(i+1)*6].iter()
        .map(|&element1| FieldElement {fieldparams: element1.fieldparams,
                                                       mont_limbs: element1.mont_limbs,
                                                      })
        .collect::<Vec<_>>()
        .try_into().unwrap(); 
            fp6::Fp6Element{content :_t, constants :None}        
    }


impl<const PARAMSIZE:usize,const N: usize> ExtField<PARAMSIZE,12,N> for Fp12Field<PARAMSIZE,N> {    
    type  ElementType   = Fp12Element<PARAMSIZE,N>;    
    type  BaseFieldType = PrimeField<N>;

    fn field_interface(&self)->PrimeField<N> 
        { self.base_field.clone() }

    fn new(base_field :&Self::BaseFieldType, consts :Option<&'static ExFieldConsts<PARAMSIZE,N>>)->  Fp12Field<PARAMSIZE,N>
        { Fp12Field {base_field : (*base_field).clone(), constants: consts.unwrap()}}   
    
    fn extconsts_interface(&self) ->Option<&'static ExFieldConsts<PARAMSIZE,N>> {
        Some(self.constants)
    }
    }

impl <const PARAMSIZE:usize,const N:usize> ExtElement<PARAMSIZE,12,N> for Fp12Element<PARAMSIZE,N>{
    fn new(content :&[FieldElement<N>; 12], consts :Option<&'static ExFieldConsts<PARAMSIZE,N>>) -> Fp12Element<PARAMSIZE,N>
        {   Fp12Element{content :content.clone(), constants :consts.unwrap()}    }

    fn content_interface(&self) -> &[super::super::super::fields::prime_fields::FieldElement<N>;12]
        {   &self.content   }
    
    fn constants_interface(&self) -> Option<&'static ExFieldConsts<PARAMSIZE,N>> {
        Some(self.constants)
    }
    
    fn sqr(&self) -> Self {
        let a = get_slice_fp6(&self, 0);
        let b = get_slice_fp6(&self, 1);
        let t0 = a.sqr();            
        let t1 = b.sqr();
        let t2 = a.addto(&b).sqr();              
        let x0 = t1.mul_by_u().addto(&t0);
        let x1 = t2.substract(&t0).substract(&t1);
        Self {  content : [x0.content[0], x0.content[1],x0.content[2],x0.content[3],x0.content[4],x0.content[5],
                           x1.content[0], x1.content[1],x1.content[2],x1.content[3],x1.content[4],x1.content[5]],
                constants : self.constants}                                                                                                                                                    
    }

    fn multiply(&self, rhs:&Self) -> Self {          
        let a0 = get_slice_fp6(&self, 0);
        let b0 = get_slice_fp6(&self, 1);
        let a1 = get_slice_fp6(&rhs, 0);
        let b1 = get_slice_fp6(&rhs, 1);
        let t0 = a0.multiply(&a1);
        let t1 = b0.multiply(&b1);
        let t3 = a0.addto(&b0).multiply(&a1.addto(&b1));
        let x0 = t1.mul_by_u().addto(&t0);
        let x1 = t3.substract(&t0).substract(&t1);
        Self {  content : [x0.content[0], x0.content[1],x0.content[2],x0.content[3],x0.content[4],x0.content[5],
                           x1.content[0], x1.content[1],x1.content[2],x1.content[3],x1.content[4],x1.content[5]],
                constants : self.constants}  
    }
}

impl <const PARAMSIZE:usize,const N:usize> Fp12Element <PARAMSIZE,N>{

    pub fn conjugate(&self) -> Self {
        Self { content : [self.content[0],self.content[1],self.content[2],self.content[3],self.content[4],self.content[5],
                          self.content[6].negate(),self.content[7].negate(),self.content[8].negate(),
                          self.content[9].negate(),self.content[10].negate(),self.content[11].negate()],
               constants :self.constants }
    }

    pub fn frobinus(&self, order :u8) -> Self {        
        match  order { 1 =>{ Self { content : [self.content[0], 
                                                self.content[1].negate(),
                                                self.content[3].multiply(&self.constants.frobinus_consts[0]),
                                                self.content[2].multiply(&self.constants.frobinus_consts[0]),
                                                self.content[4].multiply(&self.constants.frobinus_consts[0].add(1u64)),
                                                self.content[5].negate().multiply(&self.constants.frobinus_consts[0].add(1u64)),
                                                self.content[6].substract(&self.content[7]).multiply(&self.constants.frobinus_consts[1]),
                                                self.content[6].addto(&self.content[7]).negate().multiply(&self.constants.frobinus_consts[1]),                                                
                                                self.content[8].addto(&self.content[9]).multiply(&self.constants.frobinus_consts[2]),
                                                self.content[8].substract(&self.content[9]).multiply(&self.constants.frobinus_consts[2]),
                                                self.content[10].substract(&self.content[11]).multiply(&self.constants.frobinus_consts[3]),
                                                self.content[10].addto(&self.content[11]).negate().multiply(&self.constants.frobinus_consts[3])
                                                ],
                                    constants :self.constants }
                            },
                       2 =>{Self { content : [self.content[0], 
                                              self.content[1],
                                              self.content[2].negate().multiply(&self.constants.frobinus_consts[0].add(1u64)),
                                              self.content[3].negate().multiply(&self.constants.frobinus_consts[0].add(1u64)),
                                              self.content[4].multiply(&self.constants.frobinus_consts[0]),
                                              self.content[5].multiply(&self.constants.frobinus_consts[0]),
                                              self.content[6].negate().multiply(&self.constants.frobinus_consts[0]),
                                              self.content[7].negate().multiply(&self.constants.frobinus_consts[0]),                                                
                                              self.content[8].negate(),
                                              self.content[9].negate(),
                                              self.content[10].multiply(&self.constants.frobinus_consts[0].add(1u64)),
                                              self.content[11].multiply(&self.constants.frobinus_consts[0].add(1u64))
                                                ],
                                   constants :self.constants }}
                       3 =>{Self { content : [self.content[0], 
                                              self.content[1].negate(),
                                              self.content[3],
                                              self.content[2],
                                              self.content[4].negate(),
                                              self.content[5],
                                              self.content[7].substract(&self.content[6]).multiply(&self.constants.frobinus_consts[2]),
                                              self.content[7].addto(&self.content[6]).multiply(&self.constants.frobinus_consts[2]),                                                
                                              self.content[8].addto(&self.content[9]).negate().multiply(&self.constants.frobinus_consts[2]),
                                              self.content[9].substract(&self.content[8]).multiply(&self.constants.frobinus_consts[2]),
                                              self.content[10].substract(&self.content[11]).multiply(&self.constants.frobinus_consts[2]),
                                              self.content[10].addto(&self.content[11]).negate().multiply(&self.constants.frobinus_consts[2]),
                                                ],
                                  constants :self.constants }}
                       _=>{ panic!("Frobinus not implemented for order {}.",order  )}            
        }
    }

    pub fn invert(&self) -> Self {        
        let a = get_slice_fp6(&self, 0); 
        let b = get_slice_fp6(&self, 1); 
        let t = a.sqr().substract(&b.sqr().mul_by_u()).invert();
        let x0 = a.multiply(&t);
        let x1 = b.negate().multiply(&t);
        Self {  content : [x0.content[0], x0.content[1],x0.content[2],x0.content[3],x0.content[4],x0.content[5],
                           x1.content[0], x1.content[1],x1.content[2],x1.content[3],x1.content[4],x1.content[5]],
                constants : self.constants}  
    }
    
    pub fn unisqr(&self) -> Self { 
        let a0= Fp2Element{content :[self.content[0],self.content[1]]};
        let a1= Fp2Element{content :[self.content[2],self.content[3]]};
        let a2= Fp2Element{content :[self.content[4],self.content[5]]};
        let a3= Fp2Element{content :[self.content[6],self.content[7]]};
        let a4= Fp2Element{content :[self.content[8],self.content[9]]};
        let a5= Fp2Element{content :[self.content[10],self.content[11]]};
        let x0 = a0.sqr();
        let x1 = a4.sqr();
        let x2 = a0.addto(&a4);
        let y0 = a1.sqr();
        let y1 = a5.sqr();
        let y2 = a1.addto(&a5);
        let z0 = a3.sqr();
        let z1 = a2.sqr();
        let z2 = a3.addto(&a2);
        let t0 = x0.addto(&x1.mul_by_u_p_1()).mulbyu8(3u8).substract(&a0.addto(&a0));
        let t1 = z0.addto(&z1.mul_by_u_p_1()).mulbyu8(3u8).substract(&a1.addto(&a1));
        let t2 = y0.addto(&y1.mul_by_u_p_1()).mulbyu8(3u8).substract(&a2.addto(&a2));
        let t3 = y2.sqr().substract(&y0).substract(&y1).mul_by_u_p_1().mulbyu8(3u8).addto(&a3.addto(&a3));
        let t4 = x2.sqr().substract(&x0).substract(&x1).mulbyu8(3u8).addto(&a4.addto(&a4));
        let t5 = z2.sqr().substract(&z0).substract(&z1).mulbyu8(3u8).addto(&a5.addto(&a5));
        Self { content : [t0.content[0],t0.content[1],t1.content[0],t1.content[1],t2.content[0],t2.content[1],
                          t3.content[0],t3.content[1],t4.content[0],t4.content[1],t5.content[0],t5.content[1],],
               constants : self.constants }
    }
    
    pub fn cyclotomic_power(&self, e:& dyn Exponent<N>, negative:bool, naf_repre:&Option<Vec<i8>>) -> Self {
        // Efficient powering inside cyclotomic sub-group (x^(p^6+1)=1) -- Scott Approache           
        let one = FieldElement{ mont_limbs:self.content[0].fieldparams.one,
                                                 fieldparams:self.content[0].fieldparams};
        let zero = FieldElement{ mont_limbs:[0;N],
                                                  fieldparams:self.content[0].fieldparams};        
        let mut result: [FieldElement<N>; 12] =[zero;12] ;
        result [0] = one; 
        let mut result = Self::new(&result,self.constants_interface());
        if naf_repre.is_none() { if let Some(array) = e.to_u64_array() 
                                        { let limbnum= e.get_len();
                                          for i in array[0..limbnum].as_ref().iter().rev()  {
                                                    for j in (0..64).rev(){ result = result.unisqr();                
                                                                                 if (i >> j) & 1 == 1 {  result = result.multiply(&self);}                
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
             if negative {result.conjugate()} else {result}
    }

    pub fn final_exponentiation(&self, use_naf:bool) ->Self {
        // It really compute 3*(f^(p^12-1)/r) (the power of a pairings is a pairings)
        // Daiki Hayashida and Kenichiro Hayasaka and Tadanori Teruya : https://eprint.iacr.org/2020/875.pdf
        let u   = self.constants.u.abs().try_into().unwrap();        
        let um1 = (self.constants.u - 1).abs().try_into().unwrap();
        let um1_sign = (self.constants.u - 1) < 0;
        let u_sign   = self.constants.u < 0;
        let mut naf_u_repr   :Option<Vec<i8>> = None;
        let mut naf_um1_repr :Option<Vec<i8>> = None;
        if use_naf { naf_u_repr   = Some(<u128 as Exponent<N>>::to_naf(&u));
                     naf_um1_repr = Some(<u128 as Exponent<N>>::to_naf(&um1));
                   };
            // Soft exponentiation part :f^((p^6-1)*(p^2+1))      
        let mut s = self.conjugate().multiply(&self.invert());
        let mut t = s.frobinus(2);
        s = t.multiply(&s);       
        
            // Hard part of the Exponentiation  :  f^(3*(p^4-p^2+1)/r) = f^((u-1)^2*(u+p)*(u^2+p^2-1)+3)         
        let mut a = s.cyclotomic_power(&um1, um1_sign, &naf_um1_repr);
        a = a.cyclotomic_power(&um1, um1_sign, &naf_um1_repr);
        t = a.frobinus(1);
        a = t.multiply(&a.cyclotomic_power(&u, u_sign, &naf_u_repr)) ;
        t = a.frobinus(2).multiply(&a.conjugate()); 
        a = a.cyclotomic_power(&u, u_sign, &naf_u_repr);
        a = a.cyclotomic_power(&u,u_sign, &naf_u_repr);
        a = t.multiply(&a);
        a = a.multiply(&s).multiply(&s.unisqr()); 
        a
    }
   
    pub fn sparse_multiply(&self, rhs:&[&[FieldElement<N>];3]) -> Self {     
    //   Multiplication with a sparse Fp12 according to M - type twiste   
        let a0 = get_slice_fp6(&self, 0);
        let b0 = get_slice_fp6(&self, 1);        
        let t0 = a0.sparse_multiply(&rhs,0);        
        let t1 = b0.sparse_multiply(&rhs,1);                
        let t2 = a0.addto(&b0).sparse_multiply(&[ &[rhs[0][0],rhs[0][1]],
                                                                                &[rhs[1][0].addto(&rhs[2][0]),rhs[1][1].addto(&rhs[2][1])],
                                                                                &[rhs[0][0].zero(),rhs[0][0].zero()]], 0);                
        let x0 = t1.mul_by_u().addto(&t0);                                      
        let x1 = t2.substract(&t0).substract(&t1);
        Self {  content : [x0.content[0], x0.content[1],x0.content[2],x0.content[3],x0.content[4],x0.content[5],
                           x1.content[0], x1.content[1],x1.content[2],x1.content[3],x1.content[4],x1.content[5]],
                constants : self.constants}  
    }    
}


impl<const PARAMSIZE:usize,const N: usize> fmt::Display for Fp12Element<PARAMSIZE,N> {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {  write!(f, "{:}", &self.to_a_string()) }
    }

impl<const PARAMSIZE:usize,const N: usize> PartialEq for Fp12Element<PARAMSIZE,N> {
        fn eq(&self, other: &Self) -> bool {    self.equal(other) }
    }   
