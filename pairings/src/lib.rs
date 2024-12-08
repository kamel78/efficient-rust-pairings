// Code developed by FARAOUN Kamel Mohamed.
// EEDIS-Laboratory. UDL-University. Algeria
// During May 2024.

pub mod fields;
pub mod extensions;
pub mod tools;
pub mod parameters;
pub mod curves;
pub mod engines;
pub mod tests;


use curves::{curve_arithmetics::EcPoint, g1::{G1Element, G1Field}, g2::{G2Element, G2Field}, gt::{GTElement, GTField}};
use engines::{  bls12_381_engine, bls12_446_engine, bls12_461_engine, bls24_315_engine, bls24_477_engine, bls24_479_engine, bls24_509_engine, 
                bls24_509_snark_engine, bls24_559_engine, bls48_277_engine, bls48_287_engine, bls48_571_engine, bls48_573_engine, bls48_575_engine, bls48_581_engine};
use extensions::g2_extfields::ExtFieldG2Element;
use fields::prime_fields::{FieldElement, PrimeField};
use tools::{arithmetic_interface::ArithmeticOperations, exponent::Exponent};

#[derive(Debug)]
pub enum CurvesNames {
    Bls12_381, 
    Bls12_446, 
    Bls12_461, 
    Bls24_315, 
    Bls24_477, 
    Bls24_479, 
    Bls24_509, 
    Bls24_509Snark, 
    Bls24_559, 
    Bls48_277, 
    Bls48_287, 
    Bls48_571, 
    Bls48_573, 
    Bls48_575, 
    Bls48_581
}
impl CurvesNames {
    pub fn all() -> [CurvesNames; 15] {
        [
            CurvesNames::Bls12_381,
            CurvesNames::Bls12_446,
            CurvesNames::Bls12_461,
            CurvesNames::Bls24_315,
            CurvesNames::Bls24_477,
            CurvesNames::Bls24_479,
            CurvesNames::Bls24_509,
            CurvesNames::Bls24_509Snark,
            CurvesNames::Bls24_559,
            CurvesNames::Bls48_277,
            CurvesNames::Bls48_287,
            CurvesNames::Bls48_571,
            CurvesNames::Bls48_573,
            CurvesNames::Bls48_575,
            CurvesNames::Bls48_581
        ]
    }
    pub fn bls12() -> [CurvesNames; 3] {
        [
            CurvesNames::Bls12_381,
            CurvesNames::Bls12_446,
            CurvesNames::Bls12_461,            
        ]
    }
    pub fn bls24() -> [CurvesNames; 6] {
            [   CurvesNames::Bls24_315,
                CurvesNames::Bls24_477,
                CurvesNames::Bls24_479,
                CurvesNames::Bls24_509,
                CurvesNames::Bls24_509Snark,
                CurvesNames::Bls24_559,
            ]
        }
    pub fn bls48() -> [CurvesNames; 6] {
            [                
                CurvesNames::Bls48_277,
                CurvesNames::Bls48_287,
                CurvesNames::Bls48_571,
                CurvesNames::Bls48_573,
                CurvesNames::Bls48_575,
                CurvesNames::Bls48_581
            ]
        }
}


#[derive(Debug)]
pub struct Pairings <const R:usize,const N:usize,const MAX_COEFS_COUNT1:usize,const PRAMASIZE:usize, const MAX_COEFS_COUNT2 :usize>
         {  pub identifier :&'static str,
            pub curvename :CurvesNames,
            pub g1 : &'static G1Field<R,N,MAX_COEFS_COUNT1>, 
            pub g2:  &'static G2Field<PRAMASIZE,R,N,MAX_COEFS_COUNT2>,
            pub gt : &'static GTField<N,PRAMASIZE>,
            pub fr : &'static PrimeField<R>}



pub trait  PairingsEngine {
    type GT;
    type G1;
    type G2;
    fn miller_loop(&self, p:&Self::G1, q:&Self::G2) -> Self::GT;
    fn paire(&self, p:&Self::G1, q:&Self::G2) -> Self::GT;
    fn multi_paire(&self, p:&[Self::G1], q:&[Self::G2]) -> Self::GT; 
 }  
 
 
pub struct BLS12{}
pub struct BLS24{}
pub struct BLS48{}
pub trait Bls12Curves{ 
    fn _381() ->&'static Pairings<4, 6, 16, 4,  4>;    
    fn _446() ->&'static Pairings<5, 7, 3, 4,  7>;    
    fn _461() ->&'static Pairings<5, 8, 10, 4,  4>;
}

impl Bls12Curves for BLS12{
    fn _381() ->&'static Pairings<4, 6, 16, 4,  4>
    {
        bls12_381_engine()
    }

    fn _446() ->&'static Pairings<5, 7, 3, 4,  7>
    {
        bls12_446_engine()
    }

    fn _461() ->&'static Pairings<5, 8, 10, 4,  4>
    {
        bls12_461_engine()
    }
}
    
pub trait BLS24Curves{
    fn _315() ->&'static Pairings<4, 5, 3, 11,  7>;
    fn _477() ->&'static Pairings<6, 8, 7, 11,  4>;
    fn _479() ->&'static Pairings<7, 8, 4, 10,  10>;
    fn _509_snark() ->&'static Pairings<7, 8, 3, 11,  7>;
    fn _509() ->&'static Pairings<7, 8, 3, 10,  10>;
    fn _559() ->&'static Pairings<8, 9, 4, 10,  4>;
}

impl BLS24Curves for BLS24 {
    fn _315() ->&'static Pairings<4, 5, 3, 11,  7>
    {
        bls24_315_engine()
    }
    
    fn _477() ->&'static Pairings<6, 8, 7, 11,  4>
    {
        bls24_477_engine()
    }

    fn _479() ->&'static Pairings<7, 8, 4, 10,  10>
    {
        bls24_479_engine()
    }

    fn _509_snark() ->&'static Pairings<7, 8, 3, 11,  7>
    {
        bls24_509_snark_engine()
    }

    fn _509() ->&'static Pairings<7, 8, 3, 10,  10>
    {
        bls24_509_engine()
    }

    fn _559() ->&'static Pairings<8, 9, 4, 10,  4>
    {
        bls24_559_engine()
    }
}

pub trait  BLS48Curves{
    fn _575() ->&'static Pairings<9, 9, 4, 24,  7>;
    fn _573() ->&'static Pairings<8, 9, 7, 24,  4>;
    fn _571() ->&'static Pairings<8, 9, 3, 36,  4>;
    fn _581() ->&'static Pairings<9, 10, 3, 24,  4>;
    fn _287() ->&'static Pairings<4, 5, 4, 35,  19>;
    fn _277() ->&'static Pairings<4, 5, 10, 24,  4>;
}

impl  BLS48Curves for BLS48 {
    fn _575() ->&'static Pairings<9, 9, 4, 24,  7>
    {
        bls48_575_engine()
    }   

    fn _581() ->&'static Pairings<9, 10, 3, 24,  4>
    {
        bls48_581_engine()
    }

    fn _573() ->&'static Pairings<8, 9, 7, 24,  4>
    {
        bls48_573_engine()
    }  

    fn _571() ->&'static Pairings<8, 9, 3, 36,  4>
    {
        bls48_571_engine()
    }   

    fn _287() ->&'static Pairings<4, 5, 4, 35,  19>
    {
        bls48_287_engine()
    }   

    fn _277() ->&'static Pairings<4, 5, 10, 24,  4>
    {
        bls48_277_engine()
    }   
}



fn double_jacobian_for_miller<const N: usize, const PARAMSIZE :usize>
                  (q : &mut EcPoint<ExtFieldG2Element<N,PARAMSIZE>>,
                   px :&FieldElement<N>,py : &FieldElement<N> )-> [ExtFieldG2Element<N,PARAMSIZE>;3] 
{        
         //    Personal implementation of Algorithm 26 from https://eprint.iacr.org/2010/354.pdf
         let mut t0 = q.x.sqr();         
         let mut t1 = q.y.sqr();
         let mut t2 = t1.sqr();
         let mut t3 = t1.addto(&q.x).sqr().substract(&t0).substract(&t2).double();
         let t4     = t0.double().addto(&t0);
         let mut t6 = t4.addto(&q.x);
         let t5     = t4.sqr();
         let z2     = q.z.sqr();
         q.x =  t5.substract(&t3.double());
         q.z =  q.z.addto(&q.y).sqr().substract(&t1).substract(&z2);
         q.y =  t3.substract(&q.x).multiply(&t4); 
         t2  = t2.double().double().double();
         q.y = q.y.substract(&t2);
         t3  = t4.multiply(&z2).double().negate();
         t6  = t6.sqr().substract(&t0).substract(&t5);
         t1  = t1.double().double();
         t6  = t6.substract(&t1);
         t0  = q.z.multiply(&z2).double();
         [t0.mulby_fp_element(py),t3.mulby_fp_element(px),t6]
 } 

fn add_jacobian_for_miller<const N: usize, const PARAMSIZE :usize>
                    ( q :&mut EcPoint<ExtFieldG2Element<N,PARAMSIZE>>,
                      xq :&ExtFieldG2Element<N,PARAMSIZE>,
                      yq :&ExtFieldG2Element<N,PARAMSIZE>,
                      px :&FieldElement<N>, py: &FieldElement<N> )-> [ExtFieldG2Element<N,PARAMSIZE>;3]                     
{        
         //    Personal implementation of Algorithm 27 from https://eprint.iacr.org/2010/354.pdf
         let mut z2 = q.z.sqr();
         let y2     = yq.sqr();
         let mut t0 = z2.multiply(&xq);
         let mut t1 = yq.addto(&q.z).sqr().substract(&y2).substract(&z2).multiply(&z2);
         let t2 = t0.substract(&q.x);
         let t3 = t2.sqr();
         let t4 = t3.double().double();
         let t5 = t4.multiply(&t2);
         let t6 = t1.substract(&q.y.double());
         let mut t9 = t6.multiply(&xq);
         let t7     = t4.multiply(&q.x);
         q.x = t6.sqr().substract(&t5).substract(&t7.double());
         q.z = q.z.addto(&t2).sqr().substract(&z2).substract(&t3);
         let mut t10 = yq.addto(&q.z);
         let t8      = t7.substract(&q.x).multiply(&t6);
         t0  = q.y.multiply(&t5).double();
         q.y = t8.substract(&t0);
         t10 = t10.sqr().substract(&y2);
         z2  = q.z.sqr();
         t10 = t10.substract(&z2);
         t9  = t9.double().substract(&t10);
         t10 = q.z.double();
         t1  = t6.double().negate();
         [t10.mulby_fp_element(py),t1.mulby_fp_element(px),t9]
}

impl <const R:usize,const N:usize,const MAX_COEFS_COUNT1:usize,const PRAMASIZE:usize, const MAX_COEFS_COUNT2 :usize> PairingsEngine 
      for Pairings<R,N,MAX_COEFS_COUNT1,PRAMASIZE,MAX_COEFS_COUNT2> {
        type GT = GTElement<N,PRAMASIZE>;
        type G1 = G1Element<R,N,MAX_COEFS_COUNT1>;
        type G2 = G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT2>;     


      fn miller_loop(&self, p:&Self::G1, q:&Self::G2) -> Self::GT 
      {         
         let mut f = self.gt.one();
         let _p = p.to_affine().point;
         let mut _q: EcPoint<ExtFieldG2Element<N, PRAMASIZE>>  = q.to_affine().point;
         let qx = _q.x;
         let qy = _q.y;
         let _loop = &<u128 as Exponent<N>>::to_naf(&(self.g2.consts.u.abs() as u128))[1..];
         let idx = if self.g2.consts.twist_type == 'D' {(0,2,1)} else {(2,1,0)};
         for i in _loop{ f = f.sqr();                                
                              let l = double_jacobian_for_miller(&mut _q,&_p.x,&_p.y);                                                                                                                                          
                              f = f.sparse_multiply(&[ l[idx.0].content(),
                                                            l[idx.1].content(),
                                                            l[idx.2].content()],self.g2.consts.twist_type);                                                    
                              if *i != 0 { let l = 
                                                if *i==1 {add_jacobian_for_miller(&mut _q, &qx, &qy, &_p.x, &_p.y)}
                                                else {add_jacobian_for_miller(&mut _q, &qx, &qy.negate(), &_p.x,&_p.y)} ;                                                                                                                                                                                                         
                                                f = f.sparse_multiply(&[  l[idx.0].content(),
                                                                              l[idx.1].content(),
                                                                              l[idx.2].content()],self.g2.consts.twist_type);
                                    }
                                    
                            }
         f        
         }      

      fn paire(&self, p:&Self::G1, q:&Self::G2) -> Self::GT 
      {
         self.miller_loop(p, q).final_exponentiation()
      }            

      fn multi_paire(&self, p_list:&[Self::G1], q_list:&[Self::G2]) -> Self::GT 
      {
         if p_list.len()!=q_list.len() {panic!("Incompatible number of elements from G1 and G2 ... ")}
         else {      let mut f = self.gt.one();
                     let _loop = &<u128 as Exponent<N>>::to_naf(&(self.g2.consts.u.abs() as u128))[1..];
                     let mut _plist = Vec::<EcPoint<FieldElement<N>>>::new();                     
                     let mut _qlist  = Vec::<EcPoint<ExtFieldG2Element<N, PRAMASIZE>>>::new();
                     for q in q_list {_qlist.push(q.to_affine().point)};
                     for p in p_list {_plist.push(p.to_affine().point)};
                     let qfix =_qlist.clone();                    
                     let idx = if self.g2.consts.twist_type == 'D' {(0,2,1)} else {(2,1,0)};
                     for ib in _loop
                           {  f = f.sqr();                                
                              for i in 0..q_list.len() 
                                    {  let l = double_jacobian_for_miller(&mut _qlist[i],&_plist[i].x,&_plist[i].y);
                                       f = f.sparse_multiply(&[ l[idx.0].content(),l[idx.1].content(),l[idx.2].content()],self.g2.consts.twist_type);
                                    }   
                              if *ib != 0 { for i in 0..q_list.len() 
                                             {  let _qy =if *ib==1 {qfix[i].y} else {qfix[i].y.negate()};
                                                let l = 
                                                add_jacobian_for_miller(&mut _qlist[i], &qfix[i].x,&_qy, &_plist[i].x, &_plist[i].y);
                                                f = f.sparse_multiply(&[ l[idx.0].content(),l[idx.1].content(),l[idx.2].content()],self.g2.consts.twist_type);
                                             } 
                                          }
                                       }
                     f.final_exponentiation()        
              }

      }
 }