// Code developed by FARAOUN Kamel Mohamed.
// EEDIS-Laboratory. UDL-University. Algeria
// During May 2024.

use base64::engine::general_purpose;
use num_bigint::{BigUint, ToBigInt};
use num_traits::{Zero,One,Num,ToPrimitive};
use std::fmt;
use std::ops::{Add, BitAnd, BitOr, Mul, Neg, Sub};
use std::str::FromStr;
use crate::curves::curve_arithmetics::*;
use crate::fields::prime_fields::{FieldElement, PrimeField};
use crate::tools::hashs::{i2osp, i2osp_pf, os2ip};
use crate::tools::arithmetic_interface::ArithmeticOperations;
use crate::tools::recoders::{recod_scalar_glv2, WDMASK, WDSIZE, WSIZE};
use base64::{self, Engine};

#[derive(Debug)]
pub struct G1SwuIsogeniesConsts<const N:usize,const MAX_COEFS_COUNT:usize> {
        pub z:FieldElement<N>,
        pub swu_a:FieldElement<N>,
        pub swu_b:FieldElement<N>,
        pub xnum :[FieldElement<N>;MAX_COEFS_COUNT],
        pub xden :[FieldElement<N>;MAX_COEFS_COUNT],
        pub ynum :[FieldElement<N>;MAX_COEFS_COUNT],
        pub yden :[FieldElement<N>;MAX_COEFS_COUNT],
        pub inv_z:FieldElement<N>,
        pub b_div_a:FieldElement<N>,        
        }

#[derive(Debug)]
pub struct G1Consts<const R:usize,const N:usize,const MAX_COEFS_COUNT :usize>
    {   pub b:  FieldElement<N>,
        pub a : FieldElement<N>,
        pub h1: u128,
        pub w : FieldElement<N>,
        pub swu_consts :G1SwuIsogeniesConsts<N,MAX_COEFS_COUNT>,
        pub lambda: FieldElement<R>,
        pub lambda_big: BigUint,
        pub base_field_numbits:usize,
        pub security_level:usize,
        pub default_generator : EcPoint<FieldElement<N>>
    }

#[derive(Clone,Copy)]
pub struct G1Element<const R:usize ,const N:usize,const MAX_COEFS_COUNT :usize>
{   pub consts : &'static G1Consts<R,N,MAX_COEFS_COUNT>,
    pub point  : EcPoint<FieldElement<N>>,
}

#[derive(Debug)]
pub struct G1Field<const R:usize,const N:usize,const MAX_COEFS_COUNT:usize>
{   pub consts : &'static G1Consts<R,N,MAX_COEFS_COUNT>,
    pub base_field: &'static PrimeField<N>,
    pub fr_field :  &'static PrimeField<R>,
}

impl <const R:usize,const N:usize,const MAX_COEFS_COUNT:usize> G1Element<R,N,MAX_COEFS_COUNT>
        {   
            pub fn addto(&self, other: &G1Element<R,N,MAX_COEFS_COUNT>) -> G1Element<R,N,MAX_COEFS_COUNT>
            {
                G1Element { point :self.point.add_jacobian(&other.point),
                            consts : self.consts,
                            }
            }
            pub fn substract(&self, other: &G1Element<R,N,MAX_COEFS_COUNT>) -> G1Element<R,N,MAX_COEFS_COUNT>
            {
                G1Element { point :self.point.add_jacobian(&other.point.negate()),
                            consts : self.consts,
                            }
            }
            pub fn negate(&self) -> G1Element<R,N,MAX_COEFS_COUNT>
            {
                G1Element { point :self.point.negate(),
                            consts : self.consts,
                            }
            }
            pub fn multiply_by_const(&self , scalar :i128) -> Self
            {
                G1Element { point :self.point.multiply_with_const(scalar), consts :self.consts }
            }
            pub fn multiply(&self , scalar :&FieldElement<R>) -> G1Element<R,N,MAX_COEFS_COUNT>{
                G1Element { point :self.point.multiply(&scalar),
                            consts :self.consts
                          }
            }
            pub fn equal(&self, other : &G1Element<R,N,MAX_COEFS_COUNT>) -> bool 
            {
                self.point.equal(&other.point)
            }
            pub fn  glv_multiply(&self, scalar : &FieldElement<R>) -> G1Element<R,N,MAX_COEFS_COUNT>
            {   let s = scalar.to_big_uint();
                if s.bits()<(scalar.fieldparams.num_of_bits /2).try_into().unwrap() {self.multiply(&scalar)}
                else {  let _2p     = self.point.double_jacobian();
                        let ff: BigUint = BigUint::from_str("255").unwrap();
                        let mut _2phi   = EcPoint{x: _2p.x.multiply(&self.consts.w), y: _2p.y, z: _2p.z.clone() };
                        let infinit = EcPoint {x:self.point.x.one(), y : self.point.x.one(), z: self.point.x.zero() };                
                        const TABLESIZE :usize = 1 << (WDSIZE - 1);
                        const BLOCKSIZE :usize = 1 << (WSIZE - 1);
                        let mut lookup = [infinit;TABLESIZE >> 1];    
                        lookup[0] = EcPoint {   x: self.point.x.multiply(&self.consts.w.addto(&self.point.x.one())).negate(),
                                                y: self.point.y.negate(),
                                                z: self.point.z.clone()};
                        lookup[1] = _2p.add_jacobian(&lookup[0]);
                        lookup[2] = _2p.add_jacobian(&lookup[1]);
                        lookup[3] = _2p.add_jacobian(&lookup[2]);
                        for i in BLOCKSIZE..(TABLESIZE >> 1)  { lookup[i] = _2phi.add_jacobian(&lookup[i - BLOCKSIZE]) };
                        let alpha = (&s).bitor(&self.consts.lambda_big).bitand(BigUint::one());
                        let r = BigUint::from_str_radix(&scalar.fieldparams.modulo_as_strhex[2..], 16).unwrap(); 
                        let a = (&r - &s) * (BigUint::one() - &alpha) + &s * &alpha;
                        let mut code = recod_scalar_glv2(&a, &self.consts.lambda_big);                
                        let mut limb = (&code).bitand((ff).clone()).to_u8().unwrap();
                        let mut fi   = limb & 1;
                        let mut sig : i8  = (limb & 2) as i8 - 1;
                        let mut idx = ((limb & WDMASK) >> 2) as usize;
                        let mut result =  if fi == 0 { G1Element {point :lookup[idx], consts : self.consts}} 
                                                                            else { G1Element {point :lookup[idx], consts : self.consts}.phi().negate()};
                        if sig == -1 {result.point =  result.point.negate()};                                                
                        code = code >> WDSIZE;
                        while code != BigUint::one() {  limb = (&code).bitand((ff).clone()).to_u8().unwrap();
                                                        fi   = limb & 1; 
                                                        sig  = (limb & 2) as i8 - 1;
                                                        idx  = ((limb & WDMASK) >> 2) as usize;
                                                        result.point = result.point.double_jacobian(); 
                                                        result.point = result.point.double_jacobian(); 
                                                        result.point = result.point.double_jacobian(); 
                                                        let tmp =  if fi == 0 { G1Element {point :lookup[idx], consts : self.consts}} 
                                                                                                    else { G1Element {point :lookup[idx], consts : self.consts}.phi().negate()};
                                                        result.point = if sig ==1 {result.point.add_jacobian(&tmp.point)}
                                                                    else {result.point.add_jacobian(&tmp.point.negate())};
                                                        code = code >> WDSIZE;
                                                    }           
                        if alpha.is_zero() { G1Element {point : result.point.negate(), consts:self.consts}}
                        else {result}}
            }

            pub fn phi(&self) -> G1Element<R,N,MAX_COEFS_COUNT>
            {
                G1Element{ point : EcPoint{x: self.point.x.multiply(&self.consts.w), y: self.point.y, z: self.point.z.clone() }, 
                           consts : self.consts}    
            }

            pub fn is_on_curve(&self) -> bool
            {   let z6 = self.point.z.sqr().multiply(&self.point.z).sqr();
                self.point.y.sqr().equal(&self.point.x.sqr().multiply(&self.point.x).addto(&self.consts.b.multiply(&z6)))
            }
            pub fn is_torsion(&self) -> bool
            {   
                self.multiply(&self.consts.lambda).equal(&self.phi())
            }
            pub fn to_string(&self ) -> String
            {
                self.point.to_string()
            } 
            pub fn to_hex_string(&self ) -> String
            {
                self.point.to_hex_string()
            } 
            pub fn encode_to_base64(&self) ->String
            {   
                general_purpose::STANDARD.encode(self.to_compressed_bytearray())
            }
            pub fn to_compressed_bytearray(&self) -> Vec<u8>
            {   
                //  Point compression/Serialization as described by ZCach serialization format
                //  https://www.ietf.org/archive/id/draft-irtf-cfrg-pairing-friendly-curves-11.html#name-zcash-serialization-format-
                let mut p = self.point.clone();
                p.to_affine();
                let c_bit: u8 = 1;
                let i_bit: u8 = if self.point.z.is_zero() {1} else {0};
                let s_bit: i8 = if self.point.z.is_zero() {0} else {if self.point.y.sign()==1 {1} else {0}};                
                let m_byte: u8 = (c_bit << 7) | (i_bit << 6) | (((s_bit + 1) as u8 >> 1) << 5);
                let numbits = self.point.x.fieldparams.num_of_bits;                
                let sizeinbytes = (numbits >> 3) + if (numbits % 8) ==0 {0} else {1};
                let mut x_string = if self.point.z.is_zero() {i2osp(0, sizeinbytes)}
                                            else {i2osp_pf(&p.x, sizeinbytes)};
                if self.consts.base_field_numbits % 8 <=5 {x_string[0] = x_string[0] | m_byte;}
                else {x_string.insert(0, m_byte);}
                x_string
            }
            pub fn to_affine(&self) -> Self
            {   
                if ! self.point.z.is_one() {    let mut tmp = self.point.clone();
                                                tmp.to_affine();
                                                G1Element { point :tmp, consts :self.consts }  
                                            }    
                else {self.clone()}
            }
        }


impl <const R:usize,const N:usize,const MAX_COEFS_COUNT : usize> G1Field<R,N,MAX_COEFS_COUNT> 
    {
        fn random_point_using_swu(&self,seed :FieldElement<N>) -> G1Element<R,N,MAX_COEFS_COUNT>{ 
            //     Simplified Shallue-van de Woestijne-Ulas Method (Simplified SWU for AB == 0)
            //     https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-hash-to-curve-05#section-6.6.3
            //     https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-hash-to-curve-05#appendix-C.2            
            let mut u = self.base_field.random_element();
            if !seed.is_zero() {u=seed};            
            let t1 = self.consts.swu_consts.z.multiply(&u.sqr());
            let mut t2 = t1.sqr();
            let mut x1 = t1.addto(&t2).invert();
            if x1.is_zero() {x1 = self.consts.swu_consts.inv_z.clone() }
            else {x1 = x1.addto(&self.base_field.one())}
            x1 = x1.multiply(&self.consts.swu_consts.b_div_a);            
            let gx1 = x1.sqr().addto(&self.consts.swu_consts.swu_a).multiply(&x1).addto(&self.consts.swu_consts.swu_b);
            let ty = gx1.sqrt();
            let mut x: FieldElement<N>;
            let mut y: FieldElement<N>;
            if !ty.is_none() {  y = ty.unwrap();
                                if u.sign()!= ty.unwrap().sign() { y = y.negate()}
                                x = x1.clone();  
                            }
            else { let x2 = t1.multiply(&x1);
                   t2 = t2.multiply(&t1);                   
                   let gx2 = gx1.multiply(&t2);                 
                   y  = gx2.sqrt().unwrap();
                   if u.sign() != y.sign() {  y = y.negate();} 
                   x = x2.clone(); 
                }
            let mut pow = x.clone();
            let mut xnum = self.consts.swu_consts.xnum[0];
            let mut xden = self.consts.swu_consts.xden[0];
            let mut ynum = self.consts.swu_consts.ynum[0];
            let mut yden = self.consts.swu_consts.yden[0];
            for i in 1..MAX_COEFS_COUNT {
                ynum = ynum.addto(&self.consts.swu_consts.ynum[i].multiply(&pow));
                yden = yden.addto(&self.consts.swu_consts.yden[i].multiply(&pow));
                if !self.consts.swu_consts.xnum[i].is_zero() {  xnum = xnum.addto(&self.consts.swu_consts.xnum[i].multiply(&pow))};
                if !self.consts.swu_consts.xden[i].is_zero() {  xden = xden.addto(&self.consts.swu_consts.xden[i].multiply(&pow))};
                pow = pow.multiply(&x);                                                                                    
            }             
            // Convert from projective coordinates to Jacobian coordinates
            let z = xden.multiply(&yden);
            x = xnum.multiply(&z).multiply(&yden);
            y = ynum.multiply(&z.sqr()).multiply(&xden).multiply(&y);
            G1Element { point : EcPoint { x ,y ,z },
                        consts :self.consts,
                      }

        }
        pub fn map_to_curve(&self,seed1 : FieldElement<N>) ->G1Element<R,N,MAX_COEFS_COUNT> 
        {  
             self.random_point_using_swu(seed1)
             
        }
        pub fn random_point_withseed(&self,seed1 : FieldElement<N>) -> G1Element<R,N,MAX_COEFS_COUNT>
            {   
                let mut rp = self.map_to_curve(seed1);
                rp = rp.multiply_by_const(self.consts.h1 as i128);
                rp    
            }
        pub fn random_point(&self) -> G1Element<R,N,MAX_COEFS_COUNT>
            {   
                let mut rp = self.map_to_curve(self.base_field.zero());
                rp = rp.multiply_by_const(self.consts.h1 as i128);                
                rp.to_affine()    
            }
        pub fn hash_to_field(&self,id :&str,mode :u8) -> G1Element<R,N,MAX_COEFS_COUNT>
            {   
                //  mode=0: Encode to Curve (NU-encode-to-curve), mode=1: Random Oracle model (RO-hash-to-curve)
                //  https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-hash-to-curve-06#name-roadmap
                let hashs = self.base_field.hash_to_field(id,self.consts.security_level, 2);
                if mode ==0 {self.random_point_withseed(hashs[0]).to_affine()}   
                else {  let p1=self.random_point_withseed(hashs[0]);
                        let p2=self.random_point_withseed(hashs[1]);
                        p1.addto(&p2).to_affine()                        
                     } 
            }
        
        pub fn from_bytearray(&self,inbytes : &Vec<u8>) -> G1Element<R,N,MAX_COEFS_COUNT>
            {
                //  Point de-compression/de-Serialization as described by ZCach serialization format
                //  https://www.ietf.org/archive/id/draft-irtf-cfrg-pairing-friendly-curves-11.html#name-zcash-serialization-format-
                let mut input = inbytes.clone();
                let m_byte = input[0] & 0xE0;
                let numbits = self.base_field.parametres.num_of_bits;
                let sizeinbytes = (numbits >> 3) + if (numbits % 8) ==0 {0} else {1};
                if self.consts.base_field_numbits % 8 <=5 {input[0] = input[0] & 0x1F;}
                else {input.remove(0);};
                if m_byte == 0xE0 {panic!("Invalide compressed point format ...")};
                if m_byte & 0x80 !=0 { if input.len() != sizeinbytes {panic!("Invalide compressed point format ...")} 
                                     }
                else {if input.len() != (sizeinbytes*2) {panic!("Invalide compressed point format ...")}}
                if m_byte & 0x40 !=0 { if input.iter().any(|&e| e != 0) {panic!("Invalid compression of an infinity point...");} 
                                       else {G1Element {  point : EcPoint {x:self.base_field.one(), y : self.base_field.one(), z: self.base_field.zero() },
                                                          consts :self.consts}
                                            } 
                                     }
                else {  if input.len() == (sizeinbytes*2){ let x = self.base_field.from_bigint(&os2ip(&input[0..sizeinbytes]).to_bigint().unwrap());
                                                          let y = self.base_field.from_bigint(&os2ip(&input[sizeinbytes..]).to_bigint().unwrap());
                                                          G1Element {  point : EcPoint { x, y, z: self.base_field.one() },
                                                          consts :self.consts}  
                                                         }
                        else { let x = self.base_field.from_bigint(&os2ip(&input[0..sizeinbytes]).to_bigint().unwrap());
                               let y = x.sqr().multiply(&x).addto(&&self.consts.b).sqrt();
                               if y.is_none() {panic!("Invalide point: not in the curve ...")}
                               else { let y = y.unwrap();
                                      let r_sign = if m_byte & 0x20 !=0 {1} else {0}; 
                                      if (y.sign()+1) >> 1 == r_sign {G1Element {  point : EcPoint { x, y, z: self.base_field.one() },
                                                                                   consts :self.consts}  }
                                      else {G1Element {  point : EcPoint { x, y: y.negate(), z: self.base_field.one() },
                                                         consts :self.consts}  }
                                    }
                             }                                 
                     }
            }
            pub fn from_base64(&self,input :&str) -> G1Element<R,N,MAX_COEFS_COUNT>
            {
                let decoded_bytes = match general_purpose::STANDARD.decode(input) {
                    Ok(bytes) => bytes,
                    Err(_) => {
                        panic!("Failed to decode base64 string");
                    }
                };
                self.from_bytearray(&decoded_bytes)
            }

            pub fn default_generator(&self)->G1Element<R,N,MAX_COEFS_COUNT>
            {
                G1Element { point : EcPoint { x : self.consts.default_generator.x, 
                                               y : self.consts.default_generator.y, 
                                               z : self.consts.default_generator.z },
                            consts :self.consts} 
            }
    }
    

impl <const R:usize,const N:usize,const MAX_COEFS_COUNT:usize> fmt::Display for G1Element<R,N,MAX_COEFS_COUNT> {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            write!(f, "{:}", &self.to_string())
        }
    }

impl  <const R:usize,const N:usize,const MAX_COEFS_COUNT:usize> Add for G1Element<R,N,MAX_COEFS_COUNT> {
        type Output =  G1Element<R,N,MAX_COEFS_COUNT>;
            fn add(self, rhs: Self) -> Self::Output {   self.addto(&rhs) }
    }
impl  <const R:usize,const N:usize,const MAX_COEFS_COUNT:usize> Sub for G1Element<R,N,MAX_COEFS_COUNT> {
        type Output =  G1Element<R,N,MAX_COEFS_COUNT>;
            fn sub(self, rhs: Self) -> Self::Output {   self.substract(&rhs) }
    }
impl  <const R:usize,const N:usize,const MAX_COEFS_COUNT:usize> Neg for G1Element<R,N,MAX_COEFS_COUNT> {
        type Output =  G1Element<R,N,MAX_COEFS_COUNT>;
            fn neg(self) -> Self::Output { self.negate() }
        }   
impl  <const R:usize,const N:usize,const MAX_COEFS_COUNT:usize> PartialEq for G1Element<R,N,MAX_COEFS_COUNT> {
            fn eq(&self, other: &Self) -> bool {    self.equal(other) }
        }
impl <const R:usize,const N:usize,const MAX_COEFS_COUNT:usize> Mul<i128> for G1Element<R,N,MAX_COEFS_COUNT> {
        type Output = G1Element<R,N,MAX_COEFS_COUNT>;    
            fn mul(self, rhs: i128) -> Self::Output {     self.multiply_by_const(rhs)   }
        }
impl <const R:usize,const N:usize,const MAX_COEFS_COUNT:usize> Mul<G1Element<R,N,MAX_COEFS_COUNT>> for i128 {
        type Output = G1Element<R,N,MAX_COEFS_COUNT>;
        fn mul(self, rhs: G1Element<R,N,MAX_COEFS_COUNT>) -> Self::Output {  rhs.multiply_by_const(self)  }
        }
impl <const R:usize,const N:usize,const MAX_COEFS_COUNT:usize> Mul<u64> for G1Element<R,N,MAX_COEFS_COUNT> {
        type Output = G1Element<R,N,MAX_COEFS_COUNT>;    
            fn mul(self, rhs: u64) -> Self::Output {     self.multiply_by_const(rhs as i128)   }
        }
impl <const R:usize,const N:usize,const MAX_COEFS_COUNT:usize> Mul<G1Element<R,N,MAX_COEFS_COUNT>> for u64 {
        type Output = G1Element<R,N,MAX_COEFS_COUNT>;
            fn mul(self, rhs: G1Element<R,N,MAX_COEFS_COUNT>) -> Self::Output {  rhs.multiply_by_const(self as i128)  }
        }    
impl <const R:usize,const N:usize,const MAX_COEFS_COUNT:usize> Mul<i64> for G1Element<R,N,MAX_COEFS_COUNT> {
        type Output = G1Element<R,N,MAX_COEFS_COUNT>;    
            fn mul(self, rhs: i64) -> Self::Output {     self.multiply_by_const(rhs as i128)   }
        }            
impl <const R:usize,const N:usize,const MAX_COEFS_COUNT:usize> Mul<G1Element<R,N,MAX_COEFS_COUNT>> for i64 {
        type Output = G1Element<R,N,MAX_COEFS_COUNT>;
            fn mul(self, rhs: G1Element<R,N,MAX_COEFS_COUNT>) -> Self::Output {  rhs.multiply_by_const(self as i128)  }
        }    
impl <const R:usize,const N:usize,const MAX_COEFS_COUNT:usize> Mul<u8> for G1Element<R,N,MAX_COEFS_COUNT> {
        type Output = G1Element<R,N,MAX_COEFS_COUNT>;    
            fn mul(self, rhs: u8) -> Self::Output {     self.multiply_by_const(rhs as i128)   }
        }            
impl <const R:usize,const N:usize,const MAX_COEFS_COUNT:usize> Mul<G1Element<R,N,MAX_COEFS_COUNT>> for u8 {
        type Output = G1Element<R,N,MAX_COEFS_COUNT>;
            fn mul(self, rhs: G1Element<R,N,MAX_COEFS_COUNT>) -> Self::Output {  rhs.multiply_by_const(self as i128)  }
        }    
impl <const R:usize,const N:usize,const MAX_COEFS_COUNT:usize> Mul<i8> for G1Element<R,N,MAX_COEFS_COUNT> {
        type Output = G1Element<R,N,MAX_COEFS_COUNT>;    
            fn mul(self, rhs: i8) -> Self::Output {     self.multiply_by_const(rhs as i128)   }
        }            
impl <const R:usize,const N:usize,const MAX_COEFS_COUNT:usize> Mul<G1Element<R,N,MAX_COEFS_COUNT>> for i8 {
        type Output = G1Element<R,N,MAX_COEFS_COUNT>;
            fn mul(self, rhs: G1Element<R,N,MAX_COEFS_COUNT>) -> Self::Output {  rhs.multiply_by_const(self as i128)  }
        }    
impl <const R:usize,const N:usize,const MAX_COEFS_COUNT:usize> Mul<G1Element<R,N,MAX_COEFS_COUNT>> for FieldElement<R> {
            type Output = G1Element<R,N,MAX_COEFS_COUNT>;
                fn mul(self, rhs: G1Element<R,N,MAX_COEFS_COUNT>) -> Self::Output {  (&rhs).glv_multiply(&self)  }
            }    
impl <const R:usize,const N:usize,const MAX_COEFS_COUNT:usize> Mul<FieldElement<R>> for G1Element<R,N,MAX_COEFS_COUNT>  {
                type Output = G1Element<R,N,MAX_COEFS_COUNT>;
                    fn mul(self, rhs: FieldElement<R>) -> Self::Output {  self.glv_multiply(&rhs)  }
                }                