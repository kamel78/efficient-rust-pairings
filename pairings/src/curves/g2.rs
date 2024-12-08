// Code developed by FARAOUN Kamel Mohamed.
// EEDIS-Laboratory. UDL-University. Algeria
// During May 2024.

use std::{fmt, ops::{Add, Mul, Neg, Sub}};
use base64::{engine::general_purpose, Engine};
use num_bigint::BigUint;
use crate::{extensions::{ext_fields::ExFieldConsts, g2_extfields::{ExtFieldG2Element, ExtG2Field}}, 
            fields::prime_fields::{FieldElement, PrimeField}, 
            tools::{arithmetic_interface::ArithmeticOperations, hashs::i2osp}};
use super::{curve_arithmetics::EcPoint, g2_primitives::{cofactor_clean::{clean_cofactor_bls12, clean_cofactor_bls24, clean_cofactor_bls48}, 
            gls_multiplication::{gls16_multiply_bls48, gls4_multiply_bls12, gls8_multiply_bls24}, 
            phi::{ phi_bls12, phi_bls24, phi_bls48}}};


            
#[derive(Debug)]
pub struct G2SwuIsogeniesConsts<const PRAMASIZE:usize,const N:usize, const MAX_COEFS_COUNT:usize> 
        {   pub z: ExtFieldG2Element<N,PRAMASIZE>,
            pub swu_a:ExtFieldG2Element<N,PRAMASIZE>,
            pub swu_b:ExtFieldG2Element<N,PRAMASIZE>,
            pub xnum :[ExtFieldG2Element<N,PRAMASIZE>;MAX_COEFS_COUNT],
            pub xden :[ExtFieldG2Element<N,PRAMASIZE>;MAX_COEFS_COUNT],
            pub ynum :[ExtFieldG2Element<N,PRAMASIZE>;MAX_COEFS_COUNT],
            pub yden :[ExtFieldG2Element<N,PRAMASIZE>;MAX_COEFS_COUNT],
            pub inv_z:ExtFieldG2Element<N,PRAMASIZE>,
            pub j_inv_z:ExtFieldG2Element<N,PRAMASIZE>,
            pub b_div_a:ExtFieldG2Element<N,PRAMASIZE>,        
        }

#[derive(Debug)]
pub struct G2Consts<const PRAMASIZE:usize,const R:usize,const N:usize,const MAX_COEFS_COUNT :usize>
    {   pub b:  ExtFieldG2Element<N,PRAMASIZE>,
        pub a : ExtFieldG2Element<N,PRAMASIZE>,
        pub swu_consts :G2SwuIsogeniesConsts<PRAMASIZE,N,MAX_COEFS_COUNT>,
        pub extfieldparams:ExFieldConsts<PRAMASIZE,N>,
        pub u :i128,
        pub lambda: FieldElement<R>,
        pub lambda_big: BigUint,
        pub base_field_numbits:usize,
        pub security_level:usize,
        pub twist_type :char,
        pub order : usize,
        pub default_generator :EcPoint<ExtFieldG2Element<N,PRAMASIZE>>
    }

#[derive(Clone,Copy)]
pub struct G2Element<const PRAMASIZE:usize,const R:usize,const N:usize,const MAX_COEFS_COUNT :usize>    
    {   pub consts : &'static G2Consts<PRAMASIZE,R,N,MAX_COEFS_COUNT>,
        pub point  : EcPoint<ExtFieldG2Element<N,PRAMASIZE>>,
    }

#[derive(Debug)]
pub struct G2Field<const PRAMASIZE:usize, const R:usize,const N:usize,const MAX_COEFS_COUNT :usize>
{   pub consts :  &'static G2Consts<PRAMASIZE,R,N,MAX_COEFS_COUNT>,
    pub base_field: &'static ExtG2Field<N,PRAMASIZE>,
    pub fr_field :  &'static PrimeField<R>,
}


impl <const PRAMASIZE:usize,const R:usize,const N:usize,const MAX_COEFS_COUNT :usize> G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT>
        {   
            fn getorder(&self)->usize
            {   match self.point.x {    ExtFieldG2Element::Fp2_1(_) => 2,
                                        ExtFieldG2Element::Fp4_1(_) => 4,
                                        ExtFieldG2Element::Fp8_1(_) => 8,
                                        ExtFieldG2Element::Fp2_2(_) => 2,
                                        ExtFieldG2Element::Fp4_2(_) => 4,
                                        ExtFieldG2Element::Fp8_2(_) => 8,
                                        ExtFieldG2Element::Fp4_3(_) => 4,
                                        ExtFieldG2Element::Fp8_3(_) => 8,
                                    } 
            }
            pub fn addto(&self, other: &Self) -> Self
            {
                G2Element { point :self.point.add_jacobian(&other.point), consts : self.consts }
            }
            pub fn substract(&self, other: &Self) -> Self
            {
                G2Element { point :self.point.add_jacobian(&other.point.negate()), consts : self.consts }
            }
            pub fn double(&self) -> Self
            {
                G2Element { point :self.point.double_jacobian(), consts : self.consts }
            }
            pub fn negate(&self) -> Self
            {
                G2Element { point :self.point.negate(), consts : self.consts }
            }
            pub fn multiply(&self , scalar :&FieldElement<R>) -> Self
            {
                G2Element { point :self.point.multiply(&scalar),consts :self.consts }
            }
            pub fn multiply_gls(&self , scalar :&FieldElement<R>) -> Self
            {
                let _scalar: BigUint = scalar.to_big_uint();
                if _scalar.bits()>((scalar.fieldparams.num_of_bits - (scalar.fieldparams.num_of_bits / self.consts.order))).try_into().unwrap()
                {match self.point.x {   ExtFieldG2Element::Fp2_1(_) => gls4_multiply_bls12(&self, scalar),
                                        ExtFieldG2Element::Fp4_1(_) => gls8_multiply_bls24(&self, scalar),
                                        ExtFieldG2Element::Fp8_1(_) => gls16_multiply_bls48(&self, scalar),
                                        ExtFieldG2Element::Fp2_2(_) => gls4_multiply_bls12(&self, scalar),
                                        ExtFieldG2Element::Fp4_2(_) => gls8_multiply_bls24(&self, scalar),
                                        ExtFieldG2Element::Fp8_2(_) => gls16_multiply_bls48(&self, scalar),
                                        ExtFieldG2Element::Fp4_3(_) => gls8_multiply_bls24(&self, scalar),
                                        ExtFieldG2Element::Fp8_3(_) => gls16_multiply_bls48(&self, scalar),
                                   }  
                    }
                else {self.multiply( scalar)}  
            }            
            pub fn multiply_by_const(&self , scalar :i128) -> Self
            {
                G2Element { point :self.point.multiply_with_const(scalar), consts :self.consts }
            }
            pub fn equal(&self, other : &Self) -> bool 
            {
                self.point.equal(&other.point)
            }
            fn clean_cofactor(&self) -> Self
            {
                match self.point.x {    ExtFieldG2Element::Fp2_1(_) => clean_cofactor_bls12(&self),
                                        ExtFieldG2Element::Fp4_1(_) => clean_cofactor_bls24(&self),
                                        ExtFieldG2Element::Fp8_1(_) => clean_cofactor_bls48(&self),
                                        ExtFieldG2Element::Fp2_2(_) => clean_cofactor_bls12(&self),
                                        ExtFieldG2Element::Fp4_2(_) => clean_cofactor_bls24(&self),
                                        ExtFieldG2Element::Fp8_2(_) => clean_cofactor_bls48(&self),
                                        ExtFieldG2Element::Fp4_3(_) => clean_cofactor_bls24(&self),
                                        ExtFieldG2Element::Fp8_3(_) => clean_cofactor_bls48(&self),
                                   }               
            }
            pub fn phi(&self) -> Self
            {   
                match self.point.x {    ExtFieldG2Element::Fp2_1(_) => phi_bls12(&self),
                                        ExtFieldG2Element::Fp4_1(_) => phi_bls24(&self,1),
                                        ExtFieldG2Element::Fp8_1(_) => phi_bls48(&self,1),
                                        ExtFieldG2Element::Fp2_2(_) => phi_bls12(&self),
                                        ExtFieldG2Element::Fp4_2(_) => phi_bls24(&self,1),
                                        ExtFieldG2Element::Fp8_2(_) => phi_bls48(&self,1),
                                        ExtFieldG2Element::Fp4_3(_) => phi_bls24(&self,1),
                                        ExtFieldG2Element::Fp8_3(_) => phi_bls48(&self,1)
                                   }       
            }
            pub fn is_torsion(&self) -> bool
            {   
                // Check if a point P is a Torsion Point
                // According to https://eprint.iacr.org/2019/814.pdf we have to check that  u*ψ3(P)-ψ2(P)+P="Infinity"
                // Or with according to https://hackmd.io/@yelhousni/bls12_subgroup_check that :−z⋅ψ(ϕ(Q))+ϕ(Q)+Q="Infinity" when  ϕ(P) is the endomorphisme of G2 (w*x,y)
                // However, all theses are equivalents to check that u*P=ψ(P) (we can show that ψ4(P)-ψ(P)+P="Infinity" for every point on the curve)
                // Same hint pointed by Scott :https://eprint.iacr.org/2021/1130.pdf
                self.multiply_by_const(self.consts.u).equal(&self.phi())
            }
            pub fn is_on_curve(&self) -> bool
            {   
                let z6 = self.point.z.sqr().multiply(&self.point.z).sqr();
                self.point.y.sqr().equal(&self.point.x.sqr().multiply(&self.point.x).addto(&self.consts.b.multiply(&z6)))
            }
            pub fn to_string(&self) -> String
            {
                self.point.to_string()
            }             
            pub fn to_hex_string(&self) -> String
            {
                self.point.to_hex_string()
            } 
            pub fn to_affine(&self) -> Self
            {   
                if ! self.point.z.is_one() {    let mut tmp = self.point.clone();
                                                tmp.to_affine();
                                                G2Element { point :tmp, consts :self.consts }  
                                            }     
                else {  self.clone()}
            }             
            pub fn encode_to_base64(&self) ->String
            {   
                general_purpose::STANDARD.encode(self.encode_to_compressed_bytearray())
            }
            
            pub fn encode_to_compressed_bytearray(&self) -> Vec<u8>
            {   
                //  Point compression/Serialization as described by ZCach serialization format
                //  https://www.ietf.org/archive/id/draft-irtf-cfrg-pairing-friendly-curves-11.html#name-zcash-serialization-format-
                let p = self.to_affine();
                let c_bit: u8 = 1;
                let i_bit: u8 = if p.point.z.is_zero() {1} else {0};
                let s_bit: i8 = if p.point.z.is_zero() {0} else {if p.point.y.sign() == 1 {1} else {0}};   
                let m_byte: u8 = (c_bit << 7) | (i_bit << 6) | (((s_bit + 1) as u8 >> 1) << 5);
                let numbits = p.consts.base_field_numbits;                
                let sizeinbytes = (numbits >> 3) + if (numbits % 8) ==0 {0} else {1};
                let mut x_string = if self.point.z.is_zero() {i2osp(0, sizeinbytes * self.getorder())}
                                            else {p.point.x.to_i2osp_bytearray()};
                if self.consts.base_field_numbits % 8 <=5 {x_string[0] = x_string[0] | m_byte;}
                else {x_string.insert(0, m_byte);}
                x_string
            }
        }

impl <const PRAMASIZE:usize,const R:usize,const N:usize,const MAX_COEFS_COUNT : usize> G2Field<PRAMASIZE,R,N,MAX_COEFS_COUNT> 
        {   
            pub fn getorder(&self)->usize
            {   match self.base_field {    ExtG2Field::Fp2_1(_) => 2,
                                           ExtG2Field::Fp4_1(_) => 4,
                                           ExtG2Field::Fp8_1(_) => 8,
                                           ExtG2Field::Fp2_2(_) => 2,
                                           ExtG2Field::Fp4_2(_) => 4,
                                           ExtG2Field::Fp8_2(_) => 8,
                                           ExtG2Field::Fp4_3(_) => 4,
                                           ExtG2Field::Fp8_3(_) => 8                                           
                                      } 
            }
           
            // pub fn map_to_curve(&self,seed :&ExtFieldG2Element<N, PRAMASIZE>) -> G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT>
            // { 
            //     //     Simplified Shallue-van de Woestijne-Ulas Method (Simplified SWU for AB == 0)
            //     //     https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-hash-to-curve-05#section-6.6.3
            //     //     https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-hash-to-curve-05#appendix-C.2
            //     let mut u: &ExtFieldG2Element<N, PRAMASIZE> = &self.base_field.random_element();
            //     if !seed.is_zero() {u = seed};
            //     let t1 = self.consts.swu_consts.z.multiply(&u.sqr());
            //     let mut t2 = t1.sqr();
            //     let mut x1 = t1.addto(&t2).invert();
            //     if x1.is_zero() {x1 = self.consts.swu_consts.inv_z.clone() }
            //     else {x1 = x1.addto(&self.base_field.one())}
            //     x1 = x1.multiply(&self.consts.swu_consts.b_div_a);            
            //     let gx1 = x1.sqr().addto(&self.consts.swu_consts.swu_a).multiply(&x1).addto(&self.consts.swu_consts.swu_b);                               
            //     let x2 = t1.multiply(&x1);
            //     t2 = t2.multiply(&t1);
            //     let gx2 = gx1.multiply(&t2);       
            //     let isqr = gx1.is_qr();
            //     let mut x: ExtFieldG2Element<N, PRAMASIZE>;
            //     let mut y: ExtFieldG2Element<N, PRAMASIZE>;
            //     if isqr { 
            //               y = gx1.sqrt().unwrap();
            //               x = x1.clone();  
            //             }
            //     else { y = gx2.sqrt().unwrap();
            //             x = x2.clone(); 
            //         }                
            //     if u.sign()!= y.sign() { y = y.negate()}   
            //     let mut xnum = x.zero();
            //     let mut xden = x.zero();
            //     let mut ynum = x.zero();
            //     let mut yden = x.zero();
            //     for i in (0..MAX_COEFS_COUNT).rev() {
            //                     ynum = ynum.addto(&self.consts.swu_consts.ynum[i]).multiply(&x);
            //                     yden = yden.addto(&self.consts.swu_consts.yden[i]).multiply(&x);
            //                     xnum = xnum.addto(&self.consts.swu_consts.xnum[i]).multiply(&x);
            //                     xden = xden.addto(&self.consts.swu_consts.xden[i]).multiply(&x);
            //                     }             
            //     // Convert from affine coordinates to Jacobian coordinates
            //     let z = xden.multiply(&yden);
            //     x = xnum.multiply(&z).multiply(&yden);              
            //     y = ynum.multiply(&z.sqr()).multiply(&xden).multiply(&y);
            //     G2Element { point : EcPoint { x ,y ,z },
            //                 consts :self.consts,
            //               }
            // }

            pub fn map_to_curve(&self,seed :&ExtFieldG2Element<N, PRAMASIZE>) -> G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT>
            { 
                 //     Simplified Shallue-van de Woestijne-Ulas Method (Simplified SWU for AB == 0)
                //     https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-hash-to-curve-05#section-6.6.3
                //     https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-hash-to-curve-05#appendix-C.2
                //      Prpposed variant that avoid inversion whenever is the targted extension field
                let mut u: &ExtFieldG2Element<N, PRAMASIZE> = &self.base_field.random_element();
                if !seed.is_zero() {u = seed};
                let t1 = self.consts.swu_consts.z.multiply(&u.sqr());
                let t2 = t1.sqr();
                let mut d = t1.addto(&t2);
                if d.is_zero() { d = self.consts.swu_consts.j_inv_z.clone()};
                let d2 = d.sqr();
                let d3 = d2.multiply(&d);
                let x0 = d2.addto(&d).multiply(&self.consts.swu_consts.b_div_a);
                let x1 = x0.multiply(&t1);
                let f0 = x0.sqr().multiply(&x0).substract(&d2.multiply(&d3).multiply(&self.consts.swu_consts.swu_b));    
                let f1 = f0.multiply(&t1).multiply(&t2);
                let mut x: ExtFieldG2Element<N, PRAMASIZE>;
                let mut y: ExtFieldG2Element<N, PRAMASIZE>;
                if f0.is_qr() { y = f0.sqrt().unwrap();
                                x = x0.clone();
                              }
                else {  y = f1.sqrt().unwrap();
                        x = x1.clone();
                     }                              
                if (u.sign()!= y.sign())&(d.sign()!=-1) { y = y.negate()}   
                let mut pow_x :[ExtFieldG2Element<N, PRAMASIZE>;MAX_COEFS_COUNT] = [self.base_field.one();MAX_COEFS_COUNT];           
                let mut pow_z2 :[ExtFieldG2Element<N, PRAMASIZE>;MAX_COEFS_COUNT] = [self.base_field.one();MAX_COEFS_COUNT];           
                let mut pow_x_z2 :[ExtFieldG2Element<N, PRAMASIZE>;MAX_COEFS_COUNT] = [self.base_field.one();MAX_COEFS_COUNT];                      
                for i in 1..MAX_COEFS_COUNT{ pow_x[i]  = pow_x[i-1].multiply(&x);
                                                    pow_z2[i] = pow_z2[i-1].multiply(&d2); 
                                                  }
                for i in 0..MAX_COEFS_COUNT{ pow_x_z2[i] = pow_x[i].multiply(&pow_z2[MAX_COEFS_COUNT - i - 1]);}
                let mut xnum = x.zero();
                let mut xden = x.zero();
                let mut ynum = x.zero();
                let mut yden = x.zero();
                for i in 0..MAX_COEFS_COUNT {
                                ynum = ynum.addto(&self.consts.swu_consts.ynum[i].multiply(&pow_x_z2[i]));
                                yden = yden.addto(&self.consts.swu_consts.yden[i].multiply(&pow_x_z2[i]));
                                xnum = xnum.addto(&self.consts.swu_consts.xnum[i].multiply(&pow_x_z2[i]));
                                xden = xden.addto(&self.consts.swu_consts.xden[i].multiply(&pow_x_z2[i]));
                                }             
                // Convert from affine coordinates to Jacobian coordinates
                yden = yden.multiply(&d3); 
                let z = xden.multiply(&yden);
                x = xnum.multiply(&z).multiply(&yden);
                y = ynum.multiply(&z.sqr()).multiply(&xden).multiply(&y);
                G2Element { point : EcPoint { x ,y ,z },
                            consts :self.consts,
                          }
            }

            pub fn random_point_withseed(&self,seed1 : &ExtFieldG2Element<N, PRAMASIZE>) -> G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT>
                {   
                    let rp = self.map_to_curve(seed1);
                    rp.clean_cofactor()
                }
            
            pub fn random_point(&self) -> G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT>
                {   
                    self.map_to_curve(&self.base_field.zero()).clean_cofactor().to_affine()                    
                }

            pub fn random_point_trys(&self) -> G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT>
            {   
                let mut found =false;
                let mut x =  self.base_field.zero();
                let mut y =  self.base_field.zero();
                let z =  self.base_field.one();
                while !found {
                    x =  self.base_field.random_element();
                    y = x.multiply(&x.sqr()).addto(&self.consts.b);  
                    found = y.is_qr();                    
                }
                y = y.sqrt().unwrap();
                G2Element { point : EcPoint { x ,y , z},
                            consts :self.consts,
                          }
            }
    
            pub fn hash_to_field(&self,id :&str, mode :u8) -> G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT>
                {                   
                    //  mode=0: Encode to Curve (NU-encode-to-curve), mode=1: Random Oracle model (RO-hash-to-curve)
                    //  https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-hash-to-curve-06#name-roadmap                    
                    if mode ==0 {   let hashs = self.base_field.hash_to_field(id,self.consts.security_level,1);
                                    self.random_point_withseed(&hashs[0]).to_affine()}
                    else {  let hashs = self.base_field.hash_to_field(id,self.consts.security_level,2);
                            let p1=self.random_point_withseed(&hashs[0]);
                            let p2=self.random_point_withseed(&hashs[1]);
                            p1.addto(&p2).to_affine()
                         } 
                }
            
            pub fn from_bytearray(&self,inbytes : &Vec<u8>) -> G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT>
                {
                    //  Point de-compression/de-Serialization as described by ZCach serialization format
                    //  https://www.ietf.org/archive/id/draft-irtf-cfrg-pairing-friendly-curves-11.html#name-zcash-serialization-format-
                    let mut input = inbytes.clone();
                    let m_byte = input[0] & 0xE0;
                    let numbits = self.consts.base_field_numbits;
                    let sizeinbytes = (numbits >> 3) + if (numbits % 8) ==0 {0} else {1};
                    if self.consts.base_field_numbits % 8 <=5 {input[0] = input[0] & 0x1F;}
                    else {input.remove(0);};
                    if m_byte == 0xE0 {panic!("Invalide compressed point format ...")};
                    if m_byte & 0x80 !=0 { if input.len() != sizeinbytes * self.getorder() {panic!("Invalide compressed point format ...")} 
                                         }
                    else {if input.len() != (sizeinbytes * 2 * self.getorder()) {panic!("Invalide compressed point format ...")}}
                    if m_byte & 0x40 !=0 { if input.iter().any(|&e| e != 0) {panic!("Invalid compression of an infinity point...");} 
                                           else {G2Element {  point : EcPoint {x:self.base_field.one(), y : self.base_field.one(), z: self.base_field.zero() },
                                                              consts :self.consts}
                                                } 
                                         }
                    else {  if input.len() == (sizeinbytes * 2 * self.getorder()){ let x = self.base_field.from_i2osp_bytearray(&input[0..sizeinbytes * self.getorder()]);
                                                               let y = self.base_field.from_i2osp_bytearray(&input[sizeinbytes * self.getorder()..]);
                                                               G2Element {  point : EcPoint { x, y, z: self.base_field.one() },
                                                               consts :self.consts}  
                                                             }
                            else { let x = self.base_field.from_i2osp_bytearray(&input[0..]);
                                   let y = x.sqr().multiply(&x).addto(&&self.consts.b).sqrt();
                                   if y.is_none() {panic!("Invalide point: not in the curve ...")}
                                   else { let y = y.unwrap();
                                          let r_sign = if m_byte & 0x20 !=0 {1} else {0}; 
                                          if (y.sign()+1) >> 1 == r_sign {G2Element {  point : EcPoint { x, y, z: self.base_field.one() },
                                                                                       consts :self.consts}  }
                                          else {G2Element {  point : EcPoint { x, y: y.negate(), z: self.base_field.one() },
                                                             consts :self.consts}  }
                                        }
                                 }                                 
                         }
                }
    
                pub fn from_base64(&self,input :&str) -> G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT>
                {
                    let decoded_bytes = match general_purpose::STANDARD.decode(input) {
                        Ok(bytes) => bytes,
                        Err(_) => {panic!("Failed to decode base64 string");}
                    };
                    self.from_bytearray(&decoded_bytes)
                }

                pub fn default_generator(&self)->G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT>
                {
                    {G2Element {  point : EcPoint { x :self.consts.default_generator.x, 
                                                    y: self.consts.default_generator.y, 
                                                    z: self.consts.default_generator.z },
                                  consts :self.consts}  } 
                }
        }
        
    
impl <const PRAMASIZE:usize,const R:usize,const N:usize,const MAX_COEFS_COUNT:usize> fmt::Display for G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT> {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            write!(f, "{:}", &self.to_string())
        }
    }
impl  <const PRAMASIZE:usize,const R:usize,const N:usize,const MAX_COEFS_COUNT:usize> Add for G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT> {
        type Output =  G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT>;
            fn add(self, rhs: Self) -> Self::Output {   self.addto(&rhs) }
    }
impl  <const PRAMASIZE:usize,const R:usize,const N:usize,const MAX_COEFS_COUNT:usize> Sub for G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT> {
        type Output =  G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT>;
            fn sub(self, rhs: Self) -> Self::Output {   self.substract(&rhs) }
    }
impl  <const PRAMASIZE:usize,const R:usize,const N:usize,const MAX_COEFS_COUNT:usize> Neg for G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT> {
        type Output =  G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT>;
            fn neg(self) -> Self::Output { self.negate() }
        }   
impl  <const PRAMASIZE:usize,const R:usize,const N:usize,const MAX_COEFS_COUNT:usize> PartialEq for G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT> {
            fn eq(&self, other: &Self) -> bool {    self.equal(other) }
        }
impl <const PRAMASIZE:usize,const R:usize,const N:usize,const MAX_COEFS_COUNT:usize> Mul<i128> for G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT> {
        type Output = G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT>;    
            fn mul(self, rhs: i128) -> Self::Output {     self.multiply_by_const(rhs)   }
        }
impl <const PRAMASIZE:usize,const R:usize,const N:usize,const MAX_COEFS_COUNT:usize> Mul<G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT>> for i128 {
        type Output = G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT>;
        fn mul(self, rhs: G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT>) -> Self::Output {  rhs.multiply_by_const(self)  }
        }
impl <const PRAMASIZE:usize,const R:usize,const N:usize,const MAX_COEFS_COUNT:usize> Mul<u64> for G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT> {
        type Output = G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT>;    
            fn mul(self, rhs: u64) -> Self::Output {     self.multiply_by_const(rhs as i128)   }
        }
impl <const PRAMASIZE:usize,const R:usize,const N:usize,const MAX_COEFS_COUNT:usize> Mul<G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT>> for u64 {
        type Output = G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT>;
            fn mul(self, rhs: G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT>) -> Self::Output {  rhs.multiply_by_const(self as i128)  }
        }    
impl <const PRAMASIZE:usize,const R:usize,const N:usize,const MAX_COEFS_COUNT:usize> Mul<i64> for G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT> {
        type Output = G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT>;    
            fn mul(self, rhs: i64) -> Self::Output {     self.multiply_by_const(rhs as i128)   }
        }            
impl <const PRAMASIZE:usize,const R:usize,const N:usize,const MAX_COEFS_COUNT:usize> Mul<G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT>> for i64 {
        type Output = G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT>;
            fn mul(self, rhs: G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT>) -> Self::Output {  rhs.multiply_by_const(self as i128)  }
        }    
impl <const PRAMASIZE:usize,const R:usize,const N:usize,const MAX_COEFS_COUNT:usize> Mul<u8> for G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT> {
        type Output = G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT>;    
            fn mul(self, rhs: u8) -> Self::Output {     self.multiply_by_const(rhs as i128)   }
        }            
impl <const PRAMASIZE:usize,const R:usize,const N:usize,const MAX_COEFS_COUNT:usize> Mul<G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT>> for u8 {
        type Output = G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT>;
            fn mul(self, rhs: G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT>) -> Self::Output {  rhs.multiply_by_const(self as i128)  }
        }    
impl <const PRAMASIZE:usize,const R:usize,const N:usize,const MAX_COEFS_COUNT:usize> Mul<i8> for G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT> {
        type Output = G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT>;    
            fn mul(self, rhs: i8) -> Self::Output {     self.multiply_by_const(rhs as i128)   }
        }            
impl <const PRAMASIZE:usize,const R:usize,const N:usize,const MAX_COEFS_COUNT:usize> Mul<G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT>> for i8 {
        type Output = G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT>;
            fn mul(self, rhs: G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT>) -> Self::Output {  rhs.multiply_by_const(self as i128)  }
        }    
impl <const PRAMASIZE:usize,const R:usize,const N:usize,const MAX_COEFS_COUNT:usize> Mul<G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT>> for FieldElement<R> {
        type Output = G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT>;
            fn mul(self, rhs: G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT>) -> Self::Output {  rhs.multiply_gls(&self)  }
        }            
impl <const PRAMASIZE:usize,const R:usize,const N:usize,const MAX_COEFS_COUNT:usize> Mul<FieldElement<R>> for G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT>   {
        type Output = G2Element<PRAMASIZE,R,N,MAX_COEFS_COUNT>;
            fn mul(self, rhs: FieldElement<R>) -> Self::Output {  self.multiply_gls(&rhs)  }
        }                        