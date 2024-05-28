// Code developped by FARAOUN Kamel Mohamed.
// EEDIS-Laboratory. UDL-University. Algeria
// During May 2024.

use std::{ops::{BitAnd, BitOr}, str::FromStr};
use num_bigint::BigUint;
use num_traits::{FromPrimitive, One, ToPrimitive};


pub const WSIZE  :u8 = 3;
pub const WDSIZE :u8 = 6;
pub const WMASK  :u8 = 7;
pub const WDMASK :u8 = 63;


// Scalar Recoding for w-sliding multiplication algorithm

// Scalar recoding for multiplication with scalar (w=3) on E(Fp):subscalars ar all in {-1,1} (zeros elimination, its encode -1 in the output) 
pub fn recod_one_scalar(scalar :&BigUint) -> BigUint
    {   let mut mu    = scalar.bits();
        let mut x = scalar.clone();
        let ff: BigUint = BigUint::from_str("255").unwrap();
        mu = mu + (WSIZE as u64) - (mu % (WSIZE as u64));
        x = x.bitor(BigUint::one() << mu);        
        let mut code = BigUint::one();        
        while x != BigUint::one(){  let mut limb = (&x).bitand((ff).clone()).to_u8().unwrap_or(0);
                                    let sign = ((limb >> (WSIZE - 1)) & 2) as i8 - 1 ;                
                                    limb = (((limb ^ (sign as u8)) & WMASK) | 1) - ((sign >>1).abs() as u8);  
                                    code = (code << WSIZE).bitor(&BigUint::from_u8(limb as u8).unwrap());
                                    x = x >> WSIZE;
                                }
        code >> 1
    }

// Scalar recoding for GLV-2 : personalized recoding according to Algorithm 12 from the paper: "Optimizing and securing GLV multiplication over BLS pairings-friendly curves".
pub fn recod_scalar_glv2(scalar:&BigUint, lambda : &BigUint) -> BigUint
    {   let mut x1 = scalar % lambda;
        let mut x2 = scalar / lambda;
        let ff: BigUint = BigUint::from_str("255").unwrap();
        let t1 = (&x1).bitand((ff).clone()).to_u8().unwrap_or(0);
        let t2 = lambda.bitand((ff).clone()).to_u8().unwrap_or(0);
        let beta = BigUint::from((!t1 & 1)*(t2 & 1));
        x1 = &x1 + &beta * lambda;
        x2 = &x2 - &beta;
        let mut mu = ((&x1).bitor(&x2)).bits();
        mu = mu + (WDSIZE as u64) - (mu % (WSIZE as u64));
        x1 = (&x1).bitor(BigUint::one() << mu); 
        let mut code = BigUint::one();    
        while x1 != BigUint::one(){ let limb_x1 = (&x1).bitand((ff).clone()).to_u8().unwrap_or(0);
                                    let limb_x2 = (&x2).bitand((ff).clone()).to_u8().unwrap_or(0);
                                    let sign = ((limb_x1 >> WSIZE) & 1) as i8 - 1 ;  
                                    let even = !limb_x2 & 1;
                                    let mut ai = ((limb_x1 ^ (sign as u8)) | 1) & WMASK;
                                    let mut bi = (((limb_x2 as i16 + sign as i16) as u8) ^ (sign as u8)) & WMASK;
                                    let di : i16 = ai as i16 - bi as i16;
                                    let tmp = (WMASK as i16 - (((sign as u8) | 1) as i16 * di)) as u8 ;
                                    let mut inc = (tmp & (bi + WMASK)) >> WSIZE;
                                    ai = (((ai  as i16 - (bi * even) as i16)) as u8 ) & WMASK;
                                    bi = (di * even as i16 + bi as i16) as u8; 
                                    inc = ((inc * even) as i8 + ((even as i8 - 1) * (sign as i8)) as i8) as u8;
                                    let limb = even + (((sign as i8 + 1) as u8) << 1) + (((ai >> 1 ) + ((bi - 1) << 1)) << 2);
                                    code = (code << WDSIZE).bitor(&BigUint::from_u8(limb).unwrap());
                                    x1 =  x1 >> WSIZE;
                                    x2 = (x2 >> WSIZE) + BigUint::from(inc);
                                 }
        code  
    }

// Scalar recoding for GLS-4 (for E(Fp2)): first subscalar in {-1,1}, remaining sub-scalars in {-1,0,1} with alignement to the first one    
pub fn recod_scalar_gls4(scalar : &BigUint, scalar_neg : &BigUint, u : u128 ) -> BigUint 
    {   let alpha = (u |((scalar % u).to_u128().unwrap())) & 1;
        let mut x = (scalar_neg) * (BigUint::one() - alpha) + (scalar * alpha);           
        let mut decs :[u128;4]= [0;4];
        for i in 0..4 {decs[i] = (&x % u).to_u128().unwrap();
                              x = &x / u; 
                             }                
        let beta = ((decs[0] ^ 1)& 1) * (u & 1);                
        decs[0] = decs[0] + beta * u;         
        if beta == 1 { let mut i = 0;
                       while decs[i] == 0 { decs[i] = u.clone();
                                            i = i + 1;}
                       decs[i] = decs[i] - beta;                                                             
                     }
        let mut mu =0;
        for i in decs{ mu = std::cmp::max(mu,128-i.leading_zeros())}        
        decs[0] = (decs[0] >> 1) | (1 << mu);        
        for i in 1..4 { let mut indic :i128 = -1;
                               while indic != 0 { indic = ((decs[i] & (!decs[0])) & (indic as u128)) as i128;
                                                  indic = indic ^ (-indic);  
                                                  decs[i] = decs[i] + (indic.abs() as u128) ;
                                                }   
                                }
        let mut code = BigUint::one();    
        while decs[0]!=0 { let mut limb = (decs[0] & 1) as u8;
                                   decs[0] = decs[0] >> 1;      
                                   for i in 1..4 { let first = (decs[i] & 1) as u8;
                                                          limb = limb | (first << i);
                                                          decs[i] = decs[i] >> 1}
                                   code = (code << 4u8).bitor(&BigUint::from_u8(limb).unwrap());
                                  }
        code = (code << 1u8).bitor(BigUint::one() - &alpha);                          
        code                          
    }

// Scalar recoding for GLS-8 (For E(Fp4)): first and fourth subscalar in {-1,1}, remaining sub-scalars in {-1,0,1} with alignement to the first/fourth  one
// Ensuring both first and  fourth sub-scalars are odd numbers using combination with "u"
pub fn recod_scalar_gls8(scalar : &BigUint, u:u128,r:&BigUint) -> BigUint 
    {   let alpha = (u |((scalar % u).to_u128().unwrap())) & 1;
        let mut x = (r - scalar) * (BigUint::one() - alpha) + (scalar * alpha);                   
        let mut decs :[u128;8]= [0;8];        
        for i in 0..8 {decs[i] = (&x % u).to_u128().unwrap();
                              x = &x / u; 
                             }               
        let beta = ((decs[0] ^ 1)& 1) * (u & 1);                
        decs[0] = decs[0] + beta * u;         
        let mut i = 1;
        while decs[i] == 0 {decs[i] = u;
                            i = i + 1;}
        decs[i] = decs[i] - beta;                                                                         
        let beta =  !decs[4]  & 1;                                 
        decs[3] = decs[3] + beta * u;         
        let mut i = 4;
        while decs[i] == 0 {decs[i] = u;
                            i = i + 1;}                                    
        decs[i] = decs[i] - beta;                                                             
        let mut mu =0;
        for i in &decs{ mu = std::cmp::max(mu,128-i.leading_zeros())}                                    
        decs[0] = (decs[0] >> 1) | (1 << mu);
        decs[4] = (decs[4] >> 1) | (1 << mu);                                                          
        for i in 1..8 { if i !=4 
                                { let mut indic :i128 = -1;
                                  while indic != 0 { indic = ((decs[i] & (!decs[i & 4])) & (indic as u128)) as i128;
                                                     indic = indic ^ (-indic);  
                                                     decs[i] = decs[i] + (indic.abs()) as u128;
                                                    }   
                                }
                             }                                    
        let mut code = BigUint::one();    
        while decs[0]!=0 {  let mut limb = (decs[0] & 1) as u8;
                            decs[0] = decs[0] >> 1;      
                            for i in 1..8 { let first = (decs[i] & 1) as u8;
                                                   limb = limb | (first << i);
                                                   decs[i] = decs[i] >> 1
                                                 }
                            code = (code << 8u8).bitor(&BigUint::from_u8(limb).unwrap());
                         }
        code = (code << 1u8).bitor(BigUint::one() - &alpha);                          
        code                                  
    }

// Scalar recoding for GLS-16 (for E(Fp8)) : 1st,4th,8th and 12th in {-1,1}, remaining sub-scalars in {-1,0,1} with alignement to the recoded ones
// Ensuring 1st,4th,8th and 12th sub-scalars are odd numbers using combination with "u"
pub fn recod_scalar_gls16(scalar : &BigUint, u:u64, r:&BigUint) -> BigUint 
    {   let alpha = (u |((scalar % u).to_u64().unwrap())) & 1;
        let mut x = (r - scalar) * (BigUint::one() - alpha) + (scalar * alpha);           
        let mut decs :[u64;16]= [0;16];
        for i in 0..16 { decs[i] = (&x % u).to_u64().unwrap();
                                x = &x / u; 
                              }        
        let beta = ((decs[0] ^ 1) & 1) * (u & 1);      
        decs[0] = decs[0] + beta * u;         
        let mut i = 1;
        while decs[i] == 0 {decs[i] = u;
                            i = i + 1;}
        decs[i] = decs[i] - beta;                      
        for j in (4..=12).step_by(4){  let beta =  ((decs[j] ^ 1) & 1) * (u & 1);
                                                    decs[j-1] = decs[j-1] + (beta * u);
                                                    i = j;
                                                    while decs[i] == 0 { decs[i] = u;
                                                                         i = i + 1;}                                    
                                                    decs[i] = decs[i] - beta;       
                                                }                                               
        let mut mu =0;
        for i in &decs{ mu = std::cmp::max(mu,64-i.leading_zeros())}                                    
        decs[0] = (decs[0] >> 1) | (1 << mu);
        decs[4] = (decs[4] >> 1) | (1 << mu);
        decs[8] = (decs[8] >> 1) | (1 << mu);
        decs[12] = (decs[12] >> 1) | (1 << mu);
        for i in 1..16 { if (i % 4) != 0 
                                { let mut indic :i64 = -1;
                                  while indic != 0 { indic = ((decs[i] & (!decs[i - (i % 4)])) & (indic as u64)) as i64;
                                                     indic = indic ^ (-indic);  
                                                     decs[i] = decs[i] + (indic.abs()) as u64;
                                                    }   
                                }
                             }                                    
        let mut code = BigUint::one();    
        while decs[0]!=0 {  let mut limb = (decs[0] & 1) as u16;
                            decs[0] = decs[0] >> 1;      
                            for i in 1..16 { let first = (decs[i] & 1) as u16;
                                                    limb = limb | (first << i);
                                                    decs[i] = decs[i] >> 1
                                                  }
                            code = (code << 16u8).bitor(&BigUint::from_u16(limb).unwrap());
                         }
        code = (code << 1u8).bitor(BigUint::one() - &alpha);                          
        code                                  
    }
