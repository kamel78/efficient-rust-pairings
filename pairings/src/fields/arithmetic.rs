// // Code developed by FARAOUN Kamel Mohamed.
// // EEDIS-Laboratory. UDL-University. Algeria
// // During May 2024.

// pub fn equal <const N: usize> (a : &[u64; N], b : &[u64; N]) -> bool
//         // Constant-time equality comparison of two [u64; N]
//     {   let mut res=0;
//         for (i,j) in a.iter().zip(b){res += (i==j) as u8} ;
//         res == N as u8
//     }

// pub fn add <const N: usize> (a : &[u64; N], b : &[u64; N], _params :&super::prime_fields::FieldParams<N>) -> [u64; N] 
//         //  Addition with reduction on montgomery form
//     {   let mut result: [u64; N]= [0; N];
//         let mut carry :bool = false;
//         let mut bigger:i32=0;
//         for i in 0..N { (result[i],carry) = (a[i]).overflowing_add((b[i]) + (carry as u64));
//                         bigger += (((result[i] > _params.modulo[i]) as i32) - ((result[i] < _params.modulo[i]) as i32)) << i;
//                         };
//         if  bigger >= 0 {   let mut borrow:bool=false; 
//                             for i in 0..N {
//                                 (result[i],borrow) =  
//                                             (result[i]).overflowing_sub((&_params.modulo[i]) + (borrow as u64))}
//                         }
//     result
//     }

// pub fn sub <const N: usize> (a : &[u64; N], b : &[u64; N], _params :&super::prime_fields::FieldParams<N>)-> [u64; N] 
//         // Substration with reduction into montgomery form
//     {   let mut result: [u64; N]= [0; N];
//         let mut borrow: bool =false;
//         for i in 0..N {(result[i],borrow) = (a[i]).overflowing_sub((b[i]) + (borrow as u64))};
//         if borrow { let mut carry:bool=false;
//                     for i in 0..N { (result[i],carry) =  
//                                             (result[i]).overflowing_add((&_params.modulo[i]) + (carry as u64))} 
//                     }
//     result 
//     }

// pub fn neg <const N: usize> (a : &[u64; N],_params :&super::prime_fields::FieldParams<N>) ->[u64; N]
//         // Negate an element from Fp (MODULO - a) 
//     {   let mut result: [u64; N]= [0; N];
//         let mut borrow: bool =false;
//         for i in 0..N {(result[i],borrow) = (&_params.modulo[i]).overflowing_sub((a[i]) + (borrow as u64))};
//     result 
//     }

// pub fn mul_std<const N: usize> (a : &[u64; N], b : &[u64; N],_params :&super::prime_fields::FieldParams<N>)-> [u64; N] 
//     // CIOS implementation of montgomery multiplication (https://www.microsoft.com/en-us/research/wp-content/uploads/1998/06/97Acar.pdf)
// { 
//         let mut result: [u64; N]= [0 ;N];
//         let mut ci   :u128;
//         let mut t128 :u128;
//         let mut m    :u128;
//         let mut tn   :u64 = 0;
//         let mut tnp1 :u128;
//         for i in 0..N{   ci = 0;
//                                 for j in 0..N {  t128 = result[j] as u128 +  (a[j] as u128 * b[i] as u128) + ci ;
//                                                         (ci, result[j]) = (t128 >> 64 , t128 as u64); 
//                                                      }
//                                 t128 = tn as u128 + ci ;
//                                (tnp1, tn) = (t128 >> 64 , t128 as u64);  
//                                 m  = (result[0] as u128 * _params.qprime) & ((1 << 64) - 1);
//                                 ci = (result[0] as u128 + m * _params.modulo[0] as u128) >> 64;
//                                 for j in 1..N { t128 = result[j] as u128+  (m * _params.modulo[j] as u128) + ci;
//                                                        (ci, result[j-1]) = (t128 >> 64 , t128 as u64);
//                                                      }
//                                 t128 = tn as u128 + ci;
//                                (ci, result[N-1]) = (t128 >> 64 , t128 as u64);
//                                 tn = (ci + tnp1) as u64;
//                             }
//         let mut bigger:i32=0;
//         for i in 0..N { bigger += (((result[i] > _params.modulo[i]) as i32) - ((result[i] < _params.modulo[i]) as i32)) << i};
//         if bigger >= 0 { let mut borrow:bool=false; 
//                             for i in 0..N {(result[i],borrow) =  
//                                                   (result[i]).overflowing_sub((&_params.modulo[i]) + (borrow as u64))}
//                        }
//         result
// }

// pub fn mul<const N: usize> (a : &[u64; N], b : &[u64; N],_params :&super::prime_fields::FieldParams<N>)-> [u64; N] 
//         // CIOS implementation of montgomery multiplication (https://www.microsoft.com/en-us/research/wp-content/uploads/1998/06/97Acar.pdf)
//         // Optimized using trick from  https://hackmd.io/@gnark/modular_multiplication#fn1
//     {   
//         if  !_params.optimize_mul {mul_std(&a, &b, _params)}
//         else {  let mut result: [u64; N]= [0 ;N];
//                 let mut ai   :u128;
//                 let mut ci   :u128;
//                 let mut t128 :u128;
//                 let mut m    :u128; 
//                 for i in 0..N{   t128  = (result[0] as u128) + (a[0] as u128 * b[i] as u128);
//                                         (ai , result[0]) = (t128 >> 64 , t128 as u64);                                        
//                                         m  = (result[0] as u128 * _params.qprime) & ((1 << 64) - 1);
//                                         ci = (result[0] as u128 + m * _params.modulo[0] as u128) >> 64;
//                                         for j in 1..N {  t128 = result[j] as u128 +  (a[j] as u128 * b[i] as u128) + ai ;
//                                                             (ai, result[j]) = (t128 >> 64 , t128 as u64);
//                                                             t128 = result[j] as u128+  (m * _params.modulo[j] as u128) + ci;
//                                                             (ci, result[j-1]) = (t128 >> 64 , t128 as u64);
//                                                         }
//                                         result[N-1] = (ci + ai) as u64;                                        
//                                     }
//                 let mut bigger:i32=0;
//                 for i in 0..N { bigger += (((result[i] > _params.modulo[i]) as i32) - ((result[i] < _params.modulo[i]) as i32)) << i};
//                 if bigger >= 0 {    let mut borrow:bool=false; 
//                                     for i in 0..N {(result[i],borrow) =  
//                                                         (result[i]).overflowing_sub((&_params.modulo[i]) + (borrow as u64))}
//                                 }
//                 result
//             }
//     }

// pub fn sqr <const N: usize> (a : &[u64; N],_params :&super::prime_fields::FieldParams<N>)-> [u64; N] 
//         // CIOS implementation of montgomery squaring (https://www.microsoft.com/en-us/research/wp-content/uploads/1998/06/97Acar.pdf)
//         // Optimized using trick from  https://hackmd.io/@gnark/modular_multiplication#fn1        
//     {   if ! _params.optimize_sqr {mul(&a, &a, _params)}
//         else {
//                 let mut result: [u64; N]= [0 ;N];
//                 let mut ai   :u128;
//                 let mut ci   :u128;
//                 let mut t128 :u128;
//                 let mut m    :u128;
//                 let mut carry1    :bool;
//                 let mut carry2    :bool;
//                 let mut carry3    :bool;
//                 for i in 0..N{ t128  = (result[i] as u128) + (a[i] as u128 * a[i] as u128);
//                                         (ci , result[i]) = (t128 >> 64 , t128 as u64);
//                                         ai = 0;
//                                         for j in i+1..N {    t128  = a[j] as u128 * a[i] as u128 ;
//                                                                     (t128, carry1)  = t128.overflowing_add(t128) ;
//                                                                     (t128, carry2)  = t128.overflowing_add(ci + result[j] as u128) ;
//                                                                     (t128, carry3)  = t128.overflowing_add(((ai as u128) << 64) as u128) ;
//                                                                     ai = ((carry1 as u8) + (carry2 as u8) + (carry3 as u8)) as u128;
//                                                                     (ci, result[j]) = (t128 >> 64,t128 as u64);
//                                                                 }
//                 ai = ci;
//                 m  = (result[0] as u128 * _params.qprime) & ((1 << 64) - 1);
//                 ci = (result[0] as u128 + m * _params.modulo[0] as u128) >> 64;
//                 for j in 1..N {   t128 = result[j] as u128 + (m * _params.modulo[j] as u128) + ci;
//                                         (ci, result[j-1]) = (t128 >> 64 , t128 as u64);
//                                       }   
//                 result[N-1] = (ci + ai) as u64;
//                 }
//             let mut bigger:i32=0;
//             for i in 0..N {  bigger += (((result[i] > _params.modulo[i]) as i32) - ((result[i] < _params.modulo[i]) as i32)) << i};
//                                     if bigger >= 0 {    let mut borrow:bool=false; 
//                                                         for i in 0..N {(result[i],borrow) =  
//                                                                               (result[i]).overflowing_sub((&_params.modulo[i]) + (borrow as u64))}
//                                     }
//             result
//             }
//     } 

// pub fn invert <const N: usize> (a : &[u64; N],_params :&super::prime_fields::FieldParams<N>) -> [u64; N]
//         // Modular ivnertion of a modulo p. Implementation from "Efficient Software-Implementation of    
//         // Finite Fields with Applications to Cryptography" - Algorithm 16. https://link.springer.com/article/10.1007/s10440-006-9046-1        
//     {   #[inline(always)]
//         fn cmp <const N: usize> (a : &[u64; N], b : &[u64; N])-> i8 
//                 // Comparing tow integers based on [u64; N] representation. Proposed constant-time implementation
//             {   let mut result:i16 = 0;
//                 for (index, (ai,bi)) in (a.iter().zip(b.iter())).enumerate() {result += (((ai > bi) as i16) - ((ai < bi) as i16)) << index;}
//                 result.signum() as i8
//             }       
//         #[inline(always)]   
//         fn add_p <const N: usize> (result:&mut [u64; N], p: &[u64; N] ) 
//                 // In-site Addition with modulo (result =result + modulo)
//             {   let mut carry :bool = false;
//                 for i in 0..N { (result[i],carry) = (result[i]).overflowing_add(p[i] + (carry as u64))};
//             }        
//        // #[inline(always)]
//         fn sub_noreduction <const N: usize> (a : &[u64; N],result:&mut [u64; N]) -> bool 
//                 // In-site Substraction without reduction (result =result - a)
//             {   let mut borrow: bool =false;
//                 let mut borrow1: bool ;
//                 let mut borrow2: bool ;
//                 // for i in 0..N {(result[i],borrow) = (result[i]).overflowing_sub((a[i]) + (borrow as u64))};
//                 for i in 0..N {(result[i],borrow1) = (result[i]).overflowing_sub(a[i]);
//                                       (result[i],borrow2) = (result[i]).overflowing_sub(borrow as u64);
//                                      borrow= borrow1 | borrow2   };
//                 borrow 
//             }
//         let mut u = a.clone();
//         let mut v = _params.modulo;
//         let mut b = _params.rsquare; 
//         let mut c = _params.zero;
//         loop {  while u[0] & 1 == 0 {   r_shift(&mut u);                                    
//                                         if b[0] & 1 == 1 {add_p(&mut b,&_params.modulo)};  
//                                         r_shift(&mut b)  }
//                 while v[0] & 1 == 0 {   r_shift(&mut v);
//                                         if c[0] & 1 == 1 {add_p(&mut c,&_params.modulo)};  
//                                         r_shift(&mut c)  }                                    
//                 if cmp(&u , &v) == 1 {  sub_noreduction(&v,&mut u);
//                                         if sub_noreduction(&c,&mut b) {add_p(&mut b,&_params.modulo)}}
//                 else {  sub_noreduction(&u,&mut v);
//                         if sub_noreduction(&b,&mut c) {add_p(&mut c,&_params.modulo)}}
//                 let mut _null_u =0;                        
//                 let mut _null_v =0;                        
//                 for i in 1..N { _null_u = _null_u | u[i];
//                                        _null_v = _null_v | v[i]}
//                 if (u[0]==1) & (_null_u == 0) {return b }
//                 if (v[0]==1) & (_null_v == 0) {return c }
                
//             }
//     }

// pub fn pow <const N: usize> (a : &[u64; N], e :&[u64], useladder:bool,_params :&super::prime_fields::FieldParams<N>)-> [u64; N]
//         // Implements montgomery Ladder (secure but relatively slow with respect to square and multiply)
//     {   if useladder {  let mut r0 = _params.one;     
//                         let mut r1 = a.clone();                            
//                         for i in e.as_ref().iter().rev() { 
//                                 for j in (0..64).rev()  {   if (i >> j) & 1 == 0 {  r1 = mul(&r0, &r1, _params);
//                                                                                          r0 = sqr(&r0, _params); }                
//                                                             else {  r0 = mul(&r0, &r1,_params);
//                                                                     r1 = sqr(&r1,_params)}
//                                                         }
//                                                         }
//                         r0
//                         }
//         else {  let mut result = _params.one;     
//                 for i in e.as_ref().iter().rev()  {
//                             for j in (0..64).rev(){   result = sqr(&result,_params);                
//                                                         if (i >> j) & 1 == 1 {result = mul(&result, a,_params);}                
//                                                     }
//                                                     }
//                 result
//             }
//     }

// pub fn sqrt<const N: usize> (a:&[u64; N],_params :&super::prime_fields::FieldParams<N>) -> Option<[u64; N]>
//     {   if _params.sqrtid == -1 { // p mod 4 ==3
//                                     let result = pow(a,&_params.modplus1div4,false,_params); 
//                                     if equal(&mul(&sqr(&result,_params),&invert(a,_params),_params),&_params.one) { Some(result) }
//                                     else { None }
//                                 }
//         else {      // Tonelli-Shanks Algorithm (p mod 4 ==1)
//                     //  Adapted from  section 6 of 'Square roots from 1; 24, 51, 10 to Dan Shanks'  by Ezra Brown:
// 	                // https://www.maa.org/sites/default/files/pdf/upload_library/22/Polya/07468342.di020786.02p0470a.pdf 
//                     let w = pow(a,&_params.tonelli_params[1],false,_params);                    
//                     let mut y = mul(&a,&w,_params);
//                     let mut b = mul(&y,&w,_params);
//                     let mut g = mul(&_params.tonelli_params[0],&_params.rsquare,_params); // Convert to Montgomery repr
//                     let mut r = _params.sqrtid;
//                     let mut t = b.clone();
//                     for _ in 0..r-1 { t = sqr(&t, _params)}; 
//                     if equal(&t, &[0;N]) {return Some([0;N])}
//                     if !equal(&t, &_params.one) {return None}       
//                     loop  { let mut t= b.clone();
//                             let mut m= 0;
//                             while !equal(&t, &_params.one) {  t = sqr(&t, _params);
//                                                                    m = m + 1;
//                                                                 }
//                             if m==0 { return Some(y)}
//                             let mut ge = r - m - 1;
//                             t = g.clone();                            
//                             while ge > 0 { t = sqr(&t,_params);
//                                            ge = ge - 1}
//                             g = sqr(&t,_params);
//                             y = mul(&y, &t, _params);
//                             b = mul(&b, &g, _params);   
//                             r = m;
//                         }
//                 }
//     }

// pub fn r_shift <const N: usize> (result:&mut [u64; N])
//     // Right shift on bit string 
// {   for i in 0..N-1{ result[i] = (result[i] >> 1) | ((result[i+1] & 1) << 63)}
//                         result[N - 1] = result [N - 1] >> 1
// } 
// Code developed by FARAOUN Kamel Mohamed.
// EEDIS-Laboratory. UDL-University. Algeria
// During May 2024.

pub fn equal <const N: usize> (a : &[u64; N], b : &[u64; N]) -> bool
        // Constant-time equality comparison of two [u64; N]
    {   let mut res=0;
        for (i,j) in a.iter().zip(b){res += (i==j) as u8} ;
        res == N as u8
    }

pub fn add <const N: usize> (a : &[u64; N], b : &[u64; N], _params :&super::prime_fields::FieldParams<N>) -> [u64; N] 
        //  Addition with reduction on montgomery form
    {   let mut result: [u64; N]= [0; N];
        let mut carry :bool = false;
        let mut carry1 :bool;
        let mut bigger:i32=0;
        let mut tmp:u64;
        for i in 0..N {  (tmp,carry1)      = b[i].overflowing_add(carry as u64);
                                (result[i],carry) = (a[i]).overflowing_add(tmp);
                                carry = carry | carry1;
                                bigger += (((result[i] > _params.modulo[i]) as i32) - ((result[i] < _params.modulo[i]) as i32)) << i;
                                };
        if  (bigger >= 0) | carry {   let mut borrow:bool=false; 
                                        for i in 0..N { (tmp,carry)= (_params.modulo[i]).overflowing_add(borrow as u64);
                                                                (result[i],borrow) =  (result[i]).overflowing_sub(tmp);
                                                                borrow = borrow | carry
                                                            }
                                  }
        result
    }

pub fn sub <const N: usize> (a : &[u64; N], b : &[u64; N], _params :&super::prime_fields::FieldParams<N>)-> [u64; N] 
        // Substration with reduction into montgomery form
    {   let mut result: [u64; N]= [0; N];
        let mut borrow: bool =false;
        let mut carry: bool =false;
        let mut tmp:u64;
        for i in 0..N {  (tmp,carry)        = b[i].overflowing_add(borrow as u64);
                                (result[i],borrow) = (a[i]).overflowing_sub(tmp)};
                                borrow = borrow | carry;
        if borrow { let mut carry:bool=false;
                    for i in 0..N {  (tmp,borrow) = (_params.modulo[i]).overflowing_add(carry as u64);   
                                            (result[i],carry) =  (result[i]).overflowing_add(tmp);
                                            carry = carry | borrow;
                                         } 
                    }
    result 
    }

pub fn neg <const N: usize> (a : &[u64; N],_params :&super::prime_fields::FieldParams<N>) ->[u64; N]
        // Negate an element from Fp (MODULO - a) 
    {   let mut result: [u64; N]= [0; N];
        let mut borrow: bool =false;
        for i in 0..N {(result[i],borrow) = (&_params.modulo[i]).overflowing_sub((a[i]) + (borrow as u64))};
        result 
    }

pub fn mul_std<const N: usize> (a : &[u64; N], b : &[u64; N],_params :&super::prime_fields::FieldParams<N>)-> [u64; N] 
    // CIOS implementation of montgomery multiplication (https://www.microsoft.com/en-us/research/wp-content/uploads/1998/06/97Acar.pdf)
{ 
        let mut result: [u64; N]= [0 ;N];
        let mut ci   :u128;
        let mut t128 :u128;
        let mut m    :u128;
        let mut tn   :u64 = 0;
        let mut tnp1 :u128;
        for i in 0..N{   ci = 0;
                                for j in 0..N {  t128 = result[j] as u128 +  (a[j] as u128 * b[i] as u128) + ci  ;
                                                        (ci, result[j]) = (t128 >> 64 , t128 as u64); 
                                                     }
                                t128  = (tn as u128) +ci ;
                                (tnp1, tn) = (t128 >> 64 , t128 as u64);    
                                m  = (result[0] as u128 * _params.qprime) & ((1 << 64) - 1);
                                ci = (result[0] as u128 + m * _params.modulo[0] as u128) >> 64;
                                for j in 1..N { t128 = result[j] as u128+  (m * _params.modulo[j] as u128) + ci;
                                                       (ci, result[j-1]) = (t128 >> 64 , t128 as u64);
                                                     }
                                t128  = (tn as u128) + ci ;
                               (ci, result[N-1]) = (t128 >> 64 , t128 as u64);
                                tn = (ci + tnp1) as u64;
                            }
        let mut bigger:i32=0;
        for i in 0..N { bigger += (((result[i] > _params.modulo[i]) as i32) - ((result[i] < _params.modulo[i]) as i32)) << i};
        if (bigger >= 0)| (tn!=0) {   let mut borrow:bool=false; 
                                      let mut carry: bool;
                                      let mut tmp:u64;
                                      for i in 0..N { (tmp,carry)= (_params.modulo[i]).overflowing_add(borrow as u64);
                                                             (result[i],borrow) =  (result[i]).overflowing_sub(tmp);
                                                             borrow = borrow | carry
                                                           }
      }                     
        result
}

pub fn mul<const N: usize> (a : &[u64; N], b : &[u64; N],_params :&super::prime_fields::FieldParams<N>)-> [u64; N] 
        // CIOS implementation of montgomery multiplication (https://www.microsoft.com/en-us/research/wp-content/uploads/1998/06/97Acar.pdf)
        // Optimized using trick from  https://hackmd.io/@gnark/modular_multiplication#fn1
    {   
        if  !_params.optimize_mul {mul_std(&a, &b, _params)}
        else {  let mut result: [u64; N]= [0 ;N];
                let mut ai   :u128;
                let mut ci   :u128;
                let mut t128 :u128;
                let mut m    :u128; 
                for i in 0..N{   t128  = (result[0] as u128) + (a[0] as u128 * b[i] as u128);
                                        (ai , result[0]) = (t128 >> 64 , t128 as u64);                                        
                                        m  = (result[0] as u128 * _params.qprime) & ((1 << 64) - 1);
                                        ci = (result[0] as u128 + m * _params.modulo[0] as u128) >> 64;
                                        for j in 1..N {  t128 = result[j] as u128 +  (a[j] as u128 * b[i] as u128) + ai ;
                                                            (ai, result[j]) = (t128 >> 64 , t128 as u64);
                                                            t128 = result[j] as u128+  (m * _params.modulo[j] as u128) + ci;
                                                            (ci, result[j-1]) = (t128 >> 64 , t128 as u64);
                                                        }
                                        result[N-1] = (ci + ai) as u64;                                    
                                    }
                let mut bigger:i32=0;
                for i in 0..N { bigger += (((result[i] > _params.modulo[i]) as i32) - ((result[i] < _params.modulo[i]) as i32)) << i};
                if bigger >= 0 {    let mut borrow:bool=false; 
                                    for i in 0..N {(result[i],borrow) =  
                                                        (result[i]).overflowing_sub((&_params.modulo[i]) + (borrow as u64))}
                                }
                result
            }
    }

pub fn sqr <const N: usize> (a : &[u64; N],_params :&super::prime_fields::FieldParams<N>)-> [u64; N] 
        // CIOS implementation of montgomery squaring (https://www.microsoft.com/en-us/research/wp-content/uploads/1998/06/97Acar.pdf)
        // Optimized using trick from  https://hackmd.io/@gnark/modular_multiplication#fn1        
    {   if ! _params.optimize_sqr {mul(&a, &a, _params)}
        else {
                let mut result: [u64; N]= [0 ;N];
                let mut ai   :u128;
                let mut ci   :u128;
                let mut t128 :u128;
                let mut m    :u128;
                let mut carry1    :bool;
                let mut carry2    :bool;
                let mut carry3    :bool;
                for i in 0..N{ t128  = (result[i] as u128) + (a[i] as u128 * a[i] as u128);
                                        (ci , result[i]) = (t128 >> 64 , t128 as u64);
                                        ai = 0;
                                        for j in i+1..N {    t128  = a[j] as u128 * a[i] as u128 ;
                                                                    (t128, carry1)  = t128.overflowing_add(t128) ;
                                                                    (t128, carry2)  = t128.overflowing_add(ci + result[j] as u128) ;
                                                                    (t128, carry3)  = t128.overflowing_add(((ai as u128) << 64) as u128) ;
                                                                    ai = ((carry1 as u8) + (carry2 as u8) + (carry3 as u8)) as u128;
                                                                    (ci, result[j]) = (t128 >> 64,t128 as u64);
                                                                }
                                        ai = ci;
                                        m  = (result[0] as u128 * _params.qprime) & ((1 << 64) - 1);
                                        ci = (result[0] as u128 + m * _params.modulo[0] as u128) >> 64;
                                        for j in 1..N {   t128 = result[j] as u128 + (m * _params.modulo[j] as u128) + ci;
                                                                (ci, result[j-1]) = (t128 >> 64 , t128 as u64);
                                      }   
                result[N-1] = (ci + ai) as u64;
                }
            let mut bigger:i32=0;
            for i in 0..N {  bigger += (((result[i] > _params.modulo[i]) as i32) - ((result[i] < _params.modulo[i]) as i32)) << i};
                                    if bigger >= 0 {    let mut borrow:bool=false; 
                                                        for i in 0..N {(result[i],borrow) =  
                                                                              (result[i]).overflowing_sub((&_params.modulo[i]) + (borrow as u64))}
                                    }
            result
            }
    } 

pub fn invert <const N: usize> (a : &[u64; N],_params :&super::prime_fields::FieldParams<N>) -> [u64; N]
        // Modular ivnertion of a modulo p. Implementation from "Efficient Software-Implementation of    
        // Finite Fields with Applications to Cryptography" - Algorithm 16. https://link.springer.com/article/10.1007/s10440-006-9046-1        
    {   #[inline(always)]
        fn cmp <const N: usize> (a : &[u64; N], b : &[u64; N])-> i8 
                // Comparing tow integers based on [u64; N] representation. Proposed constant-time implementation
            {   let mut result:i16 = 0;
                for (index, (ai,bi)) in (a.iter().zip(b.iter())).enumerate() {result += (((ai > bi) as i16) - ((ai < bi) as i16)) << index;}
                result.signum() as i8
            }       
        #[inline(always)]   
        fn add_p <const N: usize> (result:&mut [u64; N], p: &[u64; N], last_carry:&mut u64 ) 
                // In-site Addition with modulo (result =result + modulo)
            {   let mut carry :bool = false;
                let mut carry1 :bool ;
                let mut tmp :u64;
                for i in 0..N { (tmp,carry1) = p[i].overflowing_add(carry as u64);
                                       (result[i],carry) = (result[i]).overflowing_add(tmp);
                                       carry = carry | carry1};               
                *last_carry = carry as u64;
            }        
       #[inline(always)]
       fn r_shift <const N: usize> (result:&mut [u64; N], carry :u64)
                // Right shift on bit string 
            {   for i in 0..N-1{ result[i] = (result[i] >> 1) | ((result[i+1] & 1) << 63)}
                if carry==0 {result[N - 1] = result [N - 1] >> 1}
                else {  result[N - 1] = result [N - 1] >> 1 | ((carry & 1) << 63)}
            } 
            
        fn sub_noreduction <const N: usize> (a : &[u64; N],result:&mut [u64; N]) -> bool 
                // In-site Substraction without reduction (result =result - a)
            {   let mut borrow: bool =false;
                let mut carry: bool ;
                let mut tmp: u64 ;
                for i in 0..N {(tmp,carry) = (a[i]).overflowing_add(borrow as u64);
                                      (result[i],borrow) = (result[i]).overflowing_sub(tmp);
                                      borrow= borrow | carry   };
                borrow 
            }
        let mut u = a.clone();
        let mut v = _params.modulo;
        let mut b = _params.rsquare; 
        let mut c = _params.zero;
        let carry :&mut u64=&mut 0;
        loop {  while u[0] & 1 == 0 {   r_shift(&mut u,0);   
                                        *carry =0;                                 
                                        if b[0] & 1 == 1 {add_p(&mut b,&_params.modulo,carry)};  
                                        r_shift(&mut b,*carry)  }
                while v[0] & 1 == 0 {   r_shift(&mut v,0);
                                        *carry =0;
                                        if c[0] & 1 == 1 {add_p(&mut c,&_params.modulo,carry)};  
                                        r_shift(&mut c,*carry) 
                                     }  
                                                                        
                if cmp(&u , &v) == 1 { sub_noreduction(&v,&mut u);
                                             if sub_noreduction(&c,&mut b) {add_p(&mut b,&_params.modulo,carry);}
                                            }
                else {  sub_noreduction(&u,&mut v);
                        if sub_noreduction(&b,&mut c) {add_p(&mut c,&_params.modulo,carry)}
                    }
                let mut _null_u =0;                        
                let mut _null_v =0;                        
                for i in 1..N { _null_u = _null_u | u[i];
                                       _null_v = _null_v | v[i]}
                if (u[0]==1) & (_null_u == 0) {return b }
                if (v[0]==1) & (_null_v == 0) {return c }
                
            }
    }

pub fn pow <const N: usize> (a : &[u64; N], e :&[u64], useladder:bool,_params :&super::prime_fields::FieldParams<N>)-> [u64; N]
        // Implements montgomery Ladder (secure but relatively slow with respect to square and multiply)
    {   if useladder {  let mut r0 = _params.one;     
                        let mut r1 = a.clone();                            
                        for i in e.as_ref().iter().rev() { 
                                for j in (0..64).rev()  {   if (i >> j) & 1 == 0 {  r1 = mul(&r0, &r1, _params);
                                                                                         r0 = sqr(&r0, _params); }                
                                                            else {  r0 = mul(&r0, &r1,_params);
                                                                    r1 = sqr(&r1,_params)}
                                                        }
                                                        }
                        r0
                        }
        else {  let mut result = _params.one;     
                for i in e.as_ref().iter().rev()  {
                            for j in (0..64).rev(){   result = sqr(&result,_params);                
                                                        if (i >> j) & 1 == 1 {result = mul(&result, a,_params);}                
                                                    }
                                                    }
                result
            }
    }

pub fn sqrt<const N: usize> (a:&[u64; N],_params :&super::prime_fields::FieldParams<N>) -> Option<[u64; N]>
    {   if _params.sqrtid == -1 { // p mod 4 ==3
                                    let result = pow(a,&_params.modplus1div4,false,_params); 
                                    if equal(&mul(&sqr(&result,_params),&invert(a,_params),_params),&_params.one) { Some(result) }
                                    else { None }
                                }
        else {      // Tonelli-Shanks Algorithm (p mod 4 ==1)
                    //  Adapted from  section 6 of 'Square roots from 1; 24, 51, 10 to Dan Shanks'  by Ezra Brown:
	                // https://www.maa.org/sites/default/files/pdf/upload_library/22/Polya/07468342.di020786.02p0470a.pdf 
                    let w = pow(a,&_params.tonelli_params[1],false,_params);                    
                    let mut y = mul(&a,&w,_params);
                    let mut b = mul(&y,&w,_params);
                    let mut g = mul(&_params.tonelli_params[0],&_params.rsquare,_params); // Convert to Montgomery repr
                    let mut r = _params.sqrtid;
                    let mut t = b.clone();
                    for _ in 0..r-1 { t = sqr(&t, _params)}; 
                    if equal(&t, &[0;N]) {return Some([0;N])}
                    if !equal(&t, &_params.one) {return None}       
                    loop  { let mut t= b.clone();
                            let mut m= 0;
                            while !equal(&t, &_params.one) {  t = sqr(&t, _params);
                                                                   m = m + 1;
                                                                }
                            if m==0 { return Some(y)}
                            let mut ge = r - m - 1;
                            t = g.clone();                            
                            while ge > 0 { t = sqr(&t,_params);
                                           ge = ge - 1}
                            g = sqr(&t,_params);
                            y = mul(&y, &t, _params);
                            b = mul(&b, &g, _params);   
                            r = m;
                        }
                }
    }

