// Code developed by FARAOUN Kamel Mohamed.
// EEDIS-Laboratory. UDL-University. Algeria
// During May 2024.

use std::ops::Add;

use crate::{curves::{curve_arithmetics::EcPoint, g2::G2Element}, 
            extensions::{ext_fields::ExtElement, fp2::Fp2Element, fp4::Fp4Element, fp8::Fp8Element, g2_extfields::ExtFieldG2Element}, 
            fields::prime_fields::FieldElement, 
            tools::arithmetic_interface::ArithmeticOperations};

pub fn phi_bls12<const PRAMASIZE:usize,const R: usize, const N: usize, const MAX_COEFS_COUNT: usize>
            (input :&G2Element<PRAMASIZE, R,N,MAX_COEFS_COUNT>) -> G2Element<PRAMASIZE, R,N,MAX_COEFS_COUNT>
    {   
        // Endomorphisme Phi(P)=u*P: Twist->Frobinus->Un-twist (in M-Type Mode for BLS12)
        // Twist : (x,y)->(x/w^2,y/w^3)=(x/(1+u)^(1/3), y/(1+u)^(1/2))
        // Frobinus :(x,y)->(x^P,y^P)
        // Un-Twist : (x,y)->(x*w^2,y*w^3)=(x*(1+u)^(1/3), y*(1+u)^(1/2))
        // Combinaison of the three maps is equivalent to one multiplication by constants : 
        //           (c1,c2)=(1/((1+u)^((prime-1) div 3)),1/((1+u)^((prime-1) div 2)))             
        let frobs = &input.consts.extfieldparams.frobinus_consts;
        let one = FieldElement{ fieldparams: frobs[0].fieldparams, mont_limbs:frobs[0].fieldparams.one };
        let zero = FieldElement{ fieldparams: frobs[0].fieldparams, mont_limbs:frobs[0].fieldparams.zero };
        let phiconst1 = Fp2Element {content :[ zero ,frobs[0].addto(&one)]};
        let phiconst2 = Fp2Element {content :[ frobs[2].negate() ,frobs[2]]};
        G2Element{  consts: input.consts,
                    point: EcPoint{ x: input.point.x.conjugate().multiply(&ExtFieldG2Element::Fp2(phiconst1)),
                                    y: input.point.y.conjugate().multiply(&ExtFieldG2Element::Fp2(phiconst2)),
                                    z: input.point.z.conjugate(),
                                },
                 }
    }
          

pub fn phi_bls24<const PRAMASIZE:usize,const R: usize, const N: usize, const MAX_COEFS_COUNT: usize>
    (input :&G2Element<PRAMASIZE, R,N,MAX_COEFS_COUNT>, order :i8) -> G2Element<PRAMASIZE, R,N,MAX_COEFS_COUNT>
    {
        // Endomorphisme Phi^n((P))))=(u^n)*P (for n =1,4 and -1): Twist->Frobinus->Un-twist (in D-Type Mode for BLS24)
        // Twist : (x,y)->(x*z^2,y*z^3)
        // Frobinus :(x,y)->(x^P,y^P)
        // Un-Twist : (x,y)->(x/z^2,y/z^3)
        // Combinaison of the three maps is equivalent to one multiplication by constants: (c1,c2)=(z^(2*(prime-1))),1/(z^(3*(prime-1))))
        let frobs = &input.consts.extfieldparams.frobinus_consts;
        let x = if let ExtFieldG2Element::Fp4(_x) = input.point.x { _x.content } 
                                      else { panic!("error handeling unsupported type ...") };
        let y = if let ExtFieldG2Element::Fp4(_y) = input.point.y { _y.content } 
                                      else { panic!("error handeling unsupported type ...") };
        let z = if let ExtFieldG2Element::Fp4(_z) = input.point.z { _z} 
                                          else { panic!("error handeling unsupported type ...") };
        let inconsts = if let ExtFieldG2Element::Fp4(_y) = input.point.y { _y.constants } 
                                                     else { panic!("error handeling unsupported type ...") };
        match order {   1 =>{   let phix = [x[0].substract(&x[1]).multiply(&frobs[7]),x[0].addto(&x[1]).multiply(&frobs[7]).negate(),
                                                                  x[2].multiply(&frobs[5].add(1u64)), x[3].multiply(&frobs[5].add(1u64)).negate()];
                                let phiy = [y[2].multiply(&frobs[1]).negate(),y[3].multiply(&frobs[1]), 
                                                                  y[1].multiply(&frobs[2]),y[0].multiply(&frobs[2])];
                                G2Element { consts: input.consts , 
                                point: EcPoint{ x: ExtFieldG2Element::Fp4(Fp4Element { content: phix, constants: inconsts }), 
                                    y: ExtFieldG2Element::Fp4(Fp4Element { content: phiy, constants: inconsts }), 
                                    z: ExtFieldG2Element::Fp4(z.frobinus()) }
                                }
                            } 
                        4 => {  let fb5plus1 = frobs[5].add(1u64);          
                                let phi4x = [x[0].negate().multiply(&fb5plus1),x[1].negate().multiply(&fb5plus1),
                                                                   x[2].negate().multiply(&fb5plus1), x[3].negate().multiply(&fb5plus1)];
                                G2Element { consts: input.consts , 
                                            point: EcPoint{ x: ExtFieldG2Element::Fp4(Fp4Element { content: phi4x, constants: inconsts }), 
                                                            y: ExtFieldG2Element::Fp4(Fp4Element { content: y, constants: inconsts }.negate()), 
                                                            z: ExtFieldG2Element::Fp4(z) }
                                            }
                              }    
                        -1 =>{  let phix = [x[1].substract(&x[0]).multiply(&frobs[6]),x[0].addto(&x[1]).multiply(&frobs[6]),
                                                                  x[2].multiply(&frobs[5]).negate(), x[3].multiply(&frobs[5])];
                                let phiy = [y[3].multiply(&frobs[1]), y[2].multiply(&frobs[1]), 
                                                                  y[0].multiply(&frobs[2]).negate(),y[1].multiply(&frobs[2])];
                                G2Element { consts: input.consts , 
                                    point: EcPoint{ x: ExtFieldG2Element::Fp4(Fp4Element { content: phix, constants: inconsts }), 
                                                    y: ExtFieldG2Element::Fp4(Fp4Element { content: phiy, constants: inconsts }), 
                                                    z: ExtFieldG2Element::Fp4(z.frobinus().conjugate()) //invert of the frobinus is its conjugate//
                                                }
                                    }
                             }
                        _ => {panic!("non-implemented order for EFp4 phi ...")}
                        }
                    }    
                

pub fn phi_bls48<const PRAMASIZE:usize,const R: usize, const N: usize, const MAX_COEFS_COUNT: usize>
            (input :&G2Element<PRAMASIZE, R,N,MAX_COEFS_COUNT> , order: i8) -> G2Element<PRAMASIZE, R,N,MAX_COEFS_COUNT>
{
    // Endomorphisme Phi^n((P))))=(u^n)*P (for n =1,4,8 and -1): Twist->Frobinus->Un-twist (in D-Type Mode for BLS48)   
    let frobs = &input.consts.extfieldparams.frobinus_consts;
    let x = if let ExtFieldG2Element::Fp8(_x) = input.point.x { _x.content } 
                                  else { panic!("error handeling unsupported type ...") };
    let y = if let ExtFieldG2Element::Fp8(_y) = input.point.y { _y.content } 
                                  else { panic!("error handeling unsupported type ...") };
    let z = if let ExtFieldG2Element::Fp8(_z) = input.point.z { _z.content} 
                                  else { panic!("error handeling unsupported type ...") };
    let inconsts = if let ExtFieldG2Element::Fp8(_y) = input.point.y { _y.constants } 
                                                else { panic!("error handeling unsupported type ...") };
    match order {   1 => {  let phix = [  x[3].double().multiply(&frobs[9]),x[2].double().multiply(&frobs[9]),
                                                                x[0].multiply(&frobs[8]).negate(),x[1].multiply(&frobs[8]),
                                                                x[4].substract(&x[5]).multiply(&frobs[7]),x[4].addto(&x[5]).multiply(&frobs[7]).negate(),
                                                                x[6].multiply(&frobs[5].add(1u64)),x[7].multiply(&frobs[5].add(1u64)).negate()];
                            let phiy = [  y[6].substract(&y[7]).multiply(&frobs[17]).negate(),y[6].addto(&y[7]).multiply(&frobs[17]), 
                                                                y[5].addto(&y[4]).multiply(&frobs[16]).negate(),y[5].substract(&y[4]).multiply(&frobs[16]), 
                                                                y[2].multiply(&frobs[14]).negate(),y[3].multiply(&frobs[14]), 
                                                                y[1].multiply(&frobs[15]).negate(),y[0].multiply(&frobs[15]).negate()];
                            G2Element { consts: input.consts , 
                                        point: EcPoint{ x: ExtFieldG2Element::Fp8(Fp8Element { content: phix, constants: inconsts }), 
                                                        y: ExtFieldG2Element::Fp8(Fp8Element { content: phiy, constants: inconsts }), 
                                                        z: ExtFieldG2Element::Fp8(Fp8Element { content: z, constants: inconsts }.frobinus()) }
                                        }
                        }        
                    4 => {  let fb5plus1 = frobs[5].add(1u64);          
                            let phi4x = [ x[0].multiply(&fb5plus1),x[1].multiply(&fb5plus1),
                                                                x[2].multiply(&fb5plus1),x[3].multiply(&fb5plus1),
                                                                x[4].multiply(&fb5plus1).negate(),x[5].multiply(&fb5plus1).negate(),
                                                                x[6].multiply(&fb5plus1).negate(),x[7].multiply(&fb5plus1).negate()];
                            let phi4y = [ y[1],y[0].negate(),y[3],y[2].negate(),
                                                                y[5].negate(),y[4],y[7].negate(),y[6]]; // y' =-y.conjugate * u                                                                 
                            G2Element { consts: input.consts , 
                                        point: EcPoint{ x: ExtFieldG2Element::Fp8(Fp8Element { content: phi4x, constants: inconsts }), 
                                                        y: ExtFieldG2Element::Fp8(Fp8Element { content: phi4y, constants: inconsts }), 
                                                        z: ExtFieldG2Element::Fp8(Fp8Element { content: z, constants: inconsts }.conjugate()) } // frobinus at order 4 for an Fp8 is simply the conjugate
                                        }
                            }
                    8 => {  let phi8x = [ x[0].multiply(&frobs[5]),x[1].multiply(&frobs[5]),
                                                                x[2].multiply(&frobs[5]),x[3].multiply(&frobs[5]),
                                                                x[4].multiply(&frobs[5]),x[5].multiply(&frobs[5]),
                                                                x[6].multiply(&frobs[5]),x[7].multiply(&frobs[5])];
                            G2Element { consts: input.consts , 
                                        point: EcPoint{ x: ExtFieldG2Element::Fp8(Fp8Element { content: phi8x, constants: inconsts }), 
                                                        y: ExtFieldG2Element::Fp8(Fp8Element { content: y, constants: inconsts }.negate()), 
                                                        z: ExtFieldG2Element::Fp8(Fp8Element { content: z, constants: inconsts }) } // frobinus at order 4 for an Fp8 is simply the identity
                                        }
                           }
                   -1 => {  let invphix = [x[2].double().multiply(&frobs[4]).negate(),x[3].multiply(&frobs[4]).double(),
                                                                x[1].multiply(&frobs[3]).negate(),x[0].multiply(&frobs[3]).negate(),
                                                                x[5].substract(&x[4]).multiply(&frobs[6]),x[5].addto(&x[4]).multiply(&frobs[6]),
                                                                x[6].multiply(&frobs[5]).negate(),x[7].multiply(&frobs[5])];
                            let invphiy = [y[7].multiply(&frobs[22]).negate(),y[6].multiply(&frobs[22]).negate(),
                                                                y[4].multiply(&frobs[17]).negate(),y[5].multiply(&frobs[17]),
                                                                y[3].addto(&y[2]).multiply(&frobs[15]),y[2].substract(&y[3]).multiply(&frobs[15]),
                                                                y[0].substract(&y[1]).multiply(&frobs[23]),y[0].addto(&y[1]).multiply(&frobs[23]).negate()];                                           
                            // Inveting the frobinus for an Fp8                                  
                            let invphiz = [z[0],z[1].negate(),z[2].addto(&z[3]).multiply(&frobs[0]).negate(),
                                                                                   z[3].substract(&z[2]).multiply(&frobs[0]),z[7].multiply(&frobs[1]),
                                                                                   z[6].multiply(&frobs[1]),z[4].multiply(&frobs[2]).negate(),z[5].multiply(&frobs[2])];                                                                                  
                            G2Element { consts: input.consts , 
                            point: EcPoint{ x: ExtFieldG2Element::Fp8(Fp8Element { content: invphix, constants: inconsts }), 
                                        y: ExtFieldG2Element::Fp8(Fp8Element { content: invphiy, constants: inconsts }), 
                                        z: ExtFieldG2Element::Fp8(Fp8Element { content: invphiz, constants: inconsts }) } 
                            }
                         }
                    _ => {panic!("non-implemented order for EFp8 phi ...")}
    }
}