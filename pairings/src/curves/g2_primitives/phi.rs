// Code developed by FARAOUN Kamel Mohamed.
// EEDIS-Laboratory. UDL-University. Algeria
// During May 2024.

use std::ops::{Add, Sub};

use crate::{curves::{curve_arithmetics::EcPoint, g2::G2Element}, 
            extensions::ext_fields::ExtElement, 
            extensions::towering1::fp2::Fp2Element as Fp2Element_1, 
            extensions::towering1::fp4::Fp4Element as Fp4Element_1, 
            extensions::towering1::fp8::Fp8Element as Fp8Element_1, 
            extensions::towering2::fp4::Fp4Element as Fp4Element_2, 
            extensions::towering2::fp8::Fp8Element as Fp8Element_2, 
            extensions::towering3::fp8::Fp8Element as Fp8Element_3, 
            extensions::g2_extfields::ExtFieldG2Element, 
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
        let phiconst1 = Fp2Element_1 {content :[ zero ,frobs[0].addto(&one)]};
        let phiconst2 = Fp2Element_1 {content :[ frobs[2].negate() ,frobs[2]]};
        G2Element{  consts: input.consts,
                    point: EcPoint{ x: input.point.x.conjugate().multiply(&ExtFieldG2Element::Fp2_1(phiconst1)),
                                    y: input.point.y.conjugate().multiply(&ExtFieldG2Element::Fp2_1(phiconst2)),
                                    z: input.point.z.conjugate(),
                                },
                 }
    }
          

pub fn phi_bls24<const PRAMASIZE:usize,const R: usize, const N: usize, const MAX_COEFS_COUNT: usize>
    (input :&G2Element<PRAMASIZE, R,N,MAX_COEFS_COUNT>, order :i8) -> G2Element<PRAMASIZE, R,N,MAX_COEFS_COUNT>
    {
        // Endomorphisme Phi^n((P))))=(u^n)*P (for n =1,4 and -1): Twist->Frobinus->Un-twist (in D-Type/M-Type Mode for BLS24-construction-1 and M-type for BLS24-construction-2 ) 
        // Example M-Type
        // Twist : (x,y)->(x*z^2,y*z^3)
        // Frobinus :(x,y)->(x^P,y^P)
        // Un-Twist : (x,y)->(x/z^2,y/z^3)
        // Combinaison of the three maps is equivalent to one multiplication by constants: (c1,c2)=(z^(2*(prime-1))),1/(z^(3*(prime-1))))
        let frobs = &input.consts.extfieldparams.frobinus_consts;
        let x = if let ExtFieldG2Element::Fp4_1(_x) = input.point.x { _x.content } 
                                      else {    if let ExtFieldG2Element::Fp4_2(_x) = input.point.x { _x.content } 
                                                else { panic!("error handeling unsupported type ...") }};
        let y = if let ExtFieldG2Element::Fp4_1(_y) = input.point.y { _y.content } 
                                      else {    if let ExtFieldG2Element::Fp4_2(_y) = input.point.y { _y.content } 
                                                else { panic!("error handeling unsupported type ...") }};
        let z = if let ExtFieldG2Element::Fp4_1(_z) = input.point.z { _z.content} 
                                      else {     if let ExtFieldG2Element::Fp4_2(_z) = input.point.z { _z.content} 
                                                 else { panic!("error handeling unsupported type ...") } };
        let inconsts = if let ExtFieldG2Element::Fp4_1(_y) = input.point.y { _y.constants } 
                                                     else {     if let ExtFieldG2Element::Fp4_2(_y) = input.point.y { _y.constants } 
                                                                else { panic!("error handeling unsupported type ...") } };
        let construction;
        match input.point.x {
            ExtFieldG2Element::Fp4_1(_) => {construction =1},
            ExtFieldG2Element::Fp4_2(_) => {construction =2},
            _ => {unimplemented!("error handeling unsupported type ...")},
        }                                                                 
        match order {   1 =>{   if construction ==1 {   let phix :[FieldElement<N>; 4];
                                                        let phiy :[FieldElement<N>; 4];
                                                        if input.consts.twist_type=='D'{  // D-Type Twiste
                                                                    phix = [x[0].substract(&x[1]).multiply(&frobs[7]),x[0].addto(&x[1]).multiply(&frobs[7]).negate(),
                                                                                                      x[2].multiply(&frobs[5].add(1u64)), x[3].multiply(&frobs[5].add(1u64)).negate()];
                                                                    phiy = [y[2].multiply(&frobs[1]).negate(),y[3].multiply(&frobs[1]), 
                                                                                                      y[1].multiply(&frobs[2]),y[0].multiply(&frobs[2])];
                                                                        }
                                                        else { // M-Type Twiste
                                                               phix = [x[0].addto(&x[1]).multiply(&frobs[6]).negate(),x[0].substract(&x[1]).multiply(&frobs[6]).negate(),
                                                                                                x[3].multiply(&frobs[5]), x[2].multiply(&frobs[5])];
                                                               phiy = [y[3].substract(&y[2]).multiply(&frobs[2]),y[3].addto(&y[2]).multiply(&frobs[2]), 
                                                                        y[0].addto(&y[1]).multiply(&frobs[10]),y[0].substract(&y[1]).multiply(&frobs[10])];
                                                                        }        
                                                        G2Element { consts: input.consts , 
                                                                    point: EcPoint{ x: ExtFieldG2Element::Fp4_1(Fp4Element_1 { content: phix, constants: inconsts }), 
                                                                                    y: ExtFieldG2Element::Fp4_1(Fp4Element_1 { content: phiy, constants: inconsts }), 
                                                                                    z: ExtFieldG2Element::Fp4_1(Fp4Element_1 { content: z, constants: inconsts }.frobinus()) }
                                                                 }}
                                else {  let phix = [x[0].multiply(&frobs[7]) ,x[1].negate().multiply(&frobs[7]),x[2].multiply(&frobs[8]),
                                                                            x[3].negate().multiply(&frobs[8])];
                                        let phiy = [y[0].multiply(&frobs[1]) ,y[1].negate().multiply(&frobs[1]),y[2].multiply(&frobs[2]),
                                                                            y[3].negate().multiply(&frobs[2])];                                   
                                        G2Element { consts: input.consts , 
                                                    point: EcPoint{ x: ExtFieldG2Element::Fp4_2(Fp4Element_2 { content: phix, constants: inconsts }), 
                                                                    y: ExtFieldG2Element::Fp4_2(Fp4Element_2 { content: phiy, constants: inconsts }), 
                                                                    z: ExtFieldG2Element::Fp4_2(Fp4Element_2 { content: z, constants: inconsts }.frobinus()) }
                                                }}
                                }
                             
                        4 => {  if construction ==1 {   let phi4x :[FieldElement<N>; 4];
                                                        if input.consts.twist_type=='D'{  // D-Type Twiste
                                                                        let fb5plus1 = frobs[5].add(1u64);          
                                                                        phi4x = [x[0].negate().multiply(&fb5plus1),x[1].negate().multiply(&fb5plus1),
                                                                                 x[2].negate().multiply(&fb5plus1), x[3].negate().multiply(&fb5plus1)];
                                                                        }
                                                        else {  // M-Type Twiste                                                     
                                                                phi4x = [x[0].multiply(&frobs[5]),x[1].multiply(&frobs[5]),
                                                                         x[2].multiply(&frobs[5]), x[3].multiply(&frobs[5])];
                                                                }                                                                        
                                                        G2Element { consts: input.consts , 
                                                                  point: EcPoint{ x: ExtFieldG2Element::Fp4_1(Fp4Element_1 { content: phi4x, constants: inconsts }), 
                                                                                  y: ExtFieldG2Element::Fp4_1(Fp4Element_1 { content: y, constants: inconsts }.negate()), 
                                                                                  z: ExtFieldG2Element::Fp4_1(Fp4Element_1 { content: z, constants: inconsts }) }
                                                                    }
                                                    }
                                else {  let phi4x = [x[0].multiply(&frobs[8]) ,x[1].multiply(&frobs[8]),x[2].multiply(&frobs[8]),
                                                                            x[3].multiply(&frobs[8])];
                                        G2Element { consts: input.consts , 
                                                    point: EcPoint{ x: ExtFieldG2Element::Fp4_2(Fp4Element_2 { content: phi4x, constants: inconsts }), 
                                                                    y: ExtFieldG2Element::Fp4_2(Fp4Element_2 { content: y, constants: inconsts }.negate()), 
                                                                    z: ExtFieldG2Element::Fp4_2(Fp4Element_2 { content: z, constants: inconsts }) }
                                                    }
                                }                                    
                              }    
                        -1 =>{  if construction ==1 {   let phix :[FieldElement<N>; 4];
                                                        let phiy :[FieldElement<N>; 4];
                                                        if input.consts.twist_type=='D'{  // D-Type Twiste
                                                                                        phix = [x[1].substract(&x[0]).multiply(&frobs[6]),x[0].addto(&x[1]).multiply(&frobs[6]),
                                                                                                                    x[2].multiply(&frobs[5]).negate(), x[3].multiply(&frobs[5])];
                                                                                        phiy = [y[3].multiply(&frobs[1]), y[2].multiply(&frobs[1]), 
                                                                                                                    y[0].multiply(&frobs[2]).negate(),y[1].multiply(&frobs[2])];
                                                                                        }
                                                        else {// M-Type Twiste                                                              
                                                                let fb5plus1 = frobs[5].add(1u64);          
                                                                phix = [x[1].addto(&x[0]).multiply(&frobs[7]),x[0].substract(&x[1]).multiply(&frobs[7]),
                                                                                            x[3].multiply(&fb5plus1).negate(), x[2].multiply(&fb5plus1).negate()];
                                                                phiy = [y[2].addto(&y[3]).multiply(&frobs[2]).negate(), y[3].substract(&y[2]).multiply(&frobs[2]), 
                                                                        y[0].substract(&y[1]).multiply(&frobs[10]),y[0].addto(&y[1]).multiply(&frobs[10]).negate()];
                                                                }                                                                                        
                                                        G2Element { consts: input.consts , 
                                                                    point: EcPoint{ x: ExtFieldG2Element::Fp4_1(Fp4Element_1 { content: phix, constants: inconsts }), 
                                                                                    y: ExtFieldG2Element::Fp4_1(Fp4Element_1 { content: phiy, constants: inconsts }), 
                                                                                    z: ExtFieldG2Element::Fp4_1(Fp4Element_1 { content: z, constants: inconsts }.frobinus().conjugate()) //invert of the frobinus is its conjugate//
                                                                                }
                                                                    }
                                                    }
                                else {  let phix = [x[0].multiply(&frobs[6]).negate(),x[1].multiply(&frobs[6]),
                                                                            x[2].multiply(&frobs[5]).negate(), x[3].multiply(&frobs[5])];
                                        let phiy = [y[0].multiply(&frobs[2]).negate(), y[1].multiply(&frobs[2]), 
                                                                            y[2].multiply(&frobs[1]).negate(),y[3].multiply(&frobs[1])];
                                        G2Element { consts: input.consts , 
                                                    point: EcPoint{ x: ExtFieldG2Element::Fp4_2(Fp4Element_2 { content: phix, constants: inconsts }), 
                                                                    y: ExtFieldG2Element::Fp4_2(Fp4Element_2 { content: phiy, constants: inconsts }), 
                                                                    z: ExtFieldG2Element::Fp4_2(Fp4Element_2 { content: z, constants: inconsts }.frobinus().conjugate()) //invert of the frobinus is its conjugate//
                                                                }
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
    let fb_id= &input.consts.extfieldparams.fb_id;
    let x = if let ExtFieldG2Element::Fp8_1(_x) = input.point.x { _x.content } 
                                  else { if let ExtFieldG2Element::Fp8_2(_x) = input.point.x { _x.content } 
                                         else {if let ExtFieldG2Element::Fp8_3(_x) = input.point.x { _x.content } 
                                                 else {panic!("error handeling unsupported type ...")} } };
    let y = if let ExtFieldG2Element::Fp8_1(_y) = input.point.y { _y.content } 
                                  else { if let ExtFieldG2Element::Fp8_2(_y) = input.point.y { _y.content } 
                                         else { if let ExtFieldG2Element::Fp8_3(_y) = input.point.y { _y.content } 
                                                else {panic!("error handeling unsupported type ...")} } };
    let z = if let ExtFieldG2Element::Fp8_1(_z) = input.point.z { _z.content} 
                                  else { if let ExtFieldG2Element::Fp8_2(_z) = input.point.z { _z.content} 
                                         else { if let ExtFieldG2Element::Fp8_3(_z) = input.point.z { _z.content } 
                                                else {panic!("error handeling unsupported type ...")} }};
    let inconsts = if let ExtFieldG2Element::Fp8_1(_y) = input.point.y { _y.constants } 
                                                 else { if let ExtFieldG2Element::Fp8_2(_y) = input.point.y { _y.constants } 
                                                        else { if let ExtFieldG2Element::Fp8_3(_y) = input.point.y { _y.constants } 
                                                               else {panic!("error handeling unsupported type ...")} } };
    let construction;
    match input.point.x {
        ExtFieldG2Element::Fp8_1(_) => {construction =1},
        ExtFieldG2Element::Fp8_2(_) => {construction =2},
        ExtFieldG2Element::Fp8_3(_) => {construction =3},
        _ => {unimplemented!("error handeling unsupported type ...")},
    }                                                        
    match order {   1 => {  if construction ==1 {   let phix = [  x[3].double().multiply(&frobs[9]),x[2].double().multiply(&frobs[9]),
                                                                                    x[0].multiply(&frobs[8]).negate(),x[1].multiply(&frobs[8]),
                                                                                    x[4].substract(&x[5]).multiply(&frobs[7]),x[4].addto(&x[5]).multiply(&frobs[7]).negate(),
                                                                                    x[6].multiply(&frobs[5].add(1u64)),x[7].multiply(&frobs[5].add(1u64)).negate()];
                                                    let phiy = [  y[6].substract(&y[7]).multiply(&frobs[17]).negate(),y[6].addto(&y[7]).multiply(&frobs[17]), 
                                                                                        y[5].addto(&y[4]).multiply(&frobs[16]).negate(),y[5].substract(&y[4]).multiply(&frobs[16]), 
                                                                                        y[2].multiply(&frobs[14]).negate(),y[3].multiply(&frobs[14]), 
                                                                                        y[1].multiply(&frobs[15]).negate(),y[0].multiply(&frobs[15]).negate()];
                                                    G2Element { consts: input.consts , 
                                                                    point: EcPoint{ x: ExtFieldG2Element::Fp8_1(Fp8Element_1 { content: phix, constants: inconsts }), 
                                                                                    y: ExtFieldG2Element::Fp8_1(Fp8Element_1 { content: phiy, constants: inconsts }), 
                                                                                    z: ExtFieldG2Element::Fp8_1(Fp8Element_1 { content: z, constants: inconsts }.frobinus()) }
                                                                    }
                                                }
                            else {  if construction ==3 {   let phix: [FieldElement<N>; 8] ;
                                                            let phiy: [FieldElement<N>; 8] ;
                                                            if input.consts.twist_type=='D'
                                                                    {// D-Type (with negation) twist 
                                                                    phix = [  x[2].addto(&x[3]).multiply(&frobs[3]).negate(),x[3].substract(&x[2]).multiply(&frobs[3]),
                                                                            x[1].substract(&x[0]).multiply(&frobs[4]),x[0].addto(&x[1]).multiply(&frobs[4]),
                                                                            x[5].multiply(&frobs[5]).negate(),x[4].multiply(&frobs[5]).negate(),
                                                                            x[6].substract(&x[7]).multiply(&frobs[6]),x[7].addto(&x[6]).multiply(&frobs[6]).negate()];
                                                                    phiy = [  y[4].multiply(&frobs[14]).negate(),y[5].multiply(&frobs[14]), 
                                                                            y[6].addto(&y[7]).multiply(&frobs[15]).negate(),y[7].substract(&y[6]).multiply(&frobs[15]), 
                                                                            y[0].addto(&y[1]).multiply(&frobs[16]).negate(),y[1].substract(&y[0]).multiply(&frobs[16]), 
                                                                            y[3].multiply(&frobs[17]),y[2].multiply(&frobs[17])];                                                                                                                         
                                                                    }
                                                                    else {// M-Type  twist 
                                                            // return toFp8([-2*fp48FrobConsts[12]*tmp[4],-2*fp48FrobConsts[12]*tmp[3],-fp48FrobConsts[11]*tmp[1],fp48FrobConsts[11]*tmp[2],
                                                            //               -fp48FrobConsts[10]*(tmp[6]-tmp[5]),-fp48FrobConsts[10]*(tmp[5]+tmp[6]),
                                                            //              -(fp48FrobConsts[8]-1)*(tmp[7]),(fp48FrobConsts[8]-1)*(tmp[8])]);
                                                                            phix = [  x[3].multiply(&frobs[9]).negate().double(),x[2].multiply(&frobs[9]).negate().double(),
                                                                                    x[0].multiply(&frobs[8]).negate(),x[1].multiply(&frobs[8]),
                                                                                    x[5].substract(&x[4]).multiply(&frobs[7]).negate(),x[4].addto(&x[5]).multiply(&frobs[7]).negate(),
                                                                                    x[6].multiply(&frobs[5].sub(1u64)).negate(),x[7].multiply(&frobs[5].sub(1u64))];
                                                                        // return toFp8([fp48FrobConsts[20]*(tmp[7]-tmp[8]),-fp48FrobConsts[20]*(tmp[7]+tmp[8]), -fp48FrobConsts[19]*(tmp[6]+tmp[5]),fp48FrobConsts[19]*(tmp[6]-tmp[5]),
                                                                        // -fp48FrobConsts[17]*(tmp[3]),fp48FrobConsts[17]*(tmp[4]),fp48FrobConsts[18]*tmp[2],fp48FrobConsts[18]*tmp[1]]);                                                                                    
                                                                            phiy = [  y[6].substract(&y[7]).multiply(&frobs[17]),y[6].addto(&y[7]).multiply(&frobs[17]).negate(), 
                                                                                      y[4].addto(&y[5]).multiply(&frobs[16]).negate(),y[5].substract(&y[4]).multiply(&frobs[16]), 
                                                                                    y[2].multiply(&frobs[14]).negate(),y[3].multiply(&frobs[14]), 
                                                                                    y[1].multiply(&frobs[15]),y[0].multiply(&frobs[15])];                                                                                                                         
                                                                        }
                                                                    G2Element { consts: input.consts , 
                                                                        point: EcPoint{ x: ExtFieldG2Element::Fp8_3(Fp8Element_3 { content: phix, constants: inconsts }), 
                                                                                        y: ExtFieldG2Element::Fp8_3(Fp8Element_3 { content: phiy, constants: inconsts }), 
                                                                                        z: ExtFieldG2Element::Fp8_3(Fp8Element_3 { content: z, constants: inconsts }.frobinus()) }
                                                                }
                                                    }
                                    else {  // construction ==2
                                            let phix: [FieldElement<N>; 8] ;
                                            let phiy: [FieldElement<N>; 8] ;
                                            if *fb_id == 1{ phix = [  x[1].multiply(&frobs[6]).negate(),x[0].multiply(&frobs[5]),
                                                                    x[3].multiply(&frobs[8]).negate(),x[2].multiply(&frobs[7]),
                                                                    x[4].multiply(&frobs[10]),x[5].multiply(&frobs[10]).negate(),
                                                                    x[6].multiply(&frobs[11]),x[7].multiply(&frobs[11]).negate()];
                                                            phiy = [  y[2].multiply(&frobs[24]),y[3].multiply(&frobs[24]).negate(), 
                                                                      y[1].multiply(&frobs[25]),y[0].multiply(&frobs[26]), 
                                                                      y[7].multiply(&frobs[27]),y[6].multiply(&frobs[28]), 
                                                                      y[4].multiply(&frobs[29]),y[5].multiply(&frobs[29]).negate()];}
                                            else {  // Phi variant for the BLS48_287  
                                                            phix = [  x[0].multiply(&frobs[3]),x[1].negate().multiply(&frobs[3]),
                                                                      x[2].multiply(&frobs[4]),x[3].negate().multiply(&frobs[4]),
                                                                      x[4].multiply(&frobs[5]),x[5].negate().multiply(&frobs[5]),
                                                                      x[6].multiply(&frobs[6]),x[7].negate().multiply(&frobs[6])];
                                                            phiy = [  y[1].multiply(&frobs[20]),y[0].multiply(&frobs[19]), 
                                                                      y[3].multiply(&frobs[22]),y[2].multiply(&frobs[21]), 
                                                                      y[5].multiply(&frobs[24]),y[4].multiply(&frobs[23]), 
                                                                      y[7].multiply(&frobs[26]),y[6].multiply(&frobs[25])];}                                                                                                                         
                                            G2Element { consts: input.consts , 
                                                        point: EcPoint{ x: ExtFieldG2Element::Fp8_2(Fp8Element_2 { content: phix, constants: inconsts }), 
                                                                        y: ExtFieldG2Element::Fp8_2(Fp8Element_2 { content: phiy, constants: inconsts }), 
                                                                        z: ExtFieldG2Element::Fp8_2(Fp8Element_2 { content: z, constants: inconsts }.frobinus()) }
                                                        }
                                        }
                            }
                        }        
                    4 => {  if construction ==1 {   let fb5plus1 = frobs[5].add(1u64);          
                                                    let phi4x = [ x[0].multiply(&fb5plus1),x[1].multiply(&fb5plus1),
                                                                                        x[2].multiply(&fb5plus1),x[3].multiply(&fb5plus1),
                                                                                        x[4].multiply(&fb5plus1).negate(),x[5].multiply(&fb5plus1).negate(),
                                                                                        x[6].multiply(&fb5plus1).negate(),x[7].multiply(&fb5plus1).negate()];
                                                    let phi4y = [ y[1],y[0].negate(),y[3],y[2].negate(),
                                                                                        y[5].negate(),y[4],y[7].negate(),y[6]]; // y' =-y.conjugate * u                                                                 
                                                    G2Element { consts: input.consts , 
                                                                    point: EcPoint{ x: ExtFieldG2Element::Fp8_1(Fp8Element_1 { content: phi4x, constants: inconsts }), 
                                                                                    y: ExtFieldG2Element::Fp8_1(Fp8Element_1 { content: phi4y, constants: inconsts }), 
                                                                                    z: ExtFieldG2Element::Fp8_1(Fp8Element_1 { content: z, constants: inconsts }.conjugate()) } // frobinus at order 4 for an Fp8 is simply the conjugate
                                                                    }
                                                }
                            else {if construction ==3 { let phix4: [FieldElement<N>; 8] ;
                                                        let phiy4: [FieldElement<N>; 8] ;
                                                        if input.consts.twist_type=='D'{// D-Type twist 
                                                                                        phix4 = [  x[0].multiply(&frobs[5]),x[1].multiply(&frobs[5]),
                                                                                                x[2].multiply(&frobs[5]),x[3].multiply(&frobs[5]),
                                                                                                x[4].multiply(&frobs[5]).negate(),x[5].multiply(&frobs[5]).negate(),
                                                                                                x[6].multiply(&frobs[5]).negate(),x[7].multiply(&frobs[5]).negate()];
                                                                                        phiy4 = [  y[1].negate(),y[0], y[3].negate(),y[2], 
                                                                                                y[5],y[4].negate(), y[7],y[6].negate()];                                
                                                                                        }
                                                        else {// M-Type twist 
                                                            // return toFp8([-(fp48FrobConsts[8]-1)*tmp[1],-(fp48FrobConsts[8]-1)*tmp[2],-(fp48FrobConsts[8]-1)*tmp[3],-(fp48FrobConsts[8]-1)*tmp[4],
                                                            // (fp48FrobConsts[8]-1)*(tmp[5]),(fp48FrobConsts[8]-1)*(tmp[6]),(fp48FrobConsts[8]-1)*(tmp[7]),(fp48FrobConsts[8]-1)*(tmp[8])]);                                                            
                                                                phix4 = [  x[0].multiply(&frobs[5].sub(1u64)).negate(),x[1].multiply(&frobs[5].sub(1u64)).negate(),
                                                                        x[2].multiply(&frobs[5].sub(1u64)).negate(),x[3].multiply(&frobs[5].sub(1u64)).negate(),
                                                                        x[4].multiply(&frobs[5].sub(1u64)),x[5].multiply(&frobs[5].sub(1u64)),
                                                                        x[6].multiply(&frobs[5].sub(1u64)),x[7].multiply(&frobs[5].sub(1u64))];
                                                                        // return toFp8([tmp[2],-tmp[1],tmp[4],-tmp[3], -tmp[6],tmp[5],-tmp[8],tmp[7]]);
                                                                phiy4 = [  y[1],y[0].negate(), y[3],y[2].negate(), 
                                                                            y[5].negate(),y[4], y[7].negate(),y[6]];                                
                                                            }
                                                    G2Element { consts: input.consts , 
                                                        point: EcPoint{ x: ExtFieldG2Element::Fp8_3(Fp8Element_3 { content: phix4, constants: inconsts }), 
                                                                        y: ExtFieldG2Element::Fp8_3(Fp8Element_3 { content: phiy4, constants: inconsts }), 
                                                                        z: ExtFieldG2Element::Fp8_3(Fp8Element_3 { content: z, constants: inconsts }.conjugate()) }
                                                        }
                                                    }                                
                                else {  //construction ==2
                                    let phix4: [FieldElement<N>; 8] ;
                                    let phiy4: [FieldElement<N>; 8] ;
                                    if *fb_id == 1{ phix4 = [  x[0].multiply(&frobs[10]),x[1].multiply(&frobs[10]),
                                                               x[2].multiply(&frobs[10]),x[3].multiply(&frobs[10]),
                                                               x[4].multiply(&frobs[10]).negate(),x[5].multiply(&frobs[10]).negate(),
                                                               x[6].multiply(&frobs[10]).negate(),x[7].multiply(&frobs[10]).negate()];
                                                    phiy4 = [  y[0].multiply(&frobs[0]).negate(),y[1].multiply(&frobs[0]).negate(), 
                                                               y[2].multiply(&frobs[0]).negate(),y[3].multiply(&frobs[0]).negate(), 
                                                               y[4].multiply(&frobs[0]),y[5].multiply(&frobs[0]), y[6].multiply(&frobs[0]),y[7].multiply(&frobs[0])];
                                                }
                                    else {  // Phi variant for the BLS48_287                           
                                            phix4 = [  x[0].multiply(&frobs[5]),x[1].multiply(&frobs[5]),
                                                       x[2].multiply(&frobs[5]),x[3].multiply(&frobs[5]),
                                                       x[4].multiply(&frobs[5]).negate(),x[5].multiply(&frobs[5]).negate(),
                                                       x[6].multiply(&frobs[5]).negate(),x[7].multiply(&frobs[5]).negate()];
                                            phiy4 = [  y[0].multiply(&frobs[0]),y[1].multiply(&frobs[0]), 
                                                       y[2].multiply(&frobs[0]),y[3].multiply(&frobs[0]), 
                                                       y[4].negate().multiply(&frobs[0]),y[5].negate().multiply(&frobs[0]), 
                                                       y[6].negate().multiply(&frobs[0]),y[7].negate().multiply(&frobs[0])];
                                         }                                                
                                    G2Element { consts: input.consts , 
                                                    point: EcPoint{ x: ExtFieldG2Element::Fp8_2(Fp8Element_2 { content: phix4, constants: inconsts }), 
                                                                    y: ExtFieldG2Element::Fp8_2(Fp8Element_2 { content: phiy4, constants: inconsts }), 
                                                                    z: ExtFieldG2Element::Fp8_2(Fp8Element_2 { content: z, constants: inconsts }.conjugate()) } // frobinus at order 4 for an Fp8 is simply the conjugate
                                            }
                                    }
                                }
                            }
                    8 => {  if construction ==1 {   let phi8x = [ x[0].multiply(&frobs[5]),x[1].multiply(&frobs[5]),
                                                                                            x[2].multiply(&frobs[5]),x[3].multiply(&frobs[5]),
                                                                                            x[4].multiply(&frobs[5]),x[5].multiply(&frobs[5]),
                                                                                            x[6].multiply(&frobs[5]),x[7].multiply(&frobs[5])];
                                                    G2Element { consts: input.consts , 
                                                                    point: EcPoint{ x: ExtFieldG2Element::Fp8_1(Fp8Element_1 { content: phi8x, constants: inconsts }), 
                                                                                    y: ExtFieldG2Element::Fp8_1(Fp8Element_1 { content: y, constants: inconsts }.negate()), 
                                                                                    z: ExtFieldG2Element::Fp8_1(Fp8Element_1 { content: z, constants: inconsts }) } // frobinus at order 4 for an Fp8 is simply the identity
                                                                    }
                                                }
                            else {if construction ==3 { let phix8: [FieldElement<N>; 8] ;
                                                        if input.consts.twist_type=='D'{//D-Type twist 
                                                                                        phix8 = [ x[0].multiply(&frobs[5].sub(1u64)),x[1].multiply(&frobs[5].sub(1u64)),
                                                                                                x[2].multiply(&frobs[5].sub(1u64)),x[3].multiply(&frobs[5].sub(1u64)),
                                                                                                x[4].multiply(&frobs[5].sub(1u64)),x[5].multiply(&frobs[5].sub(1u64)),
                                                                                                x[6].multiply(&frobs[5].sub(1u64)),x[7].multiply(&frobs[5].sub(1u64))];                      
                                                                                        }
                                                        else {//M-Type twist 
                                                            // return toFp8([-fp48FrobConsts[8]*tmp[1],-fp48FrobConsts[8]*tmp[2],-fp48FrobConsts[8]*tmp[3],-fp48FrobConsts[8]*tmp[4],
                                                            //     -fp48FrobConsts[8]*tmp[5],-fp48FrobConsts[8]*tmp[6],-fp48FrobConsts[8]*tmp[7],-fp48FrobConsts[8]*tmp[8]]);
                                                                phix8 = [ x[0].multiply(&frobs[5]).negate(),x[1].multiply(&frobs[5]).negate(),
                                                                            x[2].multiply(&frobs[5]).negate(),x[3].multiply(&frobs[5]).negate(),
                                                                            x[4].multiply(&frobs[5]).negate(),x[5].multiply(&frobs[5]).negate(),
                                                                            x[6].multiply(&frobs[5]).negate(),x[7].multiply(&frobs[5]).negate()];                      
                                                             }                                                                                        
                                                        G2Element { consts: input.consts , 
                                                            point: EcPoint{ x: ExtFieldG2Element::Fp8_3(Fp8Element_3 { content: phix8, constants: inconsts }), 
                                                                            y: ExtFieldG2Element::Fp8_3(Fp8Element_3 { content: y, constants: inconsts }.negate()), 
                                                                            z: ExtFieldG2Element::Fp8_3(Fp8Element_3 { content: z, constants: inconsts }) }
                                                            }
                                                        }
                                else {  // construction ==2  
                                        let phix8: [FieldElement<N>; 8] ;                  
                                        if *fb_id == 1{phix8 = [ x[0].multiply(&frobs[13]),x[1].multiply(&frobs[13]),
                                                                 x[2].multiply(&frobs[13]),x[3].multiply(&frobs[13]),
                                                                 x[4].multiply(&frobs[13]),x[5].multiply(&frobs[13]),
                                                                 x[6].multiply(&frobs[13]),x[7].multiply(&frobs[13])];
                                                      } 
                                        else {  // Phi variant for the BLS48_287           
                                                phix8 = [ x[0].multiply(&frobs[8]),x[1].multiply(&frobs[8]),
                                                          x[2].multiply(&frobs[8]),x[3].multiply(&frobs[8]),
                                                          x[4].multiply(&frobs[8]),x[5].multiply(&frobs[8]),
                                                          x[6].multiply(&frobs[8]),x[7].multiply(&frobs[8])];
                                        } 
                                        G2Element { consts: input.consts , 
                                                    point: EcPoint{ x: ExtFieldG2Element::Fp8_2(Fp8Element_2 { content: phix8, constants: inconsts }), 
                                                                    y: ExtFieldG2Element::Fp8_2(Fp8Element_2 { content: y, constants: inconsts }.negate()), 
                                                                    z: ExtFieldG2Element::Fp8_2(Fp8Element_2 { content: z, constants: inconsts }) } // frobinus at order 4 for an Fp8 is simply the identity
                                                    }
                                    }
                            }
                           }
                   -1 => {  if construction ==1 {   let invphix = [x[2].double().multiply(&frobs[4]).negate(),x[3].multiply(&frobs[4]).double(),
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
                                                                point: EcPoint{ x: ExtFieldG2Element::Fp8_1(Fp8Element_1 { content: invphix, constants: inconsts }), 
                                                                                y: ExtFieldG2Element::Fp8_1(Fp8Element_1 { content: invphiy, constants: inconsts }), 
                                                                                z: ExtFieldG2Element::Fp8_1(Fp8Element_1 { content: invphiz, constants: inconsts }) } 
                                                            }
                                                }
                            else {  if construction ==3 {  let invphix: [FieldElement<N>; 8] ;
                                                           let invphiy: [FieldElement<N>; 8] ;
                                                           let invphiz: [FieldElement<N>; 8] ;
                                                           if input.consts.twist_type=='D'{// D-Type twist 
                                                                                            invphix = [x[2].substract(&x[3]).multiply(&frobs[8]),x[2].addto(&x[3]).multiply(&frobs[8]).negate(),
                                                                                                        x[0].addto(&x[1]).multiply(&frobs[9]).negate(),x[1].substract(&x[0]).multiply(&frobs[9]),
                                                                                                        x[5].multiply(&frobs[5].sub(1u64)),x[4].multiply(&frobs[5].sub(1u64)),
                                                                                                        x[6].substract(&x[7]).multiply(&frobs[7]),x[6].addto(&x[7]).multiply(&frobs[7]).negate()];
                                                                                            invphiy = [y[4].addto(&y[5]).multiply(&frobs[15]),y[4].substract(&y[5]).multiply(&frobs[15]),
                                                                                                        y[7].multiply(&frobs[14]),y[6].multiply(&frobs[14]),
                                                                                                        y[0].multiply(&frobs[17]).negate(),y[1].multiply(&frobs[17]),
                                                                                                        y[2].addto(&y[3]).multiply(&frobs[16]),y[2].substract(&y[3]).multiply(&frobs[16])];                                           
                                                                                            // Inveting the frobinus for an Fp8                      
                                                                                            invphiz = [z[0],z[1].negate(),z[2].addto(&z[3]).multiply(&frobs[0]).negate(),
                                                                                                        z[3].substract(&z[2]).multiply(&frobs[0]),z[7].multiply(&frobs[1]).negate(),
                                                                                                        z[6].multiply(&frobs[1]).negate(),z[4].multiply(&frobs[2]).negate(),z[5].multiply(&frobs[2])];  
                                                                                        }
                                                            else {// M-Type twist 
                                                                // return toFp8([2*fp48FrobConsts[7]*tmp[3],-fp48FrobConsts[7]*2*tmp[4],-fp48FrobConsts[6]*tmp[2],-fp48FrobConsts[6]*tmp[1],
                                                                //     fp48FrobConsts[9]*(tmp[5]-tmp[6]),-fp48FrobConsts[9]*(tmp[5]+tmp[6]),fp48FrobConsts[8]*(tmp[7]),-fp48FrobConsts[8]*(tmp[8])]);    
                                                
                                                                invphix = [x[2].multiply(&frobs[4]).double(),x[3].multiply(&frobs[4]).double().negate(),
                                                                            x[1].multiply(&frobs[3]).negate(),x[0].multiply(&frobs[3].negate()),
                                                                            x[4].substract(&x[5]).multiply(&frobs[6]),x[4].addto(&x[5]).multiply(&frobs[6]).negate(),
                                                                            x[6].multiply(&frobs[5]),x[7].multiply(&frobs[5]).negate()];
                                                                // return toFp8([fp48FrobConsts[25]*(tmp[8]),fp48FrobConsts[25]*(tmp[7]), -fp48FrobConsts[20]*(tmp[5]),fp48FrobConsts[20]*(tmp[6]),
                                                                // fp48FrobConsts[18]*(tmp[3]+tmp[4]),fp48FrobConsts[18]*(tmp[3]-tmp[4]),fp48FrobConsts[26]*(tmp[1]-tmp[2]),-fp48FrobConsts[26]*(tmp[1]+tmp[2])]);                                                                                
                                                                invphiy = [y[7].multiply(&frobs[22]),y[6].multiply(&frobs[22]),
                                                                            y[4].multiply(&frobs[17]).negate(),y[5].multiply(&frobs[17]),
                                                                            y[2].addto(&y[3]).multiply(&frobs[15]),y[2].substract(&y[3]).multiply(&frobs[15]),
                                                                            y[0].substract(&y[1]).multiply(&frobs[23]),y[0].addto(&y[1]).multiply(&frobs[23]).negate()];                                           
                                                                // return toFp8([tmp[1],-tmp[2],-fp24FrobConsts[3]*(tmp[3]+tmp[4]),fp24FrobConsts[3]*(tmp[4]-tmp[3]),-fp24FrobConsts[4]*(tmp[8]),
                                                                // -fp24FrobConsts[4]*(tmp[7]),-fp24FrobConsts[5]*(tmp[5]),fp24FrobConsts[5]*(tmp[6])]);

                                                                // Inveting the frobinus for an Fp8                      
                                                                invphiz = [z[0],z[1].negate(),z[2].addto(&z[3]).multiply(&frobs[0]).negate(),
                                                                            z[3].substract(&z[2]).multiply(&frobs[0]),z[7].multiply(&frobs[1]).negate(),
                                                                            z[6].multiply(&frobs[1]).negate(),z[4].multiply(&frobs[2]).negate(),z[5].multiply(&frobs[2])];  
                                                                }                                                                                        
                                                        G2Element { consts: input.consts , 
                                                            point: EcPoint{ x: ExtFieldG2Element::Fp8_3(Fp8Element_3 { content: invphix, constants: inconsts }), 
                                                                            y: ExtFieldG2Element::Fp8_3(Fp8Element_3 { content: invphiy, constants: inconsts }), 
                                                                            z: ExtFieldG2Element::Fp8_3(Fp8Element_3 { content: invphiz, constants: inconsts }) } 
                                                        }                                                                                                        
                                                    }
                                    else {// construction ==2
                                        let invphix: [FieldElement<N>; 8] ;                  
                                        let invphiy: [FieldElement<N>; 8] ;                  
                                        let invphiz: [FieldElement<N>; 8] ;                  
                                        if *fb_id == 1{invphix = [x[1].multiply(&frobs[17]).negate(),x[0].multiply(&frobs[16]),x[3].multiply(&frobs[15]).negate(),
                                                                             x[2].multiply(&frobs[14]),x[4].multiply(&frobs[13]).negate(),x[5].multiply(&frobs[13]),
                                                                             x[6].multiply(&frobs[12]).negate(),x[7].multiply(&frobs[12])];
                                                       invphiy = [y[3].multiply(&frobs[27]),y[2].multiply(&frobs[28]),
                                                                             y[0].multiply(&frobs[29]).negate(),y[1].multiply(&frobs[29]),
                                                                             y[6].multiply(&frobs[24]).negate(),y[7].multiply(&frobs[24]),
                                                                             y[5].multiply(&frobs[25]),y[4].multiply(&frobs[26])];                                           
                                                        // Inveting the frobinus for an Fp8     
                                                        invphiz = [z[0],z[1].negate(),z[2].multiply(&frobs[0]).negate(),
                                                                                            z[3].multiply(&frobs[0]),z[5].multiply(&frobs[3]).negate(),
                                                                                            z[4].multiply(&frobs[4]),z[7].multiply(&frobs[1]).negate(),z[6].multiply(&frobs[2])];                                          
                                                        }                                                                                                 
                                        else {  // Phi variant for the BLS48_287           
                                                invphix = [x[0].multiply(&frobs[10]).negate(),x[1].multiply(&frobs[10]),x[2].multiply(&frobs[9]).negate(),
                                                           x[3].multiply(&frobs[9]),x[4].multiply(&frobs[8]).negate(),x[5].multiply(&frobs[8]),
                                                           x[6].multiply(&frobs[7]).negate(),x[7].multiply(&frobs[7])];
                                                invphiy = [ y[1].multiply(&frobs[26]),y[0].multiply(&frobs[25]),
                                                            y[3].multiply(&frobs[24]),y[2].multiply(&frobs[23]),
                                                            y[5].multiply(&frobs[22]),y[4].multiply(&frobs[21]),
                                                            y[7].multiply(&frobs[20]),y[6].multiply(&frobs[19])];                                           
                                                // Inveting the frobinus for an Fp8     
                                                invphiz = [z[0],z[1].negate(),z[2].multiply(&frobs[0]).negate(),z[3].multiply(&frobs[0]),
                                                           z[4].multiply(&frobs[2]).negate(),z[5].multiply(&frobs[2]),
                                                           z[6].multiply(&frobs[1]).negate(),z[7].multiply(&frobs[1])];                                          
                                             }
                                        G2Element { consts: input.consts , 
                                                    point: EcPoint{ x: ExtFieldG2Element::Fp8_2(Fp8Element_2 { content: invphix, constants: inconsts }), 
                                                                    y: ExtFieldG2Element::Fp8_2(Fp8Element_2 { content: invphiy, constants: inconsts }), 
                                                                    z: ExtFieldG2Element::Fp8_2(Fp8Element_2 { content: invphiz, constants: inconsts }) } 
                                                }
                                            }
                            }
                         }
                    _ => {panic!("non-implemented order for EFp8 phi ...")}
    }
}