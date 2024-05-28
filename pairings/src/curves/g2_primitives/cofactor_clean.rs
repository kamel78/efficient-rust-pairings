// Code developed by FARAOUN Kamel Mohamed.
// EEDIS-Laboratory. UDL-University. Algeria
// During May 2024.

use crate::curves::g2::G2Element;
use super::phi::{phi_bls24,phi_bls48};

pub fn clean_cofactor_bls12<const PRAMASIZE:usize, const R: usize, const N: usize, const MAX_COEFS_COUNT: usize>
                        (input :&G2Element<PRAMASIZE, R,N,MAX_COEFS_COUNT>) -> G2Element<PRAMASIZE, R,N,MAX_COEFS_COUNT>
{   //  Fast way to clean Cofactor on G2 for the BLS12 curves using Endomorphisme 
    //  based on 'Budroni-Pintore' approach (https://ia.cr/2017/419).
    let up = input.multiply_by_const(input.consts.u);
    let phip = input.phi();
    phip.phi().double().addto(&up.addto(&phip).multiply_by_const(input.consts.u).substract(&up).substract(&phip).substract(&input))   
}

pub fn clean_cofactor_bls24<const PRAMASIZE:usize, const R: usize, const N: usize, const MAX_COEFS_COUNT: usize>
                        (input :&G2Element<PRAMASIZE, R,N,MAX_COEFS_COUNT>) -> G2Element<PRAMASIZE, R,N,MAX_COEFS_COUNT>
{   //  Fast way to clean Cofactor on G2 for BLS24 curves using Endomorphisme 
    //  based on 'Budroni-Pintore' approach (https://ia.cr/2017/419).(personal use of inverted frobenius for optimization). 
    let phi4 = phi_bls24(input,4);
    let mut up = phi_bls24(&phi4,-1).multiply_by_const(input.consts.u - 1);    
    let mut phip = phi4.double().substract(&input).addto(&up);
    for _ in 0..3 {
        up = phi_bls24(&up,-1).multiply_by_const(input.consts.u);
        phip = phip.addto(&up);  
    }
    phip
}

pub fn clean_cofactor_bls48<const PRAMASIZE:usize, const R: usize, const N: usize, const MAX_COEFS_COUNT: usize>
                        (input :&G2Element<PRAMASIZE, R,N,MAX_COEFS_COUNT>) -> G2Element<PRAMASIZE, R,N,MAX_COEFS_COUNT>
{   
    //  Fast way to clean Cofactor on G2 for BLS48 curves using Endomorphisme
    // based on "Budroni-Pintore" approach (https://ia.cr/2017/419).
    //  Optimlized personaly using inversed endomorphisme
    let phi8 = phi_bls48(input,8);    
    let mut up = phi_bls48(&phi8,-1).multiply_by_const(input.consts.u - 1);    
    let mut phip = phi8.double().substract(&input).addto(&up);
    for _ in 0..7 {
        up = phi_bls48(&up,-1).multiply_by_const(input.consts.u);
        phip = phip.addto(&up);  
    }
    phip
}
