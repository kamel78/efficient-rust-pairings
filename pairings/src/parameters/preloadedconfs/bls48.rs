// Code developed by FARAOUN Kamel Mohamed.
// EEDIS-Laboratory. UDL-University. Algeria
// During May 2024.

use crate::curves::g1::{G1Consts, G1Field};
use crate::curves::g2::{G2Consts, G2Field};
use crate::curves::gt::GTField;
use crate::extensions::g2_extfields::ExtG2Field;
use crate::extensions::{fp48::Fp48Field,fp8::Fp8Field};
use crate::parameters::builders::{build_field_params, build_frob_params, build_g1_params, build_g2_params, ExFieldConsts, ExtField, FieldParams, PrimeField};
use crate::parameters::paramlist::BLS48_575_PARAMS;
use once_cell::sync::OnceCell;

static BLS48_575_FIELD_PARAMS       : OnceCell<FieldParams<9>>           = OnceCell::new();    
static BLS48_575_FROB_CONSTS        : OnceCell<ExFieldConsts<24,9>>      = OnceCell::new();
static BLS48_575_R_FIELD_PARAMS     : OnceCell<FieldParams<9>>           = OnceCell::new();    
static FP_BLS48_575                 : OnceCell<PrimeField<9>>            = OnceCell::new();  
static FP8_BLS48_575                : OnceCell<Fp8Field<24,9>>           = OnceCell::new();
static EXFP8_BLS48_575              : OnceCell<ExtG2Field<9,24>>         = OnceCell::new();
static BLS48_575_G1_CONSTS          : OnceCell<G1Consts<9,9,4>>          = OnceCell::new();
static BLS48_575_G2_CONSTS          : OnceCell<G2Consts<24,9,9,7>>       = OnceCell::new();
static FR_BLS48_575                 : OnceCell<PrimeField<9>>            = OnceCell::new();
static G1_BLS48_575                 : OnceCell<G1Field<9,9,4>>           = OnceCell::new();
static G2_BLS48_575                 : OnceCell<G2Field<24,9,9,7>>        = OnceCell::new();
static FP48_BLS48_575               : OnceCell<Fp48Field<24,9>>          = OnceCell::new();
static GT_BLS48_575                 : OnceCell<GTField<9,24>>           = OnceCell::new();

pub fn fr_bls48_575()-> &'static PrimeField<9>
{   if FR_BLS48_575.get().is_none() 
            { if BLS48_575_R_FIELD_PARAMS.get().is_none() 
                    {BLS48_575_R_FIELD_PARAMS.set(build_field_params(BLS48_575_PARAMS.r_as_strhex)).unwrap()}; 
              FR_BLS48_575.set(PrimeField::new(&BLS48_575_R_FIELD_PARAMS.get().unwrap())).unwrap();                                                                                                
            }
    FR_BLS48_575.get().unwrap()
}
pub fn fp8_bls48_575()-> &'static Fp8Field<24,9>
{   if FP8_BLS48_575.get().is_none() {  if BLS48_575_FIELD_PARAMS.get().is_none() 
                                                        {   BLS48_575_FIELD_PARAMS.set(build_field_params(BLS48_575_PARAMS.modulo_as_strhex)).unwrap()
                                                        };
                                         if FP_BLS48_575.get().is_none() 
                                                        {   FP_BLS48_575.set(PrimeField::new(&BLS48_575_FIELD_PARAMS.get().unwrap())).unwrap()
                                                        };
                                         if BLS48_575_FROB_CONSTS.get().is_none() 
                                                        {   BLS48_575_FROB_CONSTS.set(build_frob_params(BLS48_575_PARAMS, &FP_BLS48_575.get().unwrap())).unwrap()
                                                        };               
                                         FP8_BLS48_575.set(Fp8Field::new(&FP_BLS48_575.get().unwrap(), Some(&BLS48_575_FROB_CONSTS.get().unwrap()))).unwrap();                                                                                                 
                                      }
    FP8_BLS48_575.get().unwrap()
}

pub fn g1_bls48_575() -> &'static G1Field<9,9,4>
{   if G1_BLS48_575.get().is_none() 
            { if BLS48_575_R_FIELD_PARAMS.get().is_none()
                    {BLS48_575_R_FIELD_PARAMS.set(build_field_params(BLS48_575_PARAMS.r_as_strhex)).unwrap();}; 
              if BLS48_575_FIELD_PARAMS.get().is_none() 
                    {BLS48_575_FIELD_PARAMS.set(build_field_params(BLS48_575_PARAMS.modulo_as_strhex)).unwrap();}; 
              if FP_BLS48_575.get().is_none() 
                    {FP_BLS48_575.set(PrimeField::new(&BLS48_575_FIELD_PARAMS.get().unwrap())).unwrap();}                                                            
              if FR_BLS48_575.get().is_none() 
                    {FR_BLS48_575.set(PrimeField::new(&BLS48_575_R_FIELD_PARAMS.get().unwrap())).unwrap();}                                                                                                                           
              if BLS48_575_G1_CONSTS.get().is_none() 
                    {BLS48_575_G1_CONSTS.set(build_g1_params(BLS48_575_PARAMS,
                                                             &FP_BLS48_575.get().unwrap(),
                                                             &FR_BLS48_575.get().unwrap())).unwrap();}                                                                                                                                   
              G1_BLS48_575.set(G1Field { consts : BLS48_575_G1_CONSTS.get().unwrap(),
                                                  base_field: FP_BLS48_575.get().unwrap(),
                                                  fr_field : FR_BLS48_575.get().unwrap()}).unwrap() ;
            }
    G1_BLS48_575.get().unwrap()
}

pub fn g2_bls48_575() -> &'static G2Field<24,9,9,7>
{   if G2_BLS48_575.get().is_none() 
                {   if G1_BLS48_575.get().is_none() 
                            { let _ = g1_bls48_575();}
                    if BLS48_575_FROB_CONSTS.get().is_none() 
                            { BLS48_575_FROB_CONSTS.set(build_frob_params(BLS48_575_PARAMS, &FP_BLS48_575.get().unwrap())).unwrap()};                                         
                    if FP8_BLS48_575.get().is_none() 
                                        { let _ = fp8_bls48_575();}        
                    if BLS48_575_G2_CONSTS.get().is_none() 
                            {BLS48_575_G2_CONSTS.set(build_g2_params(BLS48_575_PARAMS,
                                                                     &ExtG2Field::Fp8(FP8_BLS48_575.get().unwrap()),
                                                                     &PrimeField::new(&BLS48_575_R_FIELD_PARAMS.get().unwrap()),
                                                                     BLS48_575_FROB_CONSTS.get().unwrap())).unwrap();}                                                                                                                                                                                                          
                    if EXFP8_BLS48_575.get().is_none()  {EXFP8_BLS48_575.set(ExtG2Field::Fp8(&FP8_BLS48_575.get().unwrap())).unwrap();}                                     
                    G2_BLS48_575.set(G2Field {  consts : BLS48_575_G2_CONSTS.get().unwrap(),
                                                base_field : &EXFP8_BLS48_575.get().unwrap() ,
                                                fr_field : FR_BLS48_575.get().unwrap()}).unwrap() ;
                }            
    G2_BLS48_575.get().unwrap()
}

pub fn fp48_bls48_575()-> &'static Fp48Field<24, 9>
{   if FP48_BLS48_575.get().is_none() {  if BLS48_575_FIELD_PARAMS.get().is_none() 
                                                        {   BLS48_575_FIELD_PARAMS.set(build_field_params(BLS48_575_PARAMS.modulo_as_strhex)).unwrap()
                                                        };
                                         if FP_BLS48_575.get().is_none() 
                                                        {   FP_BLS48_575.set(PrimeField::new(&BLS48_575_FIELD_PARAMS.get().unwrap())).unwrap()
                                                        };
                                         if BLS48_575_FROB_CONSTS.get().is_none() 
                                                        {   BLS48_575_FROB_CONSTS.set(build_frob_params(BLS48_575_PARAMS, &FP_BLS48_575.get().unwrap())).unwrap()
                                                        };               
                                         FP48_BLS48_575.set(Fp48Field::new(&FP_BLS48_575.get().unwrap(), Some(&BLS48_575_FROB_CONSTS.get().unwrap()))).unwrap();                                                                                                 
                                      }
    FP48_BLS48_575.get().unwrap()
}

pub fn gt_bls48_575() -> &'static GTField<9,24>
{  
        if FP48_BLS48_575.get().is_none() { let _ = fp48_bls48_575(); };
        GT_BLS48_575.set(GTField::Fp48(FP48_BLS48_575.get().unwrap())).unwrap();
        GT_BLS48_575.get().unwrap()       
}