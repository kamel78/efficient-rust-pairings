// Code developed by FARAOUN Kamel Mohamed.
// EEDIS-Laboratory. UDL-University. Algeria
// During May 2024.

use crate::curves::g2::*;
use crate::curves::g1::G1Field;
use crate::curves::gt::GTField;
use crate::extensions::fp2::Fp2Field;
use crate::extensions::fp12::Fp12Field;
use crate::extensions::g2_extfields::ExtG2Field;

use crate::parameters::builders::*;
use crate::parameters::paramlist::{BLS12_381_PARAMS,BLS12_461_PARAMS};
use once_cell::sync::OnceCell;

// Parameters, fields and constants for BLS12-381

static BLS12_381_FIELD_PARAMS   : OnceCell<FieldParams<6>>          = OnceCell::new();    
static BLS12_381_R_FIELD_PARAMS : OnceCell<FieldParams<4>>          = OnceCell::new();    
static BLS12_381_FROB_CONSTS    : OnceCell<ExFieldConsts<4,6>>      = OnceCell::new();
static BLS12_381_G1_CONSTS      : OnceCell<G1Consts<4,6,16>>        = OnceCell::new();
static BLS12_381_G2_CONSTS      : OnceCell<G2Consts<4,4,6,4>>       = OnceCell::new();
static FP_BLS12_381             : OnceCell<PrimeField<6>>           = OnceCell::new();
static FP2_BLS12_381            : OnceCell<Fp2Field<6>>             = OnceCell::new();
static EXFP2_BLS12_381          : OnceCell<ExtG2Field<6,4>>         = OnceCell::new();
static FR_BLS12_381             : OnceCell<PrimeField<4>>           = OnceCell::new();
static G1_BLS12_381             : OnceCell<G1Field<4,6,16>>         = OnceCell::new();
static G2_BLS12_381             : OnceCell<G2Field<4,4,6,4>>        = OnceCell::new();
static FP12_BLS12_381           : OnceCell<Fp12Field<4,6>>          = OnceCell::new();
static GT_BLS12_381             : OnceCell<GTField<6,4>>            = OnceCell::new();

pub fn g1_bls12_381() -> &'static G1Field<4,6,16>
{   if G1_BLS12_381.get().is_none() 
            { if BLS12_381_R_FIELD_PARAMS.get().is_none()
                    {BLS12_381_R_FIELD_PARAMS.set(build_field_params(BLS12_381_PARAMS.r_as_strhex)).unwrap()}; 
              if BLS12_381_FIELD_PARAMS.get().is_none() 
                    {BLS12_381_FIELD_PARAMS.set(build_field_params(BLS12_381_PARAMS.modulo_as_strhex)).unwrap()}; 
              if FP_BLS12_381.get().is_none() 
                    {FP_BLS12_381.set(PrimeField::new(&BLS12_381_FIELD_PARAMS.get().unwrap())).unwrap();}                                                            
              if FR_BLS12_381.get().is_none() 
                    {FR_BLS12_381.set(PrimeField::new(&BLS12_381_R_FIELD_PARAMS.get().unwrap())).unwrap();}                                                                                                                           
              if BLS12_381_G1_CONSTS.get().is_none() 
                    {BLS12_381_G1_CONSTS.set(build_g1_params(BLS12_381_PARAMS,
                                                        &FP_BLS12_381.get().unwrap(),
                                                           &FR_BLS12_381.get().unwrap())).unwrap();}                                                                                                                                   
              G1_BLS12_381.set(G1Field { consts : BLS12_381_G1_CONSTS.get().unwrap(),
                                                base_field: FP_BLS12_381.get().unwrap(),
                                                fr_field : FR_BLS12_381.get().unwrap()}).unwrap() ;
            }
    G1_BLS12_381.get().unwrap()
}

pub fn g2_bls12_381() -> &'static G2Field<4,4,6,4>
{   if G2_BLS12_381.get().is_none() 
                {   if G1_BLS12_381.get().is_none() 
                            { let _ =g1_bls12_381();}
                    if BLS12_381_FROB_CONSTS.get().is_none() 
                            { BLS12_381_FROB_CONSTS.set(build_frob_params(BLS12_381_PARAMS, &FP_BLS12_381.get().unwrap())).unwrap()}; 
                    if FP2_BLS12_381.get().is_none() 
                            { FP2_BLS12_381.set(Fp2Field::new(&FP_BLS12_381.get().unwrap(),None )).unwrap();}                                                            
                    if BLS12_381_G2_CONSTS.get().is_none() 
                            {BLS12_381_G2_CONSTS.set(build_g2_params(BLS12_381_PARAMS,
                                                                     &ExtG2Field::Fp2(FP2_BLS12_381.get().unwrap()),
                                                                     &PrimeField::new(&BLS12_381_R_FIELD_PARAMS.get().unwrap()),
                                                                     BLS12_381_FROB_CONSTS.get().unwrap())).unwrap();}                                                                                                                                   
                    if EXFP2_BLS12_381.get().is_none()  {EXFP2_BLS12_381.set(ExtG2Field::Fp2(&FP2_BLS12_381.get().unwrap())).unwrap();}                                     
                    G2_BLS12_381.set(G2Field {  consts : BLS12_381_G2_CONSTS.get().unwrap(),
                                                base_field : &EXFP2_BLS12_381.get().unwrap() ,
                                                        fr_field : FR_BLS12_381.get().unwrap()}).unwrap() ;
                }            
    G2_BLS12_381.get().unwrap()
}

pub fn fp12_bls12_381()-> &'static Fp12Field<4, 6>
{   if FP12_BLS12_381.get().is_none() {  if BLS12_381_FIELD_PARAMS.get().is_none() 
                                                        {   BLS12_381_FIELD_PARAMS.set(build_field_params(BLS12_381_PARAMS.modulo_as_strhex)).unwrap()
                                                        };
                                         if FP_BLS12_381.get().is_none() 
                                                        {   FP_BLS12_381.set(PrimeField::new(&BLS12_381_FIELD_PARAMS.get().unwrap())).unwrap()
                                                        };
                                         if BLS12_381_FROB_CONSTS.get().is_none() 
                                                        {   BLS12_381_FROB_CONSTS.set(build_frob_params(BLS12_381_PARAMS, &FP_BLS12_381.get().unwrap())).unwrap()
                                                        };               
                                         FP12_BLS12_381.set(Fp12Field::new(&FP_BLS12_381.get().unwrap(), Some(&BLS12_381_FROB_CONSTS.get().unwrap()))).unwrap();                                                                                                 
                                      }
    FP12_BLS12_381.get().unwrap()
}

pub fn gt_bls12_381() -> &'static GTField<6,4>
{  
        if FP12_BLS12_381.get().is_none() { let _ = fp12_bls12_381(); };
        GT_BLS12_381.set(GTField::Fp12(FP12_BLS12_381.get().unwrap())).unwrap();
        GT_BLS12_381.get().unwrap()       
}


// Parameters, fields and constants for BLS12-461

static BLS12_461_FIELD_PARAMS   : OnceCell<FieldParams<8>>              = OnceCell::new();    
static BLS12_461_R_FIELD_PARAMS : OnceCell<FieldParams<5>>              = OnceCell::new();    
static BLS12_461_FROB_CONSTS    : OnceCell<ExFieldConsts<4,8>>          = OnceCell::new();
static BLS12_461_G1_CONSTS      : OnceCell<G1Consts<5,8,10>>            = OnceCell::new();
static BLS12_461_G2_CONSTS      : OnceCell<G2Consts<4,5,8,4>>           = OnceCell::new();
static FP_BLS12_461             : OnceCell<PrimeField<8>>               = OnceCell::new();
static FP2_BLS12_461            : OnceCell<Fp2Field<8>>                 = OnceCell::new();
static EXFP2_BLS12_461          : OnceCell<ExtG2Field<8,4>>             = OnceCell::new();
static FR_BLS12_461             : OnceCell<PrimeField<5>>               = OnceCell::new();
static G1_BLS12_461             : OnceCell<G1Field<5,8,10>>             = OnceCell::new();
static G2_BLS12_461             : OnceCell<G2Field<4,5,8,4>>            = OnceCell::new();
static FP12_BLS12_461           : OnceCell<Fp12Field<4,8>>              = OnceCell::new();
static GT_BLS12_461             : OnceCell<GTField<8,4>>                = OnceCell::new();

pub fn fr_bls12_381()-> &'static PrimeField<4>
{   if FR_BLS12_381.get().is_none() 
            { if BLS12_381_R_FIELD_PARAMS.get().is_none() 
                    {BLS12_381_R_FIELD_PARAMS.set(build_field_params(BLS12_381_PARAMS.r_as_strhex)).unwrap()}; 
              FR_BLS12_381.set(PrimeField::new(&BLS12_381_R_FIELD_PARAMS.get().unwrap())).unwrap();                                                                                                
            }
    FR_BLS12_381.get().unwrap()
}

pub fn fr_bls12_461()-> &'static PrimeField<5>
{   if FR_BLS12_461.get().is_none() 
            { if BLS12_461_R_FIELD_PARAMS.get().is_none() 
                    {BLS12_461_R_FIELD_PARAMS.set(build_field_params(BLS12_461_PARAMS.r_as_strhex)).unwrap()}; 
              FR_BLS12_461.set(PrimeField::new(&BLS12_461_R_FIELD_PARAMS.get().unwrap())).unwrap();                                                                                                
            }
    FR_BLS12_461.get().unwrap()
}

pub fn g1_bls12_461() -> &'static G1Field<5,8,10>
{   if G1_BLS12_461.get().is_none() 
            { if BLS12_461_R_FIELD_PARAMS.get().is_none()
                    {BLS12_461_R_FIELD_PARAMS.set(build_field_params(BLS12_461_PARAMS.r_as_strhex)).unwrap()}; 
              if BLS12_461_FIELD_PARAMS.get().is_none() 
                    {BLS12_461_FIELD_PARAMS.set(build_field_params(BLS12_461_PARAMS.modulo_as_strhex)).unwrap()}; 
              if FP_BLS12_461.get().is_none() 
                    {FP_BLS12_461.set(PrimeField::new(&BLS12_461_FIELD_PARAMS.get().unwrap())).unwrap();}                                                            
              if FR_BLS12_461.get().is_none() 
                    {FR_BLS12_461.set(PrimeField::new(&BLS12_461_R_FIELD_PARAMS.get().unwrap())).unwrap();}                                                                                                                           
              if BLS12_461_G1_CONSTS.get().is_none() 
                    {BLS12_461_G1_CONSTS.set(build_g1_params(BLS12_461_PARAMS,
                                                        &FP_BLS12_461.get().unwrap(),
                                                           &FR_BLS12_461.get().unwrap())).unwrap();}                                                                                                                                   
              G1_BLS12_461.set(G1Field { consts : BLS12_461_G1_CONSTS.get().unwrap(),
                                                  base_field: FP_BLS12_461.get().unwrap(),
                                                  fr_field : FR_BLS12_461.get().unwrap()}).unwrap() ;
            }
    G1_BLS12_461.get().unwrap()
}

pub fn g2_bls12_461() -> &'static G2Field<4,5,8,4>
{   if G2_BLS12_461.get().is_none() 
                {   if G1_BLS12_461.get().is_none() 
                            { let _ = g1_bls12_461();}
                    if BLS12_461_FROB_CONSTS.get().is_none() 
                            { BLS12_461_FROB_CONSTS.set(build_frob_params(BLS12_461_PARAMS, &FP_BLS12_461.get().unwrap())).unwrap()}; 
                    if FP2_BLS12_461.get().is_none() 
                            { FP2_BLS12_461.set(Fp2Field::new(&FP_BLS12_461.get().unwrap(),None )).unwrap();}                                                            
                    if BLS12_461_G2_CONSTS.get().is_none() 
                            {BLS12_461_G2_CONSTS.set(build_g2_params(BLS12_461_PARAMS,
                                                                     &ExtG2Field::Fp2(FP2_BLS12_461.get().unwrap()),
                                                                     &PrimeField::new(&BLS12_461_R_FIELD_PARAMS.get().unwrap()),
                                                                     BLS12_461_FROB_CONSTS.get().unwrap())).unwrap();}                                                                                                                                   
                    if EXFP2_BLS12_461.get().is_none()  {EXFP2_BLS12_461.set(ExtG2Field::Fp2(&FP2_BLS12_461.get().unwrap())).unwrap();}                                     
                    G2_BLS12_461.set(G2Field {  consts : BLS12_461_G2_CONSTS.get().unwrap(),
                                                base_field : &EXFP2_BLS12_461.get().unwrap() ,
                                                fr_field : FR_BLS12_461.get().unwrap()}).unwrap() ;
                }            
    G2_BLS12_461.get().unwrap()
}

pub fn fp12_bls12_461()-> &'static Fp12Field<4, 8>
{   if FP12_BLS12_461.get().is_none() {  if BLS12_461_FIELD_PARAMS.get().is_none() 
                                                        {   BLS12_461_FIELD_PARAMS.set(build_field_params(BLS12_461_PARAMS.modulo_as_strhex)).unwrap()
                                                        };
                                         if FP_BLS12_461.get().is_none() 
                                                        {   FP_BLS12_381.set(PrimeField::new(&BLS12_381_FIELD_PARAMS.get().unwrap())).unwrap()
                                                        };
                                         if BLS12_461_FROB_CONSTS.get().is_none() 
                                                        {   BLS12_461_FROB_CONSTS.set(build_frob_params(BLS12_461_PARAMS, &FP_BLS12_461.get().unwrap())).unwrap()
                                                        };               
                                         FP12_BLS12_461.set(Fp12Field::new(&FP_BLS12_461.get().unwrap(), Some(&BLS12_461_FROB_CONSTS.get().unwrap()))).unwrap();                                                                                                 
                                      }
    FP12_BLS12_461.get().unwrap()
}

pub fn gt_bls12_461() -> &'static GTField<8,4>
{  
        if FP12_BLS12_461.get().is_none() { let _ = fp12_bls12_461(); };
        GT_BLS12_461.set(GTField::Fp12(FP12_BLS12_461.get().unwrap())).unwrap();
        GT_BLS12_461.get().unwrap()       
}

