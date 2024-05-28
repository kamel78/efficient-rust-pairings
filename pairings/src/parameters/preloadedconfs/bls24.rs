// Code developed by FARAOUN Kamel Mohamed.
// EEDIS-Laboratory. UDL-University. Algeria
// During May 2024.

use crate::curves::g2::*;
use crate::curves::gt::GTField;
use crate::{curves::g1::G1Field, extensions::fp24::Fp24Field};
use crate::extensions::fp4::Fp4Field;
use crate::extensions::g2_extfields::ExtG2Field;
use crate::parameters::builders::*;
use crate::parameters::paramlist::{BLS24_479_PARAMS,BLS24_559_PARAMS};
use once_cell::sync::OnceCell;


static BLS24_479_FIELD_PARAMS   : OnceCell<FieldParams<8>>              = OnceCell::new();    
static BLS24_479_FROB_CONSTS    : OnceCell<ExFieldConsts<10,8>>         = OnceCell::new();
static BLS24_479_R_FIELD_PARAMS : OnceCell<FieldParams<7>>              = OnceCell::new();    
static FP_BLS24_479             : OnceCell<PrimeField<8>>               = OnceCell::new();
static FP4_BLS24_479            : OnceCell<Fp4Field<10,8>>              = OnceCell::new();
static EXFP4_BLS24_479          : OnceCell<ExtG2Field<8,10>>            = OnceCell::new();
static BLS24_479_G1_CONSTS      : OnceCell<G1Consts<7,8,4>>             = OnceCell::new();
static BLS24_479_G2_CONSTS      : OnceCell<G2Consts<10,7,8,10>>         = OnceCell::new();
static FR_BLS24_479             : OnceCell<PrimeField<7>>               = OnceCell::new();
static G1_BLS24_479             : OnceCell<G1Field<7,8,4>>              = OnceCell::new();
static G2_BLS24_479             : OnceCell<G2Field<10,7,8,10>>          = OnceCell::new();
static FP24_BLS24_479           : OnceCell<Fp24Field<10,8>>             = OnceCell::new();
static GT_BLS24_479             : OnceCell<GTField<8,10>>               = OnceCell::new();

pub fn fr_bls24_479()-> &'static PrimeField<7>
{   if FR_BLS24_479.get().is_none() 
            { if BLS24_479_R_FIELD_PARAMS.get().is_none() 
                    {BLS24_479_R_FIELD_PARAMS.set(build_field_params(BLS24_479_PARAMS.r_as_strhex)).unwrap()}; 
              FR_BLS24_479.set(PrimeField::new(&BLS24_479_R_FIELD_PARAMS.get().unwrap())).unwrap();                                                                                                
            }
    FR_BLS24_479.get().unwrap()
}
pub fn fp4_bls24_479()-> &'static Fp4Field<10, 8>
{   if FP4_BLS24_479.get().is_none() {  if BLS24_479_FIELD_PARAMS.get().is_none() 
                                            {   BLS24_479_FIELD_PARAMS.set(build_field_params(BLS24_479_PARAMS.modulo_as_strhex)).unwrap()};
                                         if FP_BLS24_479.get().is_none() 
                                            {   FP_BLS24_479.set(PrimeField::new(&BLS24_479_FIELD_PARAMS.get().unwrap())).unwrap()};
                                         if BLS24_479_FROB_CONSTS.get().is_none() 
                                            {   BLS24_479_FROB_CONSTS.set(build_frob_params(BLS24_479_PARAMS, &FP_BLS24_479.get().unwrap())).unwrap()};               
                                         FP4_BLS24_479.set(Fp4Field::new( &FP_BLS24_479.get().unwrap(), Some(&BLS24_479_FROB_CONSTS.get().unwrap()))).unwrap();                                                                                                 
                                      }
    FP4_BLS24_479.get().unwrap()
}

pub fn g1_bls24_479() -> &'static G1Field<7,8,4>
{   if G1_BLS24_479.get().is_none() 
            { if BLS24_479_R_FIELD_PARAMS.get().is_none()
                    {BLS24_479_R_FIELD_PARAMS.set(build_field_params(BLS24_479_PARAMS.r_as_strhex)).unwrap()}; 
              if BLS24_479_FIELD_PARAMS.get().is_none() 
                    {BLS24_479_FIELD_PARAMS.set(build_field_params(BLS24_479_PARAMS.modulo_as_strhex)).unwrap()}; 
              if FP_BLS24_479.get().is_none() 
                    {FP_BLS24_479.set(PrimeField::new(&BLS24_479_FIELD_PARAMS.get().unwrap())).unwrap();}                                                            
              if FR_BLS24_479.get().is_none() 
                    {FR_BLS24_479.set(PrimeField::new(&BLS24_479_R_FIELD_PARAMS.get().unwrap())).unwrap();}                                                                                                                           
              if BLS24_479_G1_CONSTS.get().is_none() 
                    {BLS24_479_G1_CONSTS.set(build_g1_params(BLS24_479_PARAMS,
                                                        &FP_BLS24_479.get().unwrap(),
                                                           &FR_BLS24_479.get().unwrap())).unwrap();}                                                                                                                                   
              G1_BLS24_479.set(G1Field { consts : BLS24_479_G1_CONSTS.get().unwrap(),
                                                  base_field: FP_BLS24_479.get().unwrap(),
                                                  fr_field : FR_BLS24_479.get().unwrap()}).unwrap() ;
            }
    G1_BLS24_479.get().unwrap()
}

pub fn g2_bls24_479() -> &'static G2Field<10,7,8,10>
{   if G2_BLS24_479.get().is_none() 
                {   if G1_BLS24_479.get().is_none() 
                            { let _ = g1_bls24_479();}
                    if BLS24_479_FROB_CONSTS.get().is_none() 
                            { BLS24_479_FROB_CONSTS.set(build_frob_params(BLS24_479_PARAMS, &FP_BLS24_479.get().unwrap())).unwrap()}; 
                    if FP4_BLS24_479.get().is_none() 
                            { let _ = fp4_bls24_479();}                                                            
                    if BLS24_479_G2_CONSTS.get().is_none() 
                            {BLS24_479_G2_CONSTS.set(build_g2_params(BLS24_479_PARAMS,
                                                                     &ExtG2Field::Fp4(FP4_BLS24_479.get().unwrap()),
                                                                     &PrimeField::new(&BLS24_479_R_FIELD_PARAMS.get().unwrap()),
                                                                     BLS24_479_FROB_CONSTS.get().unwrap())).unwrap();}                                                                                                                                   
                    if EXFP4_BLS24_479.get().is_none()  {EXFP4_BLS24_479.set(ExtG2Field::Fp4(&FP4_BLS24_479.get().unwrap())).unwrap();}                                     
                    G2_BLS24_479.set(G2Field {  consts : BLS24_479_G2_CONSTS.get().unwrap(),
                                                base_field : &EXFP4_BLS24_479.get().unwrap() ,
                                                fr_field : FR_BLS24_479.get().unwrap()}).unwrap() ;
                }            
    G2_BLS24_479.get().unwrap()
}

pub fn fp24_bls24_479()-> &'static Fp24Field<10, 8>
{   if FP24_BLS24_479.get().is_none() {  if BLS24_479_FIELD_PARAMS.get().is_none() 
                                            {   BLS24_479_FIELD_PARAMS.set(build_field_params(BLS24_479_PARAMS.modulo_as_strhex)).unwrap()};
                                         if FP_BLS24_479.get().is_none() 
                                            {   FP_BLS24_479.set(PrimeField::new(&BLS24_479_FIELD_PARAMS.get().unwrap())).unwrap()};
                                         if BLS24_479_FROB_CONSTS.get().is_none() 
                                            {   BLS24_479_FROB_CONSTS.set(build_frob_params(BLS24_479_PARAMS, &FP_BLS24_479.get().unwrap())).unwrap()};               
                                         FP24_BLS24_479.set(Fp24Field::new(&FP_BLS24_479.get().unwrap(), Some(&BLS24_479_FROB_CONSTS.get().unwrap()))).unwrap();                                                                                                 
                                      }
    FP24_BLS24_479.get().unwrap()
}

pub fn gt_bls24_479() -> &'static GTField<8,10>
{  
        if FP24_BLS24_479.get().is_none() { let _ = fp24_bls24_479(); };
        GT_BLS24_479.set(GTField::Fp24(FP24_BLS24_479.get().unwrap())).unwrap();
        GT_BLS24_479.get().unwrap()       
}

static BLS24_559_FIELD_PARAMS       : OnceCell<FieldParams<9>>              = OnceCell::new();    
static BLS24_559_FROB_CONSTS        : OnceCell<ExFieldConsts<10,9>>         = OnceCell::new();
static BLS24_559_R_FIELD_PARAMS     : OnceCell<FieldParams<8>>              = OnceCell::new();   
static FP_BLS24_559                 : OnceCell<PrimeField<9>>               = OnceCell::new(); 
static FP4_BLS24_559                : OnceCell<Fp4Field<10,9>>              = OnceCell::new();
static EXFP4_BLS24_559              : OnceCell<ExtG2Field<9,10>>            = OnceCell::new();
static FP24_BLS24_559               : OnceCell<Fp24Field<10,9>>             = OnceCell::new();
static BLS24_559_G1_CONSTS          : OnceCell<G1Consts<8,9,4>>             = OnceCell::new();
static BLS24_559_G2_CONSTS          : OnceCell<G2Consts<10,8,9,4>>          = OnceCell::new();
static FR_BLS24_559                 : OnceCell<PrimeField<8>>               = OnceCell::new();
static G1_BLS24_559                 : OnceCell<G1Field<8,9,4>>              = OnceCell::new();
static G2_BLS24_559                 : OnceCell<G2Field<10,8,9,4>>           = OnceCell::new();
static GT_BLS24_559                 : OnceCell<GTField<9,10>>               = OnceCell::new();


pub fn fr_bls24_559()-> &'static PrimeField<8>
{   if FR_BLS24_559.get().is_none() 
            { if BLS24_559_R_FIELD_PARAMS.get().is_none() 
                    {BLS24_559_R_FIELD_PARAMS.set(build_field_params(BLS24_559_PARAMS.r_as_strhex)).unwrap()}; 
              FR_BLS24_559.set(PrimeField::new(&BLS24_559_R_FIELD_PARAMS.get().unwrap())).unwrap();                                                                                                
            }
    FR_BLS24_559.get().unwrap()
}

pub fn fp4_bls24_559()-> &'static Fp4Field<10, 9>
{   if FP4_BLS24_559.get().is_none() {  if BLS24_559_FIELD_PARAMS.get().is_none() 
                                            {   BLS24_559_FIELD_PARAMS.set(build_field_params(BLS24_559_PARAMS.modulo_as_strhex)).unwrap()};
                                         if FP_BLS24_559.get().is_none() 
                                            {   FP_BLS24_559.set(PrimeField::new(&BLS24_559_FIELD_PARAMS.get().unwrap())).unwrap()};
                                         if BLS24_559_FROB_CONSTS.get().is_none() 
                                            {   BLS24_559_FROB_CONSTS.set(build_frob_params(BLS24_559_PARAMS, &FP_BLS24_559.get().unwrap())).unwrap()};               
                                         FP4_BLS24_559.set(Fp4Field::new(&FP_BLS24_559.get().unwrap(), Some(&BLS24_559_FROB_CONSTS.get().unwrap()))).unwrap();                                                                                                 
                                      }
    FP4_BLS24_559.get().unwrap()
}

pub fn g1_bls24_559() -> &'static G1Field<8,9,4>
{   if G1_BLS24_559.get().is_none() 
            { if BLS24_559_R_FIELD_PARAMS.get().is_none()
                    {BLS24_559_R_FIELD_PARAMS.set(build_field_params(BLS24_559_PARAMS.r_as_strhex)).unwrap()}; 
              if BLS24_559_FIELD_PARAMS.get().is_none() 
                    {BLS24_559_FIELD_PARAMS.set(build_field_params(BLS24_559_PARAMS.modulo_as_strhex)).unwrap()}; 
              if FP_BLS24_559.get().is_none() 
                    {FP_BLS24_559.set(PrimeField::new(&BLS24_559_FIELD_PARAMS.get().unwrap())).unwrap();}                                                            
              if FR_BLS24_559.get().is_none() 
                    {FR_BLS24_559.set(PrimeField::new(&BLS24_559_R_FIELD_PARAMS.get().unwrap())).unwrap();}                                                                                                                           
              if BLS24_559_G1_CONSTS.get().is_none() 
                    {BLS24_559_G1_CONSTS.set(build_g1_params(BLS24_559_PARAMS,
                                                        &FP_BLS24_559.get().unwrap(),
                                                           &FR_BLS24_559.get().unwrap())).unwrap();}                                                                                                                                   
              G1_BLS24_559.set(G1Field { consts : BLS24_559_G1_CONSTS.get().unwrap(),
                                                  base_field: FP_BLS24_559.get().unwrap(),
                                                  fr_field : FR_BLS24_559.get().unwrap()}).unwrap() ;
            }
    G1_BLS24_559.get().unwrap()
}

pub fn g2_bls24_559() -> &'static G2Field<10,8,9,4>
{   if G2_BLS24_559.get().is_none() 
                {   if G1_BLS24_559.get().is_none() 
                            { let _ = g1_bls24_559();}
                    if BLS24_559_FROB_CONSTS.get().is_none() 
                            { BLS24_559_FROB_CONSTS.set(build_frob_params(BLS24_559_PARAMS, &FP_BLS24_559.get().unwrap())).unwrap()}; 
                    if FP4_BLS24_559.get().is_none() 
                            { let _ = fp4_bls24_559();}                                                            
                    if BLS24_559_G2_CONSTS.get().is_none() 
                            {BLS24_559_G2_CONSTS.set(build_g2_params(BLS24_559_PARAMS,
                                                                     &ExtG2Field::Fp4(FP4_BLS24_559.get().unwrap()),
                                                                     &PrimeField::new(&BLS24_559_R_FIELD_PARAMS.get().unwrap()),
                                                                     BLS24_559_FROB_CONSTS.get().unwrap())).unwrap();}                                                                                                                                   
                    if EXFP4_BLS24_559.get().is_none()  {EXFP4_BLS24_559.set(ExtG2Field::Fp4(&FP4_BLS24_559.get().unwrap())).unwrap();}                                     
                    G2_BLS24_559.set(G2Field {  consts : BLS24_559_G2_CONSTS.get().unwrap(),
                                                base_field : &EXFP4_BLS24_559.get().unwrap() ,
                                                        fr_field : FR_BLS24_559.get().unwrap()}).unwrap() ;
                }            
    G2_BLS24_559.get().unwrap()
}

pub fn fp24_bls24_559()-> &'static Fp24Field<10, 9>
{   if FP24_BLS24_559.get().is_none() {  if BLS24_559_FIELD_PARAMS.get().is_none() 
                                            {   BLS24_559_FIELD_PARAMS.set(build_field_params(BLS24_559_PARAMS.modulo_as_strhex)).unwrap()};
                                         if FP_BLS24_559.get().is_none() 
                                            {   FP_BLS24_559.set(PrimeField::new(&BLS24_559_FIELD_PARAMS.get().unwrap())).unwrap()};
                                         if BLS24_559_FROB_CONSTS.get().is_none() 
                                            {   BLS24_559_FROB_CONSTS.set(build_frob_params(BLS24_559_PARAMS, &FP_BLS24_559.get().unwrap())).unwrap()};               
                                         FP24_BLS24_559.set(Fp24Field::new(&FP_BLS24_559.get().unwrap(), Some(&BLS24_559_FROB_CONSTS.get().unwrap()))).unwrap();                                                                                                 
                                      }
    FP24_BLS24_559.get().unwrap()
}

pub fn gt_bls24_559() -> &'static GTField<9,10>
{  
        if FP24_BLS24_559.get().is_none() { let _ = fp24_bls24_559(); };
        GT_BLS24_559.set(GTField::Fp24(FP24_BLS24_559.get().unwrap())).unwrap();
        GT_BLS24_559.get().unwrap()       
}