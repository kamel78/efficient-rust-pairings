// Code developed by FARAOUN Kamel Mohamed.
// EEDIS-Laboratory. UDL-University. Algeria
// During May 2024.

use crate::curves::g1::{G1Consts, G1Field};
use crate::curves::g2::{G2Consts, G2Field};
use crate::curves::gt::GTField;
use crate::extensions::ext_fields::{ExFieldConsts, ExtField};
use crate::extensions::g2_extfields::ExtG2Field;
use crate::extensions::towering1::{fp48::Fp48Field as Fp48Field_1,fp8::Fp8Field as Fp8Field_1};
use crate::extensions::towering2::{fp48::Fp48Field as Fp48Field_2,fp8::Fp8Field as Fp8Field_2};
use crate::extensions::towering3::{fp48::Fp48Field as Fp48Field_3,fp8::Fp8Field as Fp8Field_3};
use crate::fields::prime_fields::{FieldParams, PrimeField};
use crate::parameters::builders::{build_field_params, build_frob_params, build_g1_params, build_g2_params};
use crate::parameters::paramlist::{BLS48_277_PARAMS, BLS48_287_PARAMS, BLS48_571_PARAMS, BLS48_573_PARAMS, BLS48_575_PARAMS, BLS48_581_PARAMS};
use once_cell::sync::OnceCell;

static BLS48_575_FIELD_PARAMS       : OnceCell<FieldParams<9>>           = OnceCell::new();    
static BLS48_575_FROB_CONSTS        : OnceCell<ExFieldConsts<24,9>>      = OnceCell::new();
static BLS48_575_R_FIELD_PARAMS     : OnceCell<FieldParams<9>>           = OnceCell::new();    
static FP_BLS48_575                 : OnceCell<PrimeField<9>>            = OnceCell::new();  
static FP8_BLS48_575                : OnceCell<Fp8Field_1<24,9>>           = OnceCell::new();
static EXFP8_BLS48_575              : OnceCell<ExtG2Field<9,24>>         = OnceCell::new();
static BLS48_575_G1_CONSTS          : OnceCell<G1Consts<9,9,4>>          = OnceCell::new();
static BLS48_575_G2_CONSTS          : OnceCell<G2Consts<24,9,9,7>>       = OnceCell::new();
static FR_BLS48_575                 : OnceCell<PrimeField<9>>            = OnceCell::new();
static G1_BLS48_575                 : OnceCell<G1Field<9,9,4>>           = OnceCell::new();
static G2_BLS48_575                 : OnceCell<G2Field<24,9,9,7>>        = OnceCell::new();
static FP48_BLS48_575               : OnceCell<Fp48Field_1<24,9>>          = OnceCell::new();
static GT_BLS48_575                 : OnceCell<GTField<9,24>>           = OnceCell::new();

pub fn fr_bls48_575()-> &'static PrimeField<9>
{   if FR_BLS48_575.get().is_none() 
            { if BLS48_575_R_FIELD_PARAMS.get().is_none() 
                    {BLS48_575_R_FIELD_PARAMS.set(build_field_params(BLS48_575_PARAMS.r_as_strhex)).unwrap()}; 
              FR_BLS48_575.set(PrimeField::new(&BLS48_575_R_FIELD_PARAMS.get().unwrap())).unwrap();                                                                                                
            }
    FR_BLS48_575.get().unwrap()
}
pub fn fp8_bls48_575()-> &'static Fp8Field_1<24,9>
{   if FP8_BLS48_575.get().is_none() {  if BLS48_575_FIELD_PARAMS.get().is_none() 
                                                        {   BLS48_575_FIELD_PARAMS.set(build_field_params(BLS48_575_PARAMS.modulo_as_strhex)).unwrap()
                                                        };
                                         if FP_BLS48_575.get().is_none() 
                                                        {   FP_BLS48_575.set(PrimeField::new(&BLS48_575_FIELD_PARAMS.get().unwrap())).unwrap()
                                                        };
                                         if BLS48_575_FROB_CONSTS.get().is_none() 
                                                        {   BLS48_575_FROB_CONSTS.set(build_frob_params(BLS48_575_PARAMS, &FP_BLS48_575.get().unwrap())).unwrap()
                                                        };               
                                         FP8_BLS48_575.set(Fp8Field_1::new(&FP_BLS48_575.get().unwrap(), Some(&BLS48_575_FROB_CONSTS.get().unwrap()))).unwrap();                                                                                                 
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
                                                                     &ExtG2Field::Fp8_1(FP8_BLS48_575.get().unwrap()),
                                                                     &PrimeField::new(&BLS48_575_R_FIELD_PARAMS.get().unwrap()),
                                                                     BLS48_575_FROB_CONSTS.get().unwrap())).unwrap();}                                                                                                                                                                                                          
                    if EXFP8_BLS48_575.get().is_none()  {EXFP8_BLS48_575.set(ExtG2Field::Fp8_1(&FP8_BLS48_575.get().unwrap())).unwrap();}                                     
                    G2_BLS48_575.set(G2Field {  consts : BLS48_575_G2_CONSTS.get().unwrap(),
                                                base_field : &EXFP8_BLS48_575.get().unwrap() ,
                                                fr_field : FR_BLS48_575.get().unwrap()}).unwrap() ;
                }            
    G2_BLS48_575.get().unwrap()
}

pub fn fp48_bls48_575()-> &'static Fp48Field_1<24, 9>
{   if FP48_BLS48_575.get().is_none() {  if BLS48_575_FIELD_PARAMS.get().is_none() 
                                                        {   BLS48_575_FIELD_PARAMS.set(build_field_params(BLS48_575_PARAMS.modulo_as_strhex)).unwrap()
                                                        };
                                         if FP_BLS48_575.get().is_none() 
                                                        {   FP_BLS48_575.set(PrimeField::new(&BLS48_575_FIELD_PARAMS.get().unwrap())).unwrap()
                                                        };
                                         if BLS48_575_FROB_CONSTS.get().is_none() 
                                                        {   BLS48_575_FROB_CONSTS.set(build_frob_params(BLS48_575_PARAMS, &FP_BLS48_575.get().unwrap())).unwrap()
                                                        };               
                                         FP48_BLS48_575.set(Fp48Field_1::new(&FP_BLS48_575.get().unwrap(), Some(&BLS48_575_FROB_CONSTS.get().unwrap()))).unwrap();                                                                                                 
                                      }
    FP48_BLS48_575.get().unwrap()
}

pub fn gt_bls48_575() -> &'static GTField<9,24>
{  
        if FP48_BLS48_575.get().is_none() { let _ = fp48_bls48_575(); };
        GT_BLS48_575.set(GTField::Fp48_1(FP48_BLS48_575.get().unwrap())).unwrap();
        GT_BLS48_575.get().unwrap()       
}


static BLS48_573_FIELD_PARAMS       : OnceCell<FieldParams<9>>           = OnceCell::new();    
static BLS48_573_FROB_CONSTS        : OnceCell<ExFieldConsts<24,9>>      = OnceCell::new();
static BLS48_573_R_FIELD_PARAMS     : OnceCell<FieldParams<8>>           = OnceCell::new();    
static FP_BLS48_573                 : OnceCell<PrimeField<9>>            = OnceCell::new();  
static FP8_BLS48_573                : OnceCell<Fp8Field_1<24,9>>           = OnceCell::new();
static EXFP8_BLS48_573              : OnceCell<ExtG2Field<9,24>>         = OnceCell::new();
static BLS48_573_G1_CONSTS          : OnceCell<G1Consts<8,9,7>>          = OnceCell::new();
static BLS48_573_G2_CONSTS          : OnceCell<G2Consts<24,8,9,4>>       = OnceCell::new();
static FR_BLS48_573                 : OnceCell<PrimeField<8>>            = OnceCell::new();
static G2_BLS48_573                 : OnceCell<G2Field<24,8,9,4>>        = OnceCell::new();
static G1_BLS48_573                 : OnceCell<G1Field<8,9,7>>           = OnceCell::new();
static FP48_BLS48_573               : OnceCell<Fp48Field_1<24,9>>          = OnceCell::new();
static GT_BLS48_573                 : OnceCell<GTField<9,24>>           = OnceCell::new();

pub fn fr_bls48_573()-> &'static PrimeField<8>
{   if FR_BLS48_573.get().is_none() 
            { if BLS48_573_R_FIELD_PARAMS.get().is_none() 
                    {BLS48_573_R_FIELD_PARAMS.set(build_field_params(BLS48_573_PARAMS.r_as_strhex)).unwrap()}; 
              FR_BLS48_573.set(PrimeField::new(&BLS48_573_R_FIELD_PARAMS.get().unwrap())).unwrap();                                                                                                
            }
    FR_BLS48_573.get().unwrap()
}
pub fn fp8_bls48_573()-> &'static Fp8Field_1<24,9>
{   if FP8_BLS48_573.get().is_none() {  if BLS48_573_FIELD_PARAMS.get().is_none() 
                                                        {   BLS48_573_FIELD_PARAMS.set(build_field_params(BLS48_573_PARAMS.modulo_as_strhex)).unwrap()
                                                        };
                                         if FP_BLS48_573.get().is_none() 
                                                        {   FP_BLS48_573.set(PrimeField::new(&BLS48_573_FIELD_PARAMS.get().unwrap())).unwrap()
                                                        };
                                         if BLS48_573_FROB_CONSTS.get().is_none() 
                                                        {   BLS48_573_FROB_CONSTS.set(build_frob_params(BLS48_573_PARAMS, &FP_BLS48_573.get().unwrap())).unwrap()
                                                        };               
                                         FP8_BLS48_573.set(Fp8Field_1::new(&FP_BLS48_573.get().unwrap(), Some(&BLS48_573_FROB_CONSTS.get().unwrap()))).unwrap();                                                                                                 
                                      }
    FP8_BLS48_573.get().unwrap()
}

pub fn g1_bls48_573() -> &'static G1Field<8,9,7>
{   if G1_BLS48_573.get().is_none() 
            { if BLS48_573_R_FIELD_PARAMS.get().is_none()
                    {BLS48_573_R_FIELD_PARAMS.set(build_field_params(BLS48_573_PARAMS.r_as_strhex)).unwrap();}; 
              if BLS48_573_FIELD_PARAMS.get().is_none() 
                    {BLS48_573_FIELD_PARAMS.set(build_field_params(BLS48_573_PARAMS.modulo_as_strhex)).unwrap();}; 
              if FP_BLS48_573.get().is_none() 
                    {FP_BLS48_573.set(PrimeField::new(&BLS48_573_FIELD_PARAMS.get().unwrap())).unwrap();}                                                            
              if FR_BLS48_573.get().is_none() 
                    {FR_BLS48_573.set(PrimeField::new(&BLS48_573_R_FIELD_PARAMS.get().unwrap())).unwrap();}                                                                                                                           
              if BLS48_573_G1_CONSTS.get().is_none() 
                    {BLS48_573_G1_CONSTS.set(build_g1_params(BLS48_573_PARAMS,
                                                             &FP_BLS48_573.get().unwrap(),
                                                             &FR_BLS48_573.get().unwrap())).unwrap();}                                                                                                                                   
              G1_BLS48_573.set(G1Field { consts : BLS48_573_G1_CONSTS.get().unwrap(),
                                                  base_field: FP_BLS48_573.get().unwrap(),
                                                  fr_field : FR_BLS48_573.get().unwrap()}).unwrap() ;
            }
    G1_BLS48_573.get().unwrap()
}

pub fn g2_bls48_573() -> &'static G2Field<24,8,9,4>
{   if G2_BLS48_573.get().is_none() 
                {   if G1_BLS48_573.get().is_none() 
                            { let _ = g1_bls48_573();}
                    if BLS48_573_FROB_CONSTS.get().is_none() 
                            { BLS48_573_FROB_CONSTS.set(build_frob_params(BLS48_573_PARAMS, &FP_BLS48_573.get().unwrap())).unwrap()};                                         
                    if FP8_BLS48_573.get().is_none() 
                                        { let _ = fp8_bls48_573();}        
                    if BLS48_573_G2_CONSTS.get().is_none() 
                            {BLS48_573_G2_CONSTS.set(build_g2_params(BLS48_573_PARAMS,
                                                                     &ExtG2Field::Fp8_1(FP8_BLS48_573.get().unwrap()),
                                                                     &PrimeField::new(&BLS48_573_R_FIELD_PARAMS.get().unwrap()),
                                                                     BLS48_573_FROB_CONSTS.get().unwrap())).unwrap();}                                                                                                                                                                                                          
                    if EXFP8_BLS48_573.get().is_none()  {EXFP8_BLS48_573.set(ExtG2Field::Fp8_1(&FP8_BLS48_573.get().unwrap())).unwrap();}                                     
                    G2_BLS48_573.set(G2Field {  consts : BLS48_573_G2_CONSTS.get().unwrap(),
                                                base_field : &EXFP8_BLS48_573.get().unwrap() ,
                                                fr_field : FR_BLS48_573.get().unwrap()}).unwrap() ;
                }            
    G2_BLS48_573.get().unwrap()
}

pub fn fp48_bls48_573()-> &'static Fp48Field_1<24, 9>
{   if FP48_BLS48_573.get().is_none() {  if BLS48_573_FIELD_PARAMS.get().is_none() 
                                                        {   BLS48_573_FIELD_PARAMS.set(build_field_params(BLS48_573_PARAMS.modulo_as_strhex)).unwrap()
                                                        };
                                         if FP_BLS48_573.get().is_none() 
                                                        {   FP_BLS48_573.set(PrimeField::new(&BLS48_573_FIELD_PARAMS.get().unwrap())).unwrap()
                                                        };
                                         if BLS48_573_FROB_CONSTS.get().is_none() 
                                                        {   BLS48_573_FROB_CONSTS.set(build_frob_params(BLS48_573_PARAMS, &FP_BLS48_573.get().unwrap())).unwrap()
                                                        };               
                                         FP48_BLS48_573.set(Fp48Field_1::new(&FP_BLS48_573.get().unwrap(), Some(&BLS48_573_FROB_CONSTS.get().unwrap()))).unwrap();                                                                                                 
                                      }
    FP48_BLS48_573.get().unwrap()
}

pub fn gt_bls48_573() -> &'static GTField<9,24>
{  
        if FP48_BLS48_573.get().is_none() { let _ = fp48_bls48_573(); };
        GT_BLS48_573.set(GTField::Fp48_1(FP48_BLS48_573.get().unwrap())).unwrap();
        GT_BLS48_573.get().unwrap()       
}

static BLS48_581_FIELD_PARAMS       : OnceCell<FieldParams<10>>           = OnceCell::new();    
static BLS48_581_FROB_CONSTS        : OnceCell<ExFieldConsts<24,10>>      = OnceCell::new();
static BLS48_581_R_FIELD_PARAMS     : OnceCell<FieldParams<9>>           = OnceCell::new();    
static FP_BLS48_581                 : OnceCell<PrimeField<10>>            = OnceCell::new();  
static FP8_BLS48_581                : OnceCell<Fp8Field_3<24,10>>           = OnceCell::new();
static EXFP8_BLS48_581              : OnceCell<ExtG2Field<10,24>>         = OnceCell::new();
static BLS48_581_G1_CONSTS          : OnceCell<G1Consts<9,10,3>>          = OnceCell::new();
static BLS48_581_G2_CONSTS          : OnceCell<G2Consts<24,9,10,4>>       = OnceCell::new();
static FR_BLS48_581                 : OnceCell<PrimeField<9>>            = OnceCell::new();
static G2_BLS48_581                 : OnceCell<G2Field<24,9,10,4>>        = OnceCell::new();
static G1_BLS48_581                 : OnceCell<G1Field<9,10,3>>           = OnceCell::new();
static FP48_BLS48_581               : OnceCell<Fp48Field_3<24,10>>          = OnceCell::new();
static GT_BLS48_581                 : OnceCell<GTField<10,24>>           = OnceCell::new();

pub fn fr_bls48_581()-> &'static PrimeField<9>
{   if FR_BLS48_581.get().is_none() 
            { if BLS48_581_R_FIELD_PARAMS.get().is_none() 
                    {BLS48_581_R_FIELD_PARAMS.set(build_field_params(BLS48_581_PARAMS.r_as_strhex)).unwrap()}; 
              FR_BLS48_581.set(PrimeField::new(&BLS48_581_R_FIELD_PARAMS.get().unwrap())).unwrap();                                                                                                
            }
    FR_BLS48_581.get().unwrap()
}
pub fn fp8_bls48_581()-> &'static Fp8Field_3<24,10>
{   if FP8_BLS48_581.get().is_none() {  if BLS48_581_FIELD_PARAMS.get().is_none() 
                                                        {   BLS48_581_FIELD_PARAMS.set(build_field_params(BLS48_581_PARAMS.modulo_as_strhex)).unwrap()
                                                        };
                                         if FP_BLS48_581.get().is_none() 
                                                        {   FP_BLS48_581.set(PrimeField::new(&BLS48_581_FIELD_PARAMS.get().unwrap())).unwrap()
                                                        };
                                         if BLS48_581_FROB_CONSTS.get().is_none() 
                                                        {   BLS48_581_FROB_CONSTS.set(build_frob_params(BLS48_581_PARAMS, &FP_BLS48_581.get().unwrap())).unwrap()
                                                        };               
                                         FP8_BLS48_581.set(Fp8Field_3::new(&FP_BLS48_581.get().unwrap(), Some(&BLS48_581_FROB_CONSTS.get().unwrap()))).unwrap();                                                                                                 
                                      }
    FP8_BLS48_581.get().unwrap()
}

pub fn g1_bls48_581() -> &'static G1Field<9,10,3>
{   if G1_BLS48_581.get().is_none() 
            { if BLS48_581_R_FIELD_PARAMS.get().is_none()
                    {BLS48_581_R_FIELD_PARAMS.set(build_field_params(BLS48_581_PARAMS.r_as_strhex)).unwrap();}; 
              if BLS48_581_FIELD_PARAMS.get().is_none() 
                    {BLS48_581_FIELD_PARAMS.set(build_field_params(BLS48_581_PARAMS.modulo_as_strhex)).unwrap();}; 
              if FP_BLS48_581.get().is_none() 
                    {FP_BLS48_581.set(PrimeField::new(&BLS48_581_FIELD_PARAMS.get().unwrap())).unwrap();}                                                            
              if FR_BLS48_581.get().is_none() 
                    {FR_BLS48_581.set(PrimeField::new(&BLS48_581_R_FIELD_PARAMS.get().unwrap())).unwrap();}                                                                                                                           
              if BLS48_581_G1_CONSTS.get().is_none() 
                    {BLS48_581_G1_CONSTS.set(build_g1_params(BLS48_581_PARAMS,
                                                             &FP_BLS48_581.get().unwrap(),
                                                             &FR_BLS48_581.get().unwrap())).unwrap();}                                                                                                                                   
              G1_BLS48_581.set(G1Field { consts : BLS48_581_G1_CONSTS.get().unwrap(),
                                                  base_field: FP_BLS48_581.get().unwrap(),
                                                  fr_field : FR_BLS48_581.get().unwrap()}).unwrap() ;
            }
    G1_BLS48_581.get().unwrap()
}

pub fn g2_bls48_581() -> &'static G2Field<24,9,10,4>
{   if G2_BLS48_581.get().is_none() 
                {   if G1_BLS48_581.get().is_none() 
                            { let _ = g1_bls48_581();}
                    if BLS48_581_FROB_CONSTS.get().is_none() 
                            { BLS48_581_FROB_CONSTS.set(build_frob_params(BLS48_581_PARAMS, &FP_BLS48_581.get().unwrap())).unwrap()};                                         
                    if FP8_BLS48_581.get().is_none() 
                                        { let _ = fp8_bls48_581();}        
                    if BLS48_581_G2_CONSTS.get().is_none() 
                            {BLS48_581_G2_CONSTS.set(build_g2_params(BLS48_581_PARAMS,
                                                                     &ExtG2Field::Fp8_3(FP8_BLS48_581.get().unwrap()),
                                                                     &PrimeField::new(&BLS48_581_R_FIELD_PARAMS.get().unwrap()),
                                                                     BLS48_581_FROB_CONSTS.get().unwrap())).unwrap();}                                                                                                                                                                                                          
                    if EXFP8_BLS48_581.get().is_none()  {EXFP8_BLS48_581.set(ExtG2Field::Fp8_3(&FP8_BLS48_581.get().unwrap())).unwrap();}                                     
                    G2_BLS48_581.set(G2Field {  consts : BLS48_581_G2_CONSTS.get().unwrap(),
                                                base_field : &EXFP8_BLS48_581.get().unwrap() ,
                                                fr_field : FR_BLS48_581.get().unwrap()}).unwrap() ;
                }            
    G2_BLS48_581.get().unwrap()
}

pub fn fp48_bls48_581()-> &'static Fp48Field_3<24, 10>
{   if FP48_BLS48_581.get().is_none() {  if BLS48_581_FIELD_PARAMS.get().is_none() 
                                                        {   BLS48_581_FIELD_PARAMS.set(build_field_params(BLS48_581_PARAMS.modulo_as_strhex)).unwrap()
                                                        };
                                         if FP_BLS48_581.get().is_none() 
                                                        {   FP_BLS48_581.set(PrimeField::new(&BLS48_581_FIELD_PARAMS.get().unwrap())).unwrap()
                                                        };
                                         if BLS48_581_FROB_CONSTS.get().is_none() 
                                                        {   BLS48_581_FROB_CONSTS.set(build_frob_params(BLS48_581_PARAMS, &FP_BLS48_581.get().unwrap())).unwrap()
                                                        };               
                                         FP48_BLS48_581.set(Fp48Field_3::new(&FP_BLS48_581.get().unwrap(), Some(&BLS48_581_FROB_CONSTS.get().unwrap()))).unwrap();                                                                                                 
                                      }
    FP48_BLS48_581.get().unwrap()
}

pub fn gt_bls48_581() -> &'static GTField<10,24>
{  
        if FP48_BLS48_581.get().is_none() { let _ = fp48_bls48_581(); };
        GT_BLS48_581.set(GTField::Fp48_3(FP48_BLS48_581.get().unwrap())).unwrap();
        GT_BLS48_581.get().unwrap()       
}


static BLS48_571_FIELD_PARAMS       : OnceCell<FieldParams<9>>           = OnceCell::new();    
static BLS48_571_FROB_CONSTS        : OnceCell<ExFieldConsts<36,9>>      = OnceCell::new();
static BLS48_571_R_FIELD_PARAMS     : OnceCell<FieldParams<8>>           = OnceCell::new();    
static FP_BLS48_571                 : OnceCell<PrimeField<9>>            = OnceCell::new();  
static FP8_BLS48_571                : OnceCell<Fp8Field_2<36,9>>           = OnceCell::new();
static EXFP8_BLS48_571              : OnceCell<ExtG2Field<9,36>>         = OnceCell::new();
static BLS48_571_G1_CONSTS          : OnceCell<G1Consts<8,9,3>>          = OnceCell::new();
static BLS48_571_G2_CONSTS          : OnceCell<G2Consts<36,8,9,4>>       = OnceCell::new();
static FR_BLS48_571                 : OnceCell<PrimeField<8>>            = OnceCell::new();
static G2_BLS48_571                 : OnceCell<G2Field<36,8,9,4>>        = OnceCell::new();
static G1_BLS48_571                 : OnceCell<G1Field<8,9,3>>           = OnceCell::new();
static FP48_BLS48_571               : OnceCell<Fp48Field_2<36,9>>          = OnceCell::new();
static GT_BLS48_571                 : OnceCell<GTField<9,36>>           = OnceCell::new();

pub fn fr_bls48_571()-> &'static PrimeField<8>
{   if FR_BLS48_571.get().is_none() 
            { if BLS48_571_R_FIELD_PARAMS.get().is_none() 
                    {BLS48_571_R_FIELD_PARAMS.set(build_field_params(BLS48_571_PARAMS.r_as_strhex)).unwrap()}; 
              FR_BLS48_571.set(PrimeField::new(&BLS48_571_R_FIELD_PARAMS.get().unwrap())).unwrap();                                                                                                
            }
    FR_BLS48_571.get().unwrap()
}
pub fn fp8_bls48_571()-> &'static Fp8Field_2<36,9>
{   if FP8_BLS48_571.get().is_none() {  if BLS48_571_FIELD_PARAMS.get().is_none() 
                                                        {   BLS48_571_FIELD_PARAMS.set(build_field_params(BLS48_571_PARAMS.modulo_as_strhex)).unwrap()
                                                        };
                                         if FP_BLS48_571.get().is_none() 
                                                        {   FP_BLS48_571.set(PrimeField::new(&BLS48_571_FIELD_PARAMS.get().unwrap())).unwrap()
                                                        };
                                         if BLS48_571_FROB_CONSTS.get().is_none() 
                                                        {   BLS48_571_FROB_CONSTS.set(build_frob_params(BLS48_571_PARAMS, &FP_BLS48_571.get().unwrap())).unwrap()
                                                        };               
                                         FP8_BLS48_571.set(Fp8Field_2::new(&FP_BLS48_571.get().unwrap(), Some(&BLS48_571_FROB_CONSTS.get().unwrap()))).unwrap();                                                                                                 
                                      }
    FP8_BLS48_571.get().unwrap()
}

pub fn g1_bls48_571() -> &'static G1Field<8,9,3>
{   if G1_BLS48_571.get().is_none() 
            { if BLS48_571_R_FIELD_PARAMS.get().is_none()
                    {BLS48_571_R_FIELD_PARAMS.set(build_field_params(BLS48_571_PARAMS.r_as_strhex)).unwrap();}; 
              if BLS48_571_FIELD_PARAMS.get().is_none() 
                    {BLS48_571_FIELD_PARAMS.set(build_field_params(BLS48_571_PARAMS.modulo_as_strhex)).unwrap();}; 
              if FP_BLS48_571.get().is_none() 
                    {FP_BLS48_571.set(PrimeField::new(&BLS48_571_FIELD_PARAMS.get().unwrap())).unwrap();}                                                            
              if FR_BLS48_571.get().is_none() 
                    {FR_BLS48_571.set(PrimeField::new(&BLS48_571_R_FIELD_PARAMS.get().unwrap())).unwrap();}                                                                                                                           
              if BLS48_571_G1_CONSTS.get().is_none() 
                    {BLS48_571_G1_CONSTS.set(build_g1_params(BLS48_571_PARAMS,
                                                             &FP_BLS48_571.get().unwrap(),
                                                             &FR_BLS48_571.get().unwrap())).unwrap();}                                                                                                                                   
              G1_BLS48_571.set(G1Field { consts : BLS48_571_G1_CONSTS.get().unwrap(),
                                                  base_field: FP_BLS48_571.get().unwrap(),
                                                  fr_field : FR_BLS48_571.get().unwrap()}).unwrap() ;
            }
    G1_BLS48_571.get().unwrap()
}

pub fn g2_bls48_571() -> &'static G2Field<36,8,9,4>
{   if G2_BLS48_571.get().is_none() 
                {   if G1_BLS48_571.get().is_none() 
                            { let _ = g1_bls48_571();}
                    if BLS48_571_FROB_CONSTS.get().is_none() 
                            { BLS48_571_FROB_CONSTS.set(build_frob_params(BLS48_571_PARAMS, &FP_BLS48_571.get().unwrap())).unwrap()};                                         
                    if FP8_BLS48_571.get().is_none() 
                                        { let _ = fp8_bls48_571();}        
                    if BLS48_571_G2_CONSTS.get().is_none() 
                            {BLS48_571_G2_CONSTS.set(build_g2_params(BLS48_571_PARAMS,
                                                                     &ExtG2Field::Fp8_2(FP8_BLS48_571.get().unwrap()),
                                                                     &PrimeField::new(&BLS48_571_R_FIELD_PARAMS.get().unwrap()),
                                                                     BLS48_571_FROB_CONSTS.get().unwrap())).unwrap();}                                                                                                                                                                                                          
                    if EXFP8_BLS48_571.get().is_none()  {EXFP8_BLS48_571.set(ExtG2Field::Fp8_2(&FP8_BLS48_571.get().unwrap())).unwrap();}                                     
                    G2_BLS48_571.set(G2Field {  consts : BLS48_571_G2_CONSTS.get().unwrap(),
                                                base_field : &EXFP8_BLS48_571.get().unwrap() ,
                                                fr_field : FR_BLS48_571.get().unwrap()}).unwrap() ;
                }            
    G2_BLS48_571.get().unwrap()
}

pub fn fp48_bls48_571()-> &'static Fp48Field_2<36, 9>
{   if FP48_BLS48_571.get().is_none() {  if BLS48_571_FIELD_PARAMS.get().is_none() 
                                                        {   BLS48_571_FIELD_PARAMS.set(build_field_params(BLS48_571_PARAMS.modulo_as_strhex)).unwrap()
                                                        };
                                         if FP_BLS48_571.get().is_none() 
                                                        {   FP_BLS48_571.set(PrimeField::new(&BLS48_571_FIELD_PARAMS.get().unwrap())).unwrap()
                                                        };
                                         if BLS48_571_FROB_CONSTS.get().is_none() 
                                                        {   BLS48_571_FROB_CONSTS.set(build_frob_params(BLS48_571_PARAMS, &FP_BLS48_571.get().unwrap())).unwrap()
                                                        };               
                                         FP48_BLS48_571.set(Fp48Field_2::new(&FP_BLS48_571.get().unwrap(), Some(&BLS48_571_FROB_CONSTS.get().unwrap()))).unwrap();                                                                                                 
                                      }
    FP48_BLS48_571.get().unwrap()
}

pub fn gt_bls48_571() -> &'static GTField<9,36>
{  
        if FP48_BLS48_571.get().is_none() { let _ = fp48_bls48_571(); };
        GT_BLS48_571.set(GTField::Fp48_2(FP48_BLS48_571.get().unwrap())).unwrap();
        GT_BLS48_571.get().unwrap()       
}

static BLS48_287_FIELD_PARAMS       : OnceCell<FieldParams<5>>                      = OnceCell::new();    
static BLS48_287_FROB_CONSTS        : OnceCell<ExFieldConsts<35,5>>      = OnceCell::new();
static BLS48_287_R_FIELD_PARAMS     : OnceCell<FieldParams<4>>                      = OnceCell::new();    
static FP_BLS48_287                 : OnceCell<PrimeField<5>>                       = OnceCell::new();  
static FP8_BLS48_287                : OnceCell<Fp8Field_2<35,5>>          = OnceCell::new();
static EXFP8_BLS48_287              : OnceCell<ExtG2Field<5,35>>          = OnceCell::new();
static BLS48_287_G1_CONSTS          : OnceCell<G1Consts<4,5,4>>  = OnceCell::new();
static BLS48_287_G2_CONSTS          : OnceCell<G2Consts<35,4,5,19>> = OnceCell::new();
static FR_BLS48_287                 : OnceCell<PrimeField<4>>                                     = OnceCell::new();
static G2_BLS48_287                 : OnceCell<G2Field<35,4,5,19>>  = OnceCell::new();
static G1_BLS48_287                 : OnceCell<G1Field<4,5,4>> = OnceCell::new();
static FP48_BLS48_287               : OnceCell<Fp48Field_2<35,5>>       = OnceCell::new();
static GT_BLS48_287                 : OnceCell<GTField<5,35>>           = OnceCell::new();

pub fn fr_bls48_287()-> &'static PrimeField<4>
{   if FR_BLS48_287.get().is_none() 
            { if BLS48_287_R_FIELD_PARAMS.get().is_none() 
                    {BLS48_287_R_FIELD_PARAMS.set(build_field_params(BLS48_287_PARAMS.r_as_strhex)).unwrap()}; 
              FR_BLS48_287.set(PrimeField::new(&BLS48_287_R_FIELD_PARAMS.get().unwrap())).unwrap();                                                                                                
            }
    FR_BLS48_287.get().unwrap()
}
pub fn fp8_bls48_287()-> &'static Fp8Field_2<35,5>
{   if FP8_BLS48_287.get().is_none() {  if BLS48_287_FIELD_PARAMS.get().is_none() 
                                                        {   BLS48_287_FIELD_PARAMS.set(build_field_params(BLS48_287_PARAMS.modulo_as_strhex)).unwrap()
                                                        };
                                         if FP_BLS48_287.get().is_none() 
                                                        {   FP_BLS48_287.set(PrimeField::new(&BLS48_287_FIELD_PARAMS.get().unwrap())).unwrap()
                                                        };
                                         if BLS48_287_FROB_CONSTS.get().is_none() 
                                                        {   BLS48_287_FROB_CONSTS.set(build_frob_params(BLS48_287_PARAMS, &FP_BLS48_287.get().unwrap())).unwrap()
                                                        };               
                                         FP8_BLS48_287.set(Fp8Field_2::new(&FP_BLS48_287.get().unwrap(), Some(&BLS48_287_FROB_CONSTS.get().unwrap()))).unwrap();                                                                                                 
                                      }
    FP8_BLS48_287.get().unwrap()
}

pub fn g1_bls48_287() -> &'static G1Field<4,5,4>
{   if G1_BLS48_287.get().is_none() 
            { if BLS48_287_R_FIELD_PARAMS.get().is_none()
                    {BLS48_287_R_FIELD_PARAMS.set(build_field_params(BLS48_287_PARAMS.r_as_strhex)).unwrap();}; 
              if BLS48_287_FIELD_PARAMS.get().is_none() 
                    {BLS48_287_FIELD_PARAMS.set(build_field_params(BLS48_287_PARAMS.modulo_as_strhex)).unwrap();}; 
              if FP_BLS48_287.get().is_none() 
                    {FP_BLS48_287.set(PrimeField::new(&BLS48_287_FIELD_PARAMS.get().unwrap())).unwrap();}                                                            
              if FR_BLS48_287.get().is_none() 
                    {FR_BLS48_287.set(PrimeField::new(&BLS48_287_R_FIELD_PARAMS.get().unwrap())).unwrap();}                                                                                                                           
              if BLS48_287_G1_CONSTS.get().is_none() 
                    {BLS48_287_G1_CONSTS.set(build_g1_params(BLS48_287_PARAMS,
                                                             &FP_BLS48_287.get().unwrap(),
                                                             &FR_BLS48_287.get().unwrap())).unwrap();}                                                                                                                                   
              G1_BLS48_287.set(G1Field { consts : BLS48_287_G1_CONSTS.get().unwrap(),
                                                  base_field: FP_BLS48_287.get().unwrap(),
                                                  fr_field : FR_BLS48_287.get().unwrap()}).unwrap() ;
            }
    G1_BLS48_287.get().unwrap()
}

pub fn g2_bls48_287() -> &'static G2Field<35,4,5,19>
{   if G2_BLS48_287.get().is_none() 
                {   if G1_BLS48_287.get().is_none() 
                            { let _ = g1_bls48_287();}
                    if BLS48_287_FROB_CONSTS.get().is_none() 
                            { BLS48_287_FROB_CONSTS.set(build_frob_params(BLS48_287_PARAMS, &FP_BLS48_287.get().unwrap())).unwrap()};                                         
                    if FP8_BLS48_287.get().is_none() 
                                        { let _ = fp8_bls48_287();}        
                    if BLS48_287_G2_CONSTS.get().is_none() 
                            {BLS48_287_G2_CONSTS.set(build_g2_params(BLS48_287_PARAMS,
                                                                     &ExtG2Field::Fp8_2(FP8_BLS48_287.get().unwrap()),
                                                                     &PrimeField::new(&BLS48_287_R_FIELD_PARAMS.get().unwrap()),
                                                                     BLS48_287_FROB_CONSTS.get().unwrap())).unwrap();}                                                                                                                                                                                                          
                    if EXFP8_BLS48_287.get().is_none()  {EXFP8_BLS48_287.set(ExtG2Field::Fp8_2(&FP8_BLS48_287.get().unwrap())).unwrap();}                                     
                    G2_BLS48_287.set(G2Field {  consts : BLS48_287_G2_CONSTS.get().unwrap(),
                                                base_field : &EXFP8_BLS48_287.get().unwrap() ,
                                                fr_field : FR_BLS48_287.get().unwrap()}).unwrap() ;
                }            
    G2_BLS48_287.get().unwrap()
}

pub fn fp48_bls48_287()-> &'static Fp48Field_2<35, 5>
{   if FP48_BLS48_287.get().is_none() {  if BLS48_287_FIELD_PARAMS.get().is_none() 
                                                        {   BLS48_287_FIELD_PARAMS.set(build_field_params(BLS48_287_PARAMS.modulo_as_strhex)).unwrap()
                                                        };
                                         if FP_BLS48_287.get().is_none() 
                                                        {   FP_BLS48_287.set(PrimeField::new(&BLS48_287_FIELD_PARAMS.get().unwrap())).unwrap()
                                                        };
                                         if BLS48_287_FROB_CONSTS.get().is_none() 
                                                        {   BLS48_287_FROB_CONSTS.set(build_frob_params(BLS48_287_PARAMS, &FP_BLS48_287.get().unwrap())).unwrap()
                                                        };               
                                         FP48_BLS48_287.set(Fp48Field_2::new(&FP_BLS48_287.get().unwrap(), Some(&BLS48_287_FROB_CONSTS.get().unwrap()))).unwrap();                                                                                                 
                                      }
    FP48_BLS48_287.get().unwrap()
}

pub fn gt_bls48_287() -> &'static GTField<5,35>
{  
        if FP48_BLS48_287.get().is_none() { let _ = fp48_bls48_287(); };
        GT_BLS48_287.set(GTField::Fp48_2(FP48_BLS48_287.get().unwrap())).unwrap();
        GT_BLS48_287.get().unwrap()       
}

static BLS48_277_FIELD_PARAMS       : OnceCell<FieldParams<5>>                      = OnceCell::new();    
static BLS48_277_FROB_CONSTS        : OnceCell<ExFieldConsts<24,5>>      = OnceCell::new();
static BLS48_277_R_FIELD_PARAMS     : OnceCell<FieldParams<4>>                      = OnceCell::new();    
static FP_BLS48_277                 : OnceCell<PrimeField<5>>                       = OnceCell::new();  
static FP8_BLS48_277                : OnceCell<Fp8Field_3<24,5>>          = OnceCell::new();
static EXFP8_BLS48_277              : OnceCell<ExtG2Field<5,24>>          = OnceCell::new();
static BLS48_277_G1_CONSTS          : OnceCell<G1Consts<4,5,10>>  = OnceCell::new();
static BLS48_277_G2_CONSTS          : OnceCell<G2Consts<24,4,5,4>> = OnceCell::new();
static FR_BLS48_277                 : OnceCell<PrimeField<4>>                                     = OnceCell::new();
static G2_BLS48_277                 : OnceCell<G2Field<24,4,5,4>>  = OnceCell::new();
static G1_BLS48_277                 : OnceCell<G1Field<4,5,10>> = OnceCell::new();
static FP48_BLS48_277               : OnceCell<Fp48Field_3<24,5>>       = OnceCell::new();
static GT_BLS48_277                 : OnceCell<GTField<5,24>>           = OnceCell::new();

pub fn fr_bls48_277()-> &'static PrimeField<4>
{   if FR_BLS48_277.get().is_none() 
            { if BLS48_277_R_FIELD_PARAMS.get().is_none() 
                    {BLS48_277_R_FIELD_PARAMS.set(build_field_params(BLS48_277_PARAMS.r_as_strhex)).unwrap()}; 
              FR_BLS48_277.set(PrimeField::new(&BLS48_277_R_FIELD_PARAMS.get().unwrap())).unwrap();                                                                                                
            }
    FR_BLS48_277.get().unwrap()
}
pub fn fp8_bls48_277()-> &'static Fp8Field_3<24,5>
{   if FP8_BLS48_277.get().is_none() {  if BLS48_277_FIELD_PARAMS.get().is_none() 
                                                        {   BLS48_277_FIELD_PARAMS.set(build_field_params(BLS48_277_PARAMS.modulo_as_strhex)).unwrap()
                                                        };
                                         if FP_BLS48_277.get().is_none() 
                                                        {   FP_BLS48_277.set(PrimeField::new(&BLS48_277_FIELD_PARAMS.get().unwrap())).unwrap()
                                                        };
                                         if BLS48_277_FROB_CONSTS.get().is_none() 
                                                        {   BLS48_277_FROB_CONSTS.set(build_frob_params(BLS48_277_PARAMS, &FP_BLS48_277.get().unwrap())).unwrap()
                                                        };               
                                         FP8_BLS48_277.set(Fp8Field_3::new(&FP_BLS48_277.get().unwrap(), Some(&BLS48_277_FROB_CONSTS.get().unwrap()))).unwrap();                                                                                                 
                                      }
    FP8_BLS48_277.get().unwrap()
}

pub fn g1_bls48_277() -> &'static G1Field<4,5,10>
{   if G1_BLS48_277.get().is_none() 
            { if BLS48_277_R_FIELD_PARAMS.get().is_none()
                    {BLS48_277_R_FIELD_PARAMS.set(build_field_params(BLS48_277_PARAMS.r_as_strhex)).unwrap();}; 
              if BLS48_277_FIELD_PARAMS.get().is_none() 
                    {BLS48_277_FIELD_PARAMS.set(build_field_params(BLS48_277_PARAMS.modulo_as_strhex)).unwrap();}; 
              if FP_BLS48_277.get().is_none() 
                    {FP_BLS48_277.set(PrimeField::new(&BLS48_277_FIELD_PARAMS.get().unwrap())).unwrap();}                                                            
              if FR_BLS48_277.get().is_none() 
                    {FR_BLS48_277.set(PrimeField::new(&BLS48_277_R_FIELD_PARAMS.get().unwrap())).unwrap();}                                                                                                                           
              if BLS48_277_G1_CONSTS.get().is_none() 
                    {BLS48_277_G1_CONSTS.set(build_g1_params(BLS48_277_PARAMS,
                                                             &FP_BLS48_277.get().unwrap(),
                                                             &FR_BLS48_277.get().unwrap())).unwrap();}                                                                                                                                   
              G1_BLS48_277.set(G1Field { consts : BLS48_277_G1_CONSTS.get().unwrap(),
                                                  base_field: FP_BLS48_277.get().unwrap(),
                                                  fr_field : FR_BLS48_277.get().unwrap()}).unwrap() ;
            }
    G1_BLS48_277.get().unwrap()
}

pub fn g2_bls48_277() -> &'static G2Field<24,4,5,4>
{   if G2_BLS48_277.get().is_none() 
                {   if G1_BLS48_277.get().is_none() 
                            { let _ = g1_bls48_277();}
                    if BLS48_277_FROB_CONSTS.get().is_none() 
                            { BLS48_277_FROB_CONSTS.set(build_frob_params(BLS48_277_PARAMS, &FP_BLS48_277.get().unwrap())).unwrap()};                                         
                    if FP8_BLS48_277.get().is_none() 
                                        { let _ = fp8_bls48_277();}        
                    if BLS48_277_G2_CONSTS.get().is_none() 
                            {BLS48_277_G2_CONSTS.set(build_g2_params(BLS48_277_PARAMS,
                                                                     &ExtG2Field::Fp8_3(FP8_BLS48_277.get().unwrap()),
                                                                     &PrimeField::new(&BLS48_277_R_FIELD_PARAMS.get().unwrap()),
                                                                     BLS48_277_FROB_CONSTS.get().unwrap())).unwrap();}                                                                                                                                                                                                          
                    if EXFP8_BLS48_277.get().is_none()  {EXFP8_BLS48_277.set(ExtG2Field::Fp8_3(&FP8_BLS48_277.get().unwrap())).unwrap();}                                     
                    G2_BLS48_277.set(G2Field {  consts : BLS48_277_G2_CONSTS.get().unwrap(),
                                                base_field : &EXFP8_BLS48_277.get().unwrap() ,
                                                fr_field : FR_BLS48_277.get().unwrap()}).unwrap() ;
                }            
    G2_BLS48_277.get().unwrap()
}

pub fn fp48_bls48_277()-> &'static Fp48Field_3<24, 5>
{   if FP48_BLS48_277.get().is_none() {  if BLS48_277_FIELD_PARAMS.get().is_none() 
                                                        {   BLS48_277_FIELD_PARAMS.set(build_field_params(BLS48_277_PARAMS.modulo_as_strhex)).unwrap()
                                                        };
                                         if FP_BLS48_277.get().is_none() 
                                                        {   FP_BLS48_277.set(PrimeField::new(&BLS48_277_FIELD_PARAMS.get().unwrap())).unwrap()
                                                        };
                                         if BLS48_277_FROB_CONSTS.get().is_none() 
                                                        {   BLS48_277_FROB_CONSTS.set(build_frob_params(BLS48_277_PARAMS, &FP_BLS48_277.get().unwrap())).unwrap()
                                                        };               
                                         FP48_BLS48_277.set(Fp48Field_3::new(&FP_BLS48_277.get().unwrap(), Some(&BLS48_277_FROB_CONSTS.get().unwrap()))).unwrap();                                                                                                 
                                      }
    FP48_BLS48_277.get().unwrap()
}

pub fn gt_bls48_277() -> &'static GTField<5,24>
{  
        if FP48_BLS48_277.get().is_none() { let _ = fp48_bls48_277(); };
        GT_BLS48_277.set(GTField::Fp48_3(FP48_BLS48_277.get().unwrap())).unwrap();
        GT_BLS48_277.get().unwrap()       }