// Code developed by FARAOUN Kamel Mohamed.
// EEDIS-Laboratory. UDL-University. Algeria
// During May 2024.

use crate::curves::g1::G1Consts;
use crate::curves::g2::*;
use crate::curves::gt::GTField;
use crate::extensions::ext_fields::{ExFieldConsts, ExtField};
use crate::fields::prime_fields::{FieldParams, PrimeField};
use crate::{curves::g1::G1Field, extensions::towering1::fp24::Fp24Field as Fp24Field_1,extensions::towering2::fp24::Fp24Field as Fp24Field_2};
use crate::extensions::towering1::fp4::Fp4Field as Fp4Field_1;
use crate::extensions::g2_extfields::ExtG2Field;
use crate::extensions::towering2::fp4::Fp4Field as Fp4Field_2;
use crate::parameters::builders::*;
use crate::parameters::paramlist::{BLS24_315_PARAMS, BLS24_479_PARAMS,BLS24_477_PARAMS, BLS24_509_PARAMS, BLS24_509_SNARK_PARAMS, BLS24_559_PARAMS};
use once_cell::sync::OnceCell;

static BLS24_477_FIELD_PARAMS   : OnceCell<FieldParams<8>> = OnceCell::new();    
static BLS24_477_FROB_CONSTS    : OnceCell<ExFieldConsts<11,8>> = OnceCell::new();
static BLS24_477_R_FIELD_PARAMS : OnceCell<FieldParams<6>> = OnceCell::new();    
static FP_BLS24_477             : OnceCell<PrimeField<8>> = OnceCell::new();
static FP4_BLS24_477            : OnceCell<Fp4Field_1<11,8>> = OnceCell::new();
static EXFP4_BLS24_477          : OnceCell<ExtG2Field<8,11>> = OnceCell::new();
static BLS24_477_G1_CONSTS      : OnceCell<G1Consts<6,8,7>> = OnceCell::new();
static BLS24_477_G2_CONSTS      : OnceCell<G2Consts<11,6,8,4>> = OnceCell::new();
static FR_BLS24_477             : OnceCell<PrimeField<6>> = OnceCell::new();
static G1_BLS24_477             : OnceCell<G1Field<6,8,7>> = OnceCell::new();
static G2_BLS24_477             : OnceCell<G2Field<11,6,8,4>> = OnceCell::new();
static FP24_BLS24_477           : OnceCell<Fp24Field_1<11,8>> = OnceCell::new();
static GT_BLS24_477             : OnceCell<GTField<8,11>> = OnceCell::new();

pub fn fr_bls24_477()-> &'static PrimeField<6>
{   if FR_BLS24_477.get().is_none() 
            { if BLS24_477_R_FIELD_PARAMS.get().is_none() 
                    {BLS24_477_R_FIELD_PARAMS.set(build_field_params(BLS24_477_PARAMS.r_as_strhex)).unwrap()}; 
              FR_BLS24_477.set(PrimeField::new(&BLS24_477_R_FIELD_PARAMS.get().unwrap())).unwrap();                                                                                                
            }
    FR_BLS24_477.get().unwrap()
}
pub fn fp4_bls24_477()-> &'static Fp4Field_1<11, 8>
{   if FP4_BLS24_477.get().is_none() {  if BLS24_477_FIELD_PARAMS.get().is_none() 
                                            {   BLS24_477_FIELD_PARAMS.set(build_field_params(BLS24_477_PARAMS.modulo_as_strhex)).unwrap()};
                                         if FP_BLS24_477.get().is_none() 
                                            {   FP_BLS24_477.set(PrimeField::new(&BLS24_477_FIELD_PARAMS.get().unwrap())).unwrap()};
                                         if BLS24_477_FROB_CONSTS.get().is_none() 
                                            {   BLS24_477_FROB_CONSTS.set(build_frob_params(BLS24_477_PARAMS, &FP_BLS24_477.get().unwrap())).unwrap()};               
                                         FP4_BLS24_477.set(Fp4Field_1::new( &FP_BLS24_477.get().unwrap(), Some(&BLS24_477_FROB_CONSTS.get().unwrap()))).unwrap();                                                                                                 
                                      }
    FP4_BLS24_477.get().unwrap()
}

pub fn g1_bls24_477() -> &'static G1Field<6,8,7>
{   if G1_BLS24_477.get().is_none() 
            { if BLS24_477_R_FIELD_PARAMS.get().is_none()
                    {BLS24_477_R_FIELD_PARAMS.set(build_field_params(BLS24_477_PARAMS.r_as_strhex)).unwrap()}; 
              if BLS24_477_FIELD_PARAMS.get().is_none() 
                    {BLS24_477_FIELD_PARAMS.set(build_field_params(BLS24_477_PARAMS.modulo_as_strhex)).unwrap()}; 
              if FP_BLS24_477.get().is_none() 
                    {FP_BLS24_477.set(PrimeField::new(&BLS24_477_FIELD_PARAMS.get().unwrap())).unwrap();}                                                            
              if FR_BLS24_477.get().is_none() 
                    {FR_BLS24_477.set(PrimeField::new(&BLS24_477_R_FIELD_PARAMS.get().unwrap())).unwrap();}                                                                                                                           
              if BLS24_477_G1_CONSTS.get().is_none() 
                    {BLS24_477_G1_CONSTS.set(build_g1_params(BLS24_477_PARAMS,
                                                        &FP_BLS24_477.get().unwrap(),
                                                           &FR_BLS24_477.get().unwrap())).unwrap();}                                                                                                                                   
              G1_BLS24_477.set(G1Field { consts : BLS24_477_G1_CONSTS.get().unwrap(),
                                                  base_field: FP_BLS24_477.get().unwrap(),
                                                  fr_field : FR_BLS24_477.get().unwrap()}).unwrap() ;
            }
    G1_BLS24_477.get().unwrap()
}

pub fn g2_bls24_477() -> &'static G2Field<11,6,8,4>
{   if G2_BLS24_477.get().is_none() 
                {   if G1_BLS24_477.get().is_none() 
                            { let _ = g1_bls24_477();}
                    if BLS24_477_FROB_CONSTS.get().is_none() 
                            { BLS24_477_FROB_CONSTS.set(build_frob_params(BLS24_477_PARAMS, &FP_BLS24_477.get().unwrap())).unwrap()}; 
                    if FP4_BLS24_477.get().is_none() 
                            { let _ = fp4_bls24_477();}                                                            
                    if BLS24_477_G2_CONSTS.get().is_none() 
                            {BLS24_477_G2_CONSTS.set(build_g2_params(BLS24_477_PARAMS,
                                                                     &ExtG2Field::Fp4_1(FP4_BLS24_477.get().unwrap()),
                                                                     &PrimeField::new(&BLS24_477_R_FIELD_PARAMS.get().unwrap()),
                                                                     BLS24_477_FROB_CONSTS.get().unwrap())).unwrap();}                                                                                                                                   
                    if EXFP4_BLS24_477.get().is_none()  {EXFP4_BLS24_477.set(ExtG2Field::Fp4_1(&FP4_BLS24_477.get().unwrap())).unwrap();}                                     
                    G2_BLS24_477.set(G2Field {  consts : BLS24_477_G2_CONSTS.get().unwrap(),
                                                base_field : &EXFP4_BLS24_477.get().unwrap() ,
                                                fr_field : FR_BLS24_477.get().unwrap()}).unwrap() ;
                }            
    G2_BLS24_477.get().unwrap()
}

pub fn fp24_bls24_477()-> &'static Fp24Field_1<11, 8>
{   if FP24_BLS24_477.get().is_none() {  if BLS24_477_FIELD_PARAMS.get().is_none() 
                                            {   BLS24_477_FIELD_PARAMS.set(build_field_params(BLS24_477_PARAMS.modulo_as_strhex)).unwrap()};
                                         if FP_BLS24_477.get().is_none() 
                                            {   FP_BLS24_477.set(PrimeField::new(&BLS24_477_FIELD_PARAMS.get().unwrap())).unwrap()};
                                         if BLS24_477_FROB_CONSTS.get().is_none() 
                                            {   BLS24_477_FROB_CONSTS.set(build_frob_params(BLS24_477_PARAMS, &FP_BLS24_477.get().unwrap())).unwrap()};               
                                         FP24_BLS24_477.set(Fp24Field_1::new(&FP_BLS24_477.get().unwrap(), Some(&BLS24_477_FROB_CONSTS.get().unwrap()))).unwrap();                                                                                                 
                                      }
    FP24_BLS24_477.get().unwrap()
}

pub fn gt_bls24_477() -> &'static GTField<8,11>
{  
        if FP24_BLS24_477.get().is_none() { let _ = fp24_bls24_477(); };
        GT_BLS24_477.set(GTField::Fp24_1(FP24_BLS24_477.get().unwrap())).unwrap();
        GT_BLS24_477.get().unwrap()       
}

static BLS24_479_FIELD_PARAMS   : OnceCell<FieldParams<8>> = OnceCell::new();    
static BLS24_479_FROB_CONSTS    : OnceCell<ExFieldConsts<10,8>> = OnceCell::new();
static BLS24_479_R_FIELD_PARAMS : OnceCell<FieldParams<7>> = OnceCell::new();    
static FP_BLS24_479             : OnceCell<PrimeField<8>> = OnceCell::new();
static FP4_BLS24_479            : OnceCell<Fp4Field_1<10,8>> = OnceCell::new();
static EXFP4_BLS24_479          : OnceCell<ExtG2Field<8,10>> = OnceCell::new();
static BLS24_479_G1_CONSTS      : OnceCell<G1Consts<7,8,4>> = OnceCell::new();
static BLS24_479_G2_CONSTS      : OnceCell<G2Consts<10,7,8,10>> = OnceCell::new();
static FR_BLS24_479             : OnceCell<PrimeField<7>> = OnceCell::new();
static G1_BLS24_479             : OnceCell<G1Field<7,8,4>> = OnceCell::new();
static G2_BLS24_479             : OnceCell<G2Field<10,7,8,10>> = OnceCell::new();
static FP24_BLS24_479           : OnceCell<Fp24Field_1<10,8>> = OnceCell::new();
static GT_BLS24_479             : OnceCell<GTField<8,10>> = OnceCell::new();

pub fn fr_bls24_479()-> &'static PrimeField<7>
{   if FR_BLS24_479.get().is_none() 
            { if BLS24_479_R_FIELD_PARAMS.get().is_none() 
                    {BLS24_479_R_FIELD_PARAMS.set(build_field_params(BLS24_479_PARAMS.r_as_strhex)).unwrap()}; 
              FR_BLS24_479.set(PrimeField::new(&BLS24_479_R_FIELD_PARAMS.get().unwrap())).unwrap();                                                                                                
            }
    FR_BLS24_479.get().unwrap()
}
pub fn fp4_bls24_479()-> &'static Fp4Field_1<10, 8>
{   if FP4_BLS24_479.get().is_none() {  if BLS24_479_FIELD_PARAMS.get().is_none() 
                                            {   BLS24_479_FIELD_PARAMS.set(build_field_params(BLS24_479_PARAMS.modulo_as_strhex)).unwrap()};
                                         if FP_BLS24_479.get().is_none() 
                                            {   FP_BLS24_479.set(PrimeField::new(&BLS24_479_FIELD_PARAMS.get().unwrap())).unwrap()};
                                         if BLS24_479_FROB_CONSTS.get().is_none() 
                                            {   BLS24_479_FROB_CONSTS.set(build_frob_params(BLS24_479_PARAMS, &FP_BLS24_479.get().unwrap())).unwrap()};               
                                         FP4_BLS24_479.set(Fp4Field_1::new( &FP_BLS24_479.get().unwrap(), Some(&BLS24_479_FROB_CONSTS.get().unwrap()))).unwrap();                                                                                                 
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
                                                                     &ExtG2Field::Fp4_1(FP4_BLS24_479.get().unwrap()),
                                                                     &PrimeField::new(&BLS24_479_R_FIELD_PARAMS.get().unwrap()),
                                                                     BLS24_479_FROB_CONSTS.get().unwrap())).unwrap();}                                                                                                                                   
                    if EXFP4_BLS24_479.get().is_none()  {EXFP4_BLS24_479.set(ExtG2Field::Fp4_1(&FP4_BLS24_479.get().unwrap())).unwrap();}                                     
                    G2_BLS24_479.set(G2Field {  consts : BLS24_479_G2_CONSTS.get().unwrap(),
                                                base_field : &EXFP4_BLS24_479.get().unwrap() ,
                                                fr_field : FR_BLS24_479.get().unwrap()}).unwrap() ;
                }            
    G2_BLS24_479.get().unwrap()
}

pub fn fp24_bls24_479()-> &'static Fp24Field_1<10, 8>
{   if FP24_BLS24_479.get().is_none() {  if BLS24_479_FIELD_PARAMS.get().is_none() 
                                            {   BLS24_479_FIELD_PARAMS.set(build_field_params(BLS24_479_PARAMS.modulo_as_strhex)).unwrap()};
                                         if FP_BLS24_479.get().is_none() 
                                            {   FP_BLS24_479.set(PrimeField::new(&BLS24_479_FIELD_PARAMS.get().unwrap())).unwrap()};
                                         if BLS24_479_FROB_CONSTS.get().is_none() 
                                            {   BLS24_479_FROB_CONSTS.set(build_frob_params(BLS24_479_PARAMS, &FP_BLS24_479.get().unwrap())).unwrap()};               
                                         FP24_BLS24_479.set(Fp24Field_1::new(&FP_BLS24_479.get().unwrap(), Some(&BLS24_479_FROB_CONSTS.get().unwrap()))).unwrap();                                                                                                 
                                      }
    FP24_BLS24_479.get().unwrap()
}

pub fn gt_bls24_479() -> &'static GTField<8,10>
{  
        if FP24_BLS24_479.get().is_none() { let _ = fp24_bls24_479(); };
        GT_BLS24_479.set(GTField::Fp24_1(FP24_BLS24_479.get().unwrap())).unwrap();
        GT_BLS24_479.get().unwrap()       
}

static BLS24_559_FIELD_PARAMS       : OnceCell<FieldParams<9>> = OnceCell::new();    
static BLS24_559_FROB_CONSTS        : OnceCell<ExFieldConsts<10,9>> = OnceCell::new();
static BLS24_559_R_FIELD_PARAMS     : OnceCell<FieldParams<8>> = OnceCell::new();   
static FP_BLS24_559                 : OnceCell<PrimeField<9>> = OnceCell::new(); 
static FP4_BLS24_559                : OnceCell<Fp4Field_1<10,9>> = OnceCell::new();
static EXFP4_BLS24_559              : OnceCell<ExtG2Field<9,10>> = OnceCell::new();
static FP24_BLS24_559               : OnceCell<Fp24Field_1<10,9>> = OnceCell::new();
static BLS24_559_G1_CONSTS          : OnceCell<G1Consts<8,9,4>> = OnceCell::new();
static BLS24_559_G2_CONSTS          : OnceCell<G2Consts<10,8,9,4>> = OnceCell::new();
static FR_BLS24_559                 : OnceCell<PrimeField<8>> = OnceCell::new();
static G1_BLS24_559                 : OnceCell<G1Field<8,9,4>> = OnceCell::new();
static G2_BLS24_559                 : OnceCell<G2Field<10,8,9,4>> = OnceCell::new();
static GT_BLS24_559                 : OnceCell<GTField<9,10>> = OnceCell::new();


pub fn fr_bls24_559()-> &'static PrimeField<8>
{   if FR_BLS24_559.get().is_none() 
            { if BLS24_559_R_FIELD_PARAMS.get().is_none() 
                    {BLS24_559_R_FIELD_PARAMS.set(build_field_params(BLS24_559_PARAMS.r_as_strhex)).unwrap()}; 
              FR_BLS24_559.set(PrimeField::new(&BLS24_559_R_FIELD_PARAMS.get().unwrap())).unwrap();                                                                                                
            }
    FR_BLS24_559.get().unwrap()
}

pub fn fp4_bls24_559()-> &'static Fp4Field_1<10, 9>
{   if FP4_BLS24_559.get().is_none() {  if BLS24_559_FIELD_PARAMS.get().is_none() 
                                            {   BLS24_559_FIELD_PARAMS.set(build_field_params(BLS24_559_PARAMS.modulo_as_strhex)).unwrap()};
                                         if FP_BLS24_559.get().is_none() 
                                            {   FP_BLS24_559.set(PrimeField::new(&BLS24_559_FIELD_PARAMS.get().unwrap())).unwrap()};
                                         if BLS24_559_FROB_CONSTS.get().is_none() 
                                            {   BLS24_559_FROB_CONSTS.set(build_frob_params(BLS24_559_PARAMS, &FP_BLS24_559.get().unwrap())).unwrap()};               
                                         FP4_BLS24_559.set(Fp4Field_1::new(&FP_BLS24_559.get().unwrap(), Some(&BLS24_559_FROB_CONSTS.get().unwrap()))).unwrap();                                                                                                 
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
                                                                     &ExtG2Field::Fp4_1(FP4_BLS24_559.get().unwrap()),
                                                                     &PrimeField::new(&BLS24_559_R_FIELD_PARAMS.get().unwrap()),
                                                                     BLS24_559_FROB_CONSTS.get().unwrap())).unwrap();}                                                                                                                                   
                    if EXFP4_BLS24_559.get().is_none()  {EXFP4_BLS24_559.set(ExtG2Field::Fp4_1(&FP4_BLS24_559.get().unwrap())).unwrap();}                                     
                    G2_BLS24_559.set(G2Field {  consts : BLS24_559_G2_CONSTS.get().unwrap(),
                                                base_field : &EXFP4_BLS24_559.get().unwrap() ,
                                                        fr_field : FR_BLS24_559.get().unwrap()}).unwrap() ;
                }            
    G2_BLS24_559.get().unwrap()
}

pub fn fp24_bls24_559()-> &'static Fp24Field_1<10, 9>
{   if FP24_BLS24_559.get().is_none() {  if BLS24_559_FIELD_PARAMS.get().is_none() 
                                            {   BLS24_559_FIELD_PARAMS.set(build_field_params(BLS24_559_PARAMS.modulo_as_strhex)).unwrap()};
                                         if FP_BLS24_559.get().is_none() 
                                            {   FP_BLS24_559.set(PrimeField::new(&BLS24_559_FIELD_PARAMS.get().unwrap())).unwrap()};
                                         if BLS24_559_FROB_CONSTS.get().is_none() 
                                            {   BLS24_559_FROB_CONSTS.set(build_frob_params(BLS24_559_PARAMS, &FP_BLS24_559.get().unwrap())).unwrap()};               
                                         FP24_BLS24_559.set(Fp24Field_1::new(&FP_BLS24_559.get().unwrap(), Some(&BLS24_559_FROB_CONSTS.get().unwrap()))).unwrap();                                                                                                 
                                      }
    FP24_BLS24_559.get().unwrap()
}

pub fn gt_bls24_559() -> &'static GTField<9,10>
{  
        if FP24_BLS24_559.get().is_none() { let _ = fp24_bls24_559(); };
        GT_BLS24_559.set(GTField::Fp24_1(FP24_BLS24_559.get().unwrap())).unwrap();
        GT_BLS24_559.get().unwrap()       
}

static BLS24_315_FIELD_PARAMS       : OnceCell<FieldParams<5>> = OnceCell::new();    
static BLS24_315_FROB_CONSTS        : OnceCell<ExFieldConsts<11,5>> = OnceCell::new();
static BLS24_315_R_FIELD_PARAMS     : OnceCell<FieldParams<4>> = OnceCell::new();   
static FP_BLS24_315                 : OnceCell<PrimeField<5>> = OnceCell::new(); 
static FP4_BLS24_315                : OnceCell<Fp4Field_2<11,5>> = OnceCell::new();
static EXFP4_BLS24_315              : OnceCell<ExtG2Field<5,11>> = OnceCell::new();
static FP24_BLS24_315               : OnceCell<Fp24Field_2<11,5>> = OnceCell::new();
static BLS24_315_G1_CONSTS          : OnceCell<G1Consts<4,5,3>> = OnceCell::new();
static BLS24_315_G2_CONSTS          : OnceCell<G2Consts<11,4,5,7>> = OnceCell::new();
static FR_BLS24_315                 : OnceCell<PrimeField<4>> = OnceCell::new();
static G1_BLS24_315                 : OnceCell<G1Field<4,5,3>> = OnceCell::new();
static G2_BLS24_315                 : OnceCell<G2Field<11,4,5,7>> = OnceCell::new();
static GT_BLS24_315                 : OnceCell<GTField<5,11>> = OnceCell::new();

pub fn fr_bls24_315()-> &'static PrimeField<4>
{   if FR_BLS24_315.get().is_none() 
            { if BLS24_315_R_FIELD_PARAMS.get().is_none() 
                    {BLS24_315_R_FIELD_PARAMS.set(build_field_params(BLS24_315_PARAMS.r_as_strhex)).unwrap()}; 
              FR_BLS24_315.set(PrimeField::new(&BLS24_315_R_FIELD_PARAMS.get().unwrap())).unwrap();                                                                                                
            }
    FR_BLS24_315.get().unwrap()
}
pub fn fp4_bls24_315()-> &'static Fp4Field_2<11, 5>
{   if FP4_BLS24_315.get().is_none() {  if BLS24_315_FIELD_PARAMS.get().is_none() 
                                            {   BLS24_315_FIELD_PARAMS.set(build_field_params(BLS24_315_PARAMS.modulo_as_strhex)).unwrap()};
                                         if FP_BLS24_315.get().is_none() 
                                            {   FP_BLS24_315.set(PrimeField::new(&BLS24_315_FIELD_PARAMS.get().unwrap())).unwrap()};
                                         if BLS24_315_FROB_CONSTS.get().is_none() 
                                            {   BLS24_315_FROB_CONSTS.set(build_frob_params(BLS24_315_PARAMS, &FP_BLS24_315.get().unwrap())).unwrap()};               
                                         FP4_BLS24_315.set(Fp4Field_2::new( &FP_BLS24_315.get().unwrap(), Some(&BLS24_315_FROB_CONSTS.get().unwrap()))).unwrap();                                                                                                 
                                      }
    FP4_BLS24_315.get().unwrap()
}

pub fn g1_bls24_315() -> &'static G1Field<4,5,3>
{   if G1_BLS24_315.get().is_none() 
            { 
              if BLS24_315_FIELD_PARAMS.get().is_none() 
                    {BLS24_315_FIELD_PARAMS.set(build_field_params(BLS24_315_PARAMS.modulo_as_strhex)).unwrap()}; 
              if BLS24_315_R_FIELD_PARAMS.get().is_none()
                    {BLS24_315_R_FIELD_PARAMS.set(build_field_params(BLS24_315_PARAMS.r_as_strhex)).unwrap()};                     
              if FP_BLS24_315.get().is_none() 
                    {FP_BLS24_315.set(PrimeField::new(&BLS24_315_FIELD_PARAMS.get().unwrap())).unwrap();}                                                                                               
              if FR_BLS24_315.get().is_none() 
                    {FR_BLS24_315.set(PrimeField::new(&BLS24_315_R_FIELD_PARAMS.get().unwrap())).unwrap();}                                                                                                                           
              if BLS24_315_G1_CONSTS.get().is_none() 
                    {BLS24_315_G1_CONSTS.set(build_g1_params(BLS24_315_PARAMS,
                                                        &FP_BLS24_315.get().unwrap(),
                                                           &FR_BLS24_315.get().unwrap())).unwrap();}                                                                                                                                   
              G1_BLS24_315.set(G1Field { consts : BLS24_315_G1_CONSTS.get().unwrap(),
                                                  base_field: FP_BLS24_315.get().unwrap(),
                                                  fr_field : FR_BLS24_315.get().unwrap()}).unwrap() ;                                                  
            }
    G1_BLS24_315.get().unwrap()
}

pub fn g2_bls24_315() -> &'static G2Field<11,4,5,7>
{   if G2_BLS24_315.get().is_none() 
                {   if G1_BLS24_315.get().is_none() 
                            { let _ = g1_bls24_315();}
                    if BLS24_315_FROB_CONSTS.get().is_none() 
                            { BLS24_315_FROB_CONSTS.set(build_frob_params(BLS24_315_PARAMS, &FP_BLS24_315.get().unwrap())).unwrap()}; 
                    if FP4_BLS24_315.get().is_none() 
                            { let _ = fp4_bls24_315();}                                                                                           
                    if BLS24_315_G2_CONSTS.get().is_none() 
                            {BLS24_315_G2_CONSTS.set(build_g2_params(BLS24_315_PARAMS,
                                                                     &ExtG2Field::Fp4_2(FP4_BLS24_315.get().unwrap()),
                                                                     &PrimeField::new(&BLS24_315_R_FIELD_PARAMS.get().unwrap()),
                                                                     BLS24_315_FROB_CONSTS.get().unwrap())).unwrap();}                                                                                                                                                                                                                           
                    if EXFP4_BLS24_315.get().is_none()  {EXFP4_BLS24_315.set(ExtG2Field::Fp4_2(&FP4_BLS24_315.get().unwrap())).unwrap();}                                     
                    G2_BLS24_315.set(G2Field {  consts : BLS24_315_G2_CONSTS.get().unwrap(),
                                                base_field : &EXFP4_BLS24_315.get().unwrap() ,
                                                fr_field : FR_BLS24_315.get().unwrap()}).unwrap() ;
                }            
    G2_BLS24_315.get().unwrap()
}

pub fn fp24_bls24_315()-> &'static Fp24Field_2<11, 5>
{  if FP24_BLS24_315.get().is_none() {  if BLS24_315_FIELD_PARAMS.get().is_none() 
                                            {   BLS24_315_FIELD_PARAMS.set(build_field_params(BLS24_315_PARAMS.modulo_as_strhex)).unwrap()};
                                         if FP_BLS24_315.get().is_none() 
                                            {   FP_BLS24_315.set(PrimeField::new(&BLS24_315_FIELD_PARAMS.get().unwrap())).unwrap()};
                                         if BLS24_315_FROB_CONSTS.get().is_none() 
                                            {   BLS24_315_FROB_CONSTS.set(build_frob_params(BLS24_315_PARAMS, &FP_BLS24_315.get().unwrap())).unwrap()};               
                                         FP24_BLS24_315.set(Fp24Field_2::new(&FP_BLS24_315.get().unwrap(), Some(&BLS24_315_FROB_CONSTS.get().unwrap()))).unwrap();                                                                                                 
                                      }
    FP24_BLS24_315.get().unwrap()
}

pub fn gt_bls24_315() -> &'static GTField<5,11>
{  
        if FP24_BLS24_315.get().is_none() { let _ = fp24_bls24_315(); };
        GT_BLS24_315.set(GTField::Fp24_2(FP24_BLS24_315.get().unwrap())).unwrap();
        GT_BLS24_315.get().unwrap()       
}

static BLS24_509_SNARK_FIELD_PARAMS       : OnceCell<FieldParams<8>> = OnceCell::new();    
static BLS24_509_SNARK_FROB_CONSTS        : OnceCell<ExFieldConsts<11,8>> = OnceCell::new();
static BLS24_509_SNARK_R_FIELD_PARAMS     : OnceCell<FieldParams<7>> = OnceCell::new();   
static FP_BLS24_509_SNARK                 : OnceCell<PrimeField<8>> = OnceCell::new(); 
static FP4_BLS24_509_SNARK                : OnceCell<Fp4Field_2<11,8>> = OnceCell::new();
static EXFP4_BLS24_509_SNARK              : OnceCell<ExtG2Field<8,11>> = OnceCell::new();
static FP24_BLS24_509_SNARK               : OnceCell<Fp24Field_2<11,8>> = OnceCell::new();
static BLS24_509_SNARK_G1_CONSTS          : OnceCell<G1Consts<7,8,3>> = OnceCell::new();
static BLS24_509_SNARK_G2_CONSTS          : OnceCell<G2Consts<11,7,8,7>> = OnceCell::new();
static FR_BLS24_509_SNARK                 : OnceCell<PrimeField<7>> = OnceCell::new();
static G1_BLS24_509_SNARK                 : OnceCell<G1Field<7,8,3>> = OnceCell::new();
static G2_BLS24_509_SNARK                 : OnceCell<G2Field<11,7,8,7>> = OnceCell::new();
static GT_BLS24_509_SNARK                 : OnceCell<GTField<8,11>> = OnceCell::new();


pub fn fr_bls24_509_snark()-> &'static PrimeField<7>
{   if FR_BLS24_509_SNARK.get().is_none() 
            { if BLS24_509_SNARK_R_FIELD_PARAMS.get().is_none() 
                    {BLS24_509_SNARK_R_FIELD_PARAMS.set(build_field_params(BLS24_509_SNARK_PARAMS.r_as_strhex)).unwrap()}; 
              FR_BLS24_509_SNARK.set(PrimeField::new(&BLS24_509_SNARK_R_FIELD_PARAMS.get().unwrap())).unwrap();                                                                                                
            }
    FR_BLS24_509_SNARK.get().unwrap()
}

pub fn fp4_bls24_509_snark()-> &'static Fp4Field_2<11, 8>
{   if FP4_BLS24_509_SNARK.get().is_none() {  if BLS24_509_SNARK_FIELD_PARAMS.get().is_none() 
                                            {   BLS24_509_SNARK_FIELD_PARAMS.set(build_field_params(BLS24_509_SNARK_PARAMS.modulo_as_strhex)).unwrap()};
                                         if FP_BLS24_509_SNARK.get().is_none() 
                                            {   FP_BLS24_509_SNARK.set(PrimeField::new(&BLS24_509_SNARK_FIELD_PARAMS.get().unwrap())).unwrap()};
                                         if BLS24_509_SNARK_FROB_CONSTS.get().is_none() 
                                            {   BLS24_509_SNARK_FROB_CONSTS.set(build_frob_params(BLS24_509_SNARK_PARAMS, &FP_BLS24_509_SNARK.get().unwrap())).unwrap()};               
                                         FP4_BLS24_509_SNARK.set(Fp4Field_2::new(&FP_BLS24_509_SNARK.get().unwrap(), Some(&BLS24_509_SNARK_FROB_CONSTS.get().unwrap()))).unwrap();                                                                                                 
                                      }
    FP4_BLS24_509_SNARK.get().unwrap()
}

pub fn g1_bls24_509_snark() -> &'static G1Field<7,8,3>
{   if G1_BLS24_509_SNARK.get().is_none() 
            { 
              if BLS24_509_SNARK_FIELD_PARAMS.get().is_none() 
                    {BLS24_509_SNARK_FIELD_PARAMS.set(build_field_params(BLS24_509_SNARK_PARAMS.modulo_as_strhex)).unwrap()}; 
              if BLS24_509_SNARK_R_FIELD_PARAMS.get().is_none()
                    {BLS24_509_SNARK_R_FIELD_PARAMS.set(build_field_params(BLS24_509_SNARK_PARAMS.r_as_strhex)).unwrap()};                     
              if FP_BLS24_509_SNARK.get().is_none() 
                    {FP_BLS24_509_SNARK.set(PrimeField::new(&BLS24_509_SNARK_FIELD_PARAMS.get().unwrap())).unwrap();}                                                                                               
              if FR_BLS24_509_SNARK.get().is_none() 
                    {FR_BLS24_509_SNARK.set(PrimeField::new(&BLS24_509_SNARK_R_FIELD_PARAMS.get().unwrap())).unwrap();}                                                                                                                           
              if BLS24_509_SNARK_G1_CONSTS.get().is_none() 
                    {BLS24_509_SNARK_G1_CONSTS.set(build_g1_params(BLS24_509_SNARK_PARAMS,
                                                        &FP_BLS24_509_SNARK.get().unwrap(),
                                                           &FR_BLS24_509_SNARK.get().unwrap())).unwrap();}                                                                                                                                   
              G1_BLS24_509_SNARK.set(G1Field { consts : BLS24_509_SNARK_G1_CONSTS.get().unwrap(),
                                                  base_field: FP_BLS24_509_SNARK.get().unwrap(),
                                                  fr_field : FR_BLS24_509_SNARK.get().unwrap()}).unwrap() ;                                                  
            }
    G1_BLS24_509_SNARK.get().unwrap()
}

pub fn g2_bls24_509_snark() -> &'static G2Field<11,7,8,7>
{   if G2_BLS24_509_SNARK.get().is_none() 
                {   if G1_BLS24_509_SNARK.get().is_none() 
                            { let _ = g1_bls24_509_snark();}
                    if BLS24_509_SNARK_FROB_CONSTS.get().is_none() 
                            { BLS24_509_SNARK_FROB_CONSTS.set(build_frob_params(BLS24_509_SNARK_PARAMS, &FP_BLS24_509_SNARK.get().unwrap())).unwrap()}; 
                    if FP4_BLS24_509_SNARK.get().is_none() 
                            { let _ = fp4_bls24_509_snark();}                                                                                           
                    if BLS24_509_SNARK_G2_CONSTS.get().is_none() 
                            {BLS24_509_SNARK_G2_CONSTS.set(build_g2_params(BLS24_509_SNARK_PARAMS,
                                                                     &ExtG2Field::Fp4_2(FP4_BLS24_509_SNARK.get().unwrap()),
                                                                     &PrimeField::new(&BLS24_509_SNARK_R_FIELD_PARAMS.get().unwrap()),
                                                                     BLS24_509_SNARK_FROB_CONSTS.get().unwrap())).unwrap();}                                                                                                                                                                                                                           
                    if EXFP4_BLS24_509_SNARK.get().is_none()  {EXFP4_BLS24_509_SNARK.set(ExtG2Field::Fp4_2(&FP4_BLS24_509_SNARK.get().unwrap())).unwrap();}                                     
                    G2_BLS24_509_SNARK.set(G2Field {  consts : BLS24_509_SNARK_G2_CONSTS.get().unwrap(),
                                                base_field : &EXFP4_BLS24_509_SNARK.get().unwrap() ,
                                                fr_field : FR_BLS24_509_SNARK.get().unwrap()}).unwrap() ;
                }            
    G2_BLS24_509_SNARK.get().unwrap()
}

pub fn fp24_bls24_509_snark()-> &'static Fp24Field_2<11, 8>
{  if FP24_BLS24_509_SNARK.get().is_none() {  if BLS24_509_SNARK_FIELD_PARAMS.get().is_none() 
                                            {   BLS24_509_SNARK_FIELD_PARAMS.set(build_field_params(BLS24_509_SNARK_PARAMS.modulo_as_strhex)).unwrap()};
                                         if FP_BLS24_509_SNARK.get().is_none() 
                                            {   FP_BLS24_509_SNARK.set(PrimeField::new(&BLS24_509_SNARK_FIELD_PARAMS.get().unwrap())).unwrap()};
                                         if BLS24_509_SNARK_FROB_CONSTS.get().is_none() 
                                            {   BLS24_509_SNARK_FROB_CONSTS.set(build_frob_params(BLS24_509_SNARK_PARAMS, &FP_BLS24_509_SNARK.get().unwrap())).unwrap()};               
                                         FP24_BLS24_509_SNARK.set(Fp24Field_2::new(&FP_BLS24_509_SNARK.get().unwrap(), Some(&BLS24_509_SNARK_FROB_CONSTS.get().unwrap()))).unwrap();                                                                                                 
                                      }
    FP24_BLS24_509_SNARK.get().unwrap()
}

pub fn gt_bls24_509_snark() -> &'static GTField<8,11>
{  
        if FP24_BLS24_509_SNARK.get().is_none() { let _ = fp24_bls24_509_snark(); };
        GT_BLS24_509_SNARK.set(GTField::Fp24_2(FP24_BLS24_509_SNARK.get().unwrap())).unwrap();
        GT_BLS24_509_SNARK.get().unwrap()       
}

static BLS24_509_FIELD_PARAMS       : OnceCell<FieldParams<8>> = OnceCell::new();    
static BLS24_509_FROB_CONSTS        : OnceCell<ExFieldConsts<10,8>> = OnceCell::new();
static BLS24_509_R_FIELD_PARAMS     : OnceCell<FieldParams<7>> = OnceCell::new();   
static FP_BLS24_509                 : OnceCell<PrimeField<8>> = OnceCell::new(); 
static FP4_BLS24_509                : OnceCell<Fp4Field_1<10,8>> = OnceCell::new();
static EXFP4_BLS24_509              : OnceCell<ExtG2Field<8,10>> = OnceCell::new();
static FP24_BLS24_509               : OnceCell<Fp24Field_1<10,8>> = OnceCell::new();
static BLS24_509_G1_CONSTS          : OnceCell<G1Consts<7,8,3>> = OnceCell::new();
static BLS24_509_G2_CONSTS          : OnceCell<G2Consts<10,7,8,10>> = OnceCell::new();
static FR_BLS24_509                 : OnceCell<PrimeField<7>> = OnceCell::new();
static G1_BLS24_509                 : OnceCell<G1Field<7,8,3>> = OnceCell::new();
static G2_BLS24_509                 : OnceCell<G2Field<10,7,8,10>> = OnceCell::new();
static GT_BLS24_509                 : OnceCell<GTField<8,10>> = OnceCell::new();

pub fn fr_bls24_509()-> &'static PrimeField<7>
{   if FR_BLS24_509.get().is_none() 
            { if BLS24_509_R_FIELD_PARAMS.get().is_none() 
                    {BLS24_509_R_FIELD_PARAMS.set(build_field_params(BLS24_509_PARAMS.r_as_strhex)).unwrap()}; 
              FR_BLS24_509.set(PrimeField::new(&BLS24_509_R_FIELD_PARAMS.get().unwrap())).unwrap();                                                                                                
            }
    FR_BLS24_509.get().unwrap()
}

pub fn fp4_bls24_509()-> &'static Fp4Field_1<10, 8>
{   if FP4_BLS24_509.get().is_none() {  if BLS24_509_FIELD_PARAMS.get().is_none() 
                                            {   BLS24_509_FIELD_PARAMS.set(build_field_params(BLS24_509_PARAMS.modulo_as_strhex)).unwrap()};
                                         if FP_BLS24_509.get().is_none() 
                                            {   FP_BLS24_509.set(PrimeField::new(&BLS24_509_FIELD_PARAMS.get().unwrap())).unwrap()};
                                         if BLS24_509_FROB_CONSTS.get().is_none() 
                                            {   BLS24_509_FROB_CONSTS.set(build_frob_params(BLS24_509_PARAMS, &FP_BLS24_509.get().unwrap())).unwrap()};               
                                         FP4_BLS24_509.set(Fp4Field_1::new(&FP_BLS24_509.get().unwrap(), Some(&BLS24_509_FROB_CONSTS.get().unwrap()))).unwrap();                                                                                                 
                                      }
    FP4_BLS24_509.get().unwrap()
}

pub fn g1_bls24_509() -> &'static G1Field<7,8,3>
{   if G1_BLS24_509.get().is_none() 
            { 
              if BLS24_509_FIELD_PARAMS.get().is_none() 
                    {BLS24_509_FIELD_PARAMS.set(build_field_params(BLS24_509_PARAMS.modulo_as_strhex)).unwrap()}; 
              if BLS24_509_R_FIELD_PARAMS.get().is_none()
                    {BLS24_509_R_FIELD_PARAMS.set(build_field_params(BLS24_509_PARAMS.r_as_strhex)).unwrap()};                     
              if FP_BLS24_509.get().is_none() 
                    {FP_BLS24_509.set(PrimeField::new(&BLS24_509_FIELD_PARAMS.get().unwrap())).unwrap();}                                                                                               
              if FR_BLS24_509.get().is_none() 
                    {FR_BLS24_509.set(PrimeField::new(&BLS24_509_R_FIELD_PARAMS.get().unwrap())).unwrap();}                                                                                                                           
              if BLS24_509_G1_CONSTS.get().is_none() 
                    {BLS24_509_G1_CONSTS.set(build_g1_params(BLS24_509_PARAMS,
                                                        &FP_BLS24_509.get().unwrap(),
                                                           &FR_BLS24_509.get().unwrap())).unwrap();}                                                                                                                                   
              G1_BLS24_509.set(G1Field { consts : BLS24_509_G1_CONSTS.get().unwrap(),
                                                  base_field: FP_BLS24_509.get().unwrap(),
                                                  fr_field : FR_BLS24_509.get().unwrap()}).unwrap() ;                                                  
            }
    G1_BLS24_509.get().unwrap()
}

pub fn g2_bls24_509() -> &'static G2Field<10,7,8,10>
{   if G2_BLS24_509.get().is_none() 
                {   if G1_BLS24_509.get().is_none() 
                            { let _ = g1_bls24_509();}
                    if BLS24_509_FROB_CONSTS.get().is_none() 
                            { BLS24_509_FROB_CONSTS.set(build_frob_params(BLS24_509_PARAMS, &FP_BLS24_509.get().unwrap())).unwrap()}; 
                    if FP4_BLS24_509.get().is_none() 
                            { let _ = fp4_bls24_509();}                                                                                           
                    if BLS24_509_G2_CONSTS.get().is_none() 
                            {BLS24_509_G2_CONSTS.set(build_g2_params(BLS24_509_PARAMS,
                                                                     &ExtG2Field::Fp4_1(FP4_BLS24_509.get().unwrap()),
                                                                     &PrimeField::new(&BLS24_509_R_FIELD_PARAMS.get().unwrap()),
                                                                     BLS24_509_FROB_CONSTS.get().unwrap())).unwrap();}                                                                                                                                                                                                                           
                    if EXFP4_BLS24_509.get().is_none()  {EXFP4_BLS24_509.set(ExtG2Field::Fp4_1(&FP4_BLS24_509.get().unwrap())).unwrap();}                                     
                    G2_BLS24_509.set(G2Field {  consts : BLS24_509_G2_CONSTS.get().unwrap(),
                                                base_field : &EXFP4_BLS24_509.get().unwrap() ,
                                                fr_field : FR_BLS24_509.get().unwrap()}).unwrap() ;
                }            
    G2_BLS24_509.get().unwrap()
}

pub fn fp24_bls24_509()-> &'static Fp24Field_1<10, 8>
{  if FP24_BLS24_509.get().is_none() {  if BLS24_509_FIELD_PARAMS.get().is_none() 
                                            {   BLS24_509_FIELD_PARAMS.set(build_field_params(BLS24_509_PARAMS.modulo_as_strhex)).unwrap()};
                                         if FP_BLS24_509.get().is_none() 
                                            {   FP_BLS24_509.set(PrimeField::new(&BLS24_509_FIELD_PARAMS.get().unwrap())).unwrap()};
                                         if BLS24_509_FROB_CONSTS.get().is_none() 
                                            {   BLS24_509_FROB_CONSTS.set(build_frob_params(BLS24_509_PARAMS, &FP_BLS24_509.get().unwrap())).unwrap()};               
                                         FP24_BLS24_509.set(Fp24Field_1::new(&FP_BLS24_509.get().unwrap(), Some(&BLS24_509_FROB_CONSTS.get().unwrap()))).unwrap();                                                                                                 
                                      }
    FP24_BLS24_509.get().unwrap()
}

pub fn gt_bls24_509() -> &'static GTField<8,10>
{  
        if FP24_BLS24_509.get().is_none() { let _ = fp24_bls24_509(); };
        GT_BLS24_509.set(GTField::Fp24_1(FP24_BLS24_509.get().unwrap())).unwrap();
        GT_BLS24_509.get().unwrap()       
}