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

// Parameters, fields and constants for BLS48-575

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
{
    FR_BLS48_575.get_or_init(|| {
        let r_params = BLS48_575_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS48_575_PARAMS.r_as_strhex));
        PrimeField::new(r_params)
    })
}

pub fn fp8_bls48_575()-> &'static Fp8Field_1<24,9>
{
    FP8_BLS48_575.get_or_init(|| {
        let field_params = BLS48_575_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS48_575_PARAMS.modulo_as_strhex));
        let fp = FP_BLS48_575.get_or_init(|| PrimeField::new(field_params));
        let frob = BLS48_575_FROB_CONSTS
            .get_or_init(|| build_frob_params(BLS48_575_PARAMS, fp));
        Fp8Field_1::new(fp, Some(frob))
    })
}

pub fn g1_bls48_575() -> &'static G1Field<9,9,4>
{
    G1_BLS48_575.get_or_init(|| {
        let r_params = BLS48_575_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS48_575_PARAMS.r_as_strhex));
        let field_params = BLS48_575_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS48_575_PARAMS.modulo_as_strhex));
        let fp = FP_BLS48_575.get_or_init(|| PrimeField::new(field_params));
        let fr = FR_BLS48_575.get_or_init(|| PrimeField::new(r_params));
        let g1_consts = BLS48_575_G1_CONSTS
            .get_or_init(|| build_g1_params(BLS48_575_PARAMS, fp, fr));
        G1Field { consts: g1_consts, base_field: fp, fr_field: fr }
    })
}

pub fn g2_bls48_575() -> &'static G2Field<24,9,9,7>
{
    G2_BLS48_575.get_or_init(|| {
        let g1 = g1_bls48_575();
        let fp = g1.base_field;
        let fr = g1.fr_field;
        let frob = BLS48_575_FROB_CONSTS
            .get_or_init(|| build_frob_params(BLS48_575_PARAMS, fp));
        let fp8 = fp8_bls48_575();
        let r_params = BLS48_575_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS48_575_PARAMS.r_as_strhex));
        let g2_consts = BLS48_575_G2_CONSTS.get_or_init(|| build_g2_params(
            BLS48_575_PARAMS,
            &ExtG2Field::Fp8_1(fp8),
            &PrimeField::new(r_params),
            frob,
        ));
        let exfp8 = EXFP8_BLS48_575.get_or_init(|| ExtG2Field::Fp8_1(fp8));
        G2Field { consts: g2_consts, base_field: exfp8, fr_field: fr }
    })
}

pub fn fp48_bls48_575()-> &'static Fp48Field_1<24, 9>
{
    FP48_BLS48_575.get_or_init(|| {
        let field_params = BLS48_575_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS48_575_PARAMS.modulo_as_strhex));
        let fp = FP_BLS48_575.get_or_init(|| PrimeField::new(field_params));
        let frob = BLS48_575_FROB_CONSTS
            .get_or_init(|| build_frob_params(BLS48_575_PARAMS, fp));
        Fp48Field_1::new(fp, Some(frob))
    })
}

pub fn gt_bls48_575() -> &'static GTField<9,24>
{
    GT_BLS48_575.get_or_init(|| GTField::Fp48_1(fp48_bls48_575()))
}

// Parameters, fields and constants for BLS48-573

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
{
    FR_BLS48_573.get_or_init(|| {
        let r_params = BLS48_573_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS48_573_PARAMS.r_as_strhex));
        PrimeField::new(r_params)
    })
}

pub fn fp8_bls48_573()-> &'static Fp8Field_1<24,9>
{
    FP8_BLS48_573.get_or_init(|| {
        let field_params = BLS48_573_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS48_573_PARAMS.modulo_as_strhex));
        let fp = FP_BLS48_573.get_or_init(|| PrimeField::new(field_params));
        let frob = BLS48_573_FROB_CONSTS
            .get_or_init(|| build_frob_params(BLS48_573_PARAMS, fp));
        Fp8Field_1::new(fp, Some(frob))
    })
}

pub fn g1_bls48_573() -> &'static G1Field<8,9,7>
{
    G1_BLS48_573.get_or_init(|| {
        let r_params = BLS48_573_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS48_573_PARAMS.r_as_strhex));
        let field_params = BLS48_573_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS48_573_PARAMS.modulo_as_strhex));
        let fp = FP_BLS48_573.get_or_init(|| PrimeField::new(field_params));
        let fr = FR_BLS48_573.get_or_init(|| PrimeField::new(r_params));
        let g1_consts = BLS48_573_G1_CONSTS
            .get_or_init(|| build_g1_params(BLS48_573_PARAMS, fp, fr));
        G1Field { consts: g1_consts, base_field: fp, fr_field: fr }
    })
}

pub fn g2_bls48_573() -> &'static G2Field<24,8,9,4>
{
    G2_BLS48_573.get_or_init(|| {
        let g1 = g1_bls48_573();
        let fp = g1.base_field;
        let fr = g1.fr_field;
        let frob = BLS48_573_FROB_CONSTS
            .get_or_init(|| build_frob_params(BLS48_573_PARAMS, fp));
        let fp8 = fp8_bls48_573();
        let r_params = BLS48_573_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS48_573_PARAMS.r_as_strhex));
        let g2_consts = BLS48_573_G2_CONSTS.get_or_init(|| build_g2_params(
            BLS48_573_PARAMS,
            &ExtG2Field::Fp8_1(fp8),
            &PrimeField::new(r_params),
            frob,
        ));
        let exfp8 = EXFP8_BLS48_573.get_or_init(|| ExtG2Field::Fp8_1(fp8));
        G2Field { consts: g2_consts, base_field: exfp8, fr_field: fr }
    })
}

pub fn fp48_bls48_573()-> &'static Fp48Field_1<24, 9>
{
    FP48_BLS48_573.get_or_init(|| {
        let field_params = BLS48_573_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS48_573_PARAMS.modulo_as_strhex));
        let fp = FP_BLS48_573.get_or_init(|| PrimeField::new(field_params));
        let frob = BLS48_573_FROB_CONSTS
            .get_or_init(|| build_frob_params(BLS48_573_PARAMS, fp));
        Fp48Field_1::new(fp, Some(frob))
    })
}

pub fn gt_bls48_573() -> &'static GTField<9,24>
{
    GT_BLS48_573.get_or_init(|| GTField::Fp48_1(fp48_bls48_573()))
}

// Parameters, fields and constants for BLS48-581

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
{
    FR_BLS48_581.get_or_init(|| {
        let r_params = BLS48_581_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS48_581_PARAMS.r_as_strhex));
        PrimeField::new(r_params)
    })
}

pub fn fp8_bls48_581()-> &'static Fp8Field_3<24,10>
{
    FP8_BLS48_581.get_or_init(|| {
        let field_params = BLS48_581_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS48_581_PARAMS.modulo_as_strhex));
        let fp = FP_BLS48_581.get_or_init(|| PrimeField::new(field_params));
        let frob = BLS48_581_FROB_CONSTS
            .get_or_init(|| build_frob_params(BLS48_581_PARAMS, fp));
        Fp8Field_3::new(fp, Some(frob))
    })
}

pub fn g1_bls48_581() -> &'static G1Field<9,10,3>
{
    G1_BLS48_581.get_or_init(|| {
        let r_params = BLS48_581_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS48_581_PARAMS.r_as_strhex));
        let field_params = BLS48_581_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS48_581_PARAMS.modulo_as_strhex));
        let fp = FP_BLS48_581.get_or_init(|| PrimeField::new(field_params));
        let fr = FR_BLS48_581.get_or_init(|| PrimeField::new(r_params));
        let g1_consts = BLS48_581_G1_CONSTS
            .get_or_init(|| build_g1_params(BLS48_581_PARAMS, fp, fr));
        G1Field { consts: g1_consts, base_field: fp, fr_field: fr }
    })
}

pub fn g2_bls48_581() -> &'static G2Field<24,9,10,4>
{
    G2_BLS48_581.get_or_init(|| {
        let g1 = g1_bls48_581();
        let fp = g1.base_field;
        let fr = g1.fr_field;
        let frob = BLS48_581_FROB_CONSTS
            .get_or_init(|| build_frob_params(BLS48_581_PARAMS, fp));
        let fp8 = fp8_bls48_581();
        let r_params = BLS48_581_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS48_581_PARAMS.r_as_strhex));
        let g2_consts = BLS48_581_G2_CONSTS.get_or_init(|| build_g2_params(
            BLS48_581_PARAMS,
            &ExtG2Field::Fp8_3(fp8),
            &PrimeField::new(r_params),
            frob,
        ));
        let exfp8 = EXFP8_BLS48_581.get_or_init(|| ExtG2Field::Fp8_3(fp8));
        G2Field { consts: g2_consts, base_field: exfp8, fr_field: fr }
    })
}

pub fn fp48_bls48_581()-> &'static Fp48Field_3<24, 10>
{
    FP48_BLS48_581.get_or_init(|| {
        let field_params = BLS48_581_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS48_581_PARAMS.modulo_as_strhex));
        let fp = FP_BLS48_581.get_or_init(|| PrimeField::new(field_params));
        let frob = BLS48_581_FROB_CONSTS
            .get_or_init(|| build_frob_params(BLS48_581_PARAMS, fp));
        Fp48Field_3::new(fp, Some(frob))
    })
}

pub fn gt_bls48_581() -> &'static GTField<10,24>
{
    GT_BLS48_581.get_or_init(|| GTField::Fp48_3(fp48_bls48_581()))
}

// Parameters, fields and constants for BLS48-571

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
{
    FR_BLS48_571.get_or_init(|| {
        let r_params = BLS48_571_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS48_571_PARAMS.r_as_strhex));
        PrimeField::new(r_params)
    })
}

pub fn fp8_bls48_571()-> &'static Fp8Field_2<36,9>
{
    FP8_BLS48_571.get_or_init(|| {
        let field_params = BLS48_571_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS48_571_PARAMS.modulo_as_strhex));
        let fp = FP_BLS48_571.get_or_init(|| PrimeField::new(field_params));
        let frob = BLS48_571_FROB_CONSTS
            .get_or_init(|| build_frob_params(BLS48_571_PARAMS, fp));
        Fp8Field_2::new(fp, Some(frob))
    })
}

pub fn g1_bls48_571() -> &'static G1Field<8,9,3>
{
    G1_BLS48_571.get_or_init(|| {
        let r_params = BLS48_571_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS48_571_PARAMS.r_as_strhex));
        let field_params = BLS48_571_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS48_571_PARAMS.modulo_as_strhex));
        let fp = FP_BLS48_571.get_or_init(|| PrimeField::new(field_params));
        let fr = FR_BLS48_571.get_or_init(|| PrimeField::new(r_params));
        let g1_consts = BLS48_571_G1_CONSTS
            .get_or_init(|| build_g1_params(BLS48_571_PARAMS, fp, fr));
        G1Field { consts: g1_consts, base_field: fp, fr_field: fr }
    })
}

pub fn g2_bls48_571() -> &'static G2Field<36,8,9,4>
{
    G2_BLS48_571.get_or_init(|| {
        let g1 = g1_bls48_571();
        let fp = g1.base_field;
        let fr = g1.fr_field;
        let frob = BLS48_571_FROB_CONSTS
            .get_or_init(|| build_frob_params(BLS48_571_PARAMS, fp));
        let fp8 = fp8_bls48_571();
        let r_params = BLS48_571_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS48_571_PARAMS.r_as_strhex));
        let g2_consts = BLS48_571_G2_CONSTS.get_or_init(|| build_g2_params(
            BLS48_571_PARAMS,
            &ExtG2Field::Fp8_2(fp8),
            &PrimeField::new(r_params),
            frob,
        ));
        let exfp8 = EXFP8_BLS48_571.get_or_init(|| ExtG2Field::Fp8_2(fp8));
        G2Field { consts: g2_consts, base_field: exfp8, fr_field: fr }
    })
}

pub fn fp48_bls48_571()-> &'static Fp48Field_2<36, 9>
{
    FP48_BLS48_571.get_or_init(|| {
        let field_params = BLS48_571_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS48_571_PARAMS.modulo_as_strhex));
        let fp = FP_BLS48_571.get_or_init(|| PrimeField::new(field_params));
        let frob = BLS48_571_FROB_CONSTS
            .get_or_init(|| build_frob_params(BLS48_571_PARAMS, fp));
        Fp48Field_2::new(fp, Some(frob))
    })
}

pub fn gt_bls48_571() -> &'static GTField<9,36>
{
    GT_BLS48_571.get_or_init(|| GTField::Fp48_2(fp48_bls48_571()))
}

// Parameters, fields and constants for BLS48-287

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
{
    FR_BLS48_287.get_or_init(|| {
        let r_params = BLS48_287_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS48_287_PARAMS.r_as_strhex));
        PrimeField::new(r_params)
    })
}

pub fn fp8_bls48_287()-> &'static Fp8Field_2<35,5>
{
    FP8_BLS48_287.get_or_init(|| {
        let field_params = BLS48_287_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS48_287_PARAMS.modulo_as_strhex));
        let fp = FP_BLS48_287.get_or_init(|| PrimeField::new(field_params));
        let frob = BLS48_287_FROB_CONSTS
            .get_or_init(|| build_frob_params(BLS48_287_PARAMS, fp));
        Fp8Field_2::new(fp, Some(frob))
    })
}

pub fn g1_bls48_287() -> &'static G1Field<4,5,4>
{
    G1_BLS48_287.get_or_init(|| {
        let r_params = BLS48_287_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS48_287_PARAMS.r_as_strhex));
        let field_params = BLS48_287_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS48_287_PARAMS.modulo_as_strhex));
        let fp = FP_BLS48_287.get_or_init(|| PrimeField::new(field_params));
        let fr = FR_BLS48_287.get_or_init(|| PrimeField::new(r_params));
        let g1_consts = BLS48_287_G1_CONSTS
            .get_or_init(|| build_g1_params(BLS48_287_PARAMS, fp, fr));
        G1Field { consts: g1_consts, base_field: fp, fr_field: fr }
    })
}

pub fn g2_bls48_287() -> &'static G2Field<35,4,5,19>
{
    G2_BLS48_287.get_or_init(|| {
        let g1 = g1_bls48_287();
        let fp = g1.base_field;
        let fr = g1.fr_field;
        let frob = BLS48_287_FROB_CONSTS
            .get_or_init(|| build_frob_params(BLS48_287_PARAMS, fp));
        let fp8 = fp8_bls48_287();
        let r_params = BLS48_287_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS48_287_PARAMS.r_as_strhex));
        let g2_consts = BLS48_287_G2_CONSTS.get_or_init(|| build_g2_params(
            BLS48_287_PARAMS,
            &ExtG2Field::Fp8_2(fp8),
            &PrimeField::new(r_params),
            frob,
        ));
        let exfp8 = EXFP8_BLS48_287.get_or_init(|| ExtG2Field::Fp8_2(fp8));
        G2Field { consts: g2_consts, base_field: exfp8, fr_field: fr }
    })
}

pub fn fp48_bls48_287()-> &'static Fp48Field_2<35, 5>
{
    FP48_BLS48_287.get_or_init(|| {
        let field_params = BLS48_287_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS48_287_PARAMS.modulo_as_strhex));
        let fp = FP_BLS48_287.get_or_init(|| PrimeField::new(field_params));
        let frob = BLS48_287_FROB_CONSTS
            .get_or_init(|| build_frob_params(BLS48_287_PARAMS, fp));
        Fp48Field_2::new(fp, Some(frob))
    })
}

pub fn gt_bls48_287() -> &'static GTField<5,35>
{
    GT_BLS48_287.get_or_init(|| GTField::Fp48_2(fp48_bls48_287()))
}

// Parameters, fields and constants for BLS48-277

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
{
    FR_BLS48_277.get_or_init(|| {
        let r_params = BLS48_277_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS48_277_PARAMS.r_as_strhex));
        PrimeField::new(r_params)
    })
}

pub fn fp8_bls48_277()-> &'static Fp8Field_3<24,5>
{
    FP8_BLS48_277.get_or_init(|| {
        let field_params = BLS48_277_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS48_277_PARAMS.modulo_as_strhex));
        let fp = FP_BLS48_277.get_or_init(|| PrimeField::new(field_params));
        let frob = BLS48_277_FROB_CONSTS
            .get_or_init(|| build_frob_params(BLS48_277_PARAMS, fp));
        Fp8Field_3::new(fp, Some(frob))
    })
}

pub fn g1_bls48_277() -> &'static G1Field<4,5,10>
{
    G1_BLS48_277.get_or_init(|| {
        let r_params = BLS48_277_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS48_277_PARAMS.r_as_strhex));
        let field_params = BLS48_277_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS48_277_PARAMS.modulo_as_strhex));
        let fp = FP_BLS48_277.get_or_init(|| PrimeField::new(field_params));
        let fr = FR_BLS48_277.get_or_init(|| PrimeField::new(r_params));
        let g1_consts = BLS48_277_G1_CONSTS
            .get_or_init(|| build_g1_params(BLS48_277_PARAMS, fp, fr));
        G1Field { consts: g1_consts, base_field: fp, fr_field: fr }
    })
}

pub fn g2_bls48_277() -> &'static G2Field<24,4,5,4>
{
    G2_BLS48_277.get_or_init(|| {
        let g1 = g1_bls48_277();
        let fp = g1.base_field;
        let fr = g1.fr_field;
        let frob = BLS48_277_FROB_CONSTS
            .get_or_init(|| build_frob_params(BLS48_277_PARAMS, fp));
        let fp8 = fp8_bls48_277();
        let r_params = BLS48_277_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS48_277_PARAMS.r_as_strhex));
        let g2_consts = BLS48_277_G2_CONSTS.get_or_init(|| build_g2_params(
            BLS48_277_PARAMS,
            &ExtG2Field::Fp8_3(fp8),
            &PrimeField::new(r_params),
            frob,
        ));
        let exfp8 = EXFP8_BLS48_277.get_or_init(|| ExtG2Field::Fp8_3(fp8));
        G2Field { consts: g2_consts, base_field: exfp8, fr_field: fr }
    })
}

pub fn fp48_bls48_277()-> &'static Fp48Field_3<24, 5>
{
    FP48_BLS48_277.get_or_init(|| {
        let field_params = BLS48_277_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS48_277_PARAMS.modulo_as_strhex));
        let fp = FP_BLS48_277.get_or_init(|| PrimeField::new(field_params));
        let frob = BLS48_277_FROB_CONSTS
            .get_or_init(|| build_frob_params(BLS48_277_PARAMS, fp));
        Fp48Field_3::new(fp, Some(frob))
    })
}

pub fn gt_bls48_277() -> &'static GTField<5,24>
{
    GT_BLS48_277.get_or_init(|| GTField::Fp48_3(fp48_bls48_277()))
}
