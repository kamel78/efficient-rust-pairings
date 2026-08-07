// Code developed by FARAOUN Kamel Mohamed.
// EEDIS-Laboratory. UDL-University. Algeria
// During May 2024.

use crate::curves::g2::*;
use crate::curves::g1::{G1Consts, G1Field};
use crate::curves::gt::GTField;
use crate::extensions::ext_fields::{ExFieldConsts, ExtField};
use crate::extensions::towering1::fp2::Fp2Field;
use crate::extensions::towering1::fp12::Fp12Field;
use crate::extensions::g2_extfields::ExtG2Field;

use crate::fields::prime_fields::{FieldParams, PrimeField};
use crate::parameters::builders::*;
use crate::parameters::paramlist::{BLS12_381_PARAMS, BLS12_446_PARAMS, BLS12_461_PARAMS};
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
{
    G1_BLS12_381.get_or_init(|| {
        let r_params = BLS12_381_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS12_381_PARAMS.r_as_strhex));
        let field_params = BLS12_381_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS12_381_PARAMS.modulo_as_strhex));
        let fp = FP_BLS12_381.get_or_init(|| PrimeField::new(field_params));
        let fr = FR_BLS12_381.get_or_init(|| PrimeField::new(r_params));
        let g1_consts = BLS12_381_G1_CONSTS
            .get_or_init(|| build_g1_params(BLS12_381_PARAMS, fp, fr));
        G1Field { consts: g1_consts, base_field: fp, fr_field: fr }
    })
}

pub fn g2_bls12_381() -> &'static G2Field<4,4,6,4>
{
    G2_BLS12_381.get_or_init(|| {
        let g1 = g1_bls12_381();
        let fp = g1.base_field;
        let fr = g1.fr_field;
        let frob = BLS12_381_FROB_CONSTS
            .get_or_init(|| build_frob_params(BLS12_381_PARAMS, fp));
        let fp2 = FP2_BLS12_381.get_or_init(|| Fp2Field::new(fp, None));
        let r_params = BLS12_381_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS12_381_PARAMS.r_as_strhex));
        let g2_consts = BLS12_381_G2_CONSTS.get_or_init(|| build_g2_params(
            BLS12_381_PARAMS,
            &ExtG2Field::Fp2_1(fp2),
            &PrimeField::new(r_params),
            frob,
        ));
        let exfp2 = EXFP2_BLS12_381.get_or_init(|| ExtG2Field::Fp2_1(fp2));
        G2Field { consts: g2_consts, base_field: exfp2, fr_field: fr }
    })
}

pub fn fp12_bls12_381()-> &'static Fp12Field<4, 6>
{
    FP12_BLS12_381.get_or_init(|| {
        let field_params = BLS12_381_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS12_381_PARAMS.modulo_as_strhex));
        let fp = FP_BLS12_381.get_or_init(|| PrimeField::new(field_params));
        let frob = BLS12_381_FROB_CONSTS
            .get_or_init(|| build_frob_params(BLS12_381_PARAMS, fp));
        Fp12Field::new(fp, Some(frob))
    })
}

pub fn gt_bls12_381() -> &'static GTField<6,4>
{
    GT_BLS12_381.get_or_init(|| GTField::Fp12_1(fp12_bls12_381()))
}

pub fn fr_bls12_381()-> &'static PrimeField<4>
{
    FR_BLS12_381.get_or_init(|| {
        let r_params = BLS12_381_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS12_381_PARAMS.r_as_strhex));
        PrimeField::new(r_params)
    })
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

pub fn fr_bls12_461()-> &'static PrimeField<5>
{
    FR_BLS12_461.get_or_init(|| {
        let r_params = BLS12_461_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS12_461_PARAMS.r_as_strhex));
        PrimeField::new(r_params)
    })
}

pub fn g1_bls12_461() -> &'static G1Field<5,8,10>
{
    G1_BLS12_461.get_or_init(|| {
        let r_params = BLS12_461_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS12_461_PARAMS.r_as_strhex));
        let field_params = BLS12_461_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS12_461_PARAMS.modulo_as_strhex));
        let fp = FP_BLS12_461.get_or_init(|| PrimeField::new(field_params));
        let fr = FR_BLS12_461.get_or_init(|| PrimeField::new(r_params));
        let g1_consts = BLS12_461_G1_CONSTS
            .get_or_init(|| build_g1_params(BLS12_461_PARAMS, fp, fr));
        G1Field { consts: g1_consts, base_field: fp, fr_field: fr }
    })
}

pub fn g2_bls12_461() -> &'static G2Field<4,5,8,4>
{
    G2_BLS12_461.get_or_init(|| {
        let g1 = g1_bls12_461();
        let fp = g1.base_field;
        let fr = g1.fr_field;
        let frob = BLS12_461_FROB_CONSTS
            .get_or_init(|| build_frob_params(BLS12_461_PARAMS, fp));
        let fp2 = FP2_BLS12_461.get_or_init(|| Fp2Field::new(fp, None));
        let r_params = BLS12_461_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS12_461_PARAMS.r_as_strhex));
        let g2_consts = BLS12_461_G2_CONSTS.get_or_init(|| build_g2_params(
            BLS12_461_PARAMS,
            &ExtG2Field::Fp2_1(fp2),
            &PrimeField::new(r_params),
            frob,
        ));
        let exfp2 = EXFP2_BLS12_461.get_or_init(|| ExtG2Field::Fp2_1(fp2));
        G2Field { consts: g2_consts, base_field: exfp2, fr_field: fr }
    })
}

pub fn fp12_bls12_461()-> &'static Fp12Field<4, 8>
{
    FP12_BLS12_461.get_or_init(|| {
        let field_params = BLS12_461_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS12_461_PARAMS.modulo_as_strhex));
        // FIXED: was erroneously populating FP_BLS12_381 here.
        let fp = FP_BLS12_461.get_or_init(|| PrimeField::new(field_params));
        let frob = BLS12_461_FROB_CONSTS
            .get_or_init(|| build_frob_params(BLS12_461_PARAMS, fp));
        Fp12Field::new(fp, Some(frob))
    })
}

pub fn gt_bls12_461() -> &'static GTField<8,4>
{
    GT_BLS12_461.get_or_init(|| GTField::Fp12_1(fp12_bls12_461()))
}

// Parameters, fields and constants for BLS12-446

static BLS12_446_FIELD_PARAMS   : OnceCell<FieldParams<7>>              = OnceCell::new();    
static BLS12_446_R_FIELD_PARAMS : OnceCell<FieldParams<5>>              = OnceCell::new();    
static BLS12_446_FROB_CONSTS    : OnceCell<ExFieldConsts<4,7>>          = OnceCell::new();
static BLS12_446_G1_CONSTS      : OnceCell<G1Consts<5,7,3>>            = OnceCell::new();
static BLS12_446_G2_CONSTS      : OnceCell<G2Consts<4,5,7,7>>           = OnceCell::new();
static FP_BLS12_446             : OnceCell<PrimeField<7>>               = OnceCell::new();
static FP2_BLS12_446            : OnceCell<Fp2Field<7>>                 = OnceCell::new();
static EXFP2_BLS12_446          : OnceCell<ExtG2Field<7,4>>             = OnceCell::new();
static FR_BLS12_446             : OnceCell<PrimeField<5>>               = OnceCell::new();
static G1_BLS12_446             : OnceCell<G1Field<5,7,3>>             = OnceCell::new();
static G2_BLS12_446             : OnceCell<G2Field<4,5,7,7>>            = OnceCell::new();
static FP12_BLS12_446           : OnceCell<Fp12Field<4,7>>              = OnceCell::new();
static GT_BLS12_446             : OnceCell<GTField<7,4>>                = OnceCell::new();

pub fn fr_bls12_446()-> &'static PrimeField<5>
{
    FR_BLS12_446.get_or_init(|| {
        let r_params = BLS12_446_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS12_446_PARAMS.r_as_strhex));
        PrimeField::new(r_params)
    })
}

pub fn g1_bls12_446() -> &'static G1Field<5,7,3>
{
    G1_BLS12_446.get_or_init(|| {
        let r_params = BLS12_446_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS12_446_PARAMS.r_as_strhex));
        let field_params = BLS12_446_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS12_446_PARAMS.modulo_as_strhex));
        let fp = FP_BLS12_446.get_or_init(|| PrimeField::new(field_params));
        let fr = FR_BLS12_446.get_or_init(|| PrimeField::new(r_params));
        let g1_consts = BLS12_446_G1_CONSTS
            .get_or_init(|| build_g1_params(BLS12_446_PARAMS, fp, fr));
        G1Field { consts: g1_consts, base_field: fp, fr_field: fr }
    })
}

pub fn g2_bls12_446() -> &'static G2Field<4,5,7,7>
{
    G2_BLS12_446.get_or_init(|| {
        let g1 = g1_bls12_446();
        let fp = g1.base_field;
        let fr = g1.fr_field;
        let frob = BLS12_446_FROB_CONSTS
            .get_or_init(|| build_frob_params(BLS12_446_PARAMS, fp));
        let fp2 = FP2_BLS12_446.get_or_init(|| Fp2Field::new(fp, None));
        let r_params = BLS12_446_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS12_446_PARAMS.r_as_strhex));
        let g2_consts = BLS12_446_G2_CONSTS.get_or_init(|| build_g2_params(
            BLS12_446_PARAMS,
            &ExtG2Field::Fp2_1(fp2),
            &PrimeField::new(r_params),
            frob,
        ));
        let exfp2 = EXFP2_BLS12_446.get_or_init(|| ExtG2Field::Fp2_1(fp2));
        G2Field { consts: g2_consts, base_field: exfp2, fr_field: fr }
    })
}

pub fn fp12_bls12_446()-> &'static Fp12Field<4, 7>
{
    FP12_BLS12_446.get_or_init(|| {
        let field_params = BLS12_446_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS12_446_PARAMS.modulo_as_strhex));
        // FIXED: was erroneously populating FP_BLS12_381 here.
        let fp = FP_BLS12_446.get_or_init(|| PrimeField::new(field_params));
        let frob = BLS12_446_FROB_CONSTS
            .get_or_init(|| build_frob_params(BLS12_446_PARAMS, fp));
        Fp12Field::new(fp, Some(frob))
    })
}

pub fn gt_bls12_446() -> &'static GTField<7,4>
{
    GT_BLS12_446.get_or_init(|| GTField::Fp12_1(fp12_bls12_446()))
}
