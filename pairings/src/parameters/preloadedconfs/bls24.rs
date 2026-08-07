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

// Parameters, fields and constants for BLS24-477

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
{
    FR_BLS24_477.get_or_init(|| {
        let r_params = BLS24_477_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS24_477_PARAMS.r_as_strhex));
        PrimeField::new(r_params)
    })
}

pub fn fp4_bls24_477()-> &'static Fp4Field_1<11, 8>
{
    FP4_BLS24_477.get_or_init(|| {
        let field_params = BLS24_477_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS24_477_PARAMS.modulo_as_strhex));
        let fp = FP_BLS24_477.get_or_init(|| PrimeField::new(field_params));
        let frob = BLS24_477_FROB_CONSTS
            .get_or_init(|| build_frob_params(BLS24_477_PARAMS, fp));
        Fp4Field_1::new(fp, Some(frob))
    })
}

pub fn g1_bls24_477() -> &'static G1Field<6,8,7>
{
    G1_BLS24_477.get_or_init(|| {
        let r_params = BLS24_477_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS24_477_PARAMS.r_as_strhex));
        let field_params = BLS24_477_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS24_477_PARAMS.modulo_as_strhex));
        let fp = FP_BLS24_477.get_or_init(|| PrimeField::new(field_params));
        let fr = FR_BLS24_477.get_or_init(|| PrimeField::new(r_params));
        let g1_consts = BLS24_477_G1_CONSTS
            .get_or_init(|| build_g1_params(BLS24_477_PARAMS, fp, fr));
        G1Field { consts: g1_consts, base_field: fp, fr_field: fr }
    })
}

pub fn g2_bls24_477() -> &'static G2Field<11,6,8,4>
{
    G2_BLS24_477.get_or_init(|| {
        let g1 = g1_bls24_477();
        let fp = g1.base_field;
        let fr = g1.fr_field;
        let frob = BLS24_477_FROB_CONSTS
            .get_or_init(|| build_frob_params(BLS24_477_PARAMS, fp));
        let fp4 = fp4_bls24_477();
        let r_params = BLS24_477_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS24_477_PARAMS.r_as_strhex));
        let g2_consts = BLS24_477_G2_CONSTS.get_or_init(|| build_g2_params(
            BLS24_477_PARAMS,
            &ExtG2Field::Fp4_1(fp4),
            &PrimeField::new(r_params),
            frob,
        ));
        let exfp4 = EXFP4_BLS24_477.get_or_init(|| ExtG2Field::Fp4_1(fp4));
        G2Field { consts: g2_consts, base_field: exfp4, fr_field: fr }
    })
}

pub fn fp24_bls24_477()-> &'static Fp24Field_1<11, 8>
{
    FP24_BLS24_477.get_or_init(|| {
        let field_params = BLS24_477_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS24_477_PARAMS.modulo_as_strhex));
        let fp = FP_BLS24_477.get_or_init(|| PrimeField::new(field_params));
        let frob = BLS24_477_FROB_CONSTS
            .get_or_init(|| build_frob_params(BLS24_477_PARAMS, fp));
        Fp24Field_1::new(fp, Some(frob))
    })
}

pub fn gt_bls24_477() -> &'static GTField<8,11>
{
    GT_BLS24_477.get_or_init(|| GTField::Fp24_1(fp24_bls24_477()))
}

// Parameters, fields and constants for BLS24-479

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
{
    FR_BLS24_479.get_or_init(|| {
        let r_params = BLS24_479_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS24_479_PARAMS.r_as_strhex));
        PrimeField::new(r_params)
    })
}

pub fn fp4_bls24_479()-> &'static Fp4Field_1<10, 8>
{
    FP4_BLS24_479.get_or_init(|| {
        let field_params = BLS24_479_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS24_479_PARAMS.modulo_as_strhex));
        let fp = FP_BLS24_479.get_or_init(|| PrimeField::new(field_params));
        let frob = BLS24_479_FROB_CONSTS
            .get_or_init(|| build_frob_params(BLS24_479_PARAMS, fp));
        Fp4Field_1::new(fp, Some(frob))
    })
}

pub fn g1_bls24_479() -> &'static G1Field<7,8,4>
{
    G1_BLS24_479.get_or_init(|| {
        let r_params = BLS24_479_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS24_479_PARAMS.r_as_strhex));
        let field_params = BLS24_479_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS24_479_PARAMS.modulo_as_strhex));
        let fp = FP_BLS24_479.get_or_init(|| PrimeField::new(field_params));
        let fr = FR_BLS24_479.get_or_init(|| PrimeField::new(r_params));
        let g1_consts = BLS24_479_G1_CONSTS
            .get_or_init(|| build_g1_params(BLS24_479_PARAMS, fp, fr));
        G1Field { consts: g1_consts, base_field: fp, fr_field: fr }
    })
}

pub fn g2_bls24_479() -> &'static G2Field<10,7,8,10>
{
    G2_BLS24_479.get_or_init(|| {
        let g1 = g1_bls24_479();
        let fp = g1.base_field;
        let fr = g1.fr_field;
        let frob = BLS24_479_FROB_CONSTS
            .get_or_init(|| build_frob_params(BLS24_479_PARAMS, fp));
        let fp4 = fp4_bls24_479();
        let r_params = BLS24_479_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS24_479_PARAMS.r_as_strhex));
        let g2_consts = BLS24_479_G2_CONSTS.get_or_init(|| build_g2_params(
            BLS24_479_PARAMS,
            &ExtG2Field::Fp4_1(fp4),
            &PrimeField::new(r_params),
            frob,
        ));
        let exfp4 = EXFP4_BLS24_479.get_or_init(|| ExtG2Field::Fp4_1(fp4));
        G2Field { consts: g2_consts, base_field: exfp4, fr_field: fr }
    })
}

pub fn fp24_bls24_479()-> &'static Fp24Field_1<10, 8>
{
    FP24_BLS24_479.get_or_init(|| {
        let field_params = BLS24_479_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS24_479_PARAMS.modulo_as_strhex));
        let fp = FP_BLS24_479.get_or_init(|| PrimeField::new(field_params));
        let frob = BLS24_479_FROB_CONSTS
            .get_or_init(|| build_frob_params(BLS24_479_PARAMS, fp));
        Fp24Field_1::new(fp, Some(frob))
    })
}

pub fn gt_bls24_479() -> &'static GTField<8,10>
{
    GT_BLS24_479.get_or_init(|| GTField::Fp24_1(fp24_bls24_479()))
}

// Parameters, fields and constants for BLS24-559

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
{
    FR_BLS24_559.get_or_init(|| {
        let r_params = BLS24_559_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS24_559_PARAMS.r_as_strhex));
        PrimeField::new(r_params)
    })
}

pub fn fp4_bls24_559()-> &'static Fp4Field_1<10, 9>
{
    FP4_BLS24_559.get_or_init(|| {
        let field_params = BLS24_559_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS24_559_PARAMS.modulo_as_strhex));
        let fp = FP_BLS24_559.get_or_init(|| PrimeField::new(field_params));
        let frob = BLS24_559_FROB_CONSTS
            .get_or_init(|| build_frob_params(BLS24_559_PARAMS, fp));
        Fp4Field_1::new(fp, Some(frob))
    })
}

pub fn g1_bls24_559() -> &'static G1Field<8,9,4>
{
    G1_BLS24_559.get_or_init(|| {
        let r_params = BLS24_559_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS24_559_PARAMS.r_as_strhex));
        let field_params = BLS24_559_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS24_559_PARAMS.modulo_as_strhex));
        let fp = FP_BLS24_559.get_or_init(|| PrimeField::new(field_params));
        let fr = FR_BLS24_559.get_or_init(|| PrimeField::new(r_params));
        let g1_consts = BLS24_559_G1_CONSTS
            .get_or_init(|| build_g1_params(BLS24_559_PARAMS, fp, fr));
        G1Field { consts: g1_consts, base_field: fp, fr_field: fr }
    })
}

pub fn g2_bls24_559() -> &'static G2Field<10,8,9,4>
{
    G2_BLS24_559.get_or_init(|| {
        let g1 = g1_bls24_559();
        let fp = g1.base_field;
        let fr = g1.fr_field;
        let frob = BLS24_559_FROB_CONSTS
            .get_or_init(|| build_frob_params(BLS24_559_PARAMS, fp));
        let fp4 = fp4_bls24_559();
        let r_params = BLS24_559_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS24_559_PARAMS.r_as_strhex));
        let g2_consts = BLS24_559_G2_CONSTS.get_or_init(|| build_g2_params(
            BLS24_559_PARAMS,
            &ExtG2Field::Fp4_1(fp4),
            &PrimeField::new(r_params),
            frob,
        ));
        let exfp4 = EXFP4_BLS24_559.get_or_init(|| ExtG2Field::Fp4_1(fp4));
        G2Field { consts: g2_consts, base_field: exfp4, fr_field: fr }
    })
}

pub fn fp24_bls24_559()-> &'static Fp24Field_1<10, 9>
{
    FP24_BLS24_559.get_or_init(|| {
        let field_params = BLS24_559_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS24_559_PARAMS.modulo_as_strhex));
        let fp = FP_BLS24_559.get_or_init(|| PrimeField::new(field_params));
        let frob = BLS24_559_FROB_CONSTS
            .get_or_init(|| build_frob_params(BLS24_559_PARAMS, fp));
        Fp24Field_1::new(fp, Some(frob))
    })
}

pub fn gt_bls24_559() -> &'static GTField<9,10>
{
    GT_BLS24_559.get_or_init(|| GTField::Fp24_1(fp24_bls24_559()))
}

// Parameters, fields and constants for BLS24-315

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
{
    FR_BLS24_315.get_or_init(|| {
        let r_params = BLS24_315_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS24_315_PARAMS.r_as_strhex));
        PrimeField::new(r_params)
    })
}

pub fn fp4_bls24_315()-> &'static Fp4Field_2<11, 5>
{
    FP4_BLS24_315.get_or_init(|| {
        let field_params = BLS24_315_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS24_315_PARAMS.modulo_as_strhex));
        let fp = FP_BLS24_315.get_or_init(|| PrimeField::new(field_params));
        let frob = BLS24_315_FROB_CONSTS
            .get_or_init(|| build_frob_params(BLS24_315_PARAMS, fp));
        Fp4Field_2::new(fp, Some(frob))
    })
}

pub fn g1_bls24_315() -> &'static G1Field<4,5,3>
{
    G1_BLS24_315.get_or_init(|| {
        let field_params = BLS24_315_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS24_315_PARAMS.modulo_as_strhex));
        let r_params = BLS24_315_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS24_315_PARAMS.r_as_strhex));
        let fp = FP_BLS24_315.get_or_init(|| PrimeField::new(field_params));
        let fr = FR_BLS24_315.get_or_init(|| PrimeField::new(r_params));
        let g1_consts = BLS24_315_G1_CONSTS
            .get_or_init(|| build_g1_params(BLS24_315_PARAMS, fp, fr));
        G1Field { consts: g1_consts, base_field: fp, fr_field: fr }
    })
}

pub fn g2_bls24_315() -> &'static G2Field<11,4,5,7>
{
    G2_BLS24_315.get_or_init(|| {
        let g1 = g1_bls24_315();
        let fp = g1.base_field;
        let fr = g1.fr_field;
        let frob = BLS24_315_FROB_CONSTS
            .get_or_init(|| build_frob_params(BLS24_315_PARAMS, fp));
        let fp4 = fp4_bls24_315();
        let r_params = BLS24_315_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS24_315_PARAMS.r_as_strhex));
        let g2_consts = BLS24_315_G2_CONSTS.get_or_init(|| build_g2_params(
            BLS24_315_PARAMS,
            &ExtG2Field::Fp4_2(fp4),
            &PrimeField::new(r_params),
            frob,
        ));
        let exfp4 = EXFP4_BLS24_315.get_or_init(|| ExtG2Field::Fp4_2(fp4));
        G2Field { consts: g2_consts, base_field: exfp4, fr_field: fr }
    })
}

pub fn fp24_bls24_315()-> &'static Fp24Field_2<11, 5>
{
    FP24_BLS24_315.get_or_init(|| {
        let field_params = BLS24_315_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS24_315_PARAMS.modulo_as_strhex));
        let fp = FP_BLS24_315.get_or_init(|| PrimeField::new(field_params));
        let frob = BLS24_315_FROB_CONSTS
            .get_or_init(|| build_frob_params(BLS24_315_PARAMS, fp));
        Fp24Field_2::new(fp, Some(frob))
    })
}

pub fn gt_bls24_315() -> &'static GTField<5,11>
{
    GT_BLS24_315.get_or_init(|| GTField::Fp24_2(fp24_bls24_315()))
}

// Parameters, fields and constants for BLS24-509-SNARK

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
{
    FR_BLS24_509_SNARK.get_or_init(|| {
        let r_params = BLS24_509_SNARK_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS24_509_SNARK_PARAMS.r_as_strhex));
        PrimeField::new(r_params)
    })
}

pub fn fp4_bls24_509_snark()-> &'static Fp4Field_2<11, 8>
{
    FP4_BLS24_509_SNARK.get_or_init(|| {
        let field_params = BLS24_509_SNARK_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS24_509_SNARK_PARAMS.modulo_as_strhex));
        let fp = FP_BLS24_509_SNARK.get_or_init(|| PrimeField::new(field_params));
        let frob = BLS24_509_SNARK_FROB_CONSTS
            .get_or_init(|| build_frob_params(BLS24_509_SNARK_PARAMS, fp));
        Fp4Field_2::new(fp, Some(frob))
    })
}

pub fn g1_bls24_509_snark() -> &'static G1Field<7,8,3>
{
    G1_BLS24_509_SNARK.get_or_init(|| {
        let field_params = BLS24_509_SNARK_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS24_509_SNARK_PARAMS.modulo_as_strhex));
        let r_params = BLS24_509_SNARK_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS24_509_SNARK_PARAMS.r_as_strhex));
        let fp = FP_BLS24_509_SNARK.get_or_init(|| PrimeField::new(field_params));
        let fr = FR_BLS24_509_SNARK.get_or_init(|| PrimeField::new(r_params));
        let g1_consts = BLS24_509_SNARK_G1_CONSTS
            .get_or_init(|| build_g1_params(BLS24_509_SNARK_PARAMS, fp, fr));
        G1Field { consts: g1_consts, base_field: fp, fr_field: fr }
    })
}

pub fn g2_bls24_509_snark() -> &'static G2Field<11,7,8,7>
{
    G2_BLS24_509_SNARK.get_or_init(|| {
        let g1 = g1_bls24_509_snark();
        let fp = g1.base_field;
        let fr = g1.fr_field;
        let frob = BLS24_509_SNARK_FROB_CONSTS
            .get_or_init(|| build_frob_params(BLS24_509_SNARK_PARAMS, fp));
        let fp4 = fp4_bls24_509_snark();
        let r_params = BLS24_509_SNARK_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS24_509_SNARK_PARAMS.r_as_strhex));
        let g2_consts = BLS24_509_SNARK_G2_CONSTS.get_or_init(|| build_g2_params(
            BLS24_509_SNARK_PARAMS,
            &ExtG2Field::Fp4_2(fp4),
            &PrimeField::new(r_params),
            frob,
        ));
        let exfp4 = EXFP4_BLS24_509_SNARK.get_or_init(|| ExtG2Field::Fp4_2(fp4));
        G2Field { consts: g2_consts, base_field: exfp4, fr_field: fr }
    })
}

pub fn fp24_bls24_509_snark()-> &'static Fp24Field_2<11, 8>
{
    FP24_BLS24_509_SNARK.get_or_init(|| {
        let field_params = BLS24_509_SNARK_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS24_509_SNARK_PARAMS.modulo_as_strhex));
        let fp = FP_BLS24_509_SNARK.get_or_init(|| PrimeField::new(field_params));
        let frob = BLS24_509_SNARK_FROB_CONSTS
            .get_or_init(|| build_frob_params(BLS24_509_SNARK_PARAMS, fp));
        Fp24Field_2::new(fp, Some(frob))
    })
}

pub fn gt_bls24_509_snark() -> &'static GTField<8,11>
{
    GT_BLS24_509_SNARK.get_or_init(|| GTField::Fp24_2(fp24_bls24_509_snark()))
}

// Parameters, fields and constants for BLS24-509

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
{
    FR_BLS24_509.get_or_init(|| {
        let r_params = BLS24_509_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS24_509_PARAMS.r_as_strhex));
        PrimeField::new(r_params)
    })
}

pub fn fp4_bls24_509()-> &'static Fp4Field_1<10, 8>
{
    FP4_BLS24_509.get_or_init(|| {
        let field_params = BLS24_509_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS24_509_PARAMS.modulo_as_strhex));
        let fp = FP_BLS24_509.get_or_init(|| PrimeField::new(field_params));
        let frob = BLS24_509_FROB_CONSTS
            .get_or_init(|| build_frob_params(BLS24_509_PARAMS, fp));
        Fp4Field_1::new(fp, Some(frob))
    })
}

pub fn g1_bls24_509() -> &'static G1Field<7,8,3>
{
    G1_BLS24_509.get_or_init(|| {
        let field_params = BLS24_509_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS24_509_PARAMS.modulo_as_strhex));
        let r_params = BLS24_509_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS24_509_PARAMS.r_as_strhex));
        let fp = FP_BLS24_509.get_or_init(|| PrimeField::new(field_params));
        let fr = FR_BLS24_509.get_or_init(|| PrimeField::new(r_params));
        let g1_consts = BLS24_509_G1_CONSTS
            .get_or_init(|| build_g1_params(BLS24_509_PARAMS, fp, fr));
        G1Field { consts: g1_consts, base_field: fp, fr_field: fr }
    })
}

pub fn g2_bls24_509() -> &'static G2Field<10,7,8,10>
{
    G2_BLS24_509.get_or_init(|| {
        let g1 = g1_bls24_509();
        let fp = g1.base_field;
        let fr = g1.fr_field;
        let frob = BLS24_509_FROB_CONSTS
            .get_or_init(|| build_frob_params(BLS24_509_PARAMS, fp));
        let fp4 = fp4_bls24_509();
        let r_params = BLS24_509_R_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS24_509_PARAMS.r_as_strhex));
        let g2_consts = BLS24_509_G2_CONSTS.get_or_init(|| build_g2_params(
            BLS24_509_PARAMS,
            &ExtG2Field::Fp4_1(fp4),
            &PrimeField::new(r_params),
            frob,
        ));
        let exfp4 = EXFP4_BLS24_509.get_or_init(|| ExtG2Field::Fp4_1(fp4));
        G2Field { consts: g2_consts, base_field: exfp4, fr_field: fr }
    })
}

pub fn fp24_bls24_509()-> &'static Fp24Field_1<10, 8>
{
    FP24_BLS24_509.get_or_init(|| {
        let field_params = BLS24_509_FIELD_PARAMS
            .get_or_init(|| build_field_params(BLS24_509_PARAMS.modulo_as_strhex));
        let fp = FP_BLS24_509.get_or_init(|| PrimeField::new(field_params));
        let frob = BLS24_509_FROB_CONSTS
            .get_or_init(|| build_frob_params(BLS24_509_PARAMS, fp));
        Fp24Field_1::new(fp, Some(frob))
    })
}

pub fn gt_bls24_509() -> &'static GTField<8,10>
{
    GT_BLS24_509.get_or_init(|| GTField::Fp24_1(fp24_bls24_509()))
}
