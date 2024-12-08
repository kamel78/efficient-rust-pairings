// Code developed by FARAOUN Kamel Mohamed.
// EEDIS-Laboratory. UDL-University. Algeria
// During May 2024.

use once_cell::sync::OnceCell;
use crate::{CurvesNames, Pairings};
use super::{parameters::preloadedconfs::bls12::*, 
            parameters::preloadedconfs::bls24::*, 
            parameters::preloadedconfs::bls48::*};

static BLS12_381_ENGINE       : OnceCell<Pairings<4, 6, 16, 4, 4>>    = OnceCell::new();    
static BLS12_446_ENGINE       : OnceCell<Pairings<5, 7, 3, 4, 7>>    = OnceCell::new();    
static BLS12_461_ENGINE       : OnceCell<Pairings<5, 8, 10, 4, 4>>    = OnceCell::new();    
static BLS24_315_ENGINE       : OnceCell<Pairings<4, 5, 3, 11, 7>>    = OnceCell::new();    
static BLS24_479_ENGINE       : OnceCell<Pairings<7, 8, 4, 10, 10>>   = OnceCell::new();    
static BLS24_477_ENGINE       : OnceCell<Pairings<6, 8, 7, 11, 4>>    = OnceCell::new();
static BLS24_509_SNARK_ENGINE : OnceCell<Pairings<7, 8, 3, 11, 7>>    = OnceCell::new();    
static BLS24_509_ENGINE       : OnceCell<Pairings<7, 8, 3, 10, 10>>   = OnceCell::new();    
static BLS24_559_ENGINE       : OnceCell<Pairings<8, 9, 4, 10, 4>>    = OnceCell::new();    
static BLS48_575_ENGINE       : OnceCell<Pairings<9, 9, 4, 24, 7>>    = OnceCell::new();    
static BLS48_581_ENGINE       : OnceCell<Pairings<9, 10, 3, 24, 4>>    = OnceCell::new();    
static BLS48_573_ENGINE       : OnceCell<Pairings<8, 9, 7, 24, 4>>     = OnceCell::new();
static BLS48_571_ENGINE       : OnceCell<Pairings<8, 9, 3, 36, 4>>     = OnceCell::new();
static BLS48_287_ENGINE       : OnceCell<Pairings<4, 5, 4, 35, 19>>    = OnceCell::new();
static BLS48_277_ENGINE       : OnceCell<Pairings<4, 5, 10, 24, 4>>    = OnceCell::new();

pub fn bls12_381_engine() -> &'static Pairings<4, 6, 16, 4, 4>
{
    if BLS12_381_ENGINE.get().is_none() { BLS12_381_ENGINE.set(Pairings{identifier :"BLS12-381",curvename :CurvesNames::Bls12_381, g1 : g1_bls12_381(),g2:g2_bls12_381(),gt:gt_bls12_381(), fr :fr_bls12_381()}).unwrap()};
    BLS12_381_ENGINE.get().unwrap()
}

pub fn bls12_446_engine() -> &'static Pairings<5, 7, 3, 4, 7>
{
    if BLS12_446_ENGINE.get().is_none() { BLS12_446_ENGINE.set(Pairings{identifier :"BLS12-446",curvename :CurvesNames::Bls12_446, g1 : g1_bls12_446(),g2:g2_bls12_446(),gt:gt_bls12_446(), fr :fr_bls12_446()}).unwrap()};
    BLS12_446_ENGINE.get().unwrap()
}

pub fn bls12_461_engine() -> &'static Pairings<5, 8, 10, 4, 4>
{
    if BLS12_461_ENGINE.get().is_none() { BLS12_461_ENGINE.set(Pairings{identifier :"BLS12-461",curvename :CurvesNames::Bls12_461,g1 : g1_bls12_461(),g2:g2_bls12_461(),gt:gt_bls12_461(), fr :fr_bls12_461()}).unwrap()};
    BLS12_461_ENGINE.get().unwrap()
}

pub fn bls24_315_engine() -> &'static Pairings<4, 5, 3, 11, 7>
{
    if BLS24_315_ENGINE.get().is_none() { BLS24_315_ENGINE.set(Pairings{identifier :"BLS24-315",curvename :CurvesNames::Bls24_315,g1 : g1_bls24_315(),g2:g2_bls24_315(),gt:gt_bls24_315(),fr :fr_bls24_315()}).unwrap()};
    BLS24_315_ENGINE.get().unwrap()
}

pub fn bls24_477_engine() -> &'static Pairings<6, 8, 7, 11, 4>
{
    if BLS24_477_ENGINE.get().is_none() { BLS24_477_ENGINE.set(Pairings{identifier :"BLS24-477",curvename :CurvesNames::Bls24_477,g1 : g1_bls24_477(),g2:g2_bls24_477(),gt:gt_bls24_477(),fr :fr_bls24_477()}).unwrap()};
    BLS24_477_ENGINE.get().unwrap()
}

pub fn bls24_479_engine() -> &'static Pairings<7, 8, 4, 10, 10>
{
    if BLS24_479_ENGINE.get().is_none() { BLS24_479_ENGINE.set(Pairings{identifier :"BLS24-479",curvename :CurvesNames::Bls24_479,g1 : g1_bls24_479(),g2:g2_bls24_479(),gt:gt_bls24_479(),fr :fr_bls24_479()}).unwrap()};
    BLS24_479_ENGINE.get().unwrap()
}

pub fn bls24_509_snark_engine() -> &'static Pairings<7, 8, 3, 11, 7>
{
    if BLS24_509_SNARK_ENGINE.get().is_none() { BLS24_509_SNARK_ENGINE.set(Pairings{identifier :"BLS24-509-SNARK",curvename :CurvesNames::Bls24_509Snark,g1 : g1_bls24_509_snark(),g2:g2_bls24_509_snark(),gt:gt_bls24_509_snark(),fr :fr_bls24_509_snark()}).unwrap()};
    BLS24_509_SNARK_ENGINE.get().unwrap()
}

pub fn bls24_509_engine() -> &'static Pairings<7, 8, 3, 10, 10>
{
    if BLS24_509_ENGINE.get().is_none() { BLS24_509_ENGINE.set(Pairings{identifier :"BLS24-509",curvename :CurvesNames::Bls24_509,g1 : g1_bls24_509(),g2:g2_bls24_509(),gt:gt_bls24_509(),fr :fr_bls24_509()}).unwrap()};
    BLS24_509_ENGINE.get().unwrap()
}

pub fn bls24_559_engine() -> &'static Pairings<8, 9, 4, 10, 4>
{
    if BLS24_559_ENGINE.get().is_none() { BLS24_559_ENGINE.set(Pairings{identifier :"BLS24-559",curvename :CurvesNames::Bls24_559,g1 : g1_bls24_559(),g2:g2_bls24_559(),gt:gt_bls24_559(), fr :fr_bls24_559()}).unwrap()};
    BLS24_559_ENGINE.get().unwrap()
}

pub fn bls48_575_engine() -> &'static Pairings<9, 9, 4, 24, 7>
{
    if BLS48_575_ENGINE.get().is_none() { BLS48_575_ENGINE.set(Pairings{identifier :"BLS48-575",curvename :CurvesNames::Bls48_575,g1 : g1_bls48_575(),g2:g2_bls48_575(),gt:gt_bls48_575(), fr :fr_bls48_575()}).unwrap()};
    BLS48_575_ENGINE.get().unwrap()
}

pub fn bls48_581_engine() -> &'static Pairings<9, 10, 3, 24, 4>
{
    if BLS48_581_ENGINE.get().is_none() { BLS48_581_ENGINE.set(Pairings{identifier :"BLS48-581",curvename :CurvesNames::Bls48_581,g1 : g1_bls48_581(),g2:g2_bls48_581(),gt:gt_bls48_581(), fr :fr_bls48_581()}).unwrap()};
    BLS48_581_ENGINE.get().unwrap()
}

pub fn bls48_573_engine() -> &'static Pairings<8, 9, 7, 24, 4>
{
    if BLS48_573_ENGINE.get().is_none() { BLS48_573_ENGINE.set(Pairings{identifier :"BLS48-573",curvename :CurvesNames::Bls48_573,g1 : g1_bls48_573(),g2:g2_bls48_573(),gt:gt_bls48_573(), fr :fr_bls48_573()}).unwrap()};
    BLS48_573_ENGINE.get().unwrap()
}

pub fn bls48_571_engine() -> &'static Pairings<8, 9, 3, 36, 4>
{
    if BLS48_571_ENGINE.get().is_none() { BLS48_571_ENGINE.set(Pairings{identifier :"BLS48-571",curvename :CurvesNames::Bls48_571,g1 : g1_bls48_571(),g2:g2_bls48_571(),gt:gt_bls48_571(), fr :fr_bls48_571()}).unwrap()};
    BLS48_571_ENGINE.get().unwrap()
}

pub fn bls48_287_engine() -> &'static Pairings<4, 5, 4, 35, 19>
{
    if BLS48_287_ENGINE.get().is_none() { BLS48_287_ENGINE.set(Pairings{identifier :"BLS48-287",curvename :CurvesNames::Bls48_287,g1 : g1_bls48_287(),g2:g2_bls48_287(),gt:gt_bls48_287(), fr :fr_bls48_287()}).unwrap()};
    BLS48_287_ENGINE.get().unwrap()
}

pub fn bls48_277_engine() -> &'static Pairings<4, 5, 10, 24, 4>
{
    if BLS48_277_ENGINE.get().is_none() { BLS48_277_ENGINE.set(Pairings{identifier :"BLS48-277",curvename :CurvesNames::Bls48_277, g1 : g1_bls48_277(),g2:g2_bls48_277(),gt:gt_bls48_277(), fr :fr_bls48_277()}).unwrap()};
    BLS48_277_ENGINE.get().unwrap()
}