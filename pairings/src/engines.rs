// Code developed by FARAOUN Kamel Mohamed.
// EEDIS-Laboratory. UDL-University. Algeria
// During May 2024.

use once_cell::sync::OnceCell;
use crate::Pairings;
use super::{parameters::preloadedconfs::bls12::*, 
            parameters::preloadedconfs::bls24::*, 
            parameters::preloadedconfs::bls48::*};

static BLS12_381_ENGINE       : OnceCell<Pairings<4, 6, 16, 4, 4>>    = OnceCell::new();    
static BLS12_461_ENGINE       : OnceCell<Pairings<5, 8, 10, 4, 4>>    = OnceCell::new();    
static BLS24_479_ENGINE       : OnceCell<Pairings<7, 8, 4, 10, 10>>   = OnceCell::new();    
static BLS24_559_ENGINE       : OnceCell<Pairings<8, 9, 4, 10, 4>>    = OnceCell::new();    
static BLS48_575_ENGINE       : OnceCell<Pairings<9, 9, 4, 24, 7>>    = OnceCell::new();    

pub fn bls12_381_engine() -> &'static Pairings<4, 6, 16, 4, 4>
{
    if BLS12_381_ENGINE.get().is_none() { BLS12_381_ENGINE.set(Pairings{g1 : g1_bls12_381(),g2:g2_bls12_381(),gt:gt_bls12_381(), fr :fr_bls12_381()}).unwrap()};
    BLS12_381_ENGINE.get().unwrap()
}

pub fn bls12_461_engine() -> &'static Pairings<5, 8, 10, 4, 4>
{
    if BLS12_461_ENGINE.get().is_none() { BLS12_461_ENGINE.set(Pairings{g1 : g1_bls12_461(),g2:g2_bls12_461(),gt:gt_bls12_461(), fr :fr_bls12_461()}).unwrap()};
    BLS12_461_ENGINE.get().unwrap()
}

pub fn bls24_479_engine() -> &'static Pairings<7, 8, 4, 10, 10>
{
    if BLS24_479_ENGINE.get().is_none() { BLS24_479_ENGINE.set(Pairings{g1 : g1_bls24_479(),g2:g2_bls24_479(),gt:gt_bls24_479(),fr :fr_bls24_479()}).unwrap()};
    BLS24_479_ENGINE.get().unwrap()
}

pub fn bls24_559_engine() -> &'static Pairings<8, 9, 4, 10, 4>
{
    if BLS24_559_ENGINE.get().is_none() { BLS24_559_ENGINE.set(Pairings{g1 : g1_bls24_559(),g2:g2_bls24_559(),gt:gt_bls24_559(), fr :fr_bls24_559()}).unwrap()};
    BLS24_559_ENGINE.get().unwrap()
}

pub fn bls48_575_engine() -> &'static Pairings<9, 9, 4, 24, 7>
{
    if BLS48_575_ENGINE.get().is_none() { BLS48_575_ENGINE.set(Pairings{g1 : g1_bls48_575(),g2:g2_bls48_575(),gt:gt_bls48_575(), fr :fr_bls48_575()}).unwrap()};
    BLS48_575_ENGINE.get().unwrap()
}

