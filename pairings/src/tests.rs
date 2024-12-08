use std::{arch::x86_64::_rdtsc, time::{Duration, Instant}}; 

use crate::{tools::arithmetic_interface::ArithmeticOperations, BLS24Curves, BLS48Curves, Bls12Curves, CurvesNames, Pairings, PairingsEngine, BLS12, BLS24, BLS48};

fn measure_time<F>(f: F) -> Duration
where
    F: Fn(),            
{   let trys_count = 1000;
    let start = Instant::now();
    for _ in 0..trys_count {f();}
    start.elapsed()/trys_count
}

fn measure_cycles<F>(f:F) -> u64
where 
  F:Fn()
{ let mut acc_cycles = 0;
  let trys_count = 100;
  for _ in 0..trys_count {  unsafe {let start_cycles =  _rdtsc(); 
                                    f();
                                    let end_cycles = _rdtsc();
                                    acc_cycles += end_cycles-start_cycles;
                                  }
                          }
  acc_cycles/(trys_count as u64)
}


pub fn check_pairing_for_curve <const R:usize,const N:usize,const MAX_COEFS_COUNT1:usize,const PRAMASIZE:usize, const MAX_COEFS_COUNT2 :usize>
      (engine :&Pairings <R, N, MAX_COEFS_COUNT1, PRAMASIZE, MAX_COEFS_COUNT2>)->bool
      {
        let s1= engine.fr.random_element();
        let s2= engine.fr.random_element();
        let p = engine.g1.random_point();
        let q = engine.g2.random_point();
        let e1 = engine.paire(&(s1*p), &(s2*q));
        let e2 =engine.paire(&(s2*p), &(s1*q));
        e1.equal(&e2) & !e1.equal(&engine.gt.one())
      }

pub fn bench_pairing_for_curve <const R:usize,const N:usize,const MAX_COEFS_COUNT1:usize,const PRAMASIZE:usize, const MAX_COEFS_COUNT2 :usize>
      (engine :&Pairings <R, N, MAX_COEFS_COUNT1, PRAMASIZE, MAX_COEFS_COUNT2>)
      {
        println!("------- {} --------------------------------",engine.identifier);        
        let seed = engine.g1.base_field.random_element();
        let mut duration = measure_time(|| {let _ = engine.g1.map_to_curve(seed);} );
        let mut cycles = measure_cycles(|| {let _ = engine.g1.map_to_curve(seed);} );
        println!("Time elapsed / CPU cycles(*10^3 Cycle) for Mapping to curve on G1  (SSWU) : {:?}/{}", duration,(cycles/1000) as u32);
        duration = measure_time(|| {let _ = engine.g1.hash_to_field("Some identity",0);} );
        cycles = measure_cycles(|| {let _ = engine.g1.hash_to_field("Some identity",0);} );
        println!("Time elapsed / CPU cycles(*10^3 Cycle) for Hashing to G1's Torsion sub-group : {:?}/{}", duration,(cycles/1000) as u32);
        let scalar = engine.fr.random_element();
        let p =engine.g1.random_point();
        duration = measure_time(|| {let _ = p.multiply(&scalar);});
        cycles = measure_cycles(|| {let _ = p.multiply(&scalar);});
        println!("Time elapsed / CPU cycles(*10^3 Cycle) for multiplication with scalar on G1 (w-sliding window (w=3)) : {:?}/{}", duration,(cycles/1000) as u32);
        duration = measure_time(|| {let _ = p.glv_multiply(&scalar);});
        cycles = measure_cycles(|| {let _ = p.glv_multiply(&scalar);});
        println!("Time elapsed / CPU cycles(*10^3 Cycle) for multiplication with scalar on G1 (GLV multiplication) : {:?}/{}", duration,(cycles/1000) as u32);
        
        let seed = engine.g2.base_field.random_element();
        duration = measure_time(|| {let _ = engine.g2.map_to_curve(&seed);} );
        cycles = measure_cycles(|| {let _ = engine.g2.map_to_curve(&seed);} );
        println!("Time elapsed / CPU cycles(*10^3 Cycle) for Mapping to curve on G2  (SSWU) : {:?}/{}", duration,(cycles/1000) as u32);
        duration = measure_time(|| {let _ = engine.g2.random_point_trys();} );
        cycles = measure_cycles(|| {let _ = engine.g2.random_point_trys();} );
        println!("Time elapsed / CPU cycles(*10^3 Cycle) for Mapping to curve on G2  (Random) : {:?}/{}", duration,(cycles/1000) as u32);
        duration = measure_time(|| {let _ = engine.g2.hash_to_field("Some identity",0);} );
        cycles = measure_cycles(|| {let _ = engine.g2.hash_to_field("Some identity",0);} );
        println!("Time elapsed / CPU cycles(*10^3 Cycle) for Hashing to G2's Torsion sub-group : {:?}/{}", duration,(cycles/1000) as u32);
        let scalar = engine.fr.random_element();
        let q =engine.g2.random_point();
        duration = measure_time(|| {let _ = q.multiply(&scalar);});
        cycles = measure_cycles(|| {let _ = q.multiply(&scalar);});
        println!("Time elapsed / CPU cycles(*10^3 Cycle) for multiplication with scalar on G2 (w-sliding window (w=3)) : {:?}/{}", duration,(cycles/1000) as u32);
        duration = measure_time(|| {let _ = q.multiply_gls(&scalar);});
        cycles = measure_cycles(|| {let _ = q.multiply_gls(&scalar);});
        println!("Time elapsed / CPU cycles(*10^3 Cycle) for multiplication with scalar on G2 (GLS multiplication) : {:?}/{}", duration,(cycles/1000) as u32);

        duration = measure_time(|| {engine.miller_loop(&p, &q);});
        cycles = measure_cycles(|| {engine.miller_loop(&p, &q);});
        println!("Time elapsed / CPU cycles(*10^3 Cycle) for Miller loop : {:?}/{}", duration,(cycles/1000) as u32);
        duration = measure_time(|| {engine.paire(&p, &q);});
        cycles = measure_cycles(|| {engine.paire(&p, &q);});
        println!("Time elapsed / CPU cycles(*10^3 Cycle) for pairings : {:?}/{}", duration,(cycles/1000) as u32);
        let e = engine.miller_loop(&p, &q);
        duration = measure_time(|| {let _ =e.final_exponentiation();});
        cycles = measure_cycles(|| {let _ =e.final_exponentiation();});
        println!("Time elapsed / CPU cycles(*10^3 Cycle) for Final Exponentiation: {:?}/{}", duration,(cycles/1000) as u32);        
      }

pub fn check_pairings(curve :&CurvesNames)->bool
{
  match curve {
    CurvesNames::Bls12_381 => { check_pairing_for_curve(BLS12::_381())},
    CurvesNames::Bls12_446 => { check_pairing_for_curve(BLS12::_446())},
    CurvesNames::Bls12_461 => { check_pairing_for_curve(BLS12::_461())},
    CurvesNames::Bls24_315 => { check_pairing_for_curve(BLS24::_315())},
    CurvesNames::Bls24_477 => { check_pairing_for_curve(BLS24::_477())},
    CurvesNames::Bls24_479 => { check_pairing_for_curve(BLS24::_479())},
    CurvesNames::Bls24_509 => { check_pairing_for_curve(BLS24::_509())},
    CurvesNames::Bls24_509Snark =>{ check_pairing_for_curve(BLS24::_509_snark())},
    CurvesNames::Bls24_559 => { check_pairing_for_curve(BLS24::_559())},
    CurvesNames::Bls48_277 => { check_pairing_for_curve(BLS48::_277())},
    CurvesNames::Bls48_287 => { check_pairing_for_curve(BLS48::_287())} ,
    CurvesNames::Bls48_571 => { check_pairing_for_curve(BLS48::_571())},
    CurvesNames::Bls48_573 => { check_pairing_for_curve(BLS48::_573())},
    CurvesNames::Bls48_575 => { check_pairing_for_curve(BLS48::_575())},
    CurvesNames::Bls48_581 => { check_pairing_for_curve(BLS48::_581())},
  } 
  
}

pub fn bench_pairings(curve :&CurvesNames)
{
  match curve {
    CurvesNames::Bls12_381 => { bench_pairing_for_curve(BLS12::_381());},
    CurvesNames::Bls12_446 => { bench_pairing_for_curve(BLS12::_446());},
    CurvesNames::Bls12_461 => { bench_pairing_for_curve(BLS12::_461());},
    CurvesNames::Bls24_315 => { bench_pairing_for_curve(BLS24::_315())},
    CurvesNames::Bls24_477 => { bench_pairing_for_curve(BLS24::_477())},
    CurvesNames::Bls24_479 => { bench_pairing_for_curve(BLS24::_479())},
    CurvesNames::Bls24_509 => { bench_pairing_for_curve(BLS24::_509())},
    CurvesNames::Bls24_509Snark =>{ bench_pairing_for_curve(BLS24::_509_snark())},
    CurvesNames::Bls24_559 => { bench_pairing_for_curve(BLS24::_559())},
    CurvesNames::Bls48_277 => { bench_pairing_for_curve(BLS48::_277())},
    CurvesNames::Bls48_287 => { bench_pairing_for_curve(BLS48::_287())} ,
    CurvesNames::Bls48_571 => { bench_pairing_for_curve(BLS48::_571())},
    CurvesNames::Bls48_573 => { bench_pairing_for_curve(BLS48::_573())},
    CurvesNames::Bls48_575 => { bench_pairing_for_curve(BLS48::_575())},
    CurvesNames::Bls48_581 => { bench_pairing_for_curve(BLS48::_581())},
  } 
  
}

