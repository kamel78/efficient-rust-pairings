mod data;
extern crate bls12_381;

use std::{arch::x86_64::_rdtsc, time::Instant};

use ark_ec::ProjectiveCurve;
use bls12_381::*;
use pairings::{fields::prime_fields::{Endianness, FieldElement}, BLS24Curves, BLS48Curves, Bls12Curves};
use data::BENCH_DATA;

extern crate ark_std;

use ark_bls12_381::{G1Projective, Fr};
use ark_ff::{Fp256, PrimeField};


fn main() {
     let en = pairings::BLS12::_381();
     let mut zk_scalars: Vec<Scalar> = vec![Scalar::zero(); 1000];
     let mut ark_scalars: Vec<Fp256<ark_bls12_381::FrParameters>> = vec![Fr::from(0); 1000];
     let mut prop_scalars: Vec<FieldElement<4>> = vec![en.fr.zero(); 1000];
    

    for i in 0..1000 {
        zk_scalars[i] = Scalar::from_bytes(&BENCH_DATA[i]).unwrap();
        ark_scalars[i] = Fr::from_le_bytes_mod_order(&BENCH_DATA[i]);
        prop_scalars[i] = en.fr.from_byte_array(&BENCH_DATA[i], Endianness::Big);
    }   


    let zk_g = G1Affine::generator();
    let prop_g = en.g1.default_generator();
    let ark_g = G1Projective::prime_subgroup_generator();


     // Benchmarking GLV multiplication on the Ark-library
    println!("{}","---------------------------------------------------------------");
    println!("{}","Benchmarking Scalar multiplication on the ArK-library (BLS12-381):");
    println!("{}","---------------------------------------------------------------");
    let start = Instant::now();
    for i in &ark_scalars { let _ = ark_g.mul(i.into_repr()); }
    let duration = start.elapsed()/1000;    
    println!("Time elapsed : {:?}", duration);
    let mut acc_cycles:u64 = 0;
    for i in &ark_scalars {   unsafe {    let start_cycles = _rdtsc();
                                                              let _ = ark_g.mul(i.into_repr()); 
                                                              let end_cycles = _rdtsc();
                                                              acc_cycles += end_cycles-start_cycles;
                                                         }
                                              }
                                              
                                              println!("CPU cycles taken: {}", (acc_cycles/(zk_scalars.len() as u64)));
    // Benchmarking GLV multiplication on the ZK-library
    println!("{}","---------------------------------------------------------------");
    println!("{}","Benchmarking GLV multiplication on the ZK-library (BLS12-381):");
    println!("{}","---------------------------------------------------------------");
    let start = Instant::now();
    for i in &zk_scalars { let _ = i * zk_g; }
    let duration = start.elapsed()/1000;    
    println!("Time elapsed : {:?}", duration);
    let mut acc_cycles:u64 = 0;
    for i in &zk_scalars {   unsafe {    let start_cycles = _rdtsc();
                                                let _ = i * zk_g; 
                                                let end_cycles = _rdtsc();
                                                acc_cycles += end_cycles-start_cycles;
                                           }
                                }
    println!("CPU cycles taken: {}", (acc_cycles/1000));

    // Benchmarking GLV multiplication on the Proposed-library (BLS12-381)
    println!("{}","---------------------------------------------------------------");
    println!("{}","Benchmarking GLV multiplication on the Proposed-Library (BLS12-381):");
    println!("{}","---------------------------------------------------------------");
    let start = Instant::now();
    for i in &prop_scalars { let _ =*i * prop_g; }
    let duration = start.elapsed()/1000;    
    println!("Time elapsed : {:?}", duration);
    acc_cycles = 0;
    for i in &prop_scalars {   unsafe {    let start_cycles = _rdtsc();
                                                            let _ = *i * prop_g; 
                                                            let end_cycles = _rdtsc();
                                                            acc_cycles += end_cycles-start_cycles;
                                                       }
                                            }
    println!("CPU cycles taken : {}", (acc_cycles/(prop_scalars.len() as u64)));
    
    // Benchmarking GLV multiplication on the Proposed-library (BLS12-461)
    println!("{}","---------------------------------------------------------------");
    println!("{}","Benchmarking GLV multiplication on the Proposed-Library (BLS12-461):");
    println!("{}","---------------------------------------------------------------");
    let prop_g = pairings::BLS12::_461().g1.default_generator();
    let mut prop_scalars:[FieldElement<5>;1000] = [pairings::BLS12::_461().fr.zero();1000];    
    for i in 0..1000 { prop_scalars[i] = pairings::BLS12::_461().fr.random_element() };
    let start = Instant::now();
    for i in prop_scalars { let _ =i * prop_g; }
    let duration = start.elapsed()/prop_scalars.len().try_into().unwrap();    
    println!("Time elapsed : {:?}", duration);
    acc_cycles = 0;
    for i in prop_scalars {   unsafe {    let start_cycles = _rdtsc();
                                                            let _ = i * prop_g; 
                                                            let end_cycles = _rdtsc();
                                                            acc_cycles += end_cycles-start_cycles;
                                                       }
                                            }
    println!("CPU cycles taken : {}", (acc_cycles/(prop_scalars.len() as u64)));

    // Benchmarking GLV multiplication on the Proposed-library (BLS12-446)
    println!("{}","---------------------------------------------------------------");
    println!("{}","Benchmarking GLV multiplication on the Proposed-Library (BLS12-446):");
    println!("{}","---------------------------------------------------------------");
    let prop_g = pairings::BLS12::_446().g1.default_generator();
    let mut prop_scalars:Vec<FieldElement<5>> =  vec![pairings::BLS12::_446().fr.zero();1000];    
    for i in 0..1000 { prop_scalars[i] = pairings::BLS12::_446().fr.random_element() };
    let start = Instant::now();
    for i in &prop_scalars { let _ =*i * prop_g; }
    let duration = start.elapsed()/1000;    
    println!("Time elapsed : {:?}", duration);
    acc_cycles = 0;
    for i in &prop_scalars {   unsafe {    let start_cycles = _rdtsc();
                                                            let _ = *i * prop_g; 
                                                            let end_cycles = _rdtsc();
                                                            acc_cycles += end_cycles-start_cycles;
                                                       }
                                            }
    println!("CPU cycles taken : {}", (acc_cycles/(prop_scalars.len() as u64)));
    
     // Benchmarking GLV multiplication on the Proposed-library (BLS24-479)
     println!("{}","---------------------------------------------------------------");
     println!("{}","Benchmarking GLV multiplication on the Proposed-Library (BLS24-479):");
     println!("{}","---------------------------------------------------------------");
     let prop_g = pairings::BLS24::_479().g1.default_generator();
     let mut prop_scalars:Vec<FieldElement<7>> = vec![pairings::BLS24::_479().fr.zero();1000];    
     for i in 0..1000 { prop_scalars[i] = pairings::BLS24::_479().fr.random_element() };
     let start = Instant::now();
     for i in &prop_scalars { let _ =*i * prop_g; }
     let duration = start.elapsed()/prop_scalars.len().try_into().unwrap();    
     println!("Time elapsed : {:?}", duration);
     acc_cycles = 0;
     for i in &prop_scalars {   unsafe {    let start_cycles = _rdtsc();
                                                             let _ = *i * prop_g; 
                                                             let end_cycles = _rdtsc();
                                                             acc_cycles += end_cycles-start_cycles;
                                                        }
                                             }
     println!("CPU cycles taken : {}", (acc_cycles/(prop_scalars.len() as u64)));

     // Benchmarking GLV multiplication on the Proposed-library (BLS24-315)
     println!("{}","---------------------------------------------------------------");
     println!("{}","Benchmarking GLV multiplication on the Proposed-Library (BLS24-315):");
     println!("{}","---------------------------------------------------------------");
     let prop_g = pairings::BLS24::_315().g1.default_generator();
     let mut prop_scalars:Vec<FieldElement<4>>= vec![pairings::BLS24::_315().fr.zero();1000];    
     for i in 0..1000 { prop_scalars[i] = pairings::BLS24::_315().fr.random_element() };
     let start = Instant::now();
     for i in &prop_scalars { let _ =*i * prop_g; }
     let duration = start.elapsed()/1000;    
     println!("Time elapsed : {:?}", duration);
     acc_cycles = 0;
     for i in &prop_scalars {   unsafe {    let start_cycles = _rdtsc();
                                                            let _ = *i * prop_g; 
                                                            let end_cycles = _rdtsc();
                                                            acc_cycles += end_cycles-start_cycles;
                                                       }
                                             }
     println!("CPU cycles taken : {}", (acc_cycles/(prop_scalars.len() as u64)));

      // Benchmarking GLV multiplication on the Proposed-library (BLS24-477)
      println!("{}","---------------------------------------------------------------");
      println!("{}","Benchmarking GLV multiplication on the Proposed-Library (BLS24-477):");
      println!("{}","---------------------------------------------------------------");
      let prop_g = pairings::BLS24::_477().g1.default_generator();
      let mut prop_scalars:Vec<FieldElement<6>> = vec![pairings::BLS24::_477().fr.zero();1000];    
      for i in 0..1000 { prop_scalars[i] = pairings::BLS24::_477().fr.random_element() };
      let start = Instant::now();
      for i in &prop_scalars { let _ =*i * prop_g; }
      let duration = start.elapsed()/1000;    
      println!("Time elapsed : {:?}", duration);
      acc_cycles = 0;
      for i in &prop_scalars {   unsafe {    let start_cycles = _rdtsc();
                                                             let _ = *i * prop_g; 
                                                             let end_cycles = _rdtsc();
                                                             acc_cycles += end_cycles-start_cycles;
                                                        }
                                              }
      println!("CPU cycles taken : {}", (acc_cycles/(prop_scalars.len() as u64)));

      // Benchmarking GLV multiplication on the Proposed-library (BLS24-509_snark)
      println!("{}","---------------------------------------------------------------");
      println!("{}","Benchmarking GLV multiplication on the Proposed-Library (BLS24-509-Snark):");
      println!("{}","---------------------------------------------------------------");
      let prop_g = pairings::BLS24::_509_snark().g1.default_generator();
      let mut prop_scalars:Vec<FieldElement<7>> = vec![pairings::BLS24::_509_snark().fr.zero();1000];    
      for i in 0..1000 { prop_scalars[i] = pairings::BLS24::_509_snark().fr.random_element() };
      let start = Instant::now();
      for i in &prop_scalars { let _ =*i * prop_g; }
      let duration = start.elapsed()/prop_scalars.len().try_into().unwrap();    
      println!("Time elapsed : {:?}", duration);
      acc_cycles = 0;
      for i in &prop_scalars {   unsafe {    let start_cycles = _rdtsc();
                                                             let _ = *i * prop_g; 
                                                             let end_cycles = _rdtsc();
                                                             acc_cycles += end_cycles-start_cycles;
                                                        }
                                              }
      println!("CPU cycles taken : {}", (acc_cycles/(prop_scalars.len() as u64)));
    
     // Benchmarking GLV multiplication on the Proposed-library (BLS24-509_)
     println!("{}","---------------------------------------------------------------");
     println!("{}","Benchmarking GLV multiplication on the Proposed-Library (BLS24-509):");
     println!("{}","---------------------------------------------------------------");
     let prop_g = pairings::BLS24::_509().g1.default_generator();
     let mut prop_scalars:Vec<FieldElement<7>> = vec![pairings::BLS24::_509().fr.zero();1000];    
     for i in 0..1000 { prop_scalars[i] = pairings::BLS24::_509().fr.random_element() };
     let start = Instant::now();
     for i in &prop_scalars { let _ =*i * prop_g; }
     let duration = start.elapsed()/1000;    
     println!("Time elapsed : {:?}", duration);
     acc_cycles = 0;
     for i in &prop_scalars {   unsafe {    let start_cycles = _rdtsc();
                                                            let _ = *i * prop_g; 
                                                            let end_cycles = _rdtsc();
                                                            acc_cycles += end_cycles-start_cycles;
                                                       }
                                             }
     println!("CPU cycles taken : {}", (acc_cycles/(prop_scalars.len() as u64)));

      // Benchmarking GLV multiplication on the Proposed-library (BLS24-559)
      println!("{}","---------------------------------------------------------------");
      println!("{}","Benchmarking GLV multiplication on the Proposed-Library (BLS24-559):");
      println!("{}","---------------------------------------------------------------");
      let prop_g = pairings::BLS24::_559().g1.default_generator();
      let mut prop_scalars:Vec<FieldElement<8>> = vec![pairings::BLS24::_559().fr.zero();1000];    
      for i in 0..1000 { prop_scalars[i] = pairings::BLS24::_559().fr.random_element() };
      let start = Instant::now();
      for i in &prop_scalars { let _ =*i * prop_g; }
      let duration = start.elapsed()/1000;    
      println!("Time elapsed : {:?}", duration);
      acc_cycles = 0;
      for i in &prop_scalars {   unsafe {    let start_cycles = _rdtsc();
                                                              let _ = *i * prop_g; 
                                                              let end_cycles = _rdtsc();
                                                              acc_cycles += end_cycles-start_cycles;
                                                         }
                                              }
      println!("CPU cycles taken : {}", (acc_cycles/(prop_scalars.len() as u64)));

      // Benchmarking GLV multiplication on the Proposed-library (BLS48-575)
      println!("{}","---------------------------------------------------------------");
      println!("{}","Benchmarking GLV multiplication on the Proposed-Library (BLS48-575):");
      println!("{}","---------------------------------------------------------------");
      let prop_g = pairings::BLS48::_575().g1.default_generator();
      let mut prop_scalars:Vec<FieldElement<9>> = vec![pairings::BLS48::_575().fr.zero();1000];    
      for i in 0..1000 { prop_scalars[i] = pairings::BLS48::_575().fr.random_element() };
      let start = Instant::now();
      for i in &prop_scalars { let _ =*i * prop_g; }
      let duration = start.elapsed()/prop_scalars.len().try_into().unwrap();    
      println!("Time elapsed : {:?}", duration);
      acc_cycles = 0;
      for i in &prop_scalars {   unsafe {    let start_cycles = _rdtsc();
                                                              let _ = *i * prop_g; 
                                                              let end_cycles = _rdtsc();
                                                              acc_cycles += end_cycles-start_cycles;
                                                         }
                                              }
      println!("CPU cycles taken : {}", (acc_cycles/(prop_scalars.len() as u64)));

      // Benchmarking GLV multiplication on the Proposed-library (BLS48-581)
      println!("{}","---------------------------------------------------------------");
      println!("{}","Benchmarking GLV multiplication on the Proposed-Library (BLS48-581):");
      println!("{}","---------------------------------------------------------------");
      let prop_g = pairings::BLS48::_581().g1.default_generator();
      let mut prop_scalars:Vec<FieldElement<9>>= vec![pairings::BLS48::_581().fr.zero();1000];    
      for i in 0..1000 { prop_scalars[i] = pairings::BLS48::_581().fr.random_element() };
      let start = Instant::now();
      for i in &prop_scalars { let _ =*i * prop_g; }
      let duration = start.elapsed()/1000;    
      println!("Time elapsed : {:?}", duration);
      acc_cycles = 0;
      for i in &prop_scalars {   unsafe {    let start_cycles = _rdtsc();
                                                              let _ = *i * prop_g; 
                                                              let end_cycles = _rdtsc();
                                                              acc_cycles += end_cycles-start_cycles;
                                                         }
                                              }
      println!("CPU cycles taken : {}", (acc_cycles/(prop_scalars.len() as u64)));

      // Benchmarking GLV multiplication on the Proposed-library (BLS48-573)
      println!("{}","---------------------------------------------------------------");
      println!("{}","Benchmarking GLV multiplication on the Proposed-Library (BLS48-573):");
      println!("{}","---------------------------------------------------------------");
      let prop_g = pairings::BLS48::_573().g1.default_generator();
      let mut prop_scalars:Vec<FieldElement<8>>= vec![pairings::BLS48::_573().fr.zero();1000];    
      for i in 0..1000 { prop_scalars[i] = pairings::BLS48::_573().fr.random_element() };
      let start = Instant::now();
      for i in &prop_scalars { let _ =*i * prop_g; }
      let duration = start.elapsed()/1000;    
      println!("Time elapsed : {:?}", duration);
      acc_cycles = 0;
      for i in &prop_scalars {   unsafe {    let start_cycles = _rdtsc();
                                                              let _ = *i * prop_g; 
                                                              let end_cycles = _rdtsc();
                                                              acc_cycles += end_cycles-start_cycles;
                                                         }
                                              }
      println!("CPU cycles taken : {}", (acc_cycles/(prop_scalars.len() as u64)));

      // Benchmarking GLV multiplication on the Proposed-library (BLS48-571)
      println!("{}","---------------------------------------------------------------");
      println!("{}","Benchmarking GLV multiplication on the Proposed-Library (BLS48-571):");
      println!("{}","---------------------------------------------------------------");
      let prop_g = pairings::BLS48::_571().g1.default_generator();
      let mut prop_scalars:Vec<FieldElement<8>> = vec![pairings::BLS48::_571().fr.zero();1000];    
      for i in 0..1000 { prop_scalars[i] = pairings::BLS48::_571().fr.random_element() };
      let start = Instant::now();
      for i in &prop_scalars { let _ =*i * prop_g; }
      let duration = start.elapsed()/1000;    
      println!("Time elapsed : {:?}", duration);
      acc_cycles = 0;
      for i in &prop_scalars {   unsafe {    let start_cycles = _rdtsc();
                                                              let _ = *i * prop_g; 
                                                              let end_cycles = _rdtsc();
                                                              acc_cycles += end_cycles-start_cycles;
                                                         }
                                              }
      println!("CPU cycles taken : {}", (acc_cycles/(prop_scalars.len() as u64)));

      // Benchmarking GLV multiplication on the Proposed-library (BLS48-287)
      println!("{}","---------------------------------------------------------------");
      println!("{}","Benchmarking GLV multiplication on the Proposed-Library (BLS48-287):");
      println!("{}","---------------------------------------------------------------");
      let prop_g = pairings::BLS48::_287().g1.default_generator();
      let mut prop_scalars:Vec<FieldElement<4>> = vec![pairings::BLS48::_287().fr.zero();1000];    
      for i in 0..1000 { prop_scalars[i] = pairings::BLS48::_287().fr.random_element() };
      let start = Instant::now();
      for i in &prop_scalars { let _ =*i * prop_g; }
      let duration = start.elapsed()/1000;    
      println!("Time elapsed : {:?}", duration);
      acc_cycles = 0;
      for i in &prop_scalars {   unsafe {    let start_cycles = _rdtsc();
                                                              let _ = *i * prop_g; 
                                                              let end_cycles = _rdtsc();
                                                              acc_cycles += end_cycles-start_cycles;
                                                         }
                                              }
      println!("CPU cycles taken : {}", (acc_cycles/(prop_scalars.len() as u64)));

     // Benchmarking GLV multiplication on the Proposed-library (BLS48-277)
     println!("{}","---------------------------------------------------------------");
     println!("{}","Benchmarking GLV multiplication on the Proposed-Library (BLS48-277):");
     println!("{}","---------------------------------------------------------------");
     let prop_g = pairings::BLS48::_277().g1.default_generator();
     let mut prop_scalars:Vec<FieldElement<4>>= vec![pairings::BLS48::_277().fr.zero();1000];    
     for i in 0..1000 { prop_scalars[i] = pairings::BLS48::_277().fr.random_element() };
     let start = Instant::now();
     for i in &prop_scalars { let _ =*i * prop_g; }
     let duration = start.elapsed()/prop_scalars.len().try_into().unwrap();    
     println!("Time elapsed : {:?}", duration);
     acc_cycles = 0;
     for i in &prop_scalars {   unsafe {    let start_cycles = _rdtsc();
                                                            let _ = *i * prop_g; 
                                                            let end_cycles = _rdtsc();
                                                            acc_cycles += end_cycles-start_cycles;
                                                       }
                                             }
     println!("CPU cycles taken : {}", (acc_cycles/(prop_scalars.len() as u64)));
}
