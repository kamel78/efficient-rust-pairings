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
    let mut zk_scalars :[Scalar;1000] = [Scalar::zero();1000];
    let mut ark_scalars :[Fp256<ark_bls12_381::FrParameters>;1000] = [Fr::from(0);1000];
    let mut prop_scalars:[FieldElement<4>;1000]=[en.fr.zero();1000];
    

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
    for i in ark_scalars { let _ = ark_g.mul(i.into_repr()); }
    let duration = start.elapsed()/ark_scalars.len().try_into().unwrap();    
    println!("Time elapsed : {:?}", duration);
    let mut acc_cycles:u64 = 0;
    for i in ark_scalars {   unsafe {    let start_cycles = _rdtsc();
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
    for i in zk_scalars { let _ = i * zk_g; }
    let duration = start.elapsed()/zk_scalars.len().try_into().unwrap();    
    println!("Time elapsed : {:?}", duration);
    let mut acc_cycles:u64 = 0;
    for i in zk_scalars {   unsafe {    let start_cycles = _rdtsc();
                                                let _ = i * zk_g; 
                                                let end_cycles = _rdtsc();
                                                acc_cycles += end_cycles-start_cycles;
                                           }
                                }
    println!("CPU cycles taken: {}", (acc_cycles/(zk_scalars.len() as u64)));

    // Benchmarking GLV multiplication on the Proposed-library (BLS12-381)
    println!("{}","---------------------------------------------------------------");
    println!("{}","Benchmarking GLV multiplication on the Proposed-Library (BLS12-381):");
    println!("{}","---------------------------------------------------------------");
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
    
     // Benchmarking GLV multiplication on the Proposed-library (BLS24-479)
     println!("{}","---------------------------------------------------------------");
     println!("{}","Benchmarking GLV multiplication on the Proposed-Library (BLS24-479):");
     println!("{}","---------------------------------------------------------------");
     let prop_g = pairings::BLS24::_479().g1.default_generator();
     let mut prop_scalars:[FieldElement<7>;1000] = [pairings::BLS24::_479().fr.zero();1000];    
     for i in 0..1000 { prop_scalars[i] = pairings::BLS24::_479().fr.random_element() };
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

      // Benchmarking GLV multiplication on the Proposed-library (BLS24-559)
      println!("{}","---------------------------------------------------------------");
      println!("{}","Benchmarking GLV multiplication on the Proposed-Library (BLS24-559):");
      println!("{}","---------------------------------------------------------------");
      let prop_g = pairings::BLS24::_559().g1.default_generator();
      let mut prop_scalars:[FieldElement<8>;1000] = [pairings::BLS24::_559().fr.zero();1000];    
      for i in 0..1000 { prop_scalars[i] = pairings::BLS24::_559().fr.random_element() };
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

      // Benchmarking GLV multiplication on the Proposed-library (BLS48-575)
      println!("{}","---------------------------------------------------------------");
      println!("{}","Benchmarking GLV multiplication on the Proposed-Library (BLS48-575):");
      println!("{}","---------------------------------------------------------------");
      let prop_g = pairings::BLS48::_575().g1.default_generator();
      let mut prop_scalars:[FieldElement<9>;1000] = [pairings::BLS48::_575().fr.zero();1000];    
      for i in 0..1000 { prop_scalars[i] = pairings::BLS48::_575().fr.random_element() };
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
}
