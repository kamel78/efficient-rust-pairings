mod data;
use std::sync::atomic::{AtomicUsize, Ordering};

use criterion::{black_box, criterion_group, criterion_main, Criterion};
extern crate bls12_381;


use pairings::{fields::prime_fields::Endianness, BLS24Curves, BLS48Curves, Bls12Curves, BLS12, BLS24, BLS48};
use data::BENCH_DATA;

extern crate ark_std;

use ark_ec::ProjectiveCurve;
use bls12_381::*;
use ark_bls12_381::{G1Projective, Fr};
use ark_ff::PrimeField;

fn bench_glv_multiply(c: &mut Criterion) {

    
    let zk_g = G1Affine::generator();    
    let ark_g = G1Projective::prime_subgroup_generator();

    // Loading engines
    let en_12_381 = BLS12::_381();
    let en_12_461 = BLS12::_461();
    let en_12_446 = BLS12::_446();
    
    let en_24_315 = BLS24::_315();
    let en_24_479 = BLS24::_479();
    let en_24_477 = BLS24::_477();
    let en_24_509_snark = BLS24::_509_snark();
    let en_24_509 = BLS24::_509();
    let en_24_559 = BLS24::_559();
    
    let en_48_581 = BLS48::_581();
    let en_48_573 = BLS48::_573();
    let en_48_575 = BLS48::_575();
    let en_48_571 = BLS48::_571();
    let en_48_287 = BLS48::_287();
    let en_48_277 = BLS48::_277();

    // Initializing generators
    let prop_g_381 = en_12_381.g1.default_generator();
    let prop_g_461 = en_12_461.g1.default_generator();
    let prop_g_446 = en_12_446.g1.default_generator();

    let prop_g_315 = en_24_315.g1.default_generator();
    let prop_g_479 = en_24_479.g1.default_generator();
    let prop_g_477 = en_24_477.g1.default_generator();
    let prop_g_509_snark = en_24_509_snark.g1.default_generator();
    let prop_g_509 = en_24_509.g1.default_generator();
    let prop_g_559 = en_24_559.g1.default_generator();
    
    let prop_g_581 = en_48_581.g1.default_generator();
    let prop_g_573 = en_48_573.g1.default_generator();
    let prop_g_571 = en_48_571.g1.default_generator();
    let prop_g_575 = en_48_575.g1.default_generator();
    let prop_g_287 = en_48_287.g1.default_generator();
    let prop_g_277 = en_48_277.g1.default_generator();

    // Initializing scalars

    let prop_scalars_zk_381: Vec<_> = (0..1000)
        .map(|i| Scalar::from_bytes(&BENCH_DATA[i]).unwrap())
        .collect();
    let prop_scalars_ark_381: Vec<_> = (0..1000)
        .map(|i| Fr::from_le_bytes_mod_order(&BENCH_DATA[i]))
        .collect();

    let prop_scalars_381: Vec<_> = (0..1000)
        .map(|i| en_12_381.fr.from_byte_array(&BENCH_DATA[i], Endianness::Big))
        .collect();
    let prop_scalars_461: Vec<_> = (0..1000)
        .map(|_| en_12_461.fr.random_element())
        .collect();
    let prop_scalars_446: Vec<_> = (0..1000)
        .map(|_| en_12_446.fr.random_element())
        .collect();
    
    let prop_scalars_315: Vec<_> = (0..1000)
        .map(|_| en_24_315.fr.random_element())
        .collect();
    let prop_scalars_479: Vec<_> = (0..1000)
        .map(|_| en_24_479.fr.random_element())
        .collect();
    let prop_scalars_477: Vec<_> = (0..1000)
        .map(|_| en_24_477.fr.random_element())
        .collect();
    let prop_scalars_509_snark: Vec<_> = (0..1000)
        .map(|_| en_24_509_snark.fr.random_element())
        .collect();
    let prop_scalars_509: Vec<_> = (0..1000)
        .map(|_| en_24_509.fr.random_element())
        .collect();
    let prop_scalars_559: Vec<_> = (0..1000)
        .map(|_| en_24_559.fr.random_element())
        .collect();
    
    let prop_scalars_581: Vec<_> = (0..1000)
        .map(|_| en_48_581.fr.random_element())
        .collect();
    let prop_scalars_573: Vec<_> = (0..1000)
        .map(|_| en_48_573.fr.random_element())
        .collect();
    let prop_scalars_575: Vec<_> = (0..1000)
        .map(|_| en_48_575.fr.random_element())
        .collect();
    let prop_scalars_571: Vec<_> = (0..1000)
        .map(|_| en_48_571.fr.random_element())
        .collect();
    let prop_scalars_277: Vec<_> = (0..1000)
        .map(|_| en_48_277.fr.random_element())
        .collect();
    let prop_scalars_287: Vec<_> = (0..1000)
        .map(|_| en_48_287.fr.random_element())
        .collect();

    let counter = AtomicUsize::new(0);

    let mut group = c.benchmark_group("glv_multiply_group");
    
    // BLS 12 
    group.sample_size(10000);
    group.noise_threshold(0.03); // More samples reduce variability    
    group.warm_up_time(std::time::Duration::from_secs(5)); 

    group.bench_function("glv_multiply_ZK_BLS_12_381_batch", |b| {
        b.iter_batched(
            || {
                // Prepare a batch of indices
                let start = counter.fetch_add(10, Ordering::Relaxed) % 1000; // Batch of 10
                (0..10).map(|j| &prop_scalars_zk_381[(start + j) % 1000]).collect::<Vec<_>>()
            },
            |batch| {
                for scalar in batch {
                    black_box( scalar * zk_g);
                }
            },
            criterion::BatchSize::SmallInput, // Adjust based on batch size
        );
    });

    group.bench_function("glv_multiply_ARK_BLS_12_381_batch", |b| {
        b.iter_batched(
            || {
                // Prepare a batch of indices
                let start = counter.fetch_add(10, Ordering::Relaxed) % 1000; // Batch of 10
                (0..10).map(|j| &prop_scalars_ark_381[(start + j) % 1000]).collect::<Vec<_>>()
            },
            |batch| {
                for scalar in batch {
                    let _ = black_box( ark_g.mul(scalar.into_repr()));
                }
            },
            criterion::BatchSize::SmallInput, // Adjust based on batch size
        );
    });
    
    group.bench_function("glv_multiply_BLS_12_381_batch", |b| {
        b.iter_batched(
            || {
                // Prepare a batch of indices
                let start = counter.fetch_add(10, Ordering::Relaxed) % 1000; // Batch of 10
                (0..10).map(|j| &prop_scalars_381[(start + j) % 1000]).collect::<Vec<_>>()
            },
            |batch| {
                for scalar in batch {
                    black_box(prop_g_381.glv_multiply(scalar));
                }
            },
            criterion::BatchSize::SmallInput, // Adjust based on batch size
        );
    });
    group.bench_function("glv_multiply_BLS_12_461_batch", |b| {
        b.iter_batched(
            || {
                // Prepare a batch of indices
                let start = counter.fetch_add(10, Ordering::Relaxed) % 1000; // Batch of 10
                (0..10).map(|j| &prop_scalars_461[(start + j) % 1000]).collect::<Vec<_>>()
            },
            |batch| {
                for scalar in batch {
                    black_box(prop_g_461.glv_multiply(scalar));
                }
            },
            criterion::BatchSize::SmallInput, // Adjust based on batch size
        );
    });
    group.bench_function("glv_multiply_BLS_12_446_batch", |b| {
        b.iter_batched(
            || {
                // Prepare a batch of indices
                let start = counter.fetch_add(10, Ordering::Relaxed) % 1000; // Batch of 10
                (0..10).map(|j| &prop_scalars_446[(start + j) % 1000]).collect::<Vec<_>>()
            },
            |batch| {
                for scalar in batch {
                    black_box(prop_g_446.glv_multiply(scalar));
                }
            },
            criterion::BatchSize::SmallInput, // Adjust based on batch size
        );
    });
    // // BLS 24 
    group.bench_function("glv_multiply_BLS_24_315_batch", |b| {
        b.iter_batched(
            || {
                // Prepare a batch of indices
                let start = counter.fetch_add(10, Ordering::Relaxed) % 1000; // Batch of 10
                (0..10).map(|j| &prop_scalars_315[(start + j) % 1000]).collect::<Vec<_>>()
            },
            |batch| {
                for scalar in batch {
                    black_box(prop_g_315.glv_multiply(scalar));
                }
            },
            criterion::BatchSize::SmallInput, // Adjust based on batch size
        );
    });
    group.bench_function("glv_multiply_BLS_24_479_batch", |b| {
        b.iter_batched(
            || {
                // Prepare a batch of indices
                let start = counter.fetch_add(10, Ordering::Relaxed) % 1000; // Batch of 10
                (0..10).map(|j| &prop_scalars_479[(start + j) % 1000]).collect::<Vec<_>>()
            },
            |batch| {
                for scalar in batch {
                    black_box(prop_g_479.glv_multiply(scalar));
                }
            },
            criterion::BatchSize::SmallInput, // Adjust based on batch size
        );
    });
    group.bench_function("glv_multiply_BLS_24_477_batch", |b| {
        b.iter_batched(
            || {
                // Prepare a batch of indices
                let start = counter.fetch_add(10, Ordering::Relaxed) % 1000; // Batch of 10
                (0..10).map(|j| &prop_scalars_477[(start + j) % 1000]).collect::<Vec<_>>()
            },
            |batch| {
                for scalar in batch {
                    black_box(prop_g_477.glv_multiply(scalar));
                }
            },
            criterion::BatchSize::SmallInput, // Adjust based on batch size
        );
    });    
    group.bench_function("glv_multiply_BLS_24_509_snark_batch", |b| {
        b.iter_batched(
            || {
                // Prepare a batch of indices
                let start = counter.fetch_add(10, Ordering::Relaxed) % 1000; // Batch of 10
                (0..10).map(|j| &prop_scalars_509_snark[(start + j) % 1000]).collect::<Vec<_>>()
            },
            |batch| {
                for scalar in batch {
                    black_box(prop_g_509_snark.glv_multiply(scalar));
                }
            },
            criterion::BatchSize::SmallInput, // Adjust based on batch size
        );
    });
    group.bench_function("glv_multiply_BLS_24_509_batch", |b| {
        b.iter_batched(
            || {
                // Prepare a batch of indices
                let start = counter.fetch_add(10, Ordering::Relaxed) % 1000; // Batch of 10
                (0..10).map(|j| &prop_scalars_509[(start + j) % 1000]).collect::<Vec<_>>()
            },
            |batch| {
                for scalar in batch {
                    black_box(prop_g_509.glv_multiply(scalar));
                }
            },
            criterion::BatchSize::SmallInput, // Adjust based on batch size
        );
    });
    group.bench_function("glv_multiply_BLS_24_559_batch", |b| {
        b.iter_batched(
            || {
                // Prepare a batch of indices
                let start = counter.fetch_add(10, Ordering::Relaxed) % 1000; // Batch of 10
                (0..10).map(|j| &prop_scalars_559[(start + j) % 1000]).collect::<Vec<_>>()
            },
            |batch| {
                for scalar in batch {
                    black_box(prop_g_559.glv_multiply(scalar));
                }
            },
            criterion::BatchSize::SmallInput, // Adjust based on batch size
        );
    });
    // // BLS 48 
    group.bench_function("glv_multiply_BLS_48_581_batch", |b| {
        b.iter_batched(
            || {
                // Prepare a batch of indices
                let start = counter.fetch_add(10, Ordering::Relaxed) % 1000; // Batch of 10
                (0..10).map(|j| &prop_scalars_581[(start + j) % 1000]).collect::<Vec<_>>()
            },
            |batch| {
                for scalar in batch {
                    black_box(prop_g_581.glv_multiply(scalar));
                }
            },
            criterion::BatchSize::SmallInput, // Adjust based on batch size
        );
    });
    group.bench_function("glv_multiply_BLS_48_575_batch", |b| {
        b.iter_batched(
            || {
                // Prepare a batch of indices
                let start = counter.fetch_add(10, Ordering::Relaxed) % 1000; // Batch of 10
                (0..10).map(|j| &prop_scalars_575[(start + j) % 1000]).collect::<Vec<_>>()
            },
            |batch| {
                for scalar in batch {
                    black_box(prop_g_575.glv_multiply(scalar));
                }
            },
            criterion::BatchSize::SmallInput, // Adjust based on batch size
        );
    });
    group.bench_function("glv_multiply_BLS_48_573_batch", |b| {
        b.iter_batched(
            || {
                // Prepare a batch of indices
                let start = counter.fetch_add(10, Ordering::Relaxed) % 1000; // Batch of 10
                (0..10).map(|j| &prop_scalars_573[(start + j) % 1000]).collect::<Vec<_>>()
            },
            |batch| {
                for scalar in batch {
                    black_box(prop_g_573.glv_multiply(scalar));
                }
            },
            criterion::BatchSize::SmallInput, // Adjust based on batch size
        );
    });
    group.bench_function("glv_multiply_BLS_48_571_batch", |b| {
        b.iter_batched(
            || {
                // Prepare a batch of indices
                let start = counter.fetch_add(10, Ordering::Relaxed) % 1000; // Batch of 10
                (0..10).map(|j| &prop_scalars_571[(start + j) % 1000]).collect::<Vec<_>>()
            },
            |batch| {
                for scalar in batch {
                    black_box(prop_g_571.glv_multiply(scalar));
                }
            },
            criterion::BatchSize::SmallInput, // Adjust based on batch size
        );
    });
    group.bench_function("glv_multiply_BLS_48_277_batch", |b| {
        b.iter_batched(
            || {
                // Prepare a batch of indices
                let start = counter.fetch_add(10, Ordering::Relaxed) % 1000; // Batch of 10
                (0..10).map(|j| &prop_scalars_277[(start + j) % 1000]).collect::<Vec<_>>()
            },
            |batch| {
                for scalar in batch {
                    black_box(prop_g_277.glv_multiply(scalar));
                }
            },
            criterion::BatchSize::SmallInput, // Adjust based on batch size
        );
    });
    group.bench_function("glv_multiply_BLS_48_287_batch", |b| {
        b.iter_batched(
            || {
                // Prepare a batch of indices
                let start = counter.fetch_add(10, Ordering::Relaxed) % 1000; // Batch of 10
                (0..10).map(|j| &prop_scalars_287[(start + j) % 1000]).collect::<Vec<_>>()
            },
            |batch| {
                for scalar in batch {
                    black_box(prop_g_287.glv_multiply(scalar));
                }
            },
            criterion::BatchSize::SmallInput, // Adjust based on batch size
        );
    });

    group.finish();
}

criterion_group!(benches, bench_glv_multiply);
criterion_main!(benches);
