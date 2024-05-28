// Code developed by FARAOUN Kamel Mohamed.
// EEDIS-Laboratory. UDL-University. Algeria
// During May 2024.

use num_bigint::{BigUint, ToBigInt};
use std::ops::BitOr;
use num_traits::Zero;
use sha2::{Digest, Sha256};
use crate::fields::{arithmetic, prime_fields::{FieldElement, PrimeField}};


//  I2OSP converts a nonnegative integer to an octet string of a specified length. RFC 3447, section 4.1 https://datatracker.ietf.org/doc/html/rfc3447#section-4.1
pub fn i2osp(x: usize, x_len: usize) -> Vec<u8> {
    let mut octets = vec![0u8; x_len];
    for i in (0..x_len).rev() {
        octets[i] = (x >> (8 * (x_len - 1 - i))) as u8;
    }
    octets
}

pub fn i2osp_pf<const N:usize>(x: &FieldElement<N>, x_len: usize) -> Vec<u8> {
    let mut one =[0u64;N];
    one[0] = 1;
    let from_mont = arithmetic::mul(&x.mont_limbs, &one,x.fieldparams);    
    let mut octets = Vec::<u8>::new();
    for num in from_mont {octets.extend(num.to_le_bytes())}
    octets.truncate(x_len);
    octets.reverse();
    octets
}


//  OS2IP converts an octet string to a nonnegative integer. RFC 3447, section 4.2 https://datatracker.ietf.org/doc/html/rfc3447#section-4.2
pub fn os2ip(octets: &[u8]) -> BigUint {
    let mut result = BigUint::zero();
    for &byte in octets {
        result = (&result << 8u8).bitor(& BigUint::from(byte as usize));
    }
    result
}





//  The expand_message_xmd function produces a pseudorandom byte string using a cryptographic hash function H that outputs b bits. 
//  https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-hash-to-curve-06#name-expand_message_xmd 
pub fn expand_message_xmd(message: &[u8], dst: &[u8], len: usize) -> Vec<u8> {
    let hash_length = 32; // Length of the hash function output in bytes (for example, SHA-256)

    // Step 1: Initialize H_0
    let mut h = Sha256::new();
    h.update(dst);
    h.update(&i2osp(0, 1)); // Append a single 0 byte to the destination string
    h.update(&i2osp(hash_length, 1)); // Append the length of the output in bytes
    h.update(&i2osp(len, 1)); // Append the desired output length in bytes
    h.update(&i2osp(0, 1)); // Append a single 0 byte to indicate the start of the message
    h.update(message);
    let h_0 = h.finalize();

    // Step 2: Expand H_0
    let mut output = vec![];
    for i in 0..(len / hash_length) + 1 {
        let mut h_i = Sha256::new();
        h_i.update(&h_0);
        h_i.update(&i2osp(i, 1));
        h_i.update(dst);
        h_i.update(&i2osp(1, 1)); // Append a single 1 byte to indicate a non-empty string
        output.extend_from_slice(&h_i.finalize());
    }
    output.truncate(len);
    output
}

//   The hash_to_field function hashes a bytes/string msg of any length into one or more elements of a field F.
//   https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-hash-to-curve-06#name-hashing-to-a-finite-field 

pub fn hash_bytes_to_field<const N:usize>(msg :&[u8],count :usize,field:&PrimeField<N>,dst:&[u8],sec_level:usize, ext_degree :usize) -> Vec<FieldElement<N>>
{
    let l = ((field.parametres.numlimbs * 64) + sec_level as usize) >> 3;    
    let byteslength = count * ext_degree * l;    
    let prngbytes = expand_message_xmd(msg, dst, byteslength);
    let mut result = Vec::<FieldElement<N>>::new();
    for i in 0..count*ext_degree { let elm_offset = l * i;
                                          let tv  = os2ip(&prngbytes[elm_offset..elm_offset+l]);
                                          result.push(field.from_bigint(&tv.to_bigint().unwrap()))
                                        }
    result
}

pub fn hash_string_to_field<const N:usize>(msg:&str,field:&PrimeField<N>,count :usize,sec_level:usize,ext_degree :usize)->Vec<FieldElement<N>>
{
    let bytes = msg.as_bytes();    
    hash_bytes_to_field(bytes, count, field, &[0], sec_level, ext_degree)
}

