// Code developed by FARAOUN Kamel Mohamed.
// EEDIS-Laboratory. UDL-University. Algeria
// During May 2024.

use std::fmt;
use num_bigint::ToBigInt;
use crate::{extensions::{ext_fields::{ExtElement, ExtField}, fp2::{Fp2Element, Fp2Field}, 
            fp4::{Fp4Element, Fp4Field}, 
            fp8::{Fp8Element, Fp8Field}}, 
            fields::prime_fields::{FieldElement,  PrimeField}, 
            tools::{arithmetic_interface::ArithmeticOperations, hashs::{hash_string_to_field, os2ip}}};


#[derive(Clone,Copy,Debug)]
pub enum ExtFieldG2Element <const N:usize, const PARAMSIZE:usize>
            {   Fp2(Fp2Element<N>),
                Fp4(Fp4Element<PARAMSIZE,N>),
                Fp8(Fp8Element<PARAMSIZE,N>)
            }

impl <const N:usize, const PARAMSIZE:usize> ArithmeticOperations for ExtFieldG2Element<N,PARAMSIZE>{
    fn addto(&self, other: &Self) -> Self {
        match (self, other) {   (ExtFieldG2Element::Fp2(x), ExtFieldG2Element::Fp2(y)) => x.addto(y).into(),
                                (ExtFieldG2Element::Fp4(x), ExtFieldG2Element::Fp4(y)) => x.addto(y).into(),
                                (ExtFieldG2Element::Fp8(x), ExtFieldG2Element::Fp8(y)) => x.addto(y).into(),
                                _ => unimplemented!("Addition not implemented for different types"),
                            }
    }
    fn double(&self) -> Self {
            match self {    ExtFieldG2Element::Fp2(x) => x.double().into(),
                            ExtFieldG2Element::Fp4(x)=> x.double().into(),
                            ExtFieldG2Element::Fp8(x)=> x.double().into(),
                        }       
    }
    fn substract(&self, other: &Self) -> Self {
        match (self, other) {   (ExtFieldG2Element::Fp2(x), ExtFieldG2Element::Fp2(y)) => x.substract(y).into(),
                                (ExtFieldG2Element::Fp4(x), ExtFieldG2Element::Fp4(y)) => x.substract(y).into(),
                                (ExtFieldG2Element::Fp8(x), ExtFieldG2Element::Fp8(y)) => x.substract(y).into(),
                                _ => unimplemented!("Substractionnot implemented for different types"),
                            }
    }
    fn multiply(&self, other: &Self) -> Self {
        match (self, other) {   (ExtFieldG2Element::Fp2(x), ExtFieldG2Element::Fp2(y)) => x.multiply(y).into(),
                                (ExtFieldG2Element::Fp4(x), ExtFieldG2Element::Fp4(y)) => x.multiply(y).into(),
                                (ExtFieldG2Element::Fp8(x), ExtFieldG2Element::Fp8(y)) => x.multiply(y).into(),
                                _ => unimplemented!("Addition not implemented for different types"),
                            }
    }
    fn sqr(&self) -> Self 
        {
            match self {    ExtFieldG2Element::Fp2(x) => x.sqr().into(),
                            ExtFieldG2Element::Fp4(x)=> x.sqr().into(),
                            ExtFieldG2Element::Fp8(x)=> x.sqr().into(),
                        }  
        }
    fn invert(&self) -> Self 
        {
            match self {    ExtFieldG2Element::Fp2(x) => x.invert().into(),
                            ExtFieldG2Element::Fp4(x)=> x.invert().into(),
                            ExtFieldG2Element::Fp8(x)=> x.invert().into(),
                    }  
        }
    fn negate(&self) -> Self 
        {
            match self {    ExtFieldG2Element::Fp2(x) => x.negate().into(),
                            ExtFieldG2Element::Fp4(x)=> x.negate().into(),
                            ExtFieldG2Element::Fp8(x)=> x.negate().into(),
                        }  
        }
    fn equal(&self, rhs :&Self) -> bool {
        match (self, rhs) {   (ExtFieldG2Element::Fp2(x), ExtFieldG2Element::Fp2(y)) => x.equal(y),
                                (ExtFieldG2Element::Fp4(x), ExtFieldG2Element::Fp4(y)) => x.equal(y),
                                (ExtFieldG2Element::Fp8(x), ExtFieldG2Element::Fp8(y)) => x.equal(y),
                                _ => unimplemented!("Addition not implemented for different types"),
                            }
    }
    fn is_zero(&self) -> bool {
        match self {    ExtFieldG2Element::Fp2(x) => x.is_zero(),
                        ExtFieldG2Element::Fp4(x)=> x.is_zero(),
                        ExtFieldG2Element::Fp8(x)=> x.is_zero(),
                    }       
    }
    fn is_one(&self) -> bool {
        match self {    ExtFieldG2Element::Fp2(x) => x.is_one(),
                        ExtFieldG2Element::Fp4(x)=> x.is_one(),
                        ExtFieldG2Element::Fp8(x)=> x.is_one(),
                    }       
    }

    fn to_dec_string(&self) -> String {
        match self {    ExtFieldG2Element::Fp2(x) => x.to_a_string(),
                        ExtFieldG2Element::Fp4(x)=> x.to_a_string(),
                        ExtFieldG2Element::Fp8(x)=> x.to_a_string(),
                    }       
    }
    fn to_hex_string(&self) -> String {
        match self {    ExtFieldG2Element::Fp2(x) => x.to_hex_string(),
                        ExtFieldG2Element::Fp4(x)=> x.to_hex_string(),
                        ExtFieldG2Element::Fp8(x)=> x.to_hex_string(),
                    }       
    }
    fn one(&self) -> Self {
        match self {    ExtFieldG2Element::Fp2(x) => x.one().into(),
                        ExtFieldG2Element::Fp4(x)=> x.one().into(),
                        ExtFieldG2Element::Fp8(x)=> x.one().into(),
                    }       
    }
    fn zero(&self) -> Self {
        match self {    ExtFieldG2Element::Fp2(x) => x.zero().into(),
                        ExtFieldG2Element::Fp4(x)=> x.zero().into(),
                        ExtFieldG2Element::Fp8(x)=> x.zero().into(),
                    }       
    }
}

impl  <const N:usize, const PARAMSIZE:usize> ExtFieldG2Element<N,PARAMSIZE>
{
    pub fn sqrt(&self) -> Option<Self> {
        match self {    ExtFieldG2Element::Fp2(x) => if x.sqrt().is_none() {None} else {Some(ExtFieldG2Element::Fp2(x.sqrt().unwrap()))},
                        ExtFieldG2Element::Fp4(x)=> if x.sqrt().is_none() {None} else {Some(ExtFieldG2Element::Fp4(x.sqrt().unwrap()))},
                        ExtFieldG2Element::Fp8(x)=> if x.sqrt().is_none() {None} else {Some(ExtFieldG2Element::Fp8(x.sqrt().unwrap()))},
                    }       
    }
    pub fn sign(&self) -> i8 {
        match self {    ExtFieldG2Element::Fp2(x) => x.sign(),
                        ExtFieldG2Element::Fp4(x)=> x.sign(),
                        ExtFieldG2Element::Fp8(x)=> x.sign(),
                    }       
    }
    pub fn conjugate(&self) -> Self {
        match self {    ExtFieldG2Element::Fp2(x) => x.conjugate().into(),
                        ExtFieldG2Element::Fp4(x)=> x.conjugate().into(),
                        ExtFieldG2Element::Fp8(x)=> x.conjugate().into(),
                    }       
    }
    pub fn to_i2osp_bytearray(&self) -> Vec<u8> 
    {
        match self {    ExtFieldG2Element::Fp2(x) => x.to_i2osp_bytearray(),
                        ExtFieldG2Element::Fp4(x)=> x.to_i2osp_bytearray(),
                        ExtFieldG2Element::Fp8(x)=> x.to_i2osp_bytearray(),
                    }       
    }
    pub fn mulby_fp_element(&self, b:&FieldElement<N>) -> Self 
    {
        match self {    ExtFieldG2Element::Fp2(x) => x.mulby_fp_element(&b).into(),
                        ExtFieldG2Element::Fp4(x)=> x.mulby_fp_element(&b).into(),
                        ExtFieldG2Element::Fp8(x)=> x.mulby_fp_element(&b).into(),
                    }       
    }
    pub fn content(&self) -> &[FieldElement<N>]
    {
        match self {    ExtFieldG2Element::Fp2(x) => &x.content,
                        ExtFieldG2Element::Fp4(x)=> &x.content,
                        ExtFieldG2Element::Fp8(x)=> &x.content,
                    }       
    }
    

}

impl<const N: usize, const PARAMSIZE: usize> From<Fp2Element<N>> for ExtFieldG2Element<N, PARAMSIZE> {
    fn from(fp2: Fp2Element<N>) -> Self {
        ExtFieldG2Element::Fp2(fp2)
    }
}

impl<const N: usize, const PARAMSIZE: usize> From<Fp4Element<PARAMSIZE, N>> for ExtFieldG2Element<N, PARAMSIZE> {
    fn from(fp4: Fp4Element<PARAMSIZE, N>) -> Self {
        ExtFieldG2Element::Fp4(fp4)
    }
}

impl<const N: usize, const PARAMSIZE: usize> From<Fp8Element<PARAMSIZE, N>> for ExtFieldG2Element<N, PARAMSIZE> {
    fn from(fp8: Fp8Element<PARAMSIZE, N>) -> Self {
        ExtFieldG2Element::Fp8(fp8)
    }
}


#[derive(Debug)]
pub enum ExtG2Field <const N:usize, const PARAMSIZE:usize>
            {   Fp2(&'static Fp2Field<N>),
                Fp4(&'static Fp4Field<PARAMSIZE,N>),
                Fp8(&'static Fp8Field<PARAMSIZE,N>)
            }

impl <const PARAMSIZE:usize,const N: usize> ExtG2Field<N,PARAMSIZE> 
{
    pub fn new(&self) -> ExtFieldG2Element<N,PARAMSIZE>
    {   
        match self {    ExtG2Field::Fp2(x) => x.random_element().into(),
                        ExtG2Field::Fp4(x)=> x.random_element().into(),
                        ExtG2Field::Fp8(x)=> x.random_element().into(),
                    }       
    }

    pub fn random_element(&self) -> ExtFieldG2Element<N,PARAMSIZE>
    {   
        match self {    ExtG2Field::Fp2(x) => x.random_element().into(),
                        ExtG2Field::Fp4(x)=> x.random_element().into(),
                        ExtG2Field::Fp8(x)=> x.random_element().into(),
                    }       
    }

    pub fn one(&self) -> ExtFieldG2Element<N,PARAMSIZE>
    {   
        match self {    ExtG2Field::Fp2(x) => x.one().clone().into(),
                        ExtG2Field::Fp4(x)=> x.one().clone().into(),
                        ExtG2Field::Fp8(x)=> x.one().clone().into(),
                    }       
    }

    pub fn zero(&self) -> ExtFieldG2Element<N,PARAMSIZE>
    {   
        match self {    ExtG2Field::Fp2(x) => x.zero().clone().into(),
                        ExtG2Field::Fp4(x)=> x.zero().clone().into(),
                        ExtG2Field::Fp8(x)=> x.zero().clone().into(),
                    }       
    }

    pub fn from_hex_strings(&self, source :&[&str]) -> ExtFieldG2Element<N,PARAMSIZE>
    {   
        match self {    ExtG2Field::Fp2(x) => x.from_hex_strings(source).into(),
                        ExtG2Field::Fp4(x)=> x.from_hex_strings(source).into(),
                        ExtG2Field::Fp8(x)=> x.from_hex_strings(source).into(),
                    }       
    }

    pub fn from_decimal_strings(&self, source :&[&str]) -> ExtFieldG2Element<N,PARAMSIZE>
    {   
        match self {    ExtG2Field::Fp2(x) => x.from_strings(source).into(),
                        ExtG2Field::Fp4(x)=> x.from_strings(source).into(),
                        ExtG2Field::Fp8(x)=> x.from_strings(source).into(),
                    }       
    }

    pub fn from_basefield_elements(&self, source :&[FieldElement<N>]) -> ExtFieldG2Element<N,PARAMSIZE>
    {   
        match self {    ExtG2Field::Fp2(x) => x.from_field_elements(source).into(),
                        ExtG2Field::Fp4(x)=> x.from_field_elements(source).into(),
                        ExtG2Field::Fp8(x)=> x.from_field_elements(source).into(),
                    }       
    }

    pub fn from_i2osp_bytearray(&self, source :&[u8]) -> ExtFieldG2Element<N,PARAMSIZE>
    {   let field =self.basefield();
        
        match self {    ExtG2Field::Fp2(x) => 
                                {   let sizeinbytes = source.len()/2;
                                    x.from_field_elements(&[field.from_bigint(&os2ip(&source[0..sizeinbytes]).to_bigint().unwrap()),
                                                                    field.from_bigint(&os2ip(&source[sizeinbytes..]).to_bigint().unwrap()) ]).into()},
                        ExtG2Field::Fp4(x)=>  
                                {   let sizeinbytes = source.len()/4;
                                    x.from_field_elements (&[  field.from_bigint(&os2ip(&source[0..sizeinbytes]).to_bigint().unwrap()),
                                                                       field.from_bigint(&os2ip(&source[sizeinbytes..2*sizeinbytes]).to_bigint().unwrap()),
                                                                       field.from_bigint(&os2ip(&source[2*sizeinbytes..3*sizeinbytes]).to_bigint().unwrap()),
                                                                       field.from_bigint(&os2ip(&source[3*sizeinbytes..]).to_bigint().unwrap()) ]).into()},
                        ExtG2Field::Fp8(x)=> 
                            {   let sizeinbytes =source.len()/8;
                                x.from_field_elements (&[  field.from_bigint(&os2ip(&source[0..sizeinbytes]).to_bigint().unwrap()),
                                                                  field.from_bigint(&os2ip(&source[sizeinbytes..2*sizeinbytes]).to_bigint().unwrap()),
                                                                  field.from_bigint(&os2ip(&source[2*sizeinbytes..3*sizeinbytes]).to_bigint().unwrap()),
                                                                  field.from_bigint(&os2ip(&source[3*sizeinbytes..4*sizeinbytes]).to_bigint().unwrap()),
                                                                  field.from_bigint(&os2ip(&source[4*sizeinbytes..5*sizeinbytes]).to_bigint().unwrap()),
                                                                  field.from_bigint(&os2ip(&source[5*sizeinbytes..6*sizeinbytes]).to_bigint().unwrap()),
                                                                  field.from_bigint(&os2ip(&source[6*sizeinbytes..7*sizeinbytes]).to_bigint().unwrap()),
                                                                  field.from_bigint(&os2ip(&source[7*sizeinbytes..8*sizeinbytes]).to_bigint().unwrap()) ]).into()},
                    }       
    }
    pub fn basefield(&self) ->PrimeField<N>
    {
        match self {    ExtG2Field::Fp2(x) => x.field_interface(),
                        ExtG2Field::Fp4(x)=> x.field_interface(),
                        ExtG2Field::Fp8(x)=> x.field_interface(),
        }       
    }
    pub fn hash_to_field(&self,id : &str, security_level:usize,count :usize) -> Vec<ExtFieldG2Element<N,PARAMSIZE>>
    {   
        let extorder;
        match self {    ExtG2Field::Fp2(_) => extorder = 2,
                        ExtG2Field::Fp4(_)=> extorder = 4,
                        ExtG2Field::Fp8(_)=> extorder = 8,
                   }  
        let hashvec = hash_string_to_field(id, &self.basefield(), count, security_level, extorder);
        let mut result = Vec::<ExtFieldG2Element<N,PARAMSIZE>>::new();
        let mut i: usize =0;
        while i< count { result.push(self.from_basefield_elements(&hashvec[i*extorder..(i+1)*extorder]));
                         i = i + 1;
                         }
        result
    }
}

impl<'a, const N: usize, const PARAMSIZE: usize> fmt::Display for ExtFieldG2Element<N,PARAMSIZE> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:}", &self.to_dec_string())
    }
    }