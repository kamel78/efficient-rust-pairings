// Code developed by FARAOUN Kamel Mohamed.
// EEDIS-Laboratory. UDL-University. Algeria
// During May 2024.

use std::fmt;
use num_bigint::ToBigInt;
use crate::{extensions::{ext_fields::{ExtElement, ExtField}, 
            towering1::fp2::{Fp2Element as Fp2Element_1,  Fp2Field as Fp2Field_1}, 
            towering2::fp2::{Fp2Element as Fp2Element_2,  Fp2Field as Fp2Field_2}, 
            towering1::fp4::{Fp4Element as Fp4Element_1, Fp4Field as Fp4Field_1}, 
            towering2::fp4::{Fp4Element as Fp4Element_2, Fp4Field as Fp4Field_2},             
            towering3::fp4::{Fp4Element as Fp4Element_3, Fp4Field as Fp4Field_3},             
            towering1::fp8::{Fp8Element as Fp8Element_1, Fp8Field as Fp8Field_1}, 
            towering2::fp8::{Fp8Element as Fp8Element_2, Fp8Field as Fp8Field_2}, 
            towering3::fp8::{Fp8Element as Fp8Element_3, Fp8Field as Fp8Field_3}}, 
            fields::prime_fields::{FieldElement,  PrimeField}, 
            tools::{arithmetic_interface::ArithmeticOperations, hashs::{hash_string_to_field, os2ip}}};


#[derive(Clone,Copy,Debug)]
pub enum ExtFieldG2Element <const N:usize, const PARAMSIZE:usize>
            {   Fp2_1(Fp2Element_1<N>),
                Fp4_1(Fp4Element_1<PARAMSIZE,N>),
                Fp8_1(Fp8Element_1<PARAMSIZE,N>),
                Fp2_2(Fp2Element_2<PARAMSIZE,N>),
                Fp4_2(Fp4Element_2<PARAMSIZE,N>),
                Fp8_2(Fp8Element_2<PARAMSIZE,N>),                
                Fp4_3(Fp4Element_3<PARAMSIZE,N>),
                Fp8_3(Fp8Element_3<PARAMSIZE,N>),                
            }

impl <const N:usize, const PARAMSIZE:usize> ArithmeticOperations for ExtFieldG2Element<N,PARAMSIZE>{
    fn addto(&self, other: &Self) -> Self {
        match (self, other) {   (ExtFieldG2Element::Fp2_1(x), ExtFieldG2Element::Fp2_1(y)) => x.addto(y).into(),
                                (ExtFieldG2Element::Fp4_1(x), ExtFieldG2Element::Fp4_1(y)) => x.addto(y).into(),
                                (ExtFieldG2Element::Fp8_1(x), ExtFieldG2Element::Fp8_1(y)) => x.addto(y).into(),
                                (ExtFieldG2Element::Fp2_2(x), ExtFieldG2Element::Fp2_2(y)) => x.addto(y).into(),
                                (ExtFieldG2Element::Fp4_2(x), ExtFieldG2Element::Fp4_2(y)) => x.addto(y).into(),
                                (ExtFieldG2Element::Fp8_2(x), ExtFieldG2Element::Fp8_2(y)) => x.addto(y).into(),
                                (ExtFieldG2Element::Fp4_3(x), ExtFieldG2Element::Fp4_3(y)) => x.addto(y).into(),
                                (ExtFieldG2Element::Fp8_3(x), ExtFieldG2Element::Fp8_3(y)) => x.addto(y).into(),
                                _ => unimplemented!("Addition not implemented for different types"),
                            }
    }
    fn double(&self) -> Self {
            match self {    ExtFieldG2Element::Fp2_1(x) => x.double().into(),
                            ExtFieldG2Element::Fp4_1(x)=> x.double().into(),
                            ExtFieldG2Element::Fp8_1(x)=> x.double().into(),
                            ExtFieldG2Element::Fp2_2(x) => x.double().into(),
                            ExtFieldG2Element::Fp4_2(x)=> x.double().into(),
                            ExtFieldG2Element::Fp8_2(x)=> x.double().into(),
                            ExtFieldG2Element::Fp4_3(x)=> x.double().into(),
                            ExtFieldG2Element::Fp8_3(x)=> x.double().into(),
                        }       
    }
    fn substract(&self, other: &Self) -> Self {
        match (self, other) {   (ExtFieldG2Element::Fp2_1(x), ExtFieldG2Element::Fp2_1(y)) => x.substract(y).into(),
                                (ExtFieldG2Element::Fp4_1(x), ExtFieldG2Element::Fp4_1(y)) => x.substract(y).into(),
                                (ExtFieldG2Element::Fp8_1(x), ExtFieldG2Element::Fp8_1(y)) => x.substract(y).into(),
                                (ExtFieldG2Element::Fp2_2(x), ExtFieldG2Element::Fp2_2(y)) => x.substract(y).into(),
                                (ExtFieldG2Element::Fp4_2(x), ExtFieldG2Element::Fp4_2(y)) => x.substract(y).into(),
                                (ExtFieldG2Element::Fp8_2(x), ExtFieldG2Element::Fp8_2(y)) => x.substract(y).into(),
                                (ExtFieldG2Element::Fp4_3(x), ExtFieldG2Element::Fp4_3(y)) => x.substract(y).into(),
                                (ExtFieldG2Element::Fp8_3(x), ExtFieldG2Element::Fp8_3(y)) => x.substract(y).into(),
                                _ => unimplemented!("Substractionnot implemented for different types"),
                            }
    }
    fn multiply(&self, other: &Self) -> Self {
        match (self, other) {   (ExtFieldG2Element::Fp2_1(x), ExtFieldG2Element::Fp2_1(y)) => x.multiply(y).into(),
                                (ExtFieldG2Element::Fp4_1(x), ExtFieldG2Element::Fp4_1(y)) => x.multiply(y).into(),
                                (ExtFieldG2Element::Fp8_1(x), ExtFieldG2Element::Fp8_1(y)) => x.multiply(y).into(),
                                (ExtFieldG2Element::Fp2_2(x), ExtFieldG2Element::Fp2_2(y)) => x.multiply(y).into(),
                                (ExtFieldG2Element::Fp4_2(x), ExtFieldG2Element::Fp4_2(y)) => x.multiply(y).into(),
                                (ExtFieldG2Element::Fp8_2(x), ExtFieldG2Element::Fp8_2(y)) => x.multiply(y).into(),
                                (ExtFieldG2Element::Fp4_3(x), ExtFieldG2Element::Fp4_3(y)) => x.multiply(y).into(),
                                (ExtFieldG2Element::Fp8_3(x), ExtFieldG2Element::Fp8_3(y)) => x.multiply(y).into(),
                                _ => unimplemented!("Addition not implemented for different types"),
                            }
    }
    fn sqr(&self) -> Self 
        {
            match self {    ExtFieldG2Element::Fp2_1(x) => x.sqr().into(),
                            ExtFieldG2Element::Fp4_1(x)=> x.sqr().into(),
                            ExtFieldG2Element::Fp8_1(x)=> x.sqr().into(),
                            ExtFieldG2Element::Fp2_2(x) => x.sqr().into(),
                            ExtFieldG2Element::Fp4_2(x)=> x.sqr().into(),
                            ExtFieldG2Element::Fp8_2(x)=> x.sqr().into(),
                            ExtFieldG2Element::Fp4_3(x)=> x.sqr().into(),
                            ExtFieldG2Element::Fp8_3(x)=> x.sqr().into(),
                        }  
        }
    fn invert(&self) -> Self 
        {
            match self {    ExtFieldG2Element::Fp2_1(x) => x.invert().into(),
                            ExtFieldG2Element::Fp4_1(x)=> x.invert().into(),
                            ExtFieldG2Element::Fp8_1(x)=> x.invert().into(),
                            ExtFieldG2Element::Fp2_2(x) => x.invert().into(),
                            ExtFieldG2Element::Fp4_2(x)=> x.invert().into(),
                            ExtFieldG2Element::Fp8_2(x)=> x.invert().into(),
                            ExtFieldG2Element::Fp4_3(x)=> x.invert().into(),
                            ExtFieldG2Element::Fp8_3(x)=> x.invert().into(),
                    }  
        }
    fn negate(&self) -> Self 
        {
            match self {    ExtFieldG2Element::Fp2_1(x) => x.negate().into(),
                            ExtFieldG2Element::Fp4_1(x)=> x.negate().into(),
                            ExtFieldG2Element::Fp8_1(x)=> x.negate().into(),
                            ExtFieldG2Element::Fp2_2(x) => x.negate().into(),
                            ExtFieldG2Element::Fp4_2(x)=> x.negate().into(),
                            ExtFieldG2Element::Fp8_2(x)=> x.negate().into(),
                            ExtFieldG2Element::Fp4_3(x)=> x.negate().into(),
                            ExtFieldG2Element::Fp8_3(x)=> x.negate().into(),
                        }  
        }
    fn equal(&self, rhs :&Self) -> bool {
        match (self, rhs) {   (ExtFieldG2Element::Fp2_1(x), ExtFieldG2Element::Fp2_1(y)) => x.equal(y),
                                (ExtFieldG2Element::Fp4_1(x), ExtFieldG2Element::Fp4_1(y)) => x.equal(y),
                                (ExtFieldG2Element::Fp8_1(x), ExtFieldG2Element::Fp8_1(y)) => x.equal(y),
                                (ExtFieldG2Element::Fp2_2(x), ExtFieldG2Element::Fp2_2(y)) => x.equal(y),
                                (ExtFieldG2Element::Fp4_2(x), ExtFieldG2Element::Fp4_2(y)) => x.equal(y),
                                (ExtFieldG2Element::Fp8_2(x), ExtFieldG2Element::Fp8_2(y)) => x.equal(y),
                                (ExtFieldG2Element::Fp4_3(x), ExtFieldG2Element::Fp4_3(y)) => x.equal(y),
                                (ExtFieldG2Element::Fp8_3(x), ExtFieldG2Element::Fp8_3(y)) => x.equal(y),
                                _ => unimplemented!("Addition not implemented for different types"),
                            }
    }
    fn is_zero(&self) -> bool {
        match self {    ExtFieldG2Element::Fp2_1(x) => x.is_zero(),
                        ExtFieldG2Element::Fp4_1(x)=> x.is_zero(),
                        ExtFieldG2Element::Fp8_1(x)=> x.is_zero(),
                        ExtFieldG2Element::Fp2_2(x) => x.is_zero(),
                        ExtFieldG2Element::Fp4_2(x)=> x.is_zero(),
                        ExtFieldG2Element::Fp8_2(x)=> x.is_zero(),
                        ExtFieldG2Element::Fp4_3(x)=> x.is_zero(),
                        ExtFieldG2Element::Fp8_3(x)=> x.is_zero(),
                    }       
    }
    fn is_one(&self) -> bool {
        match self {    ExtFieldG2Element::Fp2_1(x) => x.is_one(),
                        ExtFieldG2Element::Fp4_1(x)=> x.is_one(),
                        ExtFieldG2Element::Fp8_1(x)=> x.is_one(),
                        ExtFieldG2Element::Fp2_2(x) => x.is_one(),
                        ExtFieldG2Element::Fp4_2(x)=> x.is_one(),
                        ExtFieldG2Element::Fp8_2(x)=> x.is_one(),
                        ExtFieldG2Element::Fp4_3(x)=> x.is_one(),
                        ExtFieldG2Element::Fp8_3(x)=> x.is_one(),
                    }       
    }

    fn to_dec_string(&self) -> String {
        match self {    ExtFieldG2Element::Fp2_1(x) => x.to_a_string(),
                        ExtFieldG2Element::Fp4_1(x)=> x.to_a_string(),
                        ExtFieldG2Element::Fp8_1(x)=> x.to_a_string(),
                        ExtFieldG2Element::Fp2_2(x) => x.to_a_string(),
                        ExtFieldG2Element::Fp4_2(x)=> x.to_a_string(),
                        ExtFieldG2Element::Fp8_2(x)=> x.to_a_string(),
                        ExtFieldG2Element::Fp4_3(x)=> x.to_a_string(),
                        ExtFieldG2Element::Fp8_3(x)=> x.to_a_string(),
                    }       
    }
    fn to_hex_string(&self) -> String {
        match self {    ExtFieldG2Element::Fp2_1(x) => x.to_hex_string(),
                        ExtFieldG2Element::Fp4_1(x)=> x.to_hex_string(),
                        ExtFieldG2Element::Fp8_1(x)=> x.to_hex_string(),
                        ExtFieldG2Element::Fp2_2(x) => x.to_hex_string(),
                        ExtFieldG2Element::Fp4_2(x)=> x.to_hex_string(),
                        ExtFieldG2Element::Fp8_2(x)=> x.to_hex_string(),
                        ExtFieldG2Element::Fp4_3(x)=> x.to_hex_string(),
                        ExtFieldG2Element::Fp8_3(x)=> x.to_hex_string(),
                    }       
    }
    fn one(&self) -> Self {
        match self {    ExtFieldG2Element::Fp2_1(x) => x.one().into(),
                        ExtFieldG2Element::Fp4_1(x)=> x.one().into(),
                        ExtFieldG2Element::Fp8_1(x)=> x.one().into(),
                        ExtFieldG2Element::Fp2_2(x) => x.one().into(),
                        ExtFieldG2Element::Fp4_2(x)=> x.one().into(),
                        ExtFieldG2Element::Fp8_2(x)=> x.one().into(),
                        ExtFieldG2Element::Fp4_3(x)=> x.one().into(),
                        ExtFieldG2Element::Fp8_3(x)=> x.one().into(),
                    }       
    }
    fn zero(&self) -> Self {
        match self {    ExtFieldG2Element::Fp2_1(x) => x.zero().into(),
                        ExtFieldG2Element::Fp4_1(x)=> x.zero().into(),
                        ExtFieldG2Element::Fp8_1(x)=> x.zero().into(),
                        ExtFieldG2Element::Fp2_2(x) => x.zero().into(),
                        ExtFieldG2Element::Fp4_2(x)=> x.zero().into(),
                        ExtFieldG2Element::Fp8_2(x)=> x.zero().into(),
                        ExtFieldG2Element::Fp4_3(x)=> x.zero().into(),
                        ExtFieldG2Element::Fp8_3(x)=> x.zero().into(),
                    }       
    }
}

impl  <const N:usize, const PARAMSIZE:usize> ExtFieldG2Element<N,PARAMSIZE>
{
    pub fn sqrt(&self) -> Option<Self> {
        match self {    ExtFieldG2Element::Fp2_1(x) => if x.sqrt().is_none() {None} else {Some(ExtFieldG2Element::Fp2_1(x.sqrt().unwrap()))},
                        ExtFieldG2Element::Fp4_1(x)=> if x.sqrt().is_none() {None} else {Some(ExtFieldG2Element::Fp4_1(x.sqrt().unwrap()))},
                        ExtFieldG2Element::Fp8_1(x)=> if x.sqrt().is_none() {None} else {Some(ExtFieldG2Element::Fp8_1(x.sqrt().unwrap()))},
                        ExtFieldG2Element::Fp2_2(x) => if x.sqrt().is_none() {None} else {Some(ExtFieldG2Element::Fp2_2(x.sqrt().unwrap()))},
                        ExtFieldG2Element::Fp4_2(x)=> if x.sqrt().is_none() {None} else {Some(ExtFieldG2Element::Fp4_2(x.sqrt().unwrap()))},
                        ExtFieldG2Element::Fp8_2(x)=> if x.sqrt().is_none() {None} else {Some(ExtFieldG2Element::Fp8_2(x.sqrt().unwrap()))},
                        ExtFieldG2Element::Fp4_3(x)=> if x.sqrt().is_none() {None} else {Some(ExtFieldG2Element::Fp4_3(x.sqrt().unwrap()))},
                        ExtFieldG2Element::Fp8_3(x)=> if x.sqrt().is_none() {None} else {Some(ExtFieldG2Element::Fp8_3(x.sqrt().unwrap()))},
                    }       
    }
    pub fn is_qr(&self) -> bool {
        match self {    ExtFieldG2Element::Fp2_1(x) => {x.is_qr()},
                        ExtFieldG2Element::Fp4_1(x)=> {x.is_qr()},
                        ExtFieldG2Element::Fp8_1(x)=> {x.is_qr()},
                        ExtFieldG2Element::Fp2_2(x) => {x.is_qr()},
                        ExtFieldG2Element::Fp4_2(x)=> {x.is_qr()},
                        ExtFieldG2Element::Fp8_2(x)=> {x.is_qr()},
                        ExtFieldG2Element::Fp4_3(x)=> {x.is_qr()},
                        ExtFieldG2Element::Fp8_3(x)=> {x.is_qr()},
                    }       
    }
    pub fn sign(&self) -> i8 {
        match self {    ExtFieldG2Element::Fp2_1(x) => x.sign(),
                        ExtFieldG2Element::Fp4_1(x)=> x.sign(),
                        ExtFieldG2Element::Fp8_1(x)=> x.sign(),
                        ExtFieldG2Element::Fp2_2(x) => x.sign(),
                        ExtFieldG2Element::Fp4_2(x)=> x.sign(),
                        ExtFieldG2Element::Fp8_2(x)=> x.sign(),
                        ExtFieldG2Element::Fp4_3(x)=> x.sign(),
                        ExtFieldG2Element::Fp8_3(x)=> x.sign(),
                    }       
    }
    pub fn conjugate(&self) -> Self {
        match self {    ExtFieldG2Element::Fp2_1(x) => x.conjugate().into(),
                        ExtFieldG2Element::Fp4_1(x)=> x.conjugate().into(),
                        ExtFieldG2Element::Fp8_1(x)=> x.conjugate().into(),
                        ExtFieldG2Element::Fp2_2(x) => x.conjugate().into(),
                        ExtFieldG2Element::Fp4_2(x)=> x.conjugate().into(),
                        ExtFieldG2Element::Fp8_2(x)=> x.conjugate().into(),
                        ExtFieldG2Element::Fp4_3(x)=> x.conjugate().into(),
                        ExtFieldG2Element::Fp8_3(x)=> x.conjugate().into(),
                    }       
    }
    pub fn frobinus(&self) -> Self {
        match self {    ExtFieldG2Element::Fp2_1(x) => x.frobinus().into(),
                        ExtFieldG2Element::Fp4_1(x)=> x.frobinus().into(),
                        ExtFieldG2Element::Fp8_1(x)=> x.frobinus().into(),
                        ExtFieldG2Element::Fp2_2(x) => x.frobinus().into(),
                        ExtFieldG2Element::Fp4_2(x)=> x.frobinus().into(),
                        ExtFieldG2Element::Fp8_2(x)=> x.frobinus().into(),
                        ExtFieldG2Element::Fp4_3(x)=> x.frobinus().into(),
                        ExtFieldG2Element::Fp8_3(x)=> x.frobinus().into(),
                    }       
    }
    pub fn to_i2osp_bytearray(&self) -> Vec<u8> 
    {
        match self {    ExtFieldG2Element::Fp2_1(x) => x.to_i2osp_bytearray(),
                        ExtFieldG2Element::Fp4_1(x)=> x.to_i2osp_bytearray(),
                        ExtFieldG2Element::Fp8_1(x)=> x.to_i2osp_bytearray(),
                        ExtFieldG2Element::Fp2_2(x) => x.to_i2osp_bytearray(),
                        ExtFieldG2Element::Fp4_2(x)=> x.to_i2osp_bytearray(),
                        ExtFieldG2Element::Fp8_2(x)=> x.to_i2osp_bytearray(),
                        ExtFieldG2Element::Fp4_3(x)=> x.to_i2osp_bytearray(),
                        ExtFieldG2Element::Fp8_3(x)=> x.to_i2osp_bytearray(),
                    }       
    }
    pub fn mulby_fp_element(&self, b:&FieldElement<N>) -> Self 
    {
        match self {    ExtFieldG2Element::Fp2_1(x) => x.mulby_fp_element(&b).into(),
                        ExtFieldG2Element::Fp4_1(x)=> x.mulby_fp_element(&b).into(),
                        ExtFieldG2Element::Fp8_1(x)=> x.mulby_fp_element(&b).into(),
                        ExtFieldG2Element::Fp2_2(x) => x.mulby_fp_element(&b).into(),
                        ExtFieldG2Element::Fp4_2(x)=> x.mulby_fp_element(&b).into(),
                        ExtFieldG2Element::Fp8_2(x)=> x.mulby_fp_element(&b).into(),
                        ExtFieldG2Element::Fp4_3(x)=> x.mulby_fp_element(&b).into(),
                        ExtFieldG2Element::Fp8_3(x)=> x.mulby_fp_element(&b).into(),
                    }       
    }
    pub fn content(&self) -> &[FieldElement<N>]
    {
        match self {    ExtFieldG2Element::Fp2_1(x) => &x.content,
                        ExtFieldG2Element::Fp4_1(x)=> &x.content,
                        ExtFieldG2Element::Fp8_1(x)=> &x.content,
                        ExtFieldG2Element::Fp2_2(x) => &x.content,
                        ExtFieldG2Element::Fp4_2(x)=> &x.content,
                        ExtFieldG2Element::Fp8_2(x)=> &x.content,
                        ExtFieldG2Element::Fp4_3(x)=> &x.content,
                        ExtFieldG2Element::Fp8_3(x)=> &x.content,
                    }       
    }
    

}

impl<const N: usize, const PARAMSIZE: usize> From<Fp2Element_1<N>> for ExtFieldG2Element<N, PARAMSIZE> {
    fn from(fp2: Fp2Element_1<N>) -> Self {
        ExtFieldG2Element::Fp2_1(fp2)
    }
}

impl<const N: usize, const PARAMSIZE: usize> From<Fp2Element_2<PARAMSIZE,N>> for ExtFieldG2Element<N, PARAMSIZE> {
    fn from(fp2: Fp2Element_2<PARAMSIZE,N>) -> Self {
        ExtFieldG2Element::Fp2_2(fp2)
    }
}

impl<const N: usize, const PARAMSIZE: usize> From<Fp4Element_1<PARAMSIZE, N>> for ExtFieldG2Element<N, PARAMSIZE> {
    fn from(fp4: Fp4Element_1<PARAMSIZE, N>) -> Self {
        ExtFieldG2Element::Fp4_1(fp4)
    }
}

impl<const N: usize, const PARAMSIZE: usize> From<Fp4Element_2<PARAMSIZE, N>> for ExtFieldG2Element<N, PARAMSIZE> {
    fn from(fp4: Fp4Element_2<PARAMSIZE, N>) -> Self {
        ExtFieldG2Element::Fp4_2(fp4)
    }
}

impl<const N: usize, const PARAMSIZE: usize> From<Fp4Element_3<PARAMSIZE, N>> for ExtFieldG2Element<N, PARAMSIZE> {
    fn from(fp4: Fp4Element_3<PARAMSIZE, N>) -> Self {
        ExtFieldG2Element::Fp4_3(fp4)
    }
}

impl<const N: usize, const PARAMSIZE: usize> From<Fp8Element_1<PARAMSIZE, N>> for ExtFieldG2Element<N, PARAMSIZE> {
    fn from(fp8: Fp8Element_1<PARAMSIZE, N>) -> Self {
        ExtFieldG2Element::Fp8_1(fp8)
    }
}

impl<const N: usize, const PARAMSIZE: usize> From<Fp8Element_2<PARAMSIZE, N>> for ExtFieldG2Element<N, PARAMSIZE> {
    fn from(fp8: Fp8Element_2<PARAMSIZE, N>) -> Self {
        ExtFieldG2Element::Fp8_2(fp8)
    }
}

impl<const N: usize, const PARAMSIZE: usize> From<Fp8Element_3<PARAMSIZE, N>> for ExtFieldG2Element<N, PARAMSIZE> {
    fn from(fp8: Fp8Element_3<PARAMSIZE, N>) -> Self {
        ExtFieldG2Element::Fp8_3(fp8)
    }
}

#[derive(Debug)]
pub enum ExtG2Field <const N:usize, const PARAMSIZE:usize>
            {   Fp2_1(&'static Fp2Field_1<N>),
                Fp4_1(&'static Fp4Field_1<PARAMSIZE,N>),
                Fp8_1(&'static Fp8Field_1<PARAMSIZE,N>),
                Fp2_2(&'static Fp2Field_2<N>),
                Fp4_2(&'static Fp4Field_2<PARAMSIZE,N>),
                Fp8_2(&'static Fp8Field_2<PARAMSIZE,N>),
                Fp4_3(&'static Fp4Field_3<PARAMSIZE,N>),
                Fp8_3(&'static Fp8Field_3<PARAMSIZE,N>)
            }

impl <const PARAMSIZE:usize,const N: usize> ExtG2Field<N,PARAMSIZE> 
{
    pub fn new(&self) -> ExtFieldG2Element<N,PARAMSIZE>
    {   
        match self {    ExtG2Field::Fp2_1(x) => x.random_element().into(),
                        ExtG2Field::Fp4_1(x)=> x.random_element().into(),
                        ExtG2Field::Fp8_1(x)=> x.random_element().into(),
                        ExtG2Field::Fp2_2(x) => x.random_element().into(),
                        ExtG2Field::Fp4_2(x)=> x.random_element().into(),
                        ExtG2Field::Fp8_2(x)=> x.random_element().into(),
                        ExtG2Field::Fp4_3(x)=> x.random_element().into(),
                        ExtG2Field::Fp8_3(x)=> x.random_element().into(),
                    }       
    }

    pub fn random_element(&self) -> ExtFieldG2Element<N,PARAMSIZE>
    {   
        match self {    ExtG2Field::Fp2_1(x) => x.random_element().into(),
                        ExtG2Field::Fp4_1(x)=> x.random_element().into(),
                        ExtG2Field::Fp8_1(x)=> x.random_element().into(),
                        ExtG2Field::Fp2_2(x) => x.random_element().into(),
                        ExtG2Field::Fp4_2(x)=> x.random_element().into(),
                        ExtG2Field::Fp8_2(x)=> x.random_element().into(),
                        ExtG2Field::Fp4_3(x)=> x.random_element().into(),
                        ExtG2Field::Fp8_3(x)=> x.random_element().into(),
                    }       
    }

    pub fn one(&self) -> ExtFieldG2Element<N,PARAMSIZE>
    {   
        match self {    ExtG2Field::Fp2_1(x) => x.one().clone().into(),
                        ExtG2Field::Fp4_1(x)=> x.one().clone().into(),
                        ExtG2Field::Fp8_1(x)=> x.one().clone().into(),
                        ExtG2Field::Fp2_2(x) => x.one().clone().into(),
                        ExtG2Field::Fp4_2(x)=> x.one().clone().into(),
                        ExtG2Field::Fp8_2(x)=> x.one().clone().into(),
                        ExtG2Field::Fp4_3(x)=> x.one().clone().into(),
                        ExtG2Field::Fp8_3(x)=> x.one().clone().into(),
                    }       
    }

    pub fn zero(&self) -> ExtFieldG2Element<N,PARAMSIZE>
    {   
        match self {    ExtG2Field::Fp2_1(x) => x.zero().clone().into(),
                        ExtG2Field::Fp4_1(x)=> x.zero().clone().into(),
                        ExtG2Field::Fp8_1(x)=> x.zero().clone().into(),
                        ExtG2Field::Fp2_2(x) => x.zero().clone().into(),
                        ExtG2Field::Fp4_2(x)=> x.zero().clone().into(),
                        ExtG2Field::Fp8_2(x)=> x.zero().clone().into(),
                        ExtG2Field::Fp4_3(x)=> x.zero().clone().into(),
                        ExtG2Field::Fp8_3(x)=> x.zero().clone().into(),
                    }       
    }

    pub fn from_hex_strings(&self, source :&[&str]) -> ExtFieldG2Element<N,PARAMSIZE>
    {   
        match self {    ExtG2Field::Fp2_1(x) => x.from_hex_strings(source).into(),
                        ExtG2Field::Fp4_1(x)=> x.from_hex_strings(source).into(),
                        ExtG2Field::Fp8_1(x)=> x.from_hex_strings(source).into(),
                        ExtG2Field::Fp2_2(x) => x.from_hex_strings(source).into(),
                        ExtG2Field::Fp4_2(x)=> x.from_hex_strings(source).into(),
                        ExtG2Field::Fp8_2(x)=> x.from_hex_strings(source).into(),
                        ExtG2Field::Fp4_3(x)=> x.from_hex_strings(source).into(),
                        ExtG2Field::Fp8_3(x)=> x.from_hex_strings(source).into(),
                    }       
    }

    pub fn from_decimal_strings(&self, source :&[&str]) -> ExtFieldG2Element<N,PARAMSIZE>
    {   
        match self {    ExtG2Field::Fp2_1(x) => x.from_strings(source).into(),
                        ExtG2Field::Fp4_1(x)=> x.from_strings(source).into(),
                        ExtG2Field::Fp8_1(x)=> x.from_strings(source).into(),
                        ExtG2Field::Fp2_2(x) => x.from_strings(source).into(),
                        ExtG2Field::Fp4_2(x)=> x.from_strings(source).into(),
                        ExtG2Field::Fp8_2(x)=> x.from_strings(source).into(),
                        ExtG2Field::Fp4_3(x)=> x.from_strings(source).into(),
                        ExtG2Field::Fp8_3(x)=> x.from_strings(source).into(),
                    }       
    }

    pub fn from_basefield_elements(&self, source :&[FieldElement<N>]) -> ExtFieldG2Element<N,PARAMSIZE>
    {   
        match self {    ExtG2Field::Fp2_1(x) => x.from_field_elements(source).into(),
                        ExtG2Field::Fp4_1(x)=> x.from_field_elements(source).into(),
                        ExtG2Field::Fp8_1(x)=> x.from_field_elements(source).into(),
                        ExtG2Field::Fp2_2(x) => x.from_field_elements(source).into(),
                        ExtG2Field::Fp4_2(x)=> x.from_field_elements(source).into(),
                        ExtG2Field::Fp8_2(x)=> x.from_field_elements(source).into(),
                        ExtG2Field::Fp4_3(x)=> x.from_field_elements(source).into(),
                        ExtG2Field::Fp8_3(x)=> x.from_field_elements(source).into(),
                    }       
    }

    pub fn from_i2osp_bytearray(&self, source :&[u8]) -> ExtFieldG2Element<N,PARAMSIZE>
    {   let field =self.basefield();
        
        match self {    ExtG2Field::Fp2_1(x) => 
                                {   let sizeinbytes = source.len()/2;
                                    x.from_field_elements(&[field.from_bigint(&os2ip(&source[0..sizeinbytes]).to_bigint().unwrap()),
                                                                    field.from_bigint(&os2ip(&source[sizeinbytes..]).to_bigint().unwrap()) ]).into()},
                        ExtG2Field::Fp2_2(x) => 
                                {   let sizeinbytes = source.len()/2;
                                    x.from_field_elements(&[field.from_bigint(&os2ip(&source[0..sizeinbytes]).to_bigint().unwrap()),
                                                                    field.from_bigint(&os2ip(&source[sizeinbytes..]).to_bigint().unwrap()) ]).into()},                                                                    
                        ExtG2Field::Fp4_1(x) =>  
                                {   let sizeinbytes = source.len()/4;
                                    x.from_field_elements (&[  field.from_bigint(&os2ip(&source[0..sizeinbytes]).to_bigint().unwrap()),
                                                                       field.from_bigint(&os2ip(&source[sizeinbytes..2*sizeinbytes]).to_bigint().unwrap()),
                                                                       field.from_bigint(&os2ip(&source[2*sizeinbytes..3*sizeinbytes]).to_bigint().unwrap()),
                                                                       field.from_bigint(&os2ip(&source[3*sizeinbytes..]).to_bigint().unwrap()) ]).into()},
                        ExtG2Field::Fp4_2(x) =>  
                        {   let sizeinbytes = source.len()/4;
                            x.from_field_elements (&[  field.from_bigint(&os2ip(&source[0..sizeinbytes]).to_bigint().unwrap()),
                                                                field.from_bigint(&os2ip(&source[sizeinbytes..2*sizeinbytes]).to_bigint().unwrap()),
                                                                field.from_bigint(&os2ip(&source[2*sizeinbytes..3*sizeinbytes]).to_bigint().unwrap()),
                                                                field.from_bigint(&os2ip(&source[3*sizeinbytes..]).to_bigint().unwrap()) ]).into()},                                                                       
                        ExtG2Field::Fp4_3(x) =>  
                                                                {   let sizeinbytes = source.len()/4;
                                                                    x.from_field_elements (&[  field.from_bigint(&os2ip(&source[0..sizeinbytes]).to_bigint().unwrap()),
                                                                                                        field.from_bigint(&os2ip(&source[sizeinbytes..2*sizeinbytes]).to_bigint().unwrap()),
                                                                                                        field.from_bigint(&os2ip(&source[2*sizeinbytes..3*sizeinbytes]).to_bigint().unwrap()),
                                                                                                        field.from_bigint(&os2ip(&source[3*sizeinbytes..]).to_bigint().unwrap()) ]).into()},                                                                                                                                       
                        ExtG2Field::Fp8_1(x)=> 
                            {   let sizeinbytes =source.len()/8;
                                x.from_field_elements (&[  field.from_bigint(&os2ip(&source[0..sizeinbytes]).to_bigint().unwrap()),
                                                                  field.from_bigint(&os2ip(&source[sizeinbytes..2*sizeinbytes]).to_bigint().unwrap()),
                                                                  field.from_bigint(&os2ip(&source[2*sizeinbytes..3*sizeinbytes]).to_bigint().unwrap()),
                                                                  field.from_bigint(&os2ip(&source[3*sizeinbytes..4*sizeinbytes]).to_bigint().unwrap()),
                                                                  field.from_bigint(&os2ip(&source[4*sizeinbytes..5*sizeinbytes]).to_bigint().unwrap()),
                                                                  field.from_bigint(&os2ip(&source[5*sizeinbytes..6*sizeinbytes]).to_bigint().unwrap()),
                                                                  field.from_bigint(&os2ip(&source[6*sizeinbytes..7*sizeinbytes]).to_bigint().unwrap()),
                                                                  field.from_bigint(&os2ip(&source[7*sizeinbytes..8*sizeinbytes]).to_bigint().unwrap()) ]).into()},
                        ExtG2Field::Fp8_2(x)=> 
                        {   let sizeinbytes =source.len()/8;
                            x.from_field_elements (&[  field.from_bigint(&os2ip(&source[0..sizeinbytes]).to_bigint().unwrap()),
                                                            field.from_bigint(&os2ip(&source[sizeinbytes..2*sizeinbytes]).to_bigint().unwrap()),
                                                            field.from_bigint(&os2ip(&source[2*sizeinbytes..3*sizeinbytes]).to_bigint().unwrap()),
                                                            field.from_bigint(&os2ip(&source[3*sizeinbytes..4*sizeinbytes]).to_bigint().unwrap()),
                                                            field.from_bigint(&os2ip(&source[4*sizeinbytes..5*sizeinbytes]).to_bigint().unwrap()),
                                                            field.from_bigint(&os2ip(&source[5*sizeinbytes..6*sizeinbytes]).to_bigint().unwrap()),
                                                            field.from_bigint(&os2ip(&source[6*sizeinbytes..7*sizeinbytes]).to_bigint().unwrap()),
                                                            field.from_bigint(&os2ip(&source[7*sizeinbytes..8*sizeinbytes]).to_bigint().unwrap()) ]).into()},                                                                  
                        ExtG2Field::Fp8_3(x)=> 
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
        match self {    ExtG2Field::Fp2_1(x) => x.field_interface(),
                        ExtG2Field::Fp4_1(x)=> x.field_interface(),
                        ExtG2Field::Fp8_1(x)=> x.field_interface(),
                        ExtG2Field::Fp2_2(x) =><Fp2Field_2<N> as ExtField<PARAMSIZE, 2, N>>::field_interface(x),
                        ExtG2Field::Fp4_2(x)=> x.field_interface(),
                        ExtG2Field::Fp8_2(x)=> x.field_interface(),
                        ExtG2Field::Fp4_3(x)=> x.field_interface(),
                        ExtG2Field::Fp8_3(x)=> x.field_interface(),
        }       
    }
    pub fn hash_to_field(&self,id : &str, security_level:usize,count :usize) -> Vec<ExtFieldG2Element<N,PARAMSIZE>>
    {   
        let extorder;
        match self {    ExtG2Field::Fp2_1(_) => extorder = 2,
                        ExtG2Field::Fp4_1(_)=> extorder = 4,
                        ExtG2Field::Fp8_1(_)=> extorder = 8,
                        ExtG2Field::Fp2_2(_) => extorder = 2,
                        ExtG2Field::Fp4_2(_)=> extorder = 4,
                        ExtG2Field::Fp8_2(_)=> extorder = 8,
                        ExtG2Field::Fp4_3(_)=> extorder = 4,
                        ExtG2Field::Fp8_3(_)=> extorder = 8,
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