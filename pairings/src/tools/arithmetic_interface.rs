// Code developed by FARAOUN Kamel Mohamed.
// EEDIS-Laboratory. UDL-University. Algeria
// During May 2024.

pub trait ArithmeticOperations {
    fn addto(&self, other: &Self) -> Self;
    fn double(&self) -> Self;
    fn substract(&self, other: &Self) -> Self;
    fn multiply(&self, other: &Self) -> Self;
    fn sqr(&self) -> Self;
    fn invert(&self) -> Self;
    fn negate(&self) -> Self;
    fn equal(&self, rhs :&Self) -> bool;
    fn is_zero(&self) -> bool;
    fn is_one(&self) -> bool; 
    fn to_dec_string(&self) -> String;
    fn to_hex_string(&self) -> String;
    fn one(&self) -> Self;
    fn zero(&self) -> Self;
}