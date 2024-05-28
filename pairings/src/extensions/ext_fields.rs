// Code developed by FARAOUN Kamel Mohamed.
// EEDIS-Laboratory. UDL-University. Algeria
// During May 2024.

use base64::engine::general_purpose;
use base64::Engine;
use num_bigint::{BigInt, ToBigInt};
use crate::tools::hashs::{ i2osp_pf, os2ip};
use crate::{fields::arithmetic, tools::exponent::Exponent};
use crate::tools::arithmetic_interface::ArithmeticOperations;
use super::super::fields::prime_fields::{FieldElement,PrimeField};

#[derive(Clone,Debug)]
pub struct ExFieldConsts<const PARAMSIZE:usize,const N:usize> 
                {   pub frobinus_consts:[FieldElement<N>;PARAMSIZE],
                    pub u:i128,                    
                }

pub trait  ExtField<const PARAMSIZE:usize,const ORDER :usize, const N:usize>
    {   type ElementType : ExtElement<PARAMSIZE,ORDER,N>;
        type BaseFieldType ; 

        // Interface to the basefield and the constants of the implemented struct     
        fn field_interface(&self) -> PrimeField<N>;   
        fn extconsts_interface(&self) -> Option<&'static ExFieldConsts<PARAMSIZE,N>>;                       

        // Constructor of a new implemented extension from a base field  
        fn new(base_field :&Self::BaseFieldType, consts :Option<&'static ExFieldConsts<PARAMSIZE,N>>) -> Self;

        // Generate a random element from the extention field
        fn random_element(&self) -> Self::ElementType {
            let field = self.field_interface();
            let consts =self.extconsts_interface();
            let mut result =  [FieldElement { mont_limbs: [0; N],
                                                                        fieldparams: field.parametres};
                                                        ORDER];       
            for i in 0..ORDER{result[i] = FieldElement { mont_limbs: self.field_interface().random_element().mont_limbs,
                                                                fieldparams: field.parametres};         
                                    };
            Self::ElementType::new(&result, consts)
        }

        //  Generate an element in the extention field from a list of BigIntegets
        fn from_big_integers(&self, source :Vec<BigInt>) -> Self::ElementType{
            if source.len() != ORDER {panic!("Size of input does not correspond to the field's extension ...");}        
            let field = self.field_interface();
            let consts =self.extconsts_interface();
            let zero = FieldElement {   mont_limbs: [0; N],
                                                         fieldparams: field.parametres,
                                                     };        
            let mut result: [FieldElement<N>; ORDER] = [zero; ORDER];        
            for i in 0..source.len() {result[i] = field.from_bigint(&source[i]);}        
            Self::ElementType::new(&result, consts)
        }

        //  Generate an element in the extention field from a list of Hexadecimal strings 
        fn from_hex_strings(&self, source :&[&str]) -> Self::ElementType{
            if source.len() != ORDER {panic!("Size of input does not correspond to the field's extension ...");}        
            let field = self.field_interface();
            let consts =self.extconsts_interface();
            let zero = FieldElement {   mont_limbs: [0; N],
                                                        fieldparams: field.parametres,
                                                    };        
            let mut result: [FieldElement<N>; ORDER] = [zero; ORDER];        
            for i in 0..source.len() {result[i] = field.from_hex_str(&source[i]);}        
            Self::ElementType::new(&result,consts)
        }

        //  Generate an element in the extention field from a list of strings (Decimal format)
        fn from_strings(&self, source :&[&str]) -> Self::ElementType{
            if source.len() != ORDER {panic!("Size of input does not correspond to the field's extension ...");}        
            let field = self.field_interface();
            let consts =self.extconsts_interface();
            let zero = FieldElement {   mont_limbs: [0; N],
                                                         fieldparams: field.parametres,
                                                     };        
            let mut result: [FieldElement<N>; ORDER] = [zero; ORDER];        
            for i in 0..source.len() {result[i] = field.from_str(&source[i]);}        
            Self::ElementType::new(&result,consts)
        }

        
        //  Generate an element in the extention field from a list of baseField's elements (Fp)
        fn from_field_elements(&self, source :&[FieldElement<N>]) -> Self::ElementType{
            if source.len() != ORDER {panic!("Size of input does not correspond to the field's extension ...");}        
            let field = self.field_interface();
            let consts =self.extconsts_interface();
            let zero = FieldElement {   mont_limbs: [0; N],
                                                         fieldparams: field.parametres,
                                                     };        
            let mut result: [FieldElement<N>; ORDER] = [zero; ORDER];        
            for i in 0..source.len() {result[i] = source[i];}        
            Self::ElementType::new(&result,consts)
        }

        fn from_byte_array(&self, source :&[u8]) -> Self::ElementType{
            let field = self.field_interface();
            let numbits = field.parametres.num_of_bits;                
            let sizeinbytes = (numbits >> 3) + if (numbits % 8) ==0 {0} else {1}; 
            if sizeinbytes*ORDER != source.len() {panic!("Size of input does not correspond to the field's extension ...");}        
            let consts =self.extconsts_interface();
            let zero = FieldElement {   mont_limbs: [0; N],
                                                         fieldparams: field.parametres,
                                                     };        
            let mut result: [FieldElement<N>; ORDER] = [zero; ORDER];   
            for i in 0..ORDER {result[i] =  field.from_bigint(&os2ip(&source[i*sizeinbytes..(i+1)*sizeinbytes]).to_bigint().unwrap());}        
            Self::ElementType::new(&result,consts)
        }

        fn from_base64(&self, source :&str) -> Self::ElementType{
            let decoded_bytes = match general_purpose::STANDARD.decode(source) {
                Ok(bytes) => bytes,
                Err(_) => {
                    panic!("Failed to decode base64 string");
                }
            };
            self.from_byte_array(&decoded_bytes)
        }

        //  Generate Zero element of the extention field (identity element with respect to Addition)
        fn zero(&self) -> Self::ElementType{
            let field = self.field_interface();
            let consts =self.extconsts_interface();
            Self::ElementType::new(&[FieldElement { mont_limbs: field.zero().mont_limbs,
                                                            fieldparams: field.parametres,
                                                           }; ORDER], consts)
        }

        //  Generate One element of the extention field (identity element with respect to Multiplication)
        fn one(&self) -> Self::ElementType{
            let field = self.field_interface();
            let consts =self.extconsts_interface();
            let mut one = [FieldElement {   mont_limbs: field.zero().mont_limbs,
                                                                      fieldparams: field.parametres,
                                                                  }; ORDER];
            one[0] = FieldElement { mont_limbs: field.one().mont_limbs,
                                    fieldparams: field.parametres,
                                  };
            Self::ElementType::new(&one,consts)
        }        
    }

pub trait ExtElement<const PARAMSIZE:usize,const ORDER :usize, const N:usize>{    
    fn content_interface(&self) -> &[FieldElement<N>;ORDER];
    fn constants_interface(&self) -> Option<&'static ExFieldConsts<PARAMSIZE,N>>;
    fn multiply(&self, rhs:&Self)-> Self;
    fn sqr(&self)-> Self;
    fn new(content :&[FieldElement<N>; ORDER], consts :Option<&'static ExFieldConsts<PARAMSIZE,N>>) -> Self;
    
    fn addto (&self, b : &Self) -> Self where Self: Sized      
            {   let _a = self.content_interface();
                let _b = b.content_interface();
                let zero = FieldElement{mont_limbs:[0;N],fieldparams:_a[0].fieldparams};
                let mut result: [FieldElement<N>; ORDER] =[zero;ORDER] ;                
                for i in 0..ORDER {  result[i] = _a[i].addto(&_b[i])}
                Self::new( &result, self.constants_interface())
            }
    fn substract (&self, b : &Self) -> Self where Self: Sized      
            {   let _a = self.content_interface();
                let _b = b.content_interface();
                let zero = FieldElement{mont_limbs:[0;N],fieldparams:_a[0].fieldparams};
                let mut result: [FieldElement<N>; ORDER] =[zero;ORDER] ;                
                for i in 0..ORDER { result[i] = _a[i].substract(&_b[i])}
                Self::new( &result,self.constants_interface())
            }
    fn negate (&self) -> Self where Self: Sized      
            {   let _a = &self.content_interface();
                let zero = FieldElement{mont_limbs:[0;N],fieldparams:_a[0].fieldparams};
                let mut result: [FieldElement<N>; ORDER] =[zero;ORDER] ;                
                for i in 0..ORDER { result[i] = _a[i].negate()}
                Self::new( &result,self.constants_interface())
            }
    fn double (&self) -> Self where Self: Sized      
            {   let _a = &self.content_interface();
                let zero = FieldElement{mont_limbs:[0;N],fieldparams:_a[0].fieldparams};
                let mut result: [FieldElement<N>; ORDER] =[zero;ORDER] ;                
                for i in 0..ORDER { result[i] = _a[i].double()}
                Self::new( &result,self.constants_interface())
            }            
    fn equal(&self, other: &Self) -> bool {
                let _a = &self.content_interface();
                let _b = other.content_interface();
                let mut eq =true;
                for i in 0..ORDER { eq = eq & _a[i].equal(&_b[i])}
                eq
            }
    fn mulbyu8(&self,  rhs: u8)-> Self where Self: Sized 
            {   let _a = self.content_interface();
                let zero = FieldElement{mont_limbs:[0;N],fieldparams:_a[0].fieldparams};
                match  rhs {    0 => { Self::new(&[zero;ORDER],self.constants_interface())},
                                1 => Self::new(_a,self.constants_interface()),
                                2 => self.addto(&self),
                                3 => self.addto(&self).addto(&self),
                                4 =>{let double =self.addto(&self);
                                        double.addto(&double)},      
                                5 =>{let double =self.addto(&self);
                                        let fourth =double.addto(&double);
                                        fourth.addto(&self) },
                                _ =>{  let mut result: [FieldElement<N>; ORDER] =[zero;ORDER] ; 
                                        for i in 0..ORDER{  result[i] = rhs * _a[i]};
                                        Self::new(&result,self.constants_interface())
                                    }            
                            }        
            }
    fn mulby_fp_element (&self, b : &FieldElement<N>) -> Self where Self: Sized      
            {   let _a = self.content_interface();
                let zero = FieldElement{mont_limbs:[0;N],fieldparams:_a[0].fieldparams};
                let mut result: [FieldElement<N>; ORDER] =[zero;ORDER] ;                
                for i in 0..ORDER {  result[i] = _a[i].multiply(&b)}
                Self::new( &result, self.constants_interface())
            }
            
    fn pow(&self,e :& dyn Exponent<N>)-> Self where Self :Sized 
            {   let _a = self.content_interface();
                let out_params=  _a[0].fieldparams;
                let one = FieldElement{mont_limbs:out_params.one,fieldparams:_a[0].fieldparams};
                let zero = FieldElement{mont_limbs:[0;N],fieldparams:_a[0].fieldparams};
                let mut result: [FieldElement<N>; ORDER] =[zero;ORDER] ;
                result [0] = one;   
                let mut result = Self::new(&result,self.constants_interface());                            
                if let Some(array) = e.to_u64_array() {   let limbnum= e.get_len();
                                                                    for i in array[0..limbnum].as_ref().iter().rev()  {
                                                                        for j in (0..64).rev(){ result = result.sqr();                
                                                                                                     if (i >> j) & 1 == 1 { result = result.multiply(&self);}                
                                                                                                   }
                                                                                            }                                                        
                                                                }                                                                 
                result
            }
    
    fn to_a_string(&self) -> String
            {  let mut out= String::new();
               let _a = self.content_interface();
               out.push('(');               
               for i in 0..ORDER-1 { out.push_str(&_a[i].to_string());
                                        out.push_str(" , \n");
                                    }      
                out.push_str(&_a[ORDER - 1].to_string());                                        
                out.push_str(")");
                String::from(&out)
            }
    
    fn to_i2osp_bytearray(&self) -> Vec<u8>
            {  let mut out= Vec::<u8>::new();
               let _a = self.content_interface();       
               let numbits = _a[0].fieldparams.num_of_bits;                
               let sizeinbytes = (numbits >> 3) + if (numbits % 8) ==0 {0} else {1};     
               for i in 0..ORDER { out.extend(i2osp_pf(&_a[i], sizeinbytes)) ;}      
               out 
            }
    
    fn to_hex_string(&self) -> String
            {  let mut out= String::new();
               let _a = self.content_interface();
               out.push('(');               
               for i in 0..ORDER-1 { out.push_str(&_a[i].to_hex_string());
                                        out.push_str(" , \n");
                                    }      
                out.push_str(&_a[ORDER - 1].to_hex_string());                                        
                out.push_str(")");
                String::from(&out)
            }

    fn encode_to_base64(&self) -> String
            {  general_purpose::STANDARD.encode(self.to_i2osp_bytearray())
            }

    fn is_one(&self) -> bool
            {   let a = self.content_interface();
                let mut iszero = a[0].is_one();                
                for i in &a[1..] {iszero = iszero & i.is_zero(); }
                iszero
            }
    
    fn is_zero(&self) -> bool
            {   let a = self.content_interface();
                let mut iszero = a[0].is_zero();                
                for i in &a[1..] {iszero = iszero & i.is_zero(); }
                iszero
            }
    fn zero(&self) -> Self where Self: Sized  
            {
                let content =self.content_interface();
                Self::new(&[FieldElement { mont_limbs: content[0].fieldparams.zero,
                                                                fieldparams: content[0].fieldparams,
                                                  }; ORDER], self.constants_interface())
            }
    fn one(&self) -> Self where Self: Sized  
            {
                let content =self.content_interface();
                let mut result = [FieldElement { mont_limbs: content[0].fieldparams.zero,fieldparams: content[0].fieldparams}; ORDER];
                result[0] = FieldElement { mont_limbs: content[0].fieldparams.one,fieldparams: content[0].fieldparams};
                Self::new(&result, self.constants_interface())
            }         
    fn sign(&self) -> i8
    {
            // sgn0 "sign" of x: returns -1 if x is lexically larger than  -x and, else returns 1
            // https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-hash-to-curve-05#section-4.1
            let mut one :[u64;N] = [0;N];
            one[0] = 1;
            let content =self.content_interface();
            let mut sig:i8 = 0;
            for i in (0..ORDER).rev(){   let mut sigi:i8 ;
                                                let val = arithmetic::mul(&content[i].mont_limbs, &one, &content[0].fieldparams);                
                                                let mut j = N - 1;
                                                while (val[j] == content[0].fieldparams.sig_theshold[j]) & (j > 0) { j = j - 1}
                                                sigi = if val[j] > content[0].fieldparams.sig_theshold[j] {-1} else {1};        
                                                if content[i].is_zero() {sigi =0};
                                                if sig==0 {sig = sigi}
                                            }
            if sig==0 {sig = 1};                                                    
            sig                                            
    }
}



