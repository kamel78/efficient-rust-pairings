use std::io::{self};
use base64::{engine::general_purpose, Engine};
use pairings::{tests::{bench_pairings, check_pairings}, BLS24Curves, BLS48Curves, Bls12Curves, Pairings, PairingsEngine, BLS12, BLS24, BLS48};


fn bls_signature_scheme<const R:usize,const N:usize,const MAX_COEFS_COUNT1:usize,const PRAMASIZE:usize, const MAX_COEFS_COUNT2 :usize>
      (engine :&Pairings <R, N, MAX_COEFS_COUNT1, PRAMASIZE, MAX_COEFS_COUNT2>){
  // BLS Signature scheme : public keys are in G2 while siugnatures are in G1
  // An reverted scheme can also be implemented (public keys in G1 and signatures in G2) 
  // according to the targted application resuirements. 
  // ps: Hashing to elliptic curves support two modes :  
  //      mode = 0 stands for Non-uniform-Encoding, while mode = 1 stands for Random Oracle Model encoding 
  println!("----------------------------------------------------------------------------------------------");
  println!("Testing BLS signature scheme on {}",engine.identifier);
  // Key-paire genration :
  let sk = engine.fr.random_element();
  let pk = sk * engine.g2.default_generator();
  println!(" Secrete key (base64) = {}", sk.to_base64());
  println!(" Public Key  (base64) = {}", pk.encode_to_base64());       
  // BLS Signature : 
  let message = "This is a simple message to be signed. A message can be any arbitrary length string ....";
  let hashed_message = engine.g1.hash_to_field(&message, 0);
  let signature = sk * hashed_message;
  println!(" Signatue is (base64): {}", signature.encode_to_base64());
  // BLS Verification :
  let hashed_message = engine.g1.hash_to_field(&message, 0);
  let verification_result = engine.paire(&signature, &engine.g2.default_generator()) == engine.paire(&hashed_message, &pk);
  println!("Verification result : {}",if verification_result {"correct"} else {"incorrect"});
  // Faster way to do it using multi-pairings
  let verification_result = engine.multi_paire(&[signature,hashed_message], &[-engine.g2.default_generator(),pk]) == engine.gt.one();
  println!("Verification result : {}",if verification_result {"correct"} else {"incorrect"});       
  println!("----------------------------------------------------------------------------------------------\n");
}



fn  boneh_franklin_ibe<const R:usize,const N:usize,const MAX_COEFS_COUNT1:usize,const PRAMASIZE:usize, const MAX_COEFS_COUNT2 :usize>(
                engine :&Pairings <R, N, MAX_COEFS_COUNT1, PRAMASIZE, MAX_COEFS_COUNT2>){
  println!("----------------------------------------------------------------------------------------------");
  println!("Testing Boneh-Franklin Identity Based Encryption scheme on {}",engine.identifier);
  // Generation of Master Keys (Setup):
  let msk =  engine.fr.random_element();
  let mpk = msk * engine.g2.default_generator();
  println!("The Master secrete key : {} ",msk.to_base64());
  println!("The Master public key : {} \n",mpk.encode_to_base64());

  // Key extraction : generation of the user's secrete key for corresponding Identity :
  let user_identity ="ID-1";    
  let id_sk = msk * engine.g1.hash_to_field(&user_identity, 0);
  println!("User's secrete key for identity '{}' : {} \n",user_identity,id_sk.encode_to_base64()); 
  
  // Key confirmation : user can confirm the authenticity and corectness of the secrete key like follows: 
  let valide_secrete_key = engine.paire(&id_sk, &engine.g2.default_generator()) 
                                 == engine.paire(&engine.g1.hash_to_field(&user_identity, 0), &mpk);
  println!("User's secrete key confirmation : {} ",if valide_secrete_key {"Valid key\n"} else {"Invalid key\n"}); 

  //  Encryption of a message to the user using its Identity :
  let message ="This is a simple message to be signed. A message can be any arbitrary length string ....";
  println!("Plaintext message : {}\n",message);
  let message_as_bytes: Vec<u8> = message.as_bytes().to_vec();
  let a = engine.fr.random_element();
  let u = a * engine.g2.default_generator();
  let key_stream = engine.paire(&engine.g1.hash_to_field(&user_identity, 0),&mpk)
                            .pow(&a).derive_hkdf(8*message_as_bytes.len(), None);
  let encrypted_data: Vec<u8> = key_stream.iter().zip(message_as_bytes.iter()).map(|(&x1, &x2)| x1 ^ x2).collect();    
  let encrypted_message =[u.encode_to_base64(),general_purpose::STANDARD.encode(encrypted_data)];
  println!("Encrypted message : {:?}\n",encrypted_message);

  // Decryption of the message using the user's secrete key 
  let decoded_encryption = general_purpose::STANDARD.decode(&encrypted_message[1]).unwrap();
  let u = engine.g2.from_base64(&encrypted_message[0]);
  let key_stream = engine.paire(&id_sk, &u).derive_hkdf(8*decoded_encryption.len(), None); 
  let decrypted_message : Vec<u8> = key_stream.iter().zip(decoded_encryption.iter()).map(|(&x1, &x2)| x1 ^ x2).collect();    
  println!("Decrypted message : {}",std::str::from_utf8(&decrypted_message).unwrap());
  println!("----------------------------------------------------------------------------------------------\n");

} 

fn check_and_benchmark<const R:usize,const N:usize,const MAX_COEFS_COUNT1:usize,const PRAMASIZE:usize, const MAX_COEFS_COUNT2 :usize>(
  engine :&Pairings <R, N, MAX_COEFS_COUNT1, PRAMASIZE, MAX_COEFS_COUNT2>)
{
  println!("Checking pairins for {} : {}",engine.identifier,if check_pairings(&engine.curvename) {"correct "}  else {"Incorrect"}) ;  
  println!("Benchmark pairins for {:?} :",engine.identifier);
  bench_pairings(&engine.curvename);
}

fn run_chosen_action<const R:usize,const N:usize,const MAX_COEFS_COUNT1:usize,const PRAMASIZE:usize, const MAX_COEFS_COUNT2 :usize>(
  engine :&Pairings <R, N, MAX_COEFS_COUNT1, PRAMASIZE, MAX_COEFS_COUNT2>,choice :u32)
{
  match choice {
    1 => {println!("Runtime Bench-marking\n");
          check_and_benchmark(engine)  },
    2 => {println!("BLS Signature demonstration :\n");
          bls_signature_scheme(BLS12::_381())  },
    3 => {println!("Boneh-Franklin IBE demonstration :\n");
          boneh_franklin_ibe(BLS12::_461())   },
    4 =>{println!("Exiting the program. Goodbye!");            
        },
    _ => println!("Invalid choice. Please enter 1, 2, or 3."),
  }
}

fn main() {    

    loop {  println!("Please enter a choice (1, 2, or 3) for the following routines, or 4 to exit:");
    println!("(1)- Runtime bench-marking of several implemented functionalities (please run in '--release' mode for accurate results).");
    println!("(2)- Demonstration of the BLS signature scheme (KeyGen,signature and verification).");
    println!("(3)- Demonstration of the Boneh-Franklin Identity Based Encryption scheme.");
    println!("Enter 4 to leave ...");

    let mut input = String::new();
    io::stdin().read_line(&mut input).expect("Failed to read line");

    let choice1: u32 = match input.trim().parse() {
        Ok(num) => num,
        Err(_) => {
            println!("Invalid input. Please enter a number.");
            return;
        }
    };
    if choice1 == 4 {break;}
    println!(" Choose one of the following implemented curves :");
    println!(" --------------------------------------------------");
    println!("(1)- BLS12-381 (Security :128bit).");
    println!("(2)- BLS12-446 (Security :128bit).");
    println!("(3)- BLS12-461 (Security :128bit).");
    println!(" --------------------------------------------------");
    println!("(4)- BLS24-315 (Security :128bit).");
    println!("(5)- BLS24-477 (Security :192bit).");
    println!("(6)- BLS24-479 (Security :192bit).");
    println!("(7)- BLS24-509 (Security :192bit).");
    println!("(8)- BLS24-509-SNARK (Security :192bit).");
    println!("(9)- BLS24-559 (Security :192bit).");
    println!(" --------------------------------------------------");
    println!("(10)- BLS48-277 (Security :128bit).");
    println!("(11)- BLS48-287 (Security :128bit).");
    println!("(12)- BLS48-571 (Security :256bit).");
    println!("(13)- BLS48-573 (Security :256bit).");
    println!("(14)- BLS48-575 (Security :256bit).");
    println!("(15)- BLS48-581 (Security :256bit).");
    println!("Choose a curve :");

    let mut input = String::new();
    io::stdin().read_line(&mut input).expect("Failed to read line");
    let choice2: u32 = match input.trim().parse() {
      Ok(num) => num,
      Err(_) => {
          println!("Invalid input. Please enter a number.");
          return;
      }
  };
      match  choice2 {
              1=> run_chosen_action(BLS12::_381(), choice1),
              2=> run_chosen_action(BLS12::_446(), choice1),
              3=> run_chosen_action(BLS12::_461(), choice1),
              4=> run_chosen_action(BLS24::_315(), choice1),
              5=> run_chosen_action(BLS24::_477(), choice1),
              6=> run_chosen_action(BLS24::_479(), choice1),
              7=> run_chosen_action(BLS24::_509(), choice1),
              8=> run_chosen_action(BLS24::_509_snark(), choice1),
              9=> run_chosen_action(BLS24::_559(), choice1),
              10=> run_chosen_action(BLS48::_277(), choice1),
              11=> run_chosen_action(BLS48::_287(), choice1),
              12=> run_chosen_action(BLS48::_571(), choice1),
              13=> run_chosen_action(BLS48::_573(), choice1),
              14=> run_chosen_action(BLS48::_575(), choice1),
              15=> run_chosen_action(BLS48::_581(), choice1),
              _ => println!("Invalid choice. Please enter valid choice."),
                
            }  
}

}





        
