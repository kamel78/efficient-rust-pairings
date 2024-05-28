use std::{io, time::{Duration, Instant}};
use base64::{engine::general_purpose, Engine};
use pairings::{BLS24Curves, BLS48Curves, Bls12Curves, PairingsEngine};

fn test_timing(){
        fn measure_time<F>(f: F) -> Duration
        where
            F: Fn(),            
        {   let trys_count = 100;
            let start = Instant::now();
            for _ in 0..trys_count {f();}
            start.elapsed()/trys_count
        }
        let en = pairings::BLS12::_381();
        println!("------- BLS12-361 --------------------------------");
        let p = en.g1.hash_to_field("identity 1", 0);
        let q = en.g2.hash_to_field("identity 2", 0);
        let duration = measure_time(|| {en.paire(&p, &q);});
        println!("Time elapsed for pairings : {:?}", duration);
        let duration = measure_time(|| {en.miller_loop(&p, &q);});
        println!("Time elapsed for Miller loop : {:?}", duration);
        let e =en.miller_loop(&p, &q);
        let duration = measure_time(|| {let _ =e.final_exponentiation();});
        println!("Time elapsed for Final Exponentiation: {:?}", duration);
        let duration = measure_time(|| {let _ = en.g1.hash_to_field("identity 1", 0);});
        println!("Time elapsed for hashing to G1 : {:?}", duration);
        let duration = measure_time(|| {let _ = en.g2.hash_to_field("identity 1", 0);});
        println!("Time elapsed for hashing to G2 : {:?}", duration);
        let r = en.fr.random_element();
        let duration = measure_time(|| {let _ = r*p;});
        println!("Time elapsed multiplication on G1 (GLV): {:?}", duration);
        let duration = measure_time(|| {let _ = r*q;});
        println!("Time elapsed for multiplication on G2 (GLS): {:?}", duration);

        let en = pairings::BLS12::_461();
        println!("------- BLS12-461 --------------------------------");
        let p = en.g1.hash_to_field("identity 1", 0);
        let q = en.g2.hash_to_field("identity 2", 0);
        let duration = measure_time(|| {en.paire(&p, &q);});
        println!("Time elapsed for pairings : {:?}", duration);
        let duration = measure_time(|| {en.miller_loop(&p, &q);});
        println!("Time elapsed for Miller loop : {:?}", duration);
        let e =en.miller_loop(&p, &q);
        let duration = measure_time(|| {let _ =e.final_exponentiation();});
        println!("Time elapsed for Final Exponentiation: {:?}", duration);
        let duration = measure_time(|| {let _ = en.g1.hash_to_field("identity 1", 0);});
        println!("Time elapsed for hashing to G1 : {:?}", duration);
        let duration = measure_time(|| {let _ = en.g2.hash_to_field("identity 1", 0);});
        println!("Time elapsed for hashing to G2 : {:?}", duration);
        let r = en.fr.random_element();
        let duration = measure_time(|| {let _ = r*p;});
        println!("Time elapsed multiplication on G1 (GLV): {:?}", duration);
        let duration = measure_time(|| {let _ = r*q;});
        println!("Time elapsed for multiplication on G2 (GLS): {:?}", duration);
        
        let en = pairings::BLS24::_479();
        println!("------- BLS24-479 --------------------------------");
        let p = en.g1.hash_to_field("identity 1", 0);
        let q = en.g2.hash_to_field("identity 2", 0);
        let duration = measure_time(|| {en.paire(&p, &q);});
        println!("Time elapsed for pairings : {:?}", duration);
        let duration = measure_time(|| {en.miller_loop(&p, &q);});
        println!("Time elapsed for Miller loop : {:?}", duration);
        let e =en.miller_loop(&p, &q);
        let duration = measure_time(|| {let _ =e.final_exponentiation();});
        println!("Time elapsed for Final Exponentiation: {:?}", duration);
        let duration = measure_time(|| {let _ = en.g1.hash_to_field("identity 1", 0);});
        println!("Time elapsed for hashing to G1 : {:?}", duration);
        let duration = measure_time(|| {let _ = en.g2.hash_to_field("identity 1", 0);});
        println!("Time elapsed for hashing to G2 : {:?}", duration);
        let r = en.fr.random_element();
        let duration = measure_time(|| {let _ = r*p;});
        println!("Time elapsed multiplication on G1 (GLV): {:?}", duration);
        let duration = measure_time(|| {let _ = r*q;});
        println!("Time elapsed for multiplication on G2 (GLS): {:?}", duration);

        let en = pairings::BLS24::_559();
        println!("------- BLS24-559--------------------------------");
        let p = en.g1.hash_to_field("identity 1", 0);
        let q = en.g2.hash_to_field("identity 2", 0);
        let duration = measure_time(|| {en.paire(&p, &q);});
        println!("Time elapsed for pairings : {:?}", duration);
        let duration = measure_time(|| {en.miller_loop(&p, &q);});
        println!("Time elapsed for Miller loop : {:?}", duration);
        let e =en.miller_loop(&p, &q);
        let duration = measure_time(|| {let _ =e.final_exponentiation();});
        println!("Time elapsed for Final Exponentiation: {:?}", duration);
        let duration = measure_time(|| {let _ = en.g1.hash_to_field("identity 1", 0);});
        println!("Time elapsed for hashing to G1 : {:?}", duration);
        let duration = measure_time(|| {let _ = en.g2.hash_to_field("identity 1", 0);});
        println!("Time elapsed for hashing to G2 : {:?}", duration);
        let r = en.fr.random_element();
        let duration = measure_time(|| {let _ = r*p;});
        println!("Time elapsed multiplication on G1 (GLV): {:?}", duration);
        let duration = measure_time(|| {let _ = r*q;});
        println!("Time elapsed for multiplication on G2 (GLS): {:?}", duration);
        
        let en = pairings::BLS48::_575();
        println!("------- BLS48-575 --------------------------------");
        let p = en.g1.hash_to_field("identity 1", 0);
        let q = en.g2.hash_to_field("identity 2", 0);
        let duration = measure_time(|| {en.paire(&p, &q);});
        println!("Time elapsed for pairings : {:?}", duration);
        let duration = measure_time(|| {en.miller_loop(&p, &q);});
        println!("Time elapsed for Miller loop : {:?}", duration);
        let e =en.miller_loop(&p, &q);
        let duration = measure_time(|| {let _ =e.final_exponentiation();});
        println!("Time elapsed for Final Exponentiation: {:?}", duration);
        let duration = measure_time(|| {let _ = en.g1.hash_to_field("identity 1", 0);});
        println!("Time elapsed for hashing to G1 : {:?}", duration);
        let duration = measure_time(|| {let _ = en.g2.hash_to_field("identity 1", 0);});
        println!("Time elapsed for hashing to G2 : {:?}", duration);
        let r = en.fr.random_element();
        let duration = measure_time(|| {let _ = r*p;});
        println!("Time elapsed multiplication on G1 (GLV): {:?}", duration);
        let duration = measure_time(|| {let _ = r*q;});
        println!("Time elapsed for multiplication on G2 (GLS): {:?}", duration);
}

fn bls_signature_scheme(){
    // BLS Signature scheme : public keys are in G2 while siugnatures are in G1
    // An reverted scheme can also be implemented (public keys in G1 and signatures in G2) 
    // according to the targted application resuirements. 
    // ps: Hashing to elliptic curves support two modes :  
    //      mode = 0 stands for Non-uniform-Encoding, while mode = 1 stands for Random Oracle Model encoding 
    let engine = pairings::BLS12::_461();
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
}

fn  boneh_franklin_ibe(){
    let engine = pairings::BLS12::_461();
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
} 

fn main() {
  loop {  println!("Please enter a choice (1, 2, or 3) for the following routines, or 4 to exit:");
    println!("(1)- Runtime bench-marking of several implemented functionalities (please run in '--release' mode for accurate results).");
    println!("(2)- Demonstration of the BLS signature scheme (KeyGen,signature and verification).");
    println!("(3)- Demonstration of the Boneh-Franklin Identity Based Encryption scheme.");
    println!("Enter 4 to leave ...");

    let mut input = String::new();
    io::stdin().read_line(&mut input).expect("Failed to read line");

    let choice: u32 = match input.trim().parse() {
        Ok(num) => num,
        Err(_) => {
            println!("Invalid input. Please enter a number.");
            return;
        }
    };

    match choice {
        1 => {println!("Runtime Bench-marking\n");
              test_timing()  },
        2 => {println!("BLS Signature demonstration :\n");
              bls_signature_scheme()  },
        3 => {println!("Boneh-Franklin IBE demonstration :\n");
              boneh_franklin_ibe()   },
        4 =>{println!("Exiting the program. Goodbye!");
                break;
            },
        _ => println!("Invalid choice. Please enter 1, 2, or 3."),
    }}
    
}





        