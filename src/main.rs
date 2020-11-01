extern crate bio;
use bio::io::fasta;
use std::collections::HashMap;
use std::io;
use std::io::Write;
use std::env;
use std::process;

// Translation table 11 from https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi. 
// Indexing into this array is as follows: 
// T = 0, C = 1, A = 2, G = 3
// index = (16 * first_base) + (4 * second_base) + third_base 
// So... TTT = 0, TTC = 1, TTA = 2, ... , GGC = 61, GGA = 62, GGG = 63
const GENETIC_CODE: &str = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

// ERR_BAD_NT is an error value for an invalid nucleotide 
const ERR_BAD_NT: usize = 99; 

// enumeration of genetic sequence types
enum SeqType {
    DNA,
    Protein1, // single-letter amino acid code
    Protein3, // triple-letter amino acid code
}

// enumeration of translation types
enum Translation {
    OneLetter,
    ThreeLetter,
}

// given an input 'char', return a base equivalent
fn lookup(x: char) -> usize {
    match x {
        'T' => return 0,
        'C' => return 1,
        'A' => return 2,
        'G' => return 3,
        _ => return ERR_BAD_NT, // unknown base
    }
}

// translate a codon into its corresponding amino acid
fn translate(triplet: &str, t: Translation) -> String {
    let three_letter_code: HashMap<char, &str> = [
    ('A', "Ala"), 
    ('B', "???"), 
    ('C', "Cys"), 
    ('D', "Asp"), 
    ('E', "Glu"), 
    ('F', "Phe"), 
    ('G', "Gly"), 
    ('H', "His"), 
    ('I', "Ile"), 
    ('J', "???"), 
    ('K', "Lys"), 
    ('L', "Leu"), 
    ('M', "Met"), 
    ('N', "Asn"), 
    ('O', "Pyr"), 
    ('P', "Pro"), 
    ('Q', "Gln"), 
    ('R', "Arg"), 
    ('S', "Ser"), 
    ('T', "Thr"), 
    ('U', "Sel"), 
    ('V', "Val"), 
    ('W', "Trp"), 
    ('X', "???"), 
    ('Y', "Tyr"), 
    ('Z', "???"), 
    ('*', "***"), 
    ].iter().copied().collect();

    let mut codon = vec![ERR_BAD_NT; 3];

    for (i,base) in triplet.chars().enumerate() {
        codon[i] = lookup(base);
    }

    if codon.contains(&ERR_BAD_NT) {
        println!("Something went wrong in translation!");
        println!("{:?}", codon);
        process::exit(1);
    }

    let index: usize = (codon[0] * 16) + (codon[1] * 4) + codon[2];
    // translate the codon into single-letter code
    let c = GENETIC_CODE.chars().nth(index).unwrap();
    match t {
        Translation::OneLetter => return c.to_string(),
        Translation::ThreeLetter => return three_letter_code[&c].to_string(),
    }
}

// h/t https://medium.com/python-in-plain-english/from-python-to-rust-part-3-f035d780de1
fn reverse_complement(s: &str) -> String {
    let complements: HashMap<char, char> = [
        ('A', 'T'),
        ('C', 'G'),
        ('G', 'C'),
        ('T', 'A'),
    ].iter().copied().collect();

    let mut rev_comp = String::new();

    // iterate over the sequence in reverse, and add its complement to rev_comp
    for base in s.chars().rev() {
        rev_comp.push(complements[&base]);
    }

    return rev_comp;
}

// print a pretty sequence, 72 bases per line, plus base numbering
// s: sequence
// t: sequence type (DNA, Protein1 or Protein3)
fn print_seq(s: &str, t: SeqType) {
    let linelen = 72;

    let divisor = match t {
        // if we're printing a Protein3, count amino acids, not bases, 
        // so divide by 3
        SeqType::DNA | SeqType::Protein1 =>  1,
        SeqType::Protein3 => 3,
    };

    for i in 0..(s.len()/linelen) {
        let myline = &s[i*linelen..(i*linelen)+linelen];
        println!("{number:>0width$} {}", myline, number=(((i*linelen)/divisor))+1, width=3);
    }
    // print whatever is left
    if s.len() % linelen != 0 {
        println!("{number:>0width$} {}", &s[s.len()-s.len()%linelen..s.len()],
        number=((s.len()-s.len()%linelen))/divisor+1, width=3)
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut filename = "sequence.fasta";

    // the user can provide another FASTA file on the command line
    let args: Vec<String> = env::args().collect();
    if args.len() > 1 {
        filename = &args[1];
    }
    if !std::path::Path::new(filename).exists() {
        println!("File '{}' does not exist.", filename);
        process::exit(1);
    }

    println!("Reading FASTA records from file '{}'...", filename);
    let reader = fasta::Reader::from_file(filename)?;

    for record in reader.records() {
        let record = record.expect("Error during FASTA record parsing");
        println!("\nSequence ID: {}", record.id());
        println!("Sequence description:\n{}", record.desc().unwrap());

        if record.id() == "NC_005816.1" {
            println!("\nFrom 'https://www.ncbi.nlm.nih.gov/nuccore/NC_005816',");
            println!("we know that this piece of DNA encodes 9 genes.\n");
            println!("1) YP_RS22210: IS21-like element IS100 family transposase");
            println!("2) YP_RS22215: AAA family ATPase");
            println!("3) YP_RS22220: Rop family plasmid primer RNA-binding protein");
            println!("4) YP_RS22225: pesticin immunity protein");
            println!("5) YP_RS22230: pesticin");
            println!("6) YP_RS22235: hypothetical protein");
            println!("7) YP_RS22240: omptin family plasminogen activator Pla");
            println!("8) YP_RS22245: XRE family transcriptional regulator");
            println!("9) YP_RS22250: type II toxin-antitoxin system RelE/ParE family toxin");

            print!("\nWhich one would you like to view? > ");
            // flush buffer...
            io::stdout().flush().unwrap();

            // get input...
            let mut retval = String::new();
            io::stdin().read_line(&mut retval).expect("Failed to read from stdin");

            // use trim() to delete the trailing newline ('\n') char
            retval = retval.trim().to_string();

            let selection = match retval.parse::<usize>() {
                Ok(i) => i, // if good input, just return the number
                Err(_) => {
                    println!("Invalid input: '{}'", retval);
                    process::exit(1);
                },
            };

            let myseq = match selection {
                1 => &record.seq()[86..1109],
                2 => &record.seq()[1108..1888],
                3 => &record.seq()[2924..3119],
                4 => &record.seq()[4354..4780],
                5 => &record.seq()[4814..5888],
                6 => &record.seq()[6115..6421],
                7 => &record.seq()[6663..7602],
                8 => &record.seq()[7788..8088],
                9 => &record.seq()[8087..8429],
                _ => { 
                    println!("Invalid input: '{}'", selection);
                    process::exit(1);
                },
            };
            
            // convert the sequence record into a String
            let mut s = String::from_utf8(myseq.to_vec())?;

            // genes 5, 8 and 9 are on the 'reverse' strand, so we need their reverse complements
            if (selection == 5) | (selection == 8) | (selection == 9) {
                s = reverse_complement(&s)
            }

            println!("\nDNA sequence:");
            print_seq(&s, SeqType::DNA);
            println!("Length: {}\n", s.len());
        
            assert!(s.len() % 3 == 0, "Sequence length is not a multiple of 3!");
            let mut peptide1 = String::new();    
            let mut peptide3 = String::new();    
            let n_codons = s.len()/3;
            for i in 0..n_codons {
                let codon = &s[i*3..(i*3)+3]; // take a 3-base slice of the sequence
                peptide1.push_str(&translate(&codon, Translation::OneLetter)); // translate and add to the string
                peptide3.push_str(&translate(&codon, Translation::ThreeLetter)); // translate and add to the string
            }
            println!("One-letter code:");
            print_seq(&peptide1, SeqType::Protein1);
            println!("\nThree-letter code:");
            print_seq(&peptide3, SeqType::Protein3);
            println!("Length: {}\n", n_codons);
        }
    }

    Ok(()) 
}


// tests
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_translate_atg() {
        assert_eq!(translate("ATG", Translation::ThreeLetter), "Met");
    }

    #[test]
    fn test_translate_tag() {
        assert_eq!(translate("TAG", Translation::ThreeLetter), "***");
    }

    #[test]
    fn test_translate_ttt() {
        assert_eq!(translate("TTT", Translation::OneLetter), "F");
    }

    #[test]
    fn bad_translation_atg() {
        assert_eq!(translate("ATG", Translation::ThreeLetter), "Phe");
    }

    #[test]
    fn test_reverse_complement1() {
        assert_eq!(reverse_complement("ATG"), "CAT");
    }

    #[test]
    fn test_reverse_complement2() {
        assert_eq!(reverse_complement("AAAGGGAAATTT"), "AAATTTCCCTTT")
    }
}
