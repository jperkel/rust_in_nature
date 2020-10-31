extern crate bio;
use plotters::prelude::*;
use bio::io::fasta;
use std::collections::HashMap;

// String containing the standard genetic code. Indexing into this array is
// as follows: 
// A = 0, C = 1, G = 2, T = 3
// index = (16 * first_base) + (4 * second_base) + third_base 
// So... AAA = 0, AAC = 1, AAG = 2, ... , TTC = 61, TTG = 62, TTT = 63
const GENETIC_CODE: &str = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CYCLFLF"; 

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

// calculate the first x Fibonacci numbers
fn fibonacci(x: usize) -> Vec<usize> {
    assert!(x >= 2, "x must be at least 2!");
    let mut f1 = 0;
    let mut f2 = 1;
    
    // a vector to hold the computed values. We know the first two numbers
    let mut ar = vec![0; x]; 
    ar[0] = f1;
    ar[1] = f2;

    for i in 2..x {
        let f3 = f1 + f2;
        ar[i] = f3;
        f1 = f2;
        f2 = f3;
    }
    return ar;
}

// graph the Fibonacci sequence
fn fibonacci_plot(f: Vec<usize>) -> Result<(), Box<dyn std::error::Error>> {
    let n = f.len();
    // a vector of points to plot
    let mut v = vec![(0.0,0.0); n];
    let filename = "fibonacci.png";

    for i in 0..n {
        v[i] = (i as f32, f[i] as f32);
    }

    let root = BitMapBackend::new(filename, (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption(
            format!("The first {} Fibonacci numbers", n), 
            ("sans-serif", 30).into_font()
        )
        .margin(5)
        .x_label_area_size(30)
        .y_label_area_size(60)
        .build_cartesian_2d(0f32..27f32, 0f32..52500f32)?;

    chart.configure_mesh()
        .x_desc("Sequence no.")
        .y_desc("Fibonacci no.")
        .draw()?;

    chart
        .draw_series(v.iter()
        .map(|(x,y)| Circle::new((*x,*y), 3, BLUE.filled())
    ))?;

    chart
        .draw_series(LineSeries::new(v.iter()
            .map(|(x,y)| (*x,*y)), &RED,
    ))?;
    println!("Graph written to file '{}'.\n", filename);

    Ok(()) 
}

// given an input 'char', return a base equivalent
fn lookup(x: char) -> usize {
    match x {
        'A' => return 0,
        'C' => return 1,
        'G' => return 2,
        'T' => return 3,
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
        panic!();
    }

    let index: usize = (codon[0] * 16) + (codon[1] * 4) + codon[2];
    // translate the codon into single-letter code
    let c = GENETIC_CODE.chars().nth(index).unwrap();
    match t {
        Translation::OneLetter => return c.to_string(),
        Translation::ThreeLetter => return three_letter_code.get(&c).unwrap().to_string(),
    }
}


// tests to ensure translation is working
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
}

// print a pretty sequence, 72 bases per line, plus base numbering
// s: sequence
// t: sequence type (DNA, Protein1 or Protein3)
fn print_seq(s: &str, t: SeqType) {
    let linelen = 72;
    let divisor;
    match t {
        // if we're printing a protein, count amino acids, not bases, 
        // so divide by 3
        SeqType::DNA | SeqType::Protein1 => divisor = 1,
        SeqType::Protein3 => divisor = 3,
    }
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
    let n = 25; // number of Fibonacci numbers to compute
    let f = fibonacci(n); 
    println!("The first {} Fibonacci numbers:\n{:?}", n, f);
    fibonacci_plot(f)?;

    let filename = "sequence.fasta";
    println!("Reading FASTA record from file '{}'...", filename);
    let reader = fasta::Reader::from_file(filename)?;

    let mut nb_reads = 0;
    let mut nb_bases = 0;
    let mut s = String::new();

    for record in reader.records() {
        let record = record.expect("Error during FASTA record parsing");
        if record.id() == "NC_005816.1" {
            println!("Sequence ID: {}", record.id());
            println!("Sequence description:\n{}", record.desc().unwrap());
            // per https://github.com/jperkel/example_notebook/blob/master/My_sample_notebook.ipynb, 
            // there is a hypothetical protein at residues 3485-3857
            let myseq = &record.seq()[3485..3857];
            // myseq is a vector of ASCII codes; convert to a string
            s = String::from_utf8(myseq.to_vec())?;
        }
        nb_reads += 1;
        nb_bases += record.seq().len();
    }
    println!("Total number of reads: {}", nb_reads);
    println!("Total number of bases: {}\n", nb_bases);

    println!("Sequence of my gene [3485..3857]:");
    print_seq(&s, SeqType::DNA);
    println!("Length: {}\n", s.len());

    println!("Translation of my gene:\n");
    assert!(s.len() % 3 == 0, "Sequence length is not a multiple of 3!");
    let mut peptide1 = String::new();    
    let mut peptide3 = String::new();    
    let n_codons = s.len()/3;
    for i in 0..n_codons {
        let codon = &s[i*3..(i*3)+3]; // take a 3-base slice of the sequence
        peptide1.push_str(&translate(&codon, Translation::OneLetter)); // translate and add to the string
        peptide3.push_str(&translate(&codon, Translation::ThreeLetter)); // translate and add to the string
    }
    println!("Single-letter code:");
    print_seq(&peptide1, SeqType::Protein1);
    println!("\nTriple-letter code:");
    print_seq(&peptide3, SeqType::Protein3);
    println!("Length: {}\n", n_codons);

    Ok(()) 
}
