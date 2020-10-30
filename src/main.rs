extern crate bio;
use plotters::prelude::*;
use bio::io::fasta;

// String containing the standard genetic code. Indexing into this array is
// as follows: 
// A = 0, C = 1, G = 2, T = 3
// index = (16 * first_base) + (4 * second_base) + third_base 
// So... AAA = 0, AAC = 1, AAG = 2, ... , TTC = 61, TTG = 62, TTT = 63
const GENETIC_CODE: &str = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CYCLFLF"; 

// ERR_BAD_NT is an error value for an invalid nucleotide 
const ERR_BAD_NT: usize = 99; 


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

    for i in 0..n {
        v[i] = (i as f32, f[i] as f32);
    }

    let root = BitMapBackend::new("fibonacci.png", (640, 480)).into_drawing_area();
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
    println!("Graph written to 'fibonacci.png'.\n");

    Ok(()) 
}

// given an input 'char', return a base equivalent
fn lookup(x: char) -> usize {
    let base;

    match x {
        'A' => base = 0,
        'C' => base = 1,
        'G' => base = 2,
        'T' => base = 3,
        _ => base = ERR_BAD_NT, // unknown base
    }
    return base;
}

// translate a codon into a single-letter amino acid
fn translate(triplet: &str) -> char {
    let mut codon = vec![ERR_BAD_NT; 3];

    for (i,base) in triplet.chars().enumerate() {
        let val = lookup(base);
        if val == ERR_BAD_NT {
            panic!(); // invalid character
        }
        codon[i] = val;
    }

    let index: usize = (codon[0] * 16) + (codon[1] * 4) + codon[2];
    return GENETIC_CODE.chars().nth(index).unwrap();
}



fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Calculating Fibonacci sequence...\n");
    let n = 25; // number of Fibonacci numbers to compute
    let f = fibonacci(n); 
    println!("The first {} Fibonacci numbers:\n{:?}\n", n, f);
    fibonacci_plot(f)?;

    println!("Reading sequence data...");
    let reader = fasta::Reader::from_file("sequence.fasta")?;

    let mut nb_reads = 0;
    let mut nb_bases = 0;
    let mut s = String::new();

    for record in reader.records() {
        let record = record.expect("Error during FASTA record parsing");
        if record.id() == "NC_005816.1" {
            println!("{:?}", record.desc());
            let myseq = &record.seq()[3485..3857];
            s = String::from_utf8(myseq.to_vec())?;
        }
        nb_reads += 1;
        nb_bases += record.seq().len();
    }
    println!("Total number of reads: {}", nb_reads);
    println!("Total number of bases: {}\n", nb_bases);

    println!("Sequence of my gene [3485..3857]:");
    println!("{}", s);
    println!("Length: {}\n", s.len());

    println!("Translation of my gene:");
    assert!(s.len() % 3 == 0, "Sequence is not a multiple of 3 in length!");
    let mut peptide = String::new();    
    for i in 0..s.len()/3 {
        let codon = &s[i*3..(i*3)+3];
        peptide.push(translate(&codon));
    }
    println!("{}", peptide);
    println!("Length: {}\n", peptide.len());

    Ok(()) 
}
