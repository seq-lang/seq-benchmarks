use bio::alphabets::dna;
use bio::data_structures::bwt;
use bio::data_structures::fmindex;
use bio::data_structures::suffix_array::suffix_array;
use bio::io::{fasta, fastq};
use std::env;
use std::fs::File;
use std::io;
use std::time;

/**************************************************************************************************/
/******************************************* Problem 1   ******************************************/
/**************************************************************************************************/
pub struct Index {
    sa: Vec<usize>,
    idx: fmindex::FMDIndex<bwt::BWT, bwt::Less, bwt::Occ>,
}
impl Index {
    pub fn new(seq: &[u8]) -> Self {
        let alphabet = dna::n_alphabet();
        let sa = suffix_array(&seq);
        let bwt = bwt::bwt(&seq, &sa);
        let less = bwt::less(&bwt, &alphabet);
        let occ = bwt::Occ::new(&bwt, 3, &alphabet);
        let fm = fmindex::FMIndex::new(bwt, less, occ);
        let idx = fmindex::FMDIndex::from(fm);
        Index { sa, idx }
    }
}

fn fastmap(rec: &fastq::Record, idx: &Index, f: &mut io::BufWriter<File>) {
    use std::io::*;
    use std::str::from_utf8;

    let mems = &mut Vec::new();
    let mut i = 0;
    while i < rec.seq().len() {
        let (ni, sm) = idx.idx.smems(rec.seq(), i, 17);
        i = ni;
        mems.extend(sm);
    }
    for mem in mems {
        let _ = write!(f, "{}\tEM\t{}\t{}\t{}", rec.id(), mem.1, mem.2, mem.0.size);
        if mem.0.size <= 20 {
            for pos in mem.0.forward().occ(&idx.sa) {
                let _ = write!(f, "\tchr:+{}", pos);
            }
            for pos in mem.0.revcomp().occ(&idx.sa) {
                let _ = write!(f, "\tchr:-{}", pos);
            }
        } else {
            let _ = write!(f, "\t*\n");
        }
        let _ = write!(f, "\n");
    }
}

fn build_index(path: &str) -> (fasta::Record, Index) {
    use bio::io::fasta::FastaRead;

    print!("indexing... ");
    let t = time::Instant::now();
    let mut reader = fasta::Reader::from_file(path).expect("cannot open file");
    let mut record = fasta::Record::new();
    reader.read(&mut record).expect("Failed to parse record");
    assert!(!record.is_empty());
    let mut seq = record.seq().to_vec();
    seq.push(b'$');
    let idx = Index::new(&seq);
    print!("{}\nindex done!", t.elapsed().as_secs_f32());
    (record, idx)
}

/**************************************************************************************************/
/******************************************* Problem 2   ******************************************/
/**************************************************************************************************/

fn align(x: &[u8], y: &[u8], score_only: bool) -> (i32, i32) {
    use bio::alignment::distance::simd;
    use bio::alignment::pairwise;
    if !score_only {
        let score = |a: u8, b: u8| {
            if a == b'N' || b == b'N' {
                0
            } else if a == b {
                -1i32
            } else {
                -2i32
            }
        };
        let mut aligner = pairwise::Aligner::with_capacity(x.len(), y.len(), -4, -2, &score);
        let alignment = aligner.global(x, y);
        (alignment.score, 1)
    } else {
        let score = simd::levenshtein(x, y);
        (-(score as i32), 1)
    }
}

/**************************************************************************************************/
/******************************************* Main        ******************************************/
/**************************************************************************************************/

fn main() {
    use std::io::*;
    let args = env::args().collect::<Vec<String>>();
    if args[1] == "fastmap" {
        println!("indexing...");
        let mut t = time::Instant::now();
        let (record, idx) = build_index(&args[2]);
        println!("{}\naligning... ", t.elapsed().as_secs_f32());
        t = time::Instant::now();
        let reader = fastq::Reader::from_file(&args[3]).expect("cannot open file");
        let mut ff = BufWriter::new(File::create(&args[4]).unwrap());
        for r in reader.records() {
            fastmap(&r.unwrap(), &idx, &mut ff);
        }
        println!("{}\ndone!", t.elapsed().as_secs_f32());
    } else if args[1] == "align" {
        let times = &mut Vec::<f32>::new();
        for score_only in 0..2 {
            for i in 0..3 {
                let t = time::Instant::now();
                let src = BufReader::new(File::open(&args[2]).unwrap());
                let target = BufReader::new(File::open(&args[3]).unwrap());
                let mut total = 0;
                let mut score = 0;
                for (x, y) in src.lines().zip(target.lines()) {
                    let (s, t) = align(
                        x.unwrap().as_bytes(),
                        &y.unwrap().as_bytes(),
                        score_only != 0,
                    );
                    total += t;
                    score += s;
                }
                times.push(t.elapsed().as_secs_f32());
                // eprintln!("\r- {} {} {} {}", i, total, score, times.last().unwrap());
            }
            let mean = times.iter().sum::<f32>() / (times.len() as f32);
            let mut var = 0.0;
            for j in times.iter() {
                var += (j - mean) * (j - mean);
            }
            eprintln!(
                "[sw-time] rust {} {} {}",
                score_only,
                mean,
                (var / (times.len() as f32)).sqrt()
            );
        }
    } else {
        panic!("bad arg");
    }
}

// seq: 1:50 / 1:04 [5+58] [7g/3g]
// rust 5:30
