use clap::Parser;
use log::info;
use rust_htslib::bcf::{Format, Writer};
use rust_htslib::bcf::header::Header;
use rust_htslib::bcf::record::GenotypeAllele;
use polars::prelude::*;
use std::str;
use std::collections::HashMap;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    #[arg(long, default_value = "test_variants.parquet")]
    parquet_path: String,

    #[arg(long, default_value = "samples.txt")]
    samples_path: String,

    #[arg(long, default_value = "out.bcf")]
    output_path: String,

    #[arg(long, default_value = "2")]
    n_threads: usize
}

fn main() {
    env_logger::init();
    info!("Starting up");

    let args = Args::parse();
    let parquet_path = args.parquet_path;
    let output_path = args.output_path;
    let n_threads = args.n_threads;
    let mut file = std::fs::File::open(&parquet_path).unwrap();

    let mut df = ParquetReader::new(&mut file).finish().unwrap();

    df = df.select(["chrom", "pos", "ref", "alt", "eid"]).unwrap();
    let unique_vars = df.select(["chrom", "pos", "ref", "alt"]).
        unwrap().
        unique_stable(None, UniqueKeepStrategy::First, None).
        unwrap();


    let samples_path = args.samples_path;
    let samples: Vec<String> = read_samples(&samples_path);
    let n_samples = samples.len();
    info!("Read {} samples from {}", n_samples, samples_path);

    let mut header = Header::new();

    for sample in &samples {
        header.push_sample(sample.as_bytes());    
    }

    write_contigs(&mut header);

    let header_gt_line = r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#;
    header.push_record(header_gt_line.as_bytes());

    let mut vcf = Writer::from_path(output_path, &header, false, Format::Bcf).unwrap();
    let _ = vcf.set_threads(n_threads);

    let chrom_vec = unique_vars.column("chrom").unwrap().str().unwrap();
    let pos_vec = unique_vars.column("pos").unwrap().i32().unwrap();
    let ref_vec = unique_vars.column("ref").unwrap().str().unwrap();
    let alt_vec = unique_vars.column("alt").unwrap().str().unwrap();

    info!("Now writing to bcf");
    for i in 0..unique_vars.height() {
        let chrom = chrom_vec.get(i).unwrap();
        let pos = pos_vec.get(i).unwrap() as i64;
        let ref_ = ref_vec.get(i).unwrap();
        let alt = alt_vec.get(i).unwrap();

        let eids = df.clone().lazy().filter(
            col("chrom").eq(lit(chrom))
            .and(col("pos").eq(lit(pos)))
            .and(col("ref").eq(lit(ref_)))
            .and(col("alt").eq(lit(alt)))
            ).collect().
            unwrap().
            column("eid").
            unwrap().
            i32().
            unwrap().
            into_iter().
            map(|s| s.unwrap().to_string()).
            collect::<Vec<String>>();

        let eids_lookup = create_lookup(&eids);

        let mut record = vcf.empty_record();
        // info!("Writing variant: {}:{}:{}:{} with {} samples", chrom, pos, ref_, alt, eids.len());
        let rid = vcf.header().name2rid(chrom.as_bytes()).unwrap();
        record.set_rid(Some(rid));
        record.set_pos(pos - 1);
        let alleles: &[&[u8]] = &[ref_.as_bytes(), alt.as_bytes()];
        record.set_alleles(alleles).expect("Failed to set alleles");

        let hom_ref = &[GenotypeAllele::Phased(0), GenotypeAllele::Phased(0)];
        let het = &[GenotypeAllele::Unphased(0), GenotypeAllele::Phased(1)];

        let mut genotypes: Vec<&[GenotypeAllele]> = Vec::new();

        // for j in 0..n_samples {
        for sample in &samples {
            // if eids.contains(&samples[j]) {
            if eids_lookup.contains_key(sample.as_bytes()) {
                genotypes.push(het);
            } else {
                genotypes.push(hom_ref);
            }
        }

        record.push_genotypes(&genotypes.concat()).unwrap();

        vcf.write(&record).unwrap();
    }

    info!("Done writing variants to bcf");

}

fn write_contigs(header: &mut Header) {
    let contigs = [
        "chr1",
        "chr2",
        "chr3",
        "chr4",
        "chr5",
        "chr6",
        "chr7",
        "chr8",
        "chr9",
        "chr10",
        "chr11",
        "chr12",
        "chr13",
        "chr14",
        "chr15",
        "chr16",
        "chr17",
        "chr18",
        "chr19",
        "chr20",
        "chr21",
        "chr22",
        "chrX",
        "chrY",
        "chrM"
    ];

    let lengths = [
        248956422,
        242193529,
        198295559,
        190214555,
        181538259,
        170805979,
        159345973,
        145138636,
        138394717,
        133797422,
        135086622,
        133275309,
        114364328,
        107043718,
        101991189,
        90338345,
        83257441,
        80373285,
        58617616,
        64444167,
        46709983,
        50818468,
        156040895,
        57227415,
        16569
    ];

    for i in 0..contigs.len() {
        let contig_str = format!(r#"##contig=<ID={},length={}>"#, contigs[i], lengths[i]);
        header.push_record(contig_str.as_bytes());
    }
}

fn read_samples(samples_path: &str) -> Vec<String> {
    let samples = std::fs::read_to_string(samples_path).unwrap();
    let samples: Vec<String> = samples.lines().map(|s| s.to_string()).collect();

    samples
}

fn create_lookup(samples: &Vec<String>) -> HashMap<&[u8], ()> {

    let mut sample_lookup: HashMap<&[u8], ()> = HashMap::new();

    for sample in samples {
        sample_lookup.insert(sample.as_bytes(), ());
    }

    sample_lookup
}
