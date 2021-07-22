
use std::io::prelude::*;
use std::io::{BufReader,BufWriter};
use std::fs::File;

use chrono::{DateTime, Utc};
use intervaltree::*;
use rulinalg::utils::argmax;

struct Problem {
    delta: u32,
    isoforms: Vec<Vec<(u64, u64)>>,
    reads: Vec<Vec<(u64, u64)>>
}

fn load_problem(filename: &String) -> Problem {
    //set up file IO
    let fp: File = File::open(filename).unwrap();
    let mut reader: BufReader<File> = BufReader::new(fp);

    //first line has isoform count and delta
    let mut line = String::new();
    reader.read_line(&mut line).unwrap();
    let split_line: Vec<&str> = line[..line.len()-1].split(" ").collect();
    let num_isoforms: u32 = split_line[0].parse().unwrap();
    let delta: u32 = split_line[1].parse().unwrap();
    
    //parse all the isoform lines
    let mut isoforms: Vec<Vec<(u64, u64)>> = vec![];
    for _ in 0..num_isoforms {
        let mut isoform_exons: Vec<(u64, u64)> = vec![];
        line.clear();
        reader.read_line(&mut line).unwrap();
        let split_ranges = line[..line.len()-1].split(",");
        for s in split_ranges {
            let split_range: Vec<&str> = s.split("-").collect();
            let start: u64 = split_range[0].parse().unwrap();
            let end: u64 = split_range[1].parse::<u64>().unwrap();
            isoform_exons.push((start, end));
        }
        isoforms.push(isoform_exons);
    }

    //get the number of reads
    line.clear();
    reader.read_line(&mut line).unwrap();
    let num_reads: u32 = line[..line.len()-1].parse().unwrap();

    //now parse the reads
    let mut reads: Vec<Vec<(u64, u64)>> = vec![];
    for _ in 0..num_reads {
        let mut read_exons: Vec<(u64, u64)> = vec![];
        line.clear();
        reader.read_line(&mut line).unwrap();
        let split_ranges = line[..line.len()-1].split(",");
        for s in split_ranges {
            let split_range: Vec<&str> = s.split("-").collect();
            let start: u64 = split_range[0].parse().unwrap();
            let end: u64 = split_range[1].parse::<u64>().unwrap();
            read_exons.push((start, end));
        }
        reads.push(read_exons);
    }

    Problem{
        delta, isoforms, reads
    }
}

fn get_introns(exons: &Vec<(u64, u64)>) -> Vec<(u64, u64)> {
    let mut ret = vec![];
    for vec_index in 0..(exons.len()-1) {
        ret.push((exons[vec_index].1, exons[vec_index+1].0))
    }
    ret
}

fn solve_problem_ivtree(problem: Problem) -> Vec<u64> {
    let isoforms = problem.isoforms;
    let reads = problem.reads;

    let mut exons_regions = vec![];
    let mut intron_regions = vec![];

    for (iso_index, isoform) in isoforms.iter().enumerate() {
        let introns = get_introns(isoform);
        for (start, end) in isoform.iter() {
            exons_regions.push((start.clone()..end.clone(), iso_index as u64));
        }

        for (start, end) in introns.iter() {
            intron_regions.push((start.clone()..end.clone(), iso_index as u64));
        }
    }

    //create an interval tree for all the exonic and intronic regions
    let exon_tree: IntervalTree<u64, u64> = exons_regions.iter().cloned().collect();
    let intron_tree: IntervalTree<u64, u64> = intron_regions.iter().cloned().collect();

    let mut ret: Vec<u64> = vec![];
    let mut scores: Vec<f64>;

    for (read_index, read_exons) in reads.iter().enumerate() {
        println!("[{}]\tRead #{} / {}", Utc::now(), read_index, reads.len());
        scores = vec![0.0; isoforms.len()];
        let read_introns = get_introns(read_exons);
        let total_length = read_exons[read_exons.len()-1].1 - read_exons[0].0;
        let mut exon_length = 0;
        for (exon_start, exon_end) in read_exons {
            exon_length += exon_end - exon_start;
        }
        let intron_length = total_length - exon_length;

        for (exon_start, exon_end) in read_exons {
            let matches = exon_tree.query(*exon_start..*exon_end);
            for m in matches {
                let iso_start = m.range.start;
                let iso_end = m.range.end;
                let iso_index = m.value;

                let max_start = std::cmp::max(exon_start, &iso_start);
                let min_end = std::cmp::min(exon_end, &iso_end);
                scores[iso_index as usize] += 2.0 * (min_end - max_start) as f64 / exon_length as f64;
            }
        }

        for (intron_start, intron_end) in read_introns.iter() {
            let matches = intron_tree.query(*intron_start..*intron_end);
            for m in matches {
                let iso_start = m.range.start;
                let iso_end = m.range.end;
                let iso_index = m.value;

                let max_start = std::cmp::max(intron_start, &iso_start);
                let min_end = std::cmp::min(intron_end, &iso_end);
                scores[iso_index as usize] += (min_end - max_start) as f64 / intron_length as f64;
            }
        }

        let best_match = argmax(&scores).0;
        ret.push(best_match as u64);
    }

    ret
}

fn overlap_ranges(i_ranges: &Vec<(u64, u64)>, r_ranges: &Vec<(u64, u64)>) -> u64 {
    let mut overlap: u64 = 0;

    let mut i_ind: usize = 0;
    let mut r_ind: usize = 0;

    while r_ind < r_ranges.len() && i_ind < i_ranges.len() {
        let (r_start, r_end) = r_ranges[r_ind];
        let (i_start, i_end) = i_ranges[i_ind];
        let max_start = std::cmp::max(r_start, i_start);
        let min_end = std::cmp::min(r_end, i_end);
        if max_start < min_end {
            overlap += min_end - max_start;
        }

        if r_end < i_end {
            r_ind += 1;
        } else {
            i_ind += 1;
        }
    }

    overlap
}

fn solve_problem_list(problem: Problem) -> Vec<u64> {
    let isoforms = problem.isoforms;
    let reads = problem.reads;

    let mut ret: Vec<u64> = vec![];
    let mut scores: Vec<f64>;
    
    for (read_index, read_exons) in reads.iter().enumerate() {
        println!("[{}]\tRead #{} / {}", Utc::now(), read_index, reads.len());
        scores = vec![0.0; isoforms.len()];
        let read_introns = get_introns(read_exons);
        let total_length = read_exons[read_exons.len()-1].1 - read_exons[0].0;
        let mut exon_length = 0;
        for (exon_start, exon_end) in read_exons {
            exon_length += exon_end - exon_start;
        }
        let intron_length = total_length - exon_length;

        for (isoform_index, isoform_exons) in isoforms.iter().enumerate() {
            let isoform_introns = get_introns(isoform_exons);

            let exon_overlap = overlap_ranges(isoform_exons, read_exons);
            let intron_overlap = overlap_ranges(&isoform_introns, &read_introns);

            scores[isoform_index] = 2.0*(exon_overlap as f64 / exon_length as f64) + (intron_overlap as f64 / intron_length as f64);
        }

        let best_match = argmax(&scores).0;
        ret.push(best_match as u64);
    }

    ret
}

fn write_result(results: Vec<u64>, out_fn: String) {
    //set up file IO
    let fp: File = File::create(out_fn).unwrap();
    let mut writer: BufWriter<File> = BufWriter::new(fp);
    for result in results.iter() {
        writer.write((result.to_string()+"\n").as_bytes()).unwrap();
    }
}

fn main() {
    let filename: String = "../../data/70-welcome-approx.txt".to_string();
    let out_fn: String = "../../results/7.txt".to_string();

    //let filename: String = "../../data/80-big-approx.txt".to_string();
    //let out_fn: String = "../../results/8.txt".to_string();

    //let filename: String = "../../data/90-huge-approx.txt".to_string();
    //let out_fn: String = "../../results/9.txt".to_string();

    let problem = load_problem(&filename);
    //let results = solve_problem_list(problem);
    let results = solve_problem_ivtree(problem);
    write_result(results, out_fn);
}
