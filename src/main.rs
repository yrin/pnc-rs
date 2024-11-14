/*******************************************************************************
    pnc-rs : Copyright (c) 2024 Yrin Eldfjell : GPLv3

    Multi-threaded Rust-implementation of the PNC (Parallel Neighbourhood 
    Correlation) algorithm described in my M.Sc. thesis. 
    See https://doi.org/10.1371/journal.pcbi.1000063 for a description of
    the original algorithm and idea of using Neighborhood Correlation to
    identify homologous proteins.

    EXPERIMENTAL.

    NOTE: This was my first actual Rust-program ever, so please have patience 
    with non-idiomatic patterns and design choices. 
        The next obvious step is to write a distributed multi-node 
    implementation, with the computations and sharded hashmaps being spread out
    over several nodes and communicating with e.g. TCP. I don't know if it's
    easier to use a library for this or write the networking logic directly.
        It is not my intention to provide a fully capable implementation
    with this package; it may however work for some use cases as-is.

    NOTE: The output between two different runs with identical input may differ
    slightly for numerical reasons, as the order of the computation is non-
    deterministic.

    pnc-rs is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 3 of the License.

    pnc-rs is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program, see the LICENSES file.
    If not, see <https://www.gnu.org/licenses/>.
*******************************************************************************/
use crossbeam_channel::{Receiver, Sender};
use crossbeam_channel::bounded;
use rayon::prelude::*;
use std::collections::{HashMap};
use std::env;
use std::error::Error;
use std::fs::File;
use std::io::BufReader;
use std::io::BufRead;
use std::path::Path;
use std::iter::zip;
use std::str::FromStr;
use std::thread;
use std::time::Instant;

macro_rules! info {
    () => {
        eprint!("[INFO] \n")
    };
    ($PROGRAM_START_TIME:ident,$($arg:tt)*) => {{
        eprint!("[{:5}] ", $PROGRAM_START_TIME.elapsed().as_secs());
        eprint!("[INFO] ");
        eprintln!($($arg)*);
    }};
}
/*macro_rules! warn {
    () => {
        eprint!("[WARN] \n")
    };
    ($($arg:tt)*) => {{
        eprint!("[WARN] ");
        eprintln!($($arg)*);
    }};
}*/
macro_rules! error {
    () => {
        eprint!("[ERROR] \n")
    };
    ($($arg:tt)*) => {{
        eprint!("[ERROR] ");
        eprintln!($($arg)*);
    }};
}

type AccessionId = u32;
type AccessionStr = Vec<u8>;
type TargetId = AccessionId;
type QueryId = AccessionId;
type QueryIdPair = (QueryId, QueryId);
type CalcFloat = f32;
type QuerySumTerm = CalcFloat;
type QuerySquareSumTerm = CalcFloat;
type QuerySums = Vec<QuerySumTerm>;
type QuerySquareSums = Vec<QuerySquareSumTerm>;
type QueryScore = (QueryId, CalcFloat);
type QueryScoreVec = Vec<QueryScore>;
type NCScoreVec = Vec<(QueryIdPair, CalcFloat)>;
type CrossTermValue = CalcFloat;
type CrossTermsVec = Vec<(QueryIdPair, CrossTermValue)>;
type CrossTermsMap = HashMap<QueryIdPair, CrossTermValue>;

struct IntermediateComputeResult {
    sum: QuerySums,
    sum_sq: QuerySquareSums,
    cross_terms: Vec<CrossTermsVec>
}

/* TODO cmd line arg: */ 
// + TODO:  create Config struct for ips, settings, ... 
const NC_THRESHOLD: CalcFloat = 0.05;
const INCLUDE_NC_SCORE_VS_SELF: bool = false;
const CROSS_TERM_BUFFER_PER_THREAD: usize = 8_000_000;

struct Alignments {
    max_query_id: AccessionId,
    max_target_id: AccessionId,
    query_acc_lookup: HashMap<AccessionStr, QueryId>,
    target_acc_lookup: HashMap<AccessionStr, TargetId>,
    alignment_scores: Vec::<QueryScoreVec> /* Indexed by TargetId. */
}

fn parse_alignments(file: &Path) -> Alignments {
    /* Parse alignment scores in blast-tab format:
     * QUERY_ACC<whitespace>REF_ACC<whitespace>SCORE<newline>
     */
    let mut alignments: Alignments = Alignments {
        max_query_id: 0,
        max_target_id: 0,
        query_acc_lookup: HashMap::new(),
        target_acc_lookup: HashMap::new(),
        alignment_scores: Vec::<QueryScoreVec>::new()
    };

    let f = File::open(file).expect("Couldn't open alignments file.");
    let f = BufReader::new(f);
    for line in f.lines() {
        let line = line.expect("Couldn't read line.");
        let mut parts = line.split_ascii_whitespace();
        let query_id: QueryId;
        let target_id: TargetId;
        let query_acc_str = parts.next().expect("No query acc").as_bytes().to_vec();
        if !alignments.query_acc_lookup.contains_key(&query_acc_str) {
            alignments.query_acc_lookup.insert(query_acc_str, alignments.max_query_id);
            if alignments.max_query_id == QueryId::MAX {
                panic!("Ran out of accession IDs.");
            }
            query_id = alignments.max_query_id;
            alignments.max_query_id += 1;
        } else {
            query_id = *alignments.query_acc_lookup.get(&query_acc_str)
                .expect("Query accession not found.");
        }
        let target_acc_str = parts.next().expect("No target acc").as_bytes().to_vec();
        if !alignments.target_acc_lookup.contains_key(&target_acc_str) {
            alignments.target_acc_lookup.insert(target_acc_str, alignments.max_target_id);
            if alignments.max_target_id == TargetId::MAX {
                panic!("Ran out of accession IDs.");
            }
            target_id = alignments.max_target_id;
            alignments.max_target_id += 1;
        } else {
            target_id = *alignments.target_acc_lookup.get(&target_acc_str)
                .expect("Target accession not found.");
        }
        let score = parts.next().expect("No bitscore");
        let score: CalcFloat = CalcFloat::from_str(score)
            .expect("Score float parse error");
        let qs: QueryScore = (query_id, score);
        if target_id as usize == alignments.alignment_scores.len() {
            alignments.alignment_scores.push(Vec::new());
        }
        alignments.alignment_scores[target_id as usize].push(qs)
    }
    alignments
}


fn run(path: &Path) -> Result<(), Box<dyn Error>> { 
    let stime: Instant = Instant::now(); /* Program start time. */
    info!(stime, "START. Processing alignments {} from ", path.display());
    let num_cores: usize = thread::available_parallelism().unwrap().get();
    let num_calculation_threads: usize = usize::max(num_cores, 2);
    let num_cross_term_bins: u32 = num_cores.try_into().unwrap();
    
    info!(stime, "Number of calculation threads: {}", num_calculation_threads);
    info!(stime, "Number of cross-term hashmap storage threads: {}", num_cross_term_bins);
    info!(stime, "Cross-term buffer size (per thread): {}", CROSS_TERM_BUFFER_PER_THREAD);
    info!(stime, "NC-score threshold {}", NC_THRESHOLD);
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_calculation_threads).build_global().unwrap();

    let alignments = parse_alignments(path);
    info!(stime, "Finished parsing alignments.");
    info!(stime, "Number of query accessions found: {}", alignments.max_query_id + 1);
    info!(stime, "Number of target accessions found: {}", alignments.max_target_id + 1);
    let targets_per_chunk: usize = (((alignments.max_target_id as usize) + 1) / 
                                    num_calculation_threads) + 1;
    info!(stime, "Number of targets per cross-term chunk: {}", targets_per_chunk);
    let n = (alignments.max_query_id + 1) as usize;

    /* Setup the channels that take the intermediate compute results from the 
     * compute threads to the storage threads. */
    let (tx_sum, rx_sum): (Sender<QuerySums>, 
                   Receiver<QuerySums>) = bounded(1);
    let (tx_sum_sq, rx_sum_sq): (Sender<QuerySquareSums>, 
                   Receiver<QuerySquareSums>) = bounded(1);
    let (tx_ct_pool, rx_ct_pool): (Vec<Sender<CrossTermsVec>>,
                   Vec<Receiver<CrossTermsVec>>) = (0..num_cross_term_bins)
                        .map(|_| bounded(1))
                        .into_iter()
                        .collect::<Vec<(Sender<CrossTermsVec>, 
                                        Receiver<CrossTermsVec>)>>()
                        .into_iter()
                        .unzip();

    /* Data structures for storing the merged results. */
    let mut sum: QuerySums = vec![0.0; n];
    let mut sum_sq: QuerySquareSums = vec![0.0; n];
    let cross_terms_bins: Vec<CrossTermsMap> = 
        (0..num_cross_term_bins)
        .map(|_| CrossTermsMap::new())
        .collect();

    /* Spawn the merger threads that store intermediate results sent over the
     * channels. */
    let merger_thread_sum: thread::JoinHandle<_> = thread::spawn(move || {
        for sum_chunk in rx_sum {
            for (i, s) in sum_chunk.into_iter().enumerate() {
                sum[i] += s;
            }
        }
        sum
    });
    let merger_thread_sum_sq: thread::JoinHandle<_> = thread::spawn(move || {
        for sum_sq_chunk in rx_sum_sq {
            for (i, s) in sum_sq_chunk.into_iter().enumerate() {
                sum_sq[i] += s;
            }
        }
        sum_sq
    });
    let merger_thread_pool_ct: Vec<thread::JoinHandle<_>> = zip(cross_terms_bins, rx_ct_pool)
        .enumerate()
        .map(move |(bin_idx, (mut cross_terms, rx_ct))| {
            thread::spawn(move || {
                let mut iters_to_next_logging = 100;
                for cross_terms_chunk in rx_ct {
                    if iters_to_next_logging == 0 {
                        info!(stime, "Merging an intermediate result for bin {:3}, \
                              number of cross-terms in res: {:8}, \
                              total cross-terms in bin before: {}", 
                                 bin_idx, cross_terms_chunk.len(), cross_terms.len());
                        iters_to_next_logging = 100;
                    }
                    iters_to_next_logging -= 1;
                    for (key, s) in cross_terms_chunk.iter() {
                        cross_terms.entry(*key).and_modify(|e| *e += *s).or_insert(*s);
                    }
                }
                cross_terms
            })
        })
        .collect();

    /* Start the compute threads that process the raw alignment scores 
     * grouped by target. */
    alignments.alignment_scores
        .par_chunks(targets_per_chunk)
        .enumerate()
        .for_each_with((tx_sum, tx_sum_sq, tx_ct_pool), 
                       |(tx_local_sum, tx_local_sum_sq, tx_local_ct_pool), (_i_targets, targets)| {
            let flush_cross_terms = |cross_terms: Vec<CrossTermsVec>, 
                                     tx_ct_pool: &Vec<Sender<CrossTermsVec>>| 
                                         -> Vec<CrossTermsVec> {
                    cross_terms
                        .into_iter()
                        .enumerate()
                        .for_each(|(i, cross_terms)| {
                            (*tx_ct_pool)[i].send(cross_terms)
                                .expect("Internal thread error");
                        });
                (0..num_cross_term_bins)
                    .map(|_| CrossTermsVec::new())
                    .collect()
            };
            let mut res = IntermediateComputeResult {
                sum: vec![0.0; n],
                sum_sq: vec![0.0; n],
                cross_terms: (0..num_cross_term_bins)
                    .map(|_| CrossTermsVec::new())
                    .collect()
            };
            /* Loop over target ids. */
            let mut insert_counter: usize = 0;
            for query_scores in targets {
                let q_len = query_scores.len();
                for (i, &(query_id_x, score_x)) in query_scores.iter().enumerate() {
                    res.sum[(query_id_x) as usize] += score_x;
                    res.sum_sq[(query_id_x) as usize] += (score_x).powi(2);
                    for &(query_id_y, score_y) in query_scores[i..q_len].iter() {
                        if !INCLUDE_NC_SCORE_VS_SELF && (query_id_x == query_id_y) {
                            continue
                        }
                        let key: QueryIdPair = (query_id_x, query_id_y);
                        let score_prod = (score_x) * (score_y);
                        let ct_bin: usize = ((query_id_x) % num_cross_term_bins) as usize;
                        res.cross_terms[ct_bin].push((key, score_prod));
                        insert_counter += 1;
                        if insert_counter >= CROSS_TERM_BUFFER_PER_THREAD {
                            /*info!(stime,"[TARGET BIN {}] insert counter {}, q_len {}, \
                                         sum ct {:?}", 
                                  _i_targets, 
                                  insert_counter, 
                                  q_len, 
                                  res.cross_terms.iter().map(
                                      |ct| ct.len()).collect::<Vec<usize>>());*/
                            res.cross_terms = flush_cross_terms(res.cross_terms, &tx_local_ct_pool);
                            insert_counter = 0;
                        }
                    }
                }
            }
            tx_local_sum.send(res.sum).expect("Internal thread error");
            tx_local_sum_sq.send(res.sum_sq).expect("Internal thread error");
            /*info!(stime,"[TARGET BIN {}] insert counter {} FINAL, \
                    sum ct {:?}", _i_targets, insert_counter, 
                    res.cross_terms.iter().map(|ct| ct.len()).collect::<Vec<usize>>());*/
            flush_cross_terms(res.cross_terms, &tx_local_ct_pool);
        });
    let sum = merger_thread_sum.join().expect("Internal thread error");
    let sum_sq = merger_thread_sum_sq.join().expect("Internal thread error");
    let cross_terms_bins: Vec<CrossTermsMap> = merger_thread_pool_ct.into_iter().map(
        |t| t.join().expect("Internal thread error")).collect();
    info!(stime, "Finished cross-term calculations, now calculating NC-scores.");

    /* Finally, compute the NC scores. */
    let n_float: CalcFloat = n as CalcFloat;
    let nc_scores: NCScoreVec = cross_terms_bins
//        .par_iter() /* faster, but also uses more memory */
        .into_iter()
        .flatten()
        .map(|((x_u32, y_u32), sum_xy)| {
            let x = (x_u32) as usize;
            let y = (y_u32) as usize;
            let n = n_float;
            let avg_x  = sum[x] / n;
            let avg_y  = sum[y] / n;
            let r_xy_numerator = sum_xy - n*avg_x*avg_y;
            let mut root_term_xx: CalcFloat = sum_sq[x] - n*avg_x.powi(2);
            let mut root_term_yy: CalcFloat = sum_sq[y] - n*avg_y.powi(2);
            if root_term_xx < 0.0 { root_term_xx = 0.0 }
            if root_term_yy < 0.0 { root_term_yy = 0.0 }
            let r_xy_denominator = root_term_xx.sqrt() * root_term_yy.sqrt();
            let r_xy = r_xy_numerator / r_xy_denominator;
            ((x_u32, y_u32), r_xy)
        })
        .filter(|((_, _), r_xy)| 
                *r_xy >= NC_THRESHOLD)
        .collect();
    info!(stime, "Finished calculating NC-scores, now printing results.");

    /* Create a lookup table for the query accessions. */
    let mut accession_vec: Vec<(AccessionStr, QueryId)> = 
        alignments.query_acc_lookup.into_iter().collect();
    accession_vec.sort_by(|a, b| a.1.cmp(&b.1));
    let accession_names: Vec<&AccessionStr> = accession_vec.iter()
        .map(|(key, _)| key)
        .collect();

    /* Print the NC scores. */
    for &((x_id, y_id), r_xy) in nc_scores.iter() {
        let str_id_x = accession_names.get(x_id as usize).expect("Fail!");
        let str_id_y = accession_names.get(y_id as usize).expect("Fail!");
        println!("{:}\t{:}\t{:.3}", 
                 String::from_utf8_lossy(&str_id_x), 
                 String::from_utf8_lossy(&str_id_y), 
                 r_xy);
    }

    info!(stime, "END. Program finished.");
    Ok(())
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        error!("First arg must be valid blast-tab alignment file.");
        panic!()
    }
    let alignments_filename = &args[1];
    run(&Path::new(&alignments_filename)).expect("pnc-rs failed.");
}
