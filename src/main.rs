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
use std::io;
use std::io::BufReader;
use std::io::BufRead;
use std::io::BufWriter;
use std::io::Write;
use std::path::Path;
use std::iter::zip;
use std::str::FromStr;
use std::thread;
use std::time::Instant;

macro_rules! info {
    () => {
        eprint!("[INFO] \n")
    };
    ($PROGRAM_START_TIME:expr,$($arg:tt)*) => {{
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

struct ComputeResultVec {
    sum: QuerySums,
    sum_sq: QuerySquareSums,
    cross_terms_bins: Vec<CrossTermsVec>
}
impl ComputeResultVec {
    pub fn new(n_queries: usize, n_bins: u32) -> Self {
        ComputeResultVec {
            sum: vec![0.0; n_queries],
            sum_sq: vec![0.0; n_queries],
            cross_terms_bins: (0..n_bins)
                .map(|_| CrossTermsVec::new())
                .collect()
        }
    }
}

struct ComputeResultMap {
    sum: QuerySums,
    sum_sq: QuerySquareSums,
    cross_terms_bins: Vec<CrossTermsMap>
}
impl ComputeResultMap {
    pub fn new(n_queries: usize, n_bins: u32) -> Self {
        ComputeResultMap {
            sum: vec![0.0; n_queries],
            sum_sq: vec![0.0; n_queries],
            cross_terms_bins: (0..n_bins)
                .map(|_| CrossTermsMap::new())
                .collect()
        }
    }
}

struct Alignments {
    max_query_id: AccessionId,
    max_target_id: AccessionId,
    query_acc_lookup: HashMap<AccessionStr, QueryId>,
    target_acc_lookup: HashMap<AccessionStr, TargetId>,
    alignment_scores: Vec::<QueryScoreVec> /* Indexed by TargetId. */
}

struct Config {
    n_queries: usize,
    n_calculation_threads: usize,
    n_memory_threads: u32,
    nc_threshold: CalcFloat,
    include_nc_score_vs_self: bool,
    cross_term_buffer_per_thread: usize,
    targets_per_chunk: usize,
    stime: Instant
}
impl Config {
    pub fn new() -> Self {
        let n_cores: usize = thread::available_parallelism().unwrap().get();
        Config {
            n_queries: 0,
            n_calculation_threads: usize::max(n_cores, 2),
            n_memory_threads: n_cores.try_into().unwrap(),
            nc_threshold: 0.05,
            include_nc_score_vs_self: false,
            cross_term_buffer_per_thread: 8_000_000,
            targets_per_chunk: 0,
            stime: Instant::now() /* Program start time. */
        }
    }
    pub fn update_n_queries(&self, n_queries: usize) -> Self {
        let targets_per_chunk: usize = (n_queries / self.n_calculation_threads) + 1;
        Config {
            n_queries: n_queries,
            targets_per_chunk: targets_per_chunk,
            ..*self
        }
    }
}

/// Parse alignment scores in blast-tab format:
/// QUERY_ACC<whitespace>REF_ACC<whitespace>SCORE<newline>
fn parse_alignments(file: &Path) -> Alignments {
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

/// Start the compute threads that process the raw alignment scores grouped by target.
fn process_alignments(config: &Config, 
                      alignments: &Alignments, 
                      tx_sum: Sender<QuerySums>, 
                      tx_sum_sq: Sender<QuerySquareSums>, 
                      tx_ct_pool: Vec<Sender<CrossTermsVec>>) {
    alignments.alignment_scores
        .par_chunks(config.targets_per_chunk)
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
                (0..config.n_memory_threads)
                    .map(|_| CrossTermsVec::new())
                    .collect()
            };
            let mut res = ComputeResultVec::new(config.n_queries, config.n_memory_threads);
            /* Loop over target ids. */
            let mut insert_counter: usize = 0;
            for query_scores in targets {
                let q_len = query_scores.len();
                for (i, &(query_id_x, score_x)) in query_scores.iter().enumerate() {
                    res.sum[(query_id_x) as usize] += score_x;
                    res.sum_sq[(query_id_x) as usize] += (score_x).powi(2);
                    for &(query_id_y, score_y) in query_scores[i..q_len].iter() {
                        if !config.include_nc_score_vs_self && (query_id_x == query_id_y) {
                            continue
                        }
                        let key: QueryIdPair = (query_id_x, query_id_y);
                        let score_prod = (score_x) * (score_y);
                        let ct_bin: usize = ((query_id_x) % config.n_memory_threads) as usize;
                        res.cross_terms_bins[ct_bin].push((key, score_prod));
                        insert_counter += 1;
                        if insert_counter >= config.cross_term_buffer_per_thread {
                            /*info!(config.stime,"[TARGET BIN {}] insert counter {}, q_len {}, \
                                         sum ct {:?}", 
                                  _i_targets, 
                                  insert_counter, 
                                  q_len, 
                                  res.cross_terms.iter().map(
                                      |ct| ct.len()).collect::<Vec<usize>>());*/
                            res.cross_terms_bins = flush_cross_terms(res.cross_terms_bins, &tx_local_ct_pool);
                            insert_counter = 0;
                        }
                    }
                }
            }
            tx_local_sum.send(res.sum).expect("Internal thread error");
            tx_local_sum_sq.send(res.sum_sq).expect("Internal thread error");
            /*info!(config.stime,"[TARGET BIN {}] insert counter {} FINAL, \
                    sum ct {:?}", _i_targets, insert_counter, 
                    res.cross_terms.iter().map(|ct| ct.len()).collect::<Vec<usize>>());*/
            flush_cross_terms(res.cross_terms_bins, &tx_local_ct_pool);
        });
}

/// Spawn the merger threads that store the intermediate results sent over the channels.
fn spawn_merger_threads(config: &Config, 
                        rx_sum: Receiver<QuerySums>, 
                        rx_sum_sq: Receiver<QuerySquareSums>, 
                        rx_ct_pool: Vec<Receiver<CrossTermsVec>>) 
    -> (thread::JoinHandle<QuerySums>, 
        thread::JoinHandle<QuerySquareSums>, 
        Vec<thread::JoinHandle<CrossTermsMap>>) {
    let mut cr = ComputeResultMap::new(config.n_queries, config.n_memory_threads); 
    let merger_thread_sum: thread::JoinHandle<_> = thread::spawn(move || {
        for sum_chunk in rx_sum {
            for (i, s) in sum_chunk.into_iter().enumerate() {
                cr.sum[i] += s;
            }
        }
        cr.sum
    });
    let merger_thread_sum_sq: thread::JoinHandle<_> = thread::spawn(move || {
        for sum_sq_chunk in rx_sum_sq {
            for (i, s) in sum_sq_chunk.into_iter().enumerate() {
                cr.sum_sq[i] += s;
            }
        }
        cr.sum_sq
    });
    let merger_thread_pool_ct: Vec<thread::JoinHandle<_>> = zip(cr.cross_terms_bins, rx_ct_pool)
        .enumerate()
        .map(move |(_bin_idx, (mut cross_terms, rx_ct))| {
            thread::spawn(move || {
                let mut iters_to_next_logging = 100;
                for cross_terms_chunk in rx_ct {
                    if iters_to_next_logging == 0 {
                        /*info!(config.stime, "Merging an intermediate result for bin {:3}, \
                              number of cross-terms in res: {:8}, \
                              total cross-terms in bin before: {}", 
                                 _bin_idx, cross_terms_chunk.len(), cross_terms.len());*/
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
        (merger_thread_sum, merger_thread_sum_sq, merger_thread_pool_ct)
}

/// Compute the NC scores from the merged intermediate computations.
fn compute_nc_scores(num_queries: usize, 
                     merged_results: ComputeResultMap,
                     nc_threshold: CalcFloat) -> NCScoreVec {
    let n_float: CalcFloat = num_queries as CalcFloat;
    let sum = merged_results.sum;
    let sum_sq = merged_results.sum_sq;
    merged_results.cross_terms_bins
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
                *r_xy >= nc_threshold)
        .collect()
}

/// Print the NC scores to stdout
fn print_nc_scores(nc_scores: NCScoreVec, alignments: Alignments) {

    /* Create a lookup table for the query accessions. */
    let mut accession_vec: Vec<(AccessionStr, QueryId)> = 
        alignments.query_acc_lookup.into_iter().collect();
    accession_vec.sort_by(|a, b| a.1.cmp(&b.1));
    let accession_names: Vec<&AccessionStr> = accession_vec.iter()
        .map(|(key, _)| key)
        .collect();
    const WRITE_BUF_SIZE: usize = 1024 * 1024;
    let stdout = io::stdout().lock();
    let mut wr = BufWriter::with_capacity(WRITE_BUF_SIZE, stdout);

    /* Print the NC scores. */
    for &((x_id, y_id), nc_score) in nc_scores.iter() {
        let str_id_x = accession_names.get(x_id as usize).expect("Internal error 1");
        let str_id_y = accession_names.get(y_id as usize).expect("Internal error 2");
        wr.write(str_id_x).expect("I/O error during writing of results");
        wr.write(b"\t").expect("I/O error during writing of results");
        wr.write(str_id_y).expect("I/O error during writing of results");
        wr.write(b"\t").expect("I/O error during writing of results");
        // TODO: truncate to 3 decimal places:
        wr.write(nc_score.to_string().as_bytes()).expect("I/O error during writing of results");
        wr.write(b"\n").expect("I/O error during writing of results");
    }
}

/// Uses the alignments stored in `path` to generate NC-scores.
fn run(path: &Path) -> Result<(), Box<dyn Error>> { 
    let config = Config::new();
    info!(config.stime, "START. Processing alignments {} from ", path.display());
    info!(config.stime, "Number of calculation threads: {}", config.n_calculation_threads);
    info!(config.stime, "Number of cross-term hashmap storage threads: {}", config.n_memory_threads);
    info!(config.stime, "Cross-term buffer size (per thread): {}", config.cross_term_buffer_per_thread);
    info!(config.stime, "NC-score threshold {}", config.nc_threshold);
    rayon::ThreadPoolBuilder::new()
        .num_threads(config.n_calculation_threads).build_global().unwrap();

    let alignments = parse_alignments(path);
    let config = config.update_n_queries((alignments.max_query_id + 1) as usize);
    info!(config.stime, "Finished parsing alignments.");
    info!(config.stime, "Number of query accessions found: {}", alignments.max_query_id + 1);
    info!(config.stime, "Number of target accessions found: {}", alignments.max_target_id + 1);
    info!(config.stime, "Number of targets per cross-term chunk: {}", config.targets_per_chunk);

    /* Setup the channels that take the intermediate compute results from the 
     * compute threads to the storage threads. */
    let (tx_sum, rx_sum): (Sender<QuerySums>, 
                   Receiver<QuerySums>) = bounded(1);
    let (tx_sum_sq, rx_sum_sq): (Sender<QuerySquareSums>, 
                   Receiver<QuerySquareSums>) = bounded(1);
    let (tx_ct_pool, rx_ct_pool): (Vec<Sender<CrossTermsVec>>,
                   Vec<Receiver<CrossTermsVec>>) = (0..config.n_memory_threads)
                        .map(|_| bounded(1))
                        .into_iter()
                        .collect::<Vec<(Sender<CrossTermsVec>, 
                                        Receiver<CrossTermsVec>)>>()
                        .into_iter()
                        .unzip();
    
    let (merger_thread_sum,
         merger_thread_sum_sq,
         merger_thread_pool_ct): (thread::JoinHandle<QuerySums>, 
                                  thread::JoinHandle<QuerySquareSums>, 
                                  Vec<thread::JoinHandle<CrossTermsMap>>) = spawn_merger_threads(
        &config, rx_sum, rx_sum_sq, rx_ct_pool);
    process_alignments(&config, &alignments, tx_sum, tx_sum_sq, tx_ct_pool);

    let mut merged_results = ComputeResultMap::new(config.n_queries, config.n_memory_threads);
    merged_results.sum = merger_thread_sum.join().expect("Internal thread error");
    merged_results.sum_sq = merger_thread_sum_sq.join().expect("Internal thread error");
    merged_results.cross_terms_bins = merger_thread_pool_ct.into_iter().map(
        |t| t.join().expect("Internal thread error")).collect();
    info!(config.stime, "Finished cross-term calculations, now calculating NC-scores.");

    let nc_scores = compute_nc_scores(config.n_queries, merged_results, config.nc_threshold);
    info!(config.stime, "Finished calculating NC-scores, now printing results.");

    print_nc_scores(nc_scores, alignments);

    info!(config.stime, "END. Program finished.");
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
