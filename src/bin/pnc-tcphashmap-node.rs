/*******************************************************************************
    pnc-rs : Copyright (c) 2024 Yrin Eldfjell : GPLv3

    File pnc-tcphashmap-node.rs:
        TCP HashMap node server for distributed computation.

    CAUTION: DO NOT CONNECT THIS SERVER TO THE INTERNET WITHOUT A FIREWALL.

    EXPERIMENTAL.

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
use bincode;
use std;
use std::collections::{HashMap};
use std::env;
use std::net::{TcpListener, TcpStream};
use std::time::Instant;

use libpnc::network::*;

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

type AccessionId = u32;
type CalcFloat = f32;
type QueryId = AccessionId;
type QueryIdPair = (QueryId, QueryId);
type CrossTermValue = CalcFloat;
type CrossTermsMap = HashMap<QueryIdPair, CrossTermValue>;
type CrossTermsVec = Vec<(QueryIdPair, CrossTermValue)>;

const NUM_HASHMAP_BINS: u32 = 32;

fn handle_connection(mut stream: TcpStream,
                     stime: &Instant,
                     cross_terms_bins: &mut Vec<CrossTermsMap>) {
    loop {
        let (op, payload) = read_packet(&mut stream);
        match op {
            TCP_HASHMAP_OP_STORE => {
                let cross_terms_chunk: CrossTermsVec;
                cross_terms_chunk = bincode::deserialize(&payload)
                    .expect("Failed to deserialize payload.");
                for (key, s) in cross_terms_chunk.iter() {
                    let (query_x_id, _) = *key;
                    let ct_bin: usize = (query_x_id % NUM_HASHMAP_BINS) as usize;
                    cross_terms_bins[ct_bin].entry(*key)
                        .and_modify(|e| *e += *s).or_insert(*s);
                }
            },
            TCP_HASHMAP_OP_FETCH => {
                match cross_terms_bins.pop() {
                    Some(cross_terms) => {
                        let payload: Vec<u8> = bincode::serialize(&cross_terms)
                            .expect("Internal error: failed to serialize cross-terms");
                        send_packet(&mut stream, TCP_HASHMAP_RESPONSE_MAP_NOT_EMPTY, &payload);
                    },
                    None => send_packet(&mut stream, TCP_HASHMAP_RESPONSE_MAP_EMPTY, &vec![])
                }
            },
            TCP_HASHMAP_OP_TERM_CONN => {
                break;
            }
            TCP_HASHMAP_OP_SHUTDOWN => {
                info!(stime, "Terminating process by external request");
                std::process::exit(0);
            },
            invalid => panic!("Invalid OP code: {}", invalid),
        };
    }
}

fn main() {
    let stime: Instant = Instant::now(); /* Program start time. */
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("First argument must be ip:port.");
        panic!()
    }
    let mut cross_terms_bins = (0..NUM_HASHMAP_BINS)
        .map(|_| <CrossTermsMap>::default())
        .collect();
    let ip_addr = &args[1];
    let listener = TcpListener::bind(ip_addr).unwrap();
    info!(stime, "STARTING pnc-tcphashmap-node ");
    for stream in listener.incoming() {
        let stream = stream.unwrap();
        info!(stime, "New connection for {}", ip_addr);
        stream.set_read_timeout(None).expect("set_read_timeout call failed");
        handle_connection(stream, &stime, &mut cross_terms_bins);
    }
}
