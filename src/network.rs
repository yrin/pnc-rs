/*******************************************************************************
    pnc-rs : Copyright (c) 2024 Yrin Eldfjell : GPLv3

    File network.rs:
        Helper functions for TCP HashMap server and client.

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
use byteorder::{ByteOrder, NetworkEndian};
use std::io::Read;
use std::io::Write;
use std::net::{TcpStream};

pub const TCP_HASHMAP_OP_STORE: u8 = 1;
pub const TCP_HASHMAP_OP_FETCH: u8 = 2;
pub const TCP_HASHMAP_OP_TERM_CONN: u8 = 3;
pub const TCP_HASHMAP_OP_SHUTDOWN: u8 = 4;
pub const TCP_HASHMAP_RESPONSE_MAP_EMPTY: u8 = 1;
pub const TCP_HASHMAP_RESPONSE_MAP_NOT_EMPTY: u8 = 2;

/// Read a (op, payload) packet from stream.
pub fn read_packet(stream: &mut TcpStream) -> (u8, Vec<u8>) {
    let mut header: [u8; 9] = [0; 9];
    stream.read_exact(&mut header).unwrap();
    let payload_size: usize = NetworkEndian::read_u64(&header[1..9])
        .try_into().expect("Failed to read payload size.");
    let mut payload = vec![0_u8; payload_size];
    stream.read_exact(&mut payload)
        .expect("Network I/O error when reading from TcpStream.");
    let op = header[0];
    (op, payload)
}

/// Write a (op, payload) packet to stream.
pub fn send_packet(stream: &mut TcpStream, op: u8, data: &Vec<u8>) {
    let mut header: [u8; 9] = [0; 9];
    header[0] = op;
    NetworkEndian::write_u64(&mut header[1..9], data.len() as u64);
    stream.write_all(&header)
        .expect("Network I/O error when sending to TcpStream.");
    stream.write_all(data)
        .expect("Network I/O error when sending to TcpStream.");
}
