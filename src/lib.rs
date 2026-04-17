//! lastz-gxy: Rust reimplementation of lastz.
//!
//! Stage pipeline (mirrors upstream):
//!
//! ```text
//! target ──► pos_table ──┐
//!                        ├─► seed_search ──► diag_hash ──► hsp
//! query  ──► seed stream ┘                                  │
//!                                                           ▼
//! output ◄── tweener ◄── gapped_extend ◄── anchor ◄── chain
//! ```

pub mod anchor;
pub mod chain;
pub mod cli;
pub mod diag_hash;
pub mod dna;
pub mod driver;
pub mod edit_script;
pub mod gapped_extend;
pub mod hsp;
mod hsp_simd;
pub mod output;
pub mod pos_table;
pub mod scoring;
pub mod seed_search;
pub mod seeds;
pub mod sequences;
pub mod tweener;
