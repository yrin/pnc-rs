# pnc-rs
Fast multithreaded Rust-implementation of PNC (Parallel Neighbourhood Correlation).

Based on [rust-py](https://github.com/yrin/pnc-py). The PNC algorithm is described in my [M.Sc. thesis](https://kurser.math.su.se/pluginfile.php/105616/mod_folder/content/0/2023/Yrin_Eldfjell_MSc_datalogi_2023.pdf). This repo provides a functional (experimental) multithreaded Rust-implementation of the PNC variant of the NC (Neighbourhood Correlation) algorithm.

It's time complexity is better than [snc](https://github.com/arvestad/snc/) and [snc-cpp](https://github.com/arvestad/fast-neighborhood-correlation) by a factor of `log(m)`, where `m` is the number of alignments to the reference database per query sequence (assuming an equal number of alignments per query). See my thesis for details.


## Performance

Aligner settings (Diamond v. 2.1.3):
`diamond blastp --max-target-seqs 100 --sensitive --min-score 30`

CPU: 8 cores, 16 threads.

| Sample                           |   snc-py        |   snc-cpp |   pnc-py    |    pnc-rs | 
| :------------------------------- | --------------: | --------: | ----------: | --------: |
| S. cerevisiae (6k proteins)      |        6        |         1 |          0  |         0 |
| HSA+MMU (38k proteins)           |      701        |        69 |         32  |         3 |
| 13_proteomes (411k proteins)     |    36540 [est.] |      8864 |        880  |        76 |
| UniRef50 [subsampled] (4M prot.) |     N/A         |       N/A |        N/A  |      1436 |

(_In Seconds. Aligner running time not included. `snc-py` progresses consistently linearly during its main computation step, so the estimated time is a reliable lower bound for the actual time required._)

Peak memory usage for `pnc-rs`on the 4M proteins dataset was about 58 GB. The dominating term is the 4+4+4 bytes (2x query id u32's and a float32) used for storing the cross-terms (query_x_score * query_y_score) needed to complete the final NC calculation step.

"13 proteomes" refers to: _E. coli, S. cerevisiae, M. musculus, D. melanogaster, B. subtilis, R. norvegicus, H. sapiens, A. thaliana, D. rerio, C. elegans, D. discoideum, B. taurus,_ and _O. sativa_; downloaded from [UNIPROT](https://www.uniprot.org/proteomes?query=proteome_type%3A1) (2024).

The outputs of all 4 tools were sorted and compared (for HSA+MMU). Each entry deviated from the mean by at most 0.00075, which is reasonable due to the number of significant digits and usage of 32-bit floats during the calculations in both PNC versions.
