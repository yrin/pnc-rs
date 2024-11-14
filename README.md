# pnc-rs
PNC (Parallel Neighbourhood Correlation): Fast multithreaded Rust-implementation. 

Based on [rust-py](https://github.com/yrin/pnc-py) as described in my [M.Sc. thesis](https://kurser.math.su.se/pluginfile.php/105616/mod_folder/content/0/2023/Yrin_Eldfjell_MSc_datalogi_2023.pdf), this repo provides a functional (experimental) multithreaded Rust-implementation of the NC (Neighbourhood Correlation) algorithm.

It's faster compared to [snc](https://github.com/arvestad/snc/) and [snc-cpp](https://github.com/arvestad/fast-neighborhood-correlation) by a factor of `log(m)`, where `m` is the number of alignments to the reference database per query sequence (assuming an equal number of alignments per query).


## Performance

Aligner settings:
`diamond blastp --max-target-seqs 100 --sensitive --min-score 30`

| Sample                         |   snc-py        |   snc-cpp |   pnc-py    |    pnc-rs | 
| ------------------------------ | --------------: | --------: | ----------: | --------: |
| S. cerevisiae (6k proteins)    |        6        |         1 |          0  |         0 |
| HSA+MMU (38k proteins)         |      701        |        69 |         32  |         3 |
| 13_proteomes (411k proteins)   |    36540 [est.] |      8864 |        880  |        76 |
| UniRef50 (4M prot.) subsampled |     N/A         |       N/A |        N/A  |      1436 |
(_Seconds_)

"13 proteomes" refer to: E. coli, S. cerevisiae, M. musculus, D. melanogaster, B. subtilis, R. norvegicus, H. sapiens, A. thaliana, D. rerio, C. elegans, D. discoideum, B. taurus, O. sativa; downloaded from UNIPROT (2024).

The outputs of all 4 tools were sorted and compared (for HSA+MMU). Each entry deviated from the mean by at most 0.00075, which is reasonable due to the number of significant digits and usage of 32-bit floats during the calculations in both PNC versions.
