# MHG-EVO: _de Novo_ Tree-Based MHG Finder

We are now moving MHG-EVO to conda for an easier installation. 
\
\
MHG-EVO is an extension of our previous MHG-Finder package. MHG-EVO improves our prior work two-fold: (i) it scales to larger bacterial genome datasets, (ii) it provides intermediate MHGs involving subsets of input genomes. MHG-EVO leverages the MinHash-based alignment-free distance estimator Mash, enabling a fast calculation of pairwise genomic distances between input genomes. Subsequently, neighbor-joining is employed on the distance matrix to sketch a guide tree determining the MHG partition order. Given a guide tree, the MHG partitioning at each internal node involves five steps: 
(1) BLASTn for pairwise local alignments
(2) sequence pile-up
(3) initial MHGs and alignment graph construction
(4) alignment graph traversal
(5) representative sequence calculation

The guide tree used by MHG-EVO can be either user-provided or auto-estimated. And the final output of MHG-EVO comprises an MHG set per internal node, which encompasses all taxa located beneath the respective internal node, stored in the directory specified as ```-tg```.

\
![Algorithm Overview](https://github.com/yongze-yin/MHG-EVO/blob/main/algorithm.png)

```
MHG-EVO: Guide-Tree Based Homology Group Finder

optional arguments:
  -h, --help            show this help message and exit
  -g GENOME, --genome GENOME
                        Genome nucleotide sequence directory
  -o MHG_OUTPUT_DIR, --mhg_output_dir MHG_OUTPUT_DIR
                        Directory storing all MHG-EVO outputs
  -tg TEMP_DIR, --temp_dir TEMP_DIR
                        Directory storing a copy of each genome and
                        representative genomes
  -b BLASTN_DIR, --blastn_dir BLASTN_DIR
                        Directory storing blastn results
  --mash_tree_path MASH_TREE_PATH
                        Mash-estimated guide tree path
  --customized_tree_path CUSTOMIZED_TREE_PATH
                        Path to customized tree instead of using auto-
                        estimated tree
  -r REROOT, --reroot REROOT
                        Boolean value determining whether to reroot the guide
                        tree or not. If this is set to True, it will reroot
                        the guide tree changing the MHG output order to
                        minimize the tree height. If this is False, it will keep 
                        the MHG visiting order as it is.
  -k KMER_SIZE, --kmer_size KMER_SIZE
                        Kmer size for Mash, default 16
  -t THREAD, --thread THREAD
                        Number of threads
  -a ALIGNMENT_LENGTH_THRESHOLD, --alignment_length_threshold ALIGNMENT_LENGTH_THRESHOLD
                        Alignment length threshold in base pair. MHG-EVO does
                        not consider MHGs shorter than 60 base pairs by
                        default. If this threshold is too low, it will result
                        in an excessive amount of short MHGs and a longer
                        runtime.
```

MHG-EVO supports adding new genomes or continuing unfinished MHG partitioning using ```mhg-evo-add-genome```. It first locates the newly added genomes in the new guide tree, and compare the internal node differences between the old guide tree and the new guide tree, and only computes an MHG set for internal nodes which have not yet been computed before to avoid the redundant computation.

```
mhg-evo-add-genome: add new genomes to an already partitioned set of existing
genomes

optional arguments:
  -h, --help            show this help message and exit
  -og OLD_TEMP_GENOME_DIR, --old_temp_genome_dir OLD_TEMP_GENOME_DIR
                        Path to the old mhg-evo temp genome directory, this is
                        NOT the previous genome directory, but the mhg-evo
                        temporary genome directory
  -om OLD_MHG_OUTPUT_DIR, --old_mhg_output_dir OLD_MHG_OUTPUT_DIR
                        Path to the old MHG-EVO output directory which should
                        contain an MHG set for each internal node
  -g NEW_GENOME_DIR, --new_genome_dir NEW_GENOME_DIR
                        Directory contaning the new set of genomes
  -o MHG_OUTPUT_DIR, --mhg_output_dir MHG_OUTPUT_DIR
                        Directory storing the new MHG output
  -tg NEW_TEMP_GENOME_DIR, --new_temp_genome_dir NEW_TEMP_GENOME_DIR
                        Directory storing a copy of each genome and
                        representative genomes
  -b BLASTN_DIR, --blastn_dir BLASTN_DIR
                        Directory storing blastn results
  --mash_tree_path MASH_TREE_PATH
                        Mash-estimated guide tree path
  --customized_tree_path CUSTOMIZED_TREE_PATH
                        Path to customized tree instead of using auto-
                        estimated tree
  -r REROOT, --reroot REROOT
                        Boolean value determining whether to reroot the guide
                        tree or not. If this is set to True, it will reroot
                        the guide tree changing the MHG output order to
                        minimize the tree height. If this is False, it will keep 
                        the MHG visiting order as it is.
  -k KMER_SIZE, --kmer_size KMER_SIZE
                        Kmer size for Mash, default 16
  -t THREAD, --thread THREAD
                        Number of threads
  -a ALIGNMENT_LENGTH_THRESHOLD, --alignment_length_threshold ALIGNMENT_LENGTH_THRESHOLD
                        Alignment length threshold in base pair. MHG-EVO does
                        not consider MHGs shorter than 60 base pairs by
                        default. If this threshold is too low, it will result
                        in an excessive amount of short MHGs and a longer
                        runtime.
```


## Testing Case
```
python3 mhg-evo.py -g test_data/
```

