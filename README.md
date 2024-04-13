# tMHG-Finder: _de Novo_ Tree-Based MHG Finder

We are now moving tMHG-Finder to conda for an easier installation. 
\
\
tMHG-Finder is an extension of our previous MHG-Finder package. tMHG-Finder improves our prior work two-fold: (i) it scales to larger bacterial genome datasets, (ii) it provides intermediate MHGs involving subsets of input genomes. tMHG-Finder leverages the MinHash-based alignment-free distance estimator Mash, enabling a fast calculation of pairwise genomic distances between input genomes. Subsequently, neighbor-joining is employed on the distance matrix to sketch a guide tree determining the MHG partition order. Given a guide tree, the MHG partitioning at each internal node involves five steps: 
(1) BLASTn for pairwise local alignments
(2) sequence pile-up
(3) initial MHGs and alignment graph construction
(4) alignment graph traversal
(5) representative sequence calculation

The guide tree used by tMHG-Finder can be either user-provided or auto-estimated. And the final output of tMHG-Finder comprises an MHG set per internal node, which encompasses all taxa located beneath the respective internal node, stored in the directory specified as ```-tg```.

\
![Algorithm Overview](https://github.com/yongze-yin/tMHG-Finder/blob/main/algorithm.png)

## Installation Option: git clone 
Using git clone, please install the below dependencies manually:

The following command will install tMHG-Finder an all python dependencies.
```
python -m pip install .
```

Additional non-python dependencies are listed below:
> [BEDtools](https://bedtools.readthedocs.io/en/latest/)

> [MAFFT](https://mafft.cbrc.jp/alignment/software/source.html)

> [Mash>](https://github.com/marbl/Mash)


These python packages will be installed along with tMHG-Finder using the command above.
> [numpy](https://pypi.org/project/numpy/)

> [pandas](https://pypi.org/project/pandas/)

> [Networkx](https://networkx.org/)

> [Biopython](https://biopython.org/)

> [Pathos](https://pypi.org/project/pathos/)
```
$ tmhgf find --help
usage: tmhgf find [-h] [-b BLASTN_DIR] [--mash_tree_path MASH_TREE_PATH] [--customized_tree_path CUSTOMIZED_TREE_PATH] [-r REROOT] [-k KMER_SIZE] [-t THREAD] [-a ALIGNMENT_LENGTH_THRESHOLD] -g GENOME
                  [-o MHG_OUTPUT_DIR] [-tg TEMP_DIR]

options:
  -h, --help            show this help message and exit
  -b BLASTN_DIR, --blastn_dir BLASTN_DIR
                        Directory storing blastn results
  --mash_tree_path MASH_TREE_PATH
                        Mash-estimated guide tree path
  --customized_tree_path CUSTOMIZED_TREE_PATH
                        Path to customized tree instead of using auto-estimated tree
  -r REROOT, --reroot REROOT
                        Boolean value determining whether to reroot the guide tree or not. If this is set to True, it will reroot the guide tree changing the MHG output order to minimize the tree height. If this
                        is False, it will keep the MHG visiting order as it is.
  -k KMER_SIZE, --kmer_size KMER_SIZE
                        Kmer size for Mash, default 16
  -t THREAD, --thread THREAD
                        Number of threads
  -a ALIGNMENT_LENGTH_THRESHOLD, --alignment_length_threshold ALIGNMENT_LENGTH_THRESHOLD
                        Alignment length threshold in base pair. tMHG-Finder does not consider MHGs shorter than 60 base pairs by default. If this threshold is too low, it will result in an excessive amount of
                        short MHGs and a longer runtime.
  -g GENOME, --genome GENOME
                        Genome nucleotide sequence directory
  -o MHG_OUTPUT_DIR, --mhg_output_dir MHG_OUTPUT_DIR
                        Directory storing all tMHG-Finder outputs
  -tg TEMP_DIR, --temp_dir TEMP_DIR
                        Directory storing a copy of each genome and representative genomes
```

tMHG-Finder supports adding new genomes or continuing unfinished MHG partitioning using ```tmhgf add```. It first locates the newly added genomes in the new guide tree, and compare the internal node differences between the old guide tree and the new guide tree, and only computes an MHG set for internal nodes which have not yet been computed before to avoid the redundant computation.

```
$ tmhgf add --help
usage: tmhgf add [-h] [-b BLASTN_DIR] [--mash_tree_path MASH_TREE_PATH] [--customized_tree_path CUSTOMIZED_TREE_PATH] [-r REROOT] [-k KMER_SIZE] [-t THREAD] [-a ALIGNMENT_LENGTH_THRESHOLD] -og OLD_TEMP_GENOME_DIR
                 -om OLD_MHG_OUTPUT_DIR -g NEW_GENOME_DIR [-o MHG_OUTPUT_DIR] [-tg NEW_TEMP_GENOME_DIR]

options:
  -h, --help            show this help message and exit
  -b BLASTN_DIR, --blastn_dir BLASTN_DIR
                        Directory storing blastn results
  --mash_tree_path MASH_TREE_PATH
                        Mash-estimated guide tree path
  --customized_tree_path CUSTOMIZED_TREE_PATH
                        Path to customized tree instead of using auto-estimated tree
  -r REROOT, --reroot REROOT
                        Boolean value determining whether to reroot the guide tree or not. If this is set to True, it will reroot the guide tree changing the MHG output order to minimize the tree height. If this
                        is False, it will keep the MHG visiting order as it is.
  -k KMER_SIZE, --kmer_size KMER_SIZE
                        Kmer size for Mash, default 16
  -t THREAD, --thread THREAD
                        Number of threads
  -a ALIGNMENT_LENGTH_THRESHOLD, --alignment_length_threshold ALIGNMENT_LENGTH_THRESHOLD
                        Alignment length threshold in base pair. tMHG-Finder does not consider MHGs shorter than 60 base pairs by default. If this threshold is too low, it will result in an excessive amount of
                        short MHGs and a longer runtime.
  -og OLD_TEMP_GENOME_DIR, --old_temp_genome_dir OLD_TEMP_GENOME_DIR
                        Path to the old tMHG-Finder temp genome directory, this is NOT the previous genome directory, but the tMHG-Finder temporary genome directory
  -om OLD_MHG_OUTPUT_DIR, --old_mhg_output_dir OLD_MHG_OUTPUT_DIR
                        Path to the old tMHG-Finder output directory which should contain an MHG set for each internal node
  -g NEW_GENOME_DIR, --new_genome_dir NEW_GENOME_DIR
                        Directory contaning the new set of genomes
  -o MHG_OUTPUT_DIR, --mhg_output_dir MHG_OUTPUT_DIR
                        Directory storing the new MHG output
  -tg NEW_TEMP_GENOME_DIR, --new_temp_genome_dir NEW_TEMP_GENOME_DIR
                        Directory storing a copy of each genome and representative genomes
```


## Testing Case
```
tmhgf find -g test_data/
```

