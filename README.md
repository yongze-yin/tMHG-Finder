# MHG-EVO

![Algorithm Overview](https://github.com/yongze-yin/MHG-EVO/blob/main/algorithm.png)

## **MHG-EVO** _de Novo_ Tree-Based Maximal Homologous Group Finder
```
MHG-EVO: Guide-Tree Based Homology Group Finder

optional arguments:
  -h, --help            show this help message and exit
  -g GENOME, --genome GENOME
                        Genome nucleotide sequence directory (Required)
  --temp_dir TEMP_DIR   Directory storing a copy of each genome for blastn
  --kmer_size KMER_SIZE
                        Kmer size for Mash, default 16
  -t THREAD, --thread THREAD
                        Number of threads
  --mash_tree_path MASH_TREE_PATH
                        Mash-estimated guide tree path
  --blastn_dir BLASTN_DIR
                        Directory storing blastn results
  -o MHG_OUTPUT_DIR, --mhg_output_dir MHG_OUTPUT_DIR
                        Directory storing all MHG-EVO outputs
  --reroot REROOT       Boolean value determining whether to reroot the guide tree or not. If this is set to True, it
                        will reroot the guide tree changing the MHG output order to minimize the number of total
                        internal nodes. If this is False, it will keep the MHG visiting order as it is.
  --alignment_length_threshold ALIGNMENT_LENGTH_THRESHOLD
                        Alignment length threshold in base pair. MHG-EVO does not consider MHGs shorter than 60 base
                        pairs by default. If this threshold is too low, it will result in an excessive amount of short
                        MHGs and a longer runtime.
```


## Testing Case
```
python3 mhg-evo.py -g test_data/
```

