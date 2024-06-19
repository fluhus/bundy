# Bundy
A microbial abundance estimator.

## Requirements

- [Bowtie2](https://github.com/BenLangmead/bowtie2/releases)
  in the system PATH.
- A microbial reference genome.

## Usage

### Index reference genome (once per reference)

If you haven't yet, create a Bowtie2 index for the reference genome.

```
bowtie2-build my_genome.fa my_bowtie_index
```

Add `--threads N` to run on N threads.
You may also need to add `--large-index` if the reference genome is large.

### Preprocess with bundyx (once per reference)

This generates a small output that helps bundy with abundance normalizations.

```
bundyx -i my_genome.fa -r my_bowtie_index -l READ_LENGTH -o bundyx.out
```

1. Read length should match the read length in the future input files.
   It doesn't have to be exactly the same but the closer the better.
   For example, 100 can usually cover reads of lengths from 70 to 150 bases long.
2. Add `-t N` to run on N threads.
3. If running on a queue system like sbatch or qsub,
   you can break the process down into sub-jobs for more parallelization.
   Add `-p PART_NUMBER -np TOTAL_NUM_OF_PARTS`,
   for example `-p 5 -np 100` means that this is sub-job 5 out of 100.
   **Make sure that each sub-job uses a separate output file
   (bundy will unite them later).**

### Abundance estimation

Run this on each query (fastq) file.

```
bundy -i my_data.fq -r my_bowtie_index -x bundyx.out
```

1. Add `-t N` to run on N threads.
2. If the reference contains multiple contigs per species,
   bundy needs to know which contigs to group together.
   It uses a regular expression to extract the group name
   out of each sequence name.  
   For example if the sequence names are "species_1_contig_1",
   "species_1_contig_2", "species_2_contig_1", "species_2_contig_2"...
   provide `-n species_\\d+` to group by species number.
3. Use `-h` for help about additional options.
