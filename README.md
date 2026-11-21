# GoBi
The goal of this repo is to create an automated workflow, that completes the gene annotation of a gene family completely and reliably - automating tasks, that thus have been done manually and labour-intense.

## gobi_annotate.py: the week-1 pipeline in one command

Finds every copy of a protein family in a genome and checks each copy
automatically. Works for any organism with a genome FASTA or an NCBI assembly.

```
genome + query proteins
  -> makeblastdb -> tblastn
  -> group hits into loci (one locus = one gene copy; tandem copies are split)
  -> samtools faidx (locus + flanks) -> miniprot per locus (no chimeras across copies)
  -> checks: query coverage, start M / stop, frameshifts, in-frame stops,
     GT-AG introns, assembly gaps (N) near the locus
  -> summary table, protein FASTA, GFF for the NCBI Genome Data Viewer, report
  -> optional: FAMSA alignment + IQ-TREE tree (UFBoot)
```

### Setup (once)

```bash
mamba env create -f environment.yml      # creates the separate environment "gobi_annotate"
mamba activate gobi_annotate
```

The pipeline has its own environment, so course environments (e.g. `gobi`) are never
overwritten. If `mamba env create` asks "Overwrite?", it refers to `gobi_annotate` only.

### Run

```bash
# NCBI accession + UniProt accessions (both downloaded automatically)
python gobi_annotate.py --genome GCA_009733165.1 --query Q7T1K6,Q7ZZN8

# local files, with tree
python gobi_annotate.py --genome my_genome.fna --query my_queries.faa --tree
```

Results go to `runs/<genome>_<date>/`; start with `report.md`.
Genomes and BLAST databases are cached in `~/gobi_cache` (path without spaces),
so a second run on the same genome starts at the search. Steps whose output
already exists are skipped; `--force` redoes everything.

### Output

| File | Content |
|---|---|
| `report.md` | settings, tool versions, one line per locus with verdict |
| `summary.tsv` | all checks per locus |
| `hits_by_locus.tsv` | every tblastn hit and the locus it belongs to |
| `models_best.gff` | best miniprot model per locus (genome coordinates) |
| `models_all.gff` | all miniprot models (genome coordinates) |
| `proteins.faa` | translated best model per locus |
| `loci/` | window FASTA + raw miniprot output per locus |
| `tree/` | `alignment.faa`, `tree.treefile` (with `--tree`) |
| `commands.log` | every command that was run, with time stamps |

Verdicts: **complete** (covers ≥90 % of the query, M … stop, all introns GT-AG, no frameshift),
**partial** (part missing), **suspicious** (frameshift or in-frame stop: pseudogene or
sequencing/assembly error), **no_model**. These are automatic checks: look at the loci by
eye (NCBI viewer: Tracks > Configure Tracks > Custom Data > Local File > `models_best.gff`).

Important options: `--evalue` (1e-3), `--max-intron` (20000 bp), `--flank` (10000 bp),
`--min-cov` (0.9), `--refs` (extra proteins for the tree), `--threads`. See `--help`.
