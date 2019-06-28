# arcasHLA quantification v2 #

## Setup ##
```
git clone https://github.com/roseorenbuch/arcasHLA.git
```
## Extract ##
arcasHLA extract now extracts pairs with one or both mates aligned to an HLA locus (all HLA genes and pseudogenes). Add `--chr6` if you want to extract all reads aligned to chromosome 6.


## Mismatch Quantification ##

#### Necessary Input ####
- `--subject SUBJECT` : subject name (to locate genotype in genotype tsv file)
- `--sample SAMPLE` : sample name (for output file)
- `--genotypes FILE` : genotypes.tsv (see below for necessary format)

| subject      	| A1         	| A2         	| B1         	| B2         	| C1         	| C2         	|
|--------------	|------------	|------------	|------------	|------------	|------------	|------------	|
| subject_name 	| A*01:01:01 	| A*01:01:01 	| B*07:02:01 	| B*07:02:01 	| C*04:01:01 	| C*04:01:01 	|

#### Options ####
- `--mapping MAPPING` : single, 2-field, g-group (default: g-group)
- `--count_both_reads` : if both mates overlap the same mismatch, both reads contribute to read count.
- `-k INT` : k-mer length (default: 31)
- `-o, --outdir` : out directory


#### Example ####
```
./arcasHLA quant --subject Pt23 \
                 --sample Pt23_pre \
                 --genotypes riaz.genotypes.tsv \
                 /Volumes/riaz_rna/filtered/Pt23_pre.filtered.1.fq.gz \
                 /Volumes/riaz_rna/filtered/Pt23_pre.filtered.2.fq.gz
```

The resulting json file follows the following dictionary:

```python
'genotyped_genes' : [A, B, C, ...] # Based on genotypes.tsv
'gene_count'      : {gene : count, ...} # Count of unique reads, gene == "" denotes unaligned reads
'gene_abundance'  : {gene : abundance, ...} # Relative abundance based on unique read count
gene              : {#described below}
```

This dictionary also contains a dictionary for each genotyped gene with the following objects.
```python
'allele1' : str # allele 1
'allele2' : str # allele 2
'n_mismatches_sites' : int # number of mismatches
'mismatch_locations' : {mismatch: [allele1_idx, allele2_idx], ...} # mismatch locations
```
For each read pair that uniquely maps to one allele and overlaps the mismatch, add the gene eq class to a list then count the number of times each gene is seen.
```python
'mismatch_genes' : {mismatch: [Counter(allele1_genes), Counter(allele2_genes)], ...}
```
The counts and frequencies for each allele are stored as lists.
```python
'allele_count' : [allele1_count, allele2_count],
'allele_freq' : [allele1_freq, allele2_freq]
```