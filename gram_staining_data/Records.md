# A use case recording the procoess of analyzing gram-staining data

### Prerequist
- [gotree](https://github.com/evolbioinfo/gotree)
- [python >= 3.8.8](https://www.python.org/downloads/release/python-388/)
- Evolink

We prepared essential files in `data/` and `scripts/` folders to help users to obtain essentia inputs for Evolink. 
- `data/WOL_genomes_traits.tsv`: A tab-delimited file containing five columns. They are genome (WOL species-level genome ID), species_tax_id (the NCBI taxonomy ID of the species for the genome), data_source (which database the phenotype is extracting from), gram-stain (the genome is annotated as of gram positive/negative function or not) and motility (the genome is annotated as of motility, gram_staining function or not). 
- `data/WOL_tree.nwk`: A newick format tree file with all species-level genomes in WOL (https://biocore.github.io/wol/), downloaded from https://biocore.github.io/wol/data/trees/tree.nwk.
- `data/gram_staining_emapper.og_tsv`: A tab-delimited file containing three columns. They are genome (WOL species-level genome ID), gene (gene ID) and orthogroup (gene family annotated by eggNOG mapper).
- `scripts/og_tsv2pa_table.py`: A custom python to convert `data/gram_staining_emapper.og_tsv` into an orthogroup presence/absence tab-delimited table.

### Step.1 Prepare tree file
```
# change to work directory
cd gram_staining_data/
unzip data.zip
mkdir -p tmp

# filter out genomes without motility annotation; filter out G000715975 because of its failed gene annotation
awk -F"\t" 'NR>1&&$4!="NA"&&($4=="negative"||$4=="positive"){print $1"\t"$4}' data/WOL_genome_trait.tsv | sort | uniq -c |awk '{print $2"\t"$3"\t"$1}' | grep -vw "G000715975" > tmp/gram_staining.tsv

# exclude species with some of the genomes having gram_staining function and others not
cut -f1 tmp/gram_staining.tsv | sort | uniq -d > tmp/conflict_id.list
awk -F"\t" 'BEGIN{OFS="\t"}NR==FNR{a[$1];next}!($1 in a){print $0}' tmp/conflict_id.list tmp/gram_staining.tsv > tmp/gram_staining_clean.tsv

# extract a subtree from the whole WOL species tree by keeping annotated species
gotree prune -r -f <(cut -f1 tmp/gram_staining_clean.tsv) -i data/WOL_tree.nwk -o tree.nwk
```
Output `tree.nwk` is  a newick format tree.

### Step.2 Prepare phenotype file
```
cut -f1-2 tmp/gram_staining_clean.tsv > tmp/select_tips.list
# get formatted file by defining two columns. 'Tip' is for species ID and 'Status' is 0/1 indicating gram positive/gram negative, respectively.
awk -F"\t" 'BEGIN{print "Tip\tStatus"}$2=="positive"{print $1"\t0"}$2=="negative"{print $1"\t1"}' tmp/select_tips.list > trait.tsv
```
Output `trait.tsv` is a two-column file labeling teh phenotypic presence/absence of each species.

### Step.3 Prepare genotype file
```
# convert paris of gene IDs and orthogroups into an orthogroup presence/absence tab-delimited table.
SCRIPT=scripts/og_tsv2pa_table.py
gidlist=tmp/gram_staining_gid.list
cut -f1 tmp/select_tips.list > ${gidlist}
ogtsv=data/gram_staining_emapper.og_tsv
python ${SCRIPT} -i ${ogtsv} -g ${gidlist} -o tmp/gram_staining_emapper -b
sed "s/Orthogroup/#OTU ID/g" tmp/gram_staining_emapper.BiPA.table > gene.tsv
```
Output `gene.tsv` is an orthogroup (row) X species (column) tab-delimited table.

### Step.4 run Evolink and generate plots
```
python ../Evolink.py -g gene.tsv -t trait.tsv -n tree.nwk -o Evolink_output -f > output.log
```
Output:  
- `Evolink_output_plot/result.tsv`: A tab-delimited table containing five columns, orthoID (orthogroup ID), Prevalence_index, Evolink_index, scores (outlier scores assigned by outlier detection method) and significance (if a gene is significantly associated with flagella. `sig` means significant and `NA` means not). 