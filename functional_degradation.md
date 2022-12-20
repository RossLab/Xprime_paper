# Functional degradation

Software used:
- Orthofinder
- SnpEff and SnpSift
- Bowtie2
- Samtools
- Picardtools
- Freebayes
- VCFtools
- GenomicsGeneral tools (https://github.com/simonhmartin/genomics_general)
- SNPsplit
- bedtools
- STAR
- EdgeR

## I. IDENTIFYING HOMOLOGS

I've discussed with the other authors about whether gene copies on the X and X' should be called homologs, orthologs, paralogs, or alleles. Paralogs are gene copies that arise due to duplication, which isn't technically the case here. Orthologs arise due to speciation, which is also technically not the case, although it's similar in that there are two gene copies on distinct lineages that diverged from a common ancestral gene. A homolog more broadly refers to gene copies that are descended fron a common gene, so that certainly applies here. I'm not sure if the term allele extends to 'alleles' that are no longer recombining. Other papers that compare genes on differentiated X and Y chromosomes often seem to refer to the gene copies as 'homologs', so that's what I'll do here.

Anyway, I'll use Orthofinder to identify homologous X and X' gene copies. With the change-point analysis I identified the 'main' breakpoints as 4.1Mb and 62.9Mb along the X chromosome, so I'll conservatively say that 4.05-63.95Mb is the region of the X homologous to the nonrecombining portion of X'. I'll pull coding sequences that fall within this region for identifying X-X' homologs:

```
# Get all X-linked and Inv-linked transcripts and their start/end positions:
cat augustus.hints.gff3 | grep "X" | grep "mRNA" | cut -f4,5,9 | sed 's/ID=//g' | sed 's/;.*$//g' | awk '$1 > 4050000' | awk '$2 < 62950000' | cut -f3 > Xhom_trx.txt
cat augustus.hints.gff3 | grep "Inversion" | grep "mRNA" | cut -f4,5,9 | sed 's/ID=//g' | sed 's/;.*$//g' | cut -f3 > Inv_trx.txt

# Pull longest transcripts from the aa file:
cat augustus_hints.aa | awk '/^>/ {if(N>0) printf("\n"); printf("%s\t",$0);N++;next;} {printf("%s",$0);} END {if(N>0) printf("\n");}' | tr "." "\t" | awk -F '	'	'{printf("%s\t%d\n",$0,length($3));}' | sort -t '	' -k1,1 -k4,4nr | sort -t '	' -k1,1 -u -s | sed 's/	/./' | cut -f 1,2 | tr "\t" "\n"	| fold -w 60 > augustus.proteins.longest.fasta

# Subset longest transcripts:
seqtk subseq augustus.proteins.longest.fasta Xhom_trx.txt > Xhom_longest_proteins.fasta
seqtk subseq augustus.proteins.longest.fasta Inv_trx.txt > Inv_longest_proteins.fasta

grep ">" Inv_longest_proteins.fasta | wc -l
# 3470
grep ">" Xhom_longest_proteins.fasta | wc -l
# 3429
```

So 3470 genes on the X' nonrecombining region compared to 3429 within the homologous portion of the X. The X region actually covers slightly more bp than the X' region (58.8 versus 54.7). However gene content isn't actually higher on the X' because the genes sum to 4.2Mb compared to 5.2Mb on the homologous X portion.

Run Orthofinder:
```
mkdir genes
mv Inv_longest_proteins.fasta genes/Inversion.fasta
mv Xhom_longest_proteins.fasta genes/X_hom_region.fasta
orthofinder -f genes/
```

The relevant output files (for me) are:
- Orthogroups_SingleCopyOrthologues.txt
- Inversion__v__X_hom_region.tsv
- Orthogroups.GeneCount.tsv
- Orthogroups_UnassignedGenes.tsv

I'll analyse these outputs in RStudio:

```
# load in single-copy orthologs
genes_assigned <- read.table('Inversion__v__X_hom_region.tsv', fill=T, sep='\t', header=T, stringsAsFactors = F, na.strings=c("","NA")) # 
single_copy_orthologs <- read.table('Orthogroups_SingleCopyOrthologues.txt')
colnames(single_copy_orthologs) <- c('Orthogroup')
single_copy_assignments <- merge(genes_assigned, single_copy_orthologs, by=c('Orthogroup'))
nrow(single_copy_assignments) # 2321
write.table(single_copy_assignments, file='orthofinder_single_copy_assignments.tsv', quote=F, row.names=F, col.names=T, sep='\t')
```

i.e. 2321 single-copy homologs assigned between the X and X'. What about duplicate homologs?

```
ortholog_counts <- read.table('Orthogroups.GeneCount.tsv', header=T, stringsAsFactors=F)
duplications <- ortholog_counts[which(ortholog_counts$Inversion > 1 | ortholog_counts$X_hom_region > 1),]
duplications <- duplications[which(duplications$Inversion > 0 & duplications$X_hom_region > 0),]
nrow(duplications) # 296 - no. orthogroups
sum(duplications$Inversion) # 527
sum(duplications$X_hom_region) # 679
nrow(duplications[which(duplications$Inversion > duplications$X_hom_region),]) # 70
nrow(duplications[which(duplications$X_hom_region > duplications$Inversion),]) # 162
nrow(duplications[which(duplications$X_hom_region == duplications$Inversion),]) # 64
```

i.e. 162 genes where the X has more copies than the X', compared to 70 where the X' has more copies. In 64 cases they have the same number of copies. The X having more copies is most likely due to either X duplications or X' deletions, while the X' having more copies is most likely due to X deletions or X' duplications. This is because X' mutations are less likely to be purged, making X' deletions/duplications the more likely explanation.

What about the no-hits?

```
# load in no-hits (single-copy)
genes_unassigned <- read.table('Orthogroups_UnassignedGenes.tsv', fill=T, sep='\t', header=T, stringsAsFactors = F, na.strings=c("","NA")) # 
genes_unassigned[is.na(genes_unassigned)] <- "none"
X_unassigned <- genes_unassigned[which(genes_unassigned$X_hom_region != "none"),] # 317 (genes specific to X, potentially lost from inversion)
Inv_unassigned <- genes_unassigned[which(genes_unassigned$Inversion != "none"),] # 587 (genes specific to inversion, potentially gained on inversion)
# no hits (more than one copy)
nrow(ortholog_counts[which(ortholog_counts$Inversion == 0 & ortholog_counts$X_hom_region > 0),])# 42 groups
inv_no_hit <- ortholog_counts[which(ortholog_counts$Inversion == 0 & ortholog_counts$X_hom_region > 0),]
sum(inv_no_hit$X_hom_region) # 112
nrow(ortholog_counts[which(ortholog_counts$Inversion > 0 & ortholog_counts$X_hom_region == 0),]) # 16 groups
X_no_hit <- ortholog_counts[which(ortholog_counts$Inversion > 0 & ortholog_counts$X_hom_region == 0),]
sum(X_no_hit$Inversion) # 35
```

So basically, 429 genes across 359 OGs specific to the X; 622 genes across 603 OGs specific to the X'

```
# write out single-copy orthologs:
single_copy_assignments_inv <- single_copy_assignments[,c(1,2)]
single_copy_assignments_X <- single_copy_assignments[,c(1,3)]
single_copy_assignments_inv$chrom <- "Inversion"
single_copy_assignments_X$chrom <- "Xz" # Xz to distinguish from other X-linked genes for future analyses
colnames(single_copy_assignments_inv) <- c('Orthogroup', 'gene', 'chrom')
colnames(single_copy_assignments_X) <- c('Orthogroup', 'gene', 'chrom')
single_copy_assignments_all <- merge(single_copy_assignments_inv, single_copy_assignments_X, by=c('Orthogroup', 'gene', 'chrom'), all=T)
write.table(single_copy_assignments_all, file='X_inv_single_copy_assignments.tsv', sep='\t', quote=F, col.names=T, row.names=F)
```

## II. PSEUDOGENIZING MUTATIONS (DISRUPTED GENES)

I already have variant call files for the X'X, XX and X0 genotypes against Bcop_v2 from the previous section (the SNV analysis for looking at X-X' divergence). I can use these files to extract variant information for genes of interst - i.e. those on the X'. I'm only going to look at the single-copy homologs because:
- If I consider all genes within the inversion, there may be novel or lost X'-linked genes. I care about genes where both chromosomes have copies, and what proportion of those are degraded.
- When it comes to duplicates, there may be cases where e.g. the X has one copy and the X' has two copies due to duplication, and one duplicated copy has accumulated mutations but the other has not. This may not neccessarily be due to lack of recombination but rather due to redundancy of that extra gene copy, which I think is a different kettle of fish.

Basically, I'll assume X' reads will have force-mapped to the X and then I can identify where genes have accumulated a potentially pseudogenizing mutation - defined as a frameshift and/or gain/loss of stop and start codons. I can also use the X0 variant calls to ensure that any variants I identify are indeed specific to the X' - if variants are present in both the X'X and X0 calls, then they may be due to misassemblies (shared by X' and X reads mapped against the X reference). 

**Extracting variant information**

I'll use SnpEff and SnpSift for this. First I need to use the annotation to build a database. I'll use Bcop_v2 and pull just chromosomes II, III, IV and X from my annotation. SnpEff wants .gtf files, but for some reason it doesn't like the .gtf files produced by BRAKER, so I had to convert the augustus_hints.gff3 file to gtf with gffread:
```
gffread augustus.hints.gff3 -T -o augustus.hints.gtf
# pull relevant chromosomes
grep -w "II" augustus_gff3.gtf >> Bcop_v2.augustus.gtf && grep -w "III" augustus_gff3.gtf >> Bcop_v2.augustus.gtf && grep -w "IV" augustus_gff3.gtf >> Bcop_v2.augustus.gtf && grep -w "X" augustus_gff3.gtf >> Bcop_v2.augustus.gtf
```

Build database:
```
cd ~/.conda/envs/SnpEff/SnpEff/data/Bradysia_coprophila
cp /data/ross/flies/analyses/john_urban_bcop_v2_assembly/Bcop_v2-chromosomes.fasta .
cp /data/ross/flies/analyses/bradysia_sex_determination/004_genome_annotation/results/6_Bcop_v3_anno/braker/gffread/Bcop_v2.augustus.gtf .
mv Bcop_v2-chromosomes.fasta sequences.fa
mv Bcop_v2.augustus.gtf genes.gtf
java -Xmx20g -jar ~/.conda/envs/SnpEff/SnpEff/snpEff.jar build -v Bradysia_coprophila
```

Now extract the variant information required:
```
# SNPs
for file in $(ls *.SNPs.vcf.gz)
do
	base=$(basename $file ".SNPs.vcf.gz")
	java -Xmx4g -jar ~/.conda/envs/SnpEff/SnpEff/snpEff.jar -c ~/.conda/envs/SnpEff/SnpEff/snpEff.config -v -no-downstream -no-intron -no-upstream -no-utr -no-intergenic  -stats ${base}.filterstats.html -canon Bradysia_coprophila ${base}.SNPs.vcf.gz > ${base}.SNPs.canon.ann.vcf
	java -Xmx4g -jar ~/.conda/envs/SnpEff/SnpEff/SnpSift.jar extractFields ${base}.SNPs.canon.ann.vcf CHROM ANN[0].GENEID POS REF ALT isHom"(GEN[0])" ANN[*].EFFECT > ${base}_justgenessnps.txt
done

# indels
for file in $(ls *.indels.vcf.gz)
do
	base=$(basename $file ".indels.vcf.gz")
	java -Xmx4g -jar ~/.conda/envs/SnpEff/SnpEff/snpEff.jar -c ~/.conda/envs/SnpEff/SnpEff/snpEff.config -v -no-downstream -no-intron -no-upstream -no-utr -no-intergenic  -stats ${base}.filterstats.html -canon Bradysia_coprophila ${base}.indels.vcf.gz > ${base}.indels.canon.ann.vcf
	java -Xmx4g -jar ~/.conda/envs/SnpEff/SnpEff/SnpSift.jar extractFields ${base}.indels.canon.ann.vcf CHROM ANN[0].GENEID POS REF ALT isHom"(GEN[0])" ANN[*].EFFECT > ${base}_justgenesindels.txt
done
```

Below is the analysis of the outputs in RStudio. Basically, I'm finding consensus variant calls between the two replicates within each genotype (to identify fixed differences between X and X'), excluding common variants between X and X', and then identifying pseudogenizing SNPs and indels.

```
### 1. SNPs
XpX_SNPs_r1 <- read.table('XpX_adult_female_sciara_coprophila_rep1-bwa.g_justgenessnps.txt', sep="\t", header=T, fill=T, stringsAsFactors=F, na.strings=c("","NA")) # 642368 rows
XpX_SNPs_r1 <- XpX_SNPs_r1[complete.cases(XpX_SNPs_r1),] # # 40394 rows
XpX_SNPs_r1 <- XpX_SNPs_r1[,c(1,2,3,4,5,7)]
colnames(XpX_SNPs_r1) <- c('chrom', 'gene_id', 'pos', 'ref', 'alt', 'effect')
XpX_SNPs_r2 <- read.table('XpX_adult_female_sciara_coprophila_rep2-bwa.g_justgenessnps.txt', sep="\t", header=T, fill=T, stringsAsFactors=F, na.strings=c("","NA")) # 577367 rows
XpX_SNPs_r2 <- XpX_SNPs_r2[complete.cases(XpX_SNPs_r2),] # 36274 rows
XpX_SNPs_r2 <- XpX_SNPs_r2[,c(1,2,3,4,5,7)]
colnames(XpX_SNPs_r2) <- c('chrom', 'gene_id', 'pos', 'ref', 'alt', 'effect')
# merge
XpX_SNPs_merged <- merge(XpX_SNPs_r1, XpX_SNPs_r2, by=c(1:6)) # 28716 rows
# load in X0 variants - using as a control (to exclude common variants between X’X and X0 alignments that may have been called due to misassemblies on the X chromosome)
X0_SNPs_r1 <- read.table('XO_adult_male_sciara_coprophila_rep1-bwa.g_justgenessnps.txt', sep="\t", header=T, fill=T, stringsAsFactors=F, na.strings=c("","NA")) # 37831 rows
X0_SNPs_r1 <- X0_SNPs_r1[complete.cases(X0_SNPs_r1),] # 3486 rows
X0_SNPs_r1 <- X0_SNPs_r1[,c(1,2,3,4,5,7)]
colnames(X0_SNPs_r1) <- c('chrom', 'gene_id', 'pos', 'ref', 'alt', 'effect')
X0_SNPs_r2 <- read.table('XO_adult_male_sciara_coprophila_rep2-bwa.g_justgenessnps.txt', sep="\t", header=T, fill=T, stringsAsFactors=F, na.strings=c("","NA")) # 127313 rows
X0_SNPs_r2 <- X0_SNPs_r2[complete.cases(X0_SNPs_r2),] # 10894 rows
X0_SNPs_r2 <- X0_SNPs_r2[,c(1,2,3,4,5,7)]
colnames(X0_SNPs_r2) <- c('chrom', 'gene_id', 'pos', 'ref', 'alt', 'effect')
# merge
X0_SNPs_merged <- merge(X0_SNPs_r1, X0_SNPs_r2, by=c(1:6)) # 2618 rows
# remove XpX variants also present in X0
XpX_rem_X0_SNPs <- dplyr::anti_join(XpX_SNPs_merged, X0_SNPs_merged)
nrow(XpX_rem_X0_SNPs) # 27123 (1593 variants removed)
# Isolate desired variants ('start_lost', 'stop_gained', 'stop_lost)
XpX_rem_X0_SNPs %>% group_by(effect) %>% tally()
XpX_SNPs_pseudo <- XpX_SNPs_merged[which(XpX_SNPs_merged$effect == "start_lost" | XpX_SNPs_merged$effect == "stop_gained" | XpX_SNPs_merged$effect == "stop_gained&splice_region_variant" | XpX_SNPs_merged$effect == "stop_lost&splice_region_variant"),]
nrow(XpX_SNPs_pseudo) # 126 pseudogenising mutations
# N genes with pseudogenising point mutations
nrow(unique(XpX_SNPs_pseudo[c("gene_id")])) # 119 genes with at least 1 pseudogenising mutation
XpX_SNPs_X <- XpX_SNPs_pseudo[which(XpX_SNPs_pseudo$chrom == "X"),]
nrow(unique(XpX_SNPs_X[c("gene_id")])) # 78 genes with at least 1 pseudogenising mutation on X
XpX_SNPs_A <- XpX_SNPs_pseudo[which(XpX_SNPs_pseudo$chrom != "X"),]
nrow(unique(XpX_SNPs_A[c("gene_id")])) # 41 genes with at least 1 pseudogenising mutation on As

### 2. Indels
XpX_indels_r1 <- read.table('XpX_adult_female_sciara_coprophila_rep1-bwa.g_justgenesindels.txt', sep="\t", header=T, fill=T, stringsAsFactors=F, na.strings=c("","NA")) # 168139 rows
XpX_indels_r1 <- XpX_indels_r1[complete.cases(XpX_indels_r1),] # 1928 rows
XpX_indels_r1 <- XpX_indels_r1[,c(1,2,3,4,5,7)]
colnames(XpX_indels_r1) <- c('chrom', 'gene_id', 'pos', 'ref', 'alt', 'effect')
XpX_indels_r2 <- read.table('XpX_adult_female_sciara_coprophila_rep2-bwa.g_justgenesindels.txt', sep="\t", header=T, fill=T, stringsAsFactors=F, na.strings=c("","NA")) # 149583 rows
XpX_indels_r2 <- XpX_indels_r2[complete.cases(XpX_indels_r2),] # 1750 rows
colnames(XpX_indels_r2) <- c('chrom', 'gene_id', 'pos', 'ref', 'alt', 'effect')
# merge
XpX_indels_merged <- merge(XpX_indels_r1, XpX_indels_r2, by=c(1:6)) # 1125 rows
# load in X0 variants - using as a control (to exclude common variants between X’X and X0 alignments that may have been called due to misassemblies on the X chromosome)
X0_indels_r1 <- read.table('XO_adult_male_sciara_coprophila_rep1-bwa.g_justgenesindels.txt', sep="\t", header=T, fill=T, stringsAsFactors=F, na.strings=c("","NA")) # 11186 rows
X0_indels_r1 <- X0_indels_r1[complete.cases(X0_indels_r1),] # 172 rows
X0_indels_r1 <- X0_indels_r1[,c(1,2,3,4,5,7)]
colnames(X0_indels_r1) <- c('chrom', 'gene_id', 'pos', 'ref', 'alt', 'effect')
X0_indels_r2 <- read.table('XO_adult_male_sciara_coprophila_rep2-bwa.g_justgenesindels.txt', sep="\t", header=T, fill=T, stringsAsFactors=F, na.strings=c("","NA")) # 33920 rows
X0_indels_r2 <- X0_indels_r2[complete.cases(X0_indels_r2),] # 512 rows
X0_indels_r2 <- X0_indels_r2[,c(1,2,3,4,5,7)]
colnames(X0_indels_r2) <- c('chrom', 'gene_id', 'pos', 'ref', 'alt', 'effect')
# merge
X0_indels_merged <- merge(X0_indels_r1, X0_indels_r2, by=c(1:6)) # 88 rows
# remove XpX variants also present in X0
XpX_rem_X0_indels <- dplyr::anti_join(XpX_indels_merged, X0_indels_merged) # 1125 rows
nrow(XpX_rem_X0_indels) # 1075 (50 variants removed)
# Isolate desired variants
XpX_rem_X0_indels %>% group_by(effect) %>% tally()
XpX_indels_pseudo <- XpX_indels_merged[which(XpX_indels_merged$effect == "frameshift_variant" | XpX_indels_merged$effect == "frameshift_variant&splice_region_variant" | XpX_indels_merged$effect == "frameshift_variant&start_lost" | XpX_indels_merged$effect == "frameshift_variant&start_lost&splice_region_variant" | XpX_indels_merged$effect == "frameshift_variant&stop_gained" | XpX_indels_merged$effect == "frameshift_variant&stop_lost&splice_region_variant" | XpX_indels_merged$effect == "start_lost&disruptive_inframe_deletion" | XpX_indels_merged$effect == "stop_gained&disruptive_inframe_deletion" | XpX_indels_merged$effect == "stop_gained&disruptive_inframe_insertion"),]
nrow(XpX_indels_pseudo) # 388 pseudogenising mutations
# N genes with pseudogenising indels
nrow(unique(XpX_indels_pseudo[c("gene_id")])) # 257 genes with at least 1 pseudogenising mutation
XpX_indels_X <- XpX_indels_pseudo[which(XpX_indels_pseudo$chrom == "X"),]
nrow(unique(XpX_indels_X[c("gene_id")])) # 152 genes with at least 1 pseudogenising mutation on X
XpX_indels_A <- XpX_indels_pseudo[which(XpX_indels_pseudo$chrom != "X"),]
nrow(unique(XpX_indels_A[c("gene_id")])) # 105 genes with at least 1 pseudogenising mutation on As
XpX_indels_X %>% group_by(effect) %>% tally()

### Combine SNPs + indels + add Xp orthos
genes_N_SNPs <- XpX_SNPs_pseudo %>% group_by(gene_id) %>% tally()
genes_N_indels <- XpX_indels_pseudo %>% group_by(gene_id) %>% tally()
colnames(genes_N_SNPs) <- c('X_copy', 'N_pseudo_points')
colnames(genes_N_indels) <- c('X_copy', 'N_pseudo_indels')
genes_N_pseudo_mutations <- merge(genes_N_SNPs, genes_N_indels, by=c('X_copy'), all=T)
genes_N_pseudo_mutations$N_pseudo_points <- genes_N_pseudo_mutations$N_pseudo_points %>% replace_na(0)
genes_N_pseudo_mutations$N_pseudo_indels <- genes_N_pseudo_mutations$N_pseudo_indels %>% replace_na(0)
orthologs <- read.table('orthofinder_single_copy_assignments_no_t.tsv', header=T)
orthologs <- orthologs[,c(2,3)]
colnames(orthologs) <- c('Xp_copy', 'X_copy')
orthologs_N_pseudo_mutations <- merge(genes_N_pseudo_mutations, orthologs, by=c('X_copy'))
orthologs_N_pseudo_mutations <- orthologs_N_pseudo_mutations[,c(1,4,2,3)]
nrow(orthologs_N_pseudo_mutations) # 123
# write out list
orthologs_pseudo_mutations_Xp_list <- as.data.table(orthologs_N_pseudo_mutations[,c(2)])
write.table(orthologs_pseudo_mutations_Xp_list, file='list_of_Xp_genes_with_psuedo_mutations.tsv', row.names=F, col.names=F, quote=F, sep='\t')
```

Summary:
| Mutation type | N genes |
| - | - |
| Stop codon gained | 37 |
| Stop codon lost | 10 |
| Start codon lost | 10 |
| Frameshift | 113 |
| Frameshift and stop codon gained | 15 |
| Frameshift and stop codon lost | 6 |
| Frameshift and start codon lost | 8 |
| Frameshift and stop codon gained and stop codon lost | 2 |
| Frameshift and stop codon gained and start codon lost | 1 |

Overall, 123 genes contained at least one pseudogenizing mutation - 5.3% of the single-copy homologs.

## III. EXPRESSION (SILENCED GENES)

I can examine expression of X'- versus X-linked genes to identify genes where the X' copy is silenced. I'll do this with the stage-specific data generated for the Bcop_v2 genome paper (Urban et al. 2021): embryos (x3 reps), larva (x2), pupa (x2) and adult (x2). Females are pooled and are not separated by female type, so there should be roughly three X chromosomes for every one X' in the samples.

I did originally try to examine expression by simply mapping the RNAseq data to the Bcop_v3 CDS (i.e. including X and X' gene sequences), and this suggested that about 20% of genes on the X' are silenced. However, I found this difficult to trust because I know from other analyses I've done that there's a high level of mismapping between X and X', and X' reads should be more likely to mismap to X than vice-versa because the X' is more pooly assembled and contains more gaps compared to the X. For this reason, I'm applying an allele-specific expression pipeline (see https://github.com/MooHoll/Parent_of_Origin_Expression_Bumblebee) so that I can bin X and X' RNAseq reads prior to quantification. 

**Binning RNAseq reads**

1. Map X'X reads to Bcop_v2 (autosomes + X, X' reads get force-mapped to X). This is to identify X-X' variant positions along the X. I have done this before with John's data, but I'm repeating it with my data, which is higher coverage and may therefore allow for more sensitive detection of SNPs. I'm also using bowtie2 rather than BWA because the former has an option for high sensitivity alignments.

```
cat LR43_EDSW200011440-1a_HJ5JVDSXY_L2_1.trimmed_paired.fq.gz LR43_EDSW200011440-1a_HJ5FNDSXY_L2_1.trimmed_paired.fq.gz > XpX_heads_1.fq.gz
cat LR43_EDSW200011440-1a_HJ5JVDSXY_L2_2.trimmed_paired.fq.gz LR43_EDSW200011440-1a_HJ5FNDSXY_L2_2.trimmed_paired.fq.gz > XpX_heads_2.fq.gz

genome=Bcop_v2-chromosomes.fasta
genomename=$(basename $genome ".fasta")

# build bowtie index
bowtie2-build ${genome} ${genomename}.idx

for file in $(ls *_1.fq.gz)
do
	base=$(basename $file "_1.fq.gz")
    bowtie2 --very-sensitive-local --threads 20 \
    -x ${genomename}.idx \
    -1 ${base}_1.fq.gz -2 ${base}_2.fq.gz \
    | samtools view -bS - > ${base}.vs.${genomename}.bam && samtools sort ${base}.vs.${genomename}.bam > ${base}.vs.${genomename}.sorted.bam
done
```

2. Process bamfiles prior to variant calling

```
# index bams
for file in $(ls *sorted.bam)
do
samtools index ${file}
done

echo "adding read groups"
for file in $(ls *sorted.bam)
do
  	base=$(basename ${file} "_sorted.bam")
    picard \
    AddOrReplaceReadGroups \
    I=${file} \
    O=${base}_sorted_RG.bam \
    RGID=0001${base} \
    RGLB=lib1 \
    RGPL=illumina \
    RGPU=NA \
    RGSM=${base}
done

echo "deduplicating"
for file in $(ls *_sorted_RG.bam)
do
  	base=$(basename ${file} "_sorted_RG.bam")
    java -Xmx20g -jar /ceph/users/rbaird/software/picard/build/libs/picard.jar \
    MarkDuplicates \
    I=${file} \
    O=${base}_sorted_RG_deduplicated.bam \
    M=${base}_duplicate_metrics.txt \
    REMOVE_DUPLICATES=true
done
```

3. Call SNPs

```
genome=Bcop_v2-chromosomes.fasta

# index bams
for file in $(ls *_sorted_RG_deduplicated.bam)
do
    samtools index ${file}
done

# call SNPs
# min count 2 reads for alternative alleles, 
# min 5 reads per SNP, ignore complex events, indels and mnps
for file in $(ls *_sorted_RG_deduplicated.bam)
do
  	base=$(basename ${file} "_sorted_RG_deduplicated.bam")
    freebayes \
        -f ${genome} \
        -C 2 \
        --min-coverage 5 \
        -u \
        -i \
        -X \
        -b ${file} \
        > ${base}_vcf.gz
done

# filter SNPs
for file in $(ls *gz)
do
	base=$(basename ${file} "vcf.gz")
    vcftools --gzvcf ${file} --max-alleles 2 --minQ 20 --min-meanDP 10 \
    --recode --recode-INFO-all --out ${base}
done
```

4. Filtering informative SNPs. Let's think about this...
- The desired result is (i) to identify sites along the X where the X' position is variant so I can N-mask the genome for RNAseq mapping, and (ii) know what the X and X' bases are at all these positions.
- So, I want the heterozygous (0/1) SNP calls - and the alternative allele is the X' base whereas the reference allele is the X base.

```
# # filter the vcf for het sites only:
grep -e ^# -e "0/1" *.recode.vcf > XpX_alt_het_snps.vcf
# parse vcf file - require >=5 supporting reads:
python ~/software/genomics_general/VCF_processing/parseVCF.py -i XpX_alt_het_snps.vcf --minQual 30 --gtf flag=DP min=5 max=50 -o XpX_alt_het_snps.geno
```

SNPsplit will be a useful tool to bin RNAseq reads - it requires a file in the following format:
| ID | Chr | Position| SNP value | Ref/SNP |
| - | - | - | - | - |
| 18819008 | 5 | 48794752 | 1 | C/T |
| 40491905 | 11 | 63643453 | 1 | A/G |
| 44326884 | 12 | 96627819 | 1 | T/A |

The ID and SNP value fields aren't actually used for anything by SNPsplit, so they can be arbitrary values. SNPsplit will take alignment files and split them into alignments containing the ref and alt SNPs.

The output from parseVCF.py has the following format:
| #CHROM | POS | XpX_heads.vs.Bcop_v2-chromosomes.sorted.bam |
| - | - | - |
| X | 49104 | N/N |
| X | 49143 | N/N |
| X | 120800 | A/T |

The following R code will take the parseVCF.py output and generate a table that SNPsplit can use:
```
snp_file <- read.table('XpX_alt_het_snps.geno', stringsAsFactors=F, header=F)
snp_file$ID <- as.numeric(row.names(snp_file))
snp_file$SNP_value <- 1
colnames(snp_file) <- c('Chr', 'Position', 'Ref/SNP', 'SNP-ID', 'Strand')
snp_file <- snp_file[,c(4,1,2,5,3)]
write.table(snp_file, file='Xp_X_snp_file.txt', row.names=F, col.names=T, quote=F, sep='\t')
snp_file_X <- snp_file[which(snp_file$Chr == "X"),]
write.table(snp_file, file='Xp_X_snp_file_X_only.txt', row.names=F, col.names=T, quote=F, sep='\t')
```

5. N-masking the genome

Use the XpX_alt_het_snps.vcf file from above - i.e. variant positions along the X - to N-mask the Bcop_v2 genome:

```
bedtools maskfasta -fi Bcop_v2-chromosomes.fasta -bed XpX_alt_het_snps.vcf -fo Bcop_v2-chromosomes_Nmasked.fasta
```

6. Mapping RNAseq data to the N-masked genome

```
# run genomeGenerate
STAR \
--runThreadN 16 \
--runMode genomeGenerate \
--genomeSAindexNbases 12 \
--outFileNamePrefix Bcop_v2-chromosomes_Nmasked \
--genomeDir Bcop_v2-chromosomes_Nmasked.STAR \
--genomeFastaFiles Bcop_v2-chromosomes_Nmasked.fasta

# map RNAseq reads to N-masked genome.
echo "aligning RNAseq reads"
for file in $(ls *_1.trimmed.fq.gz)
do 
	base=$(basename $file "_1.trimmed.fq.gz")
	STAR \
	--runThreadN 16 \
	--alignEndsType EndToEnd \
	--outSAMattributes NH HI NM MD \
	--outSAMtype BAM Unsorted \
	--readFilesIn ${base}_1.trimmed.fq.gz ${base}_2.trimmed.fq.gz \
	--readFilesCommand zcat \
	--outTmpDir ${base}.out \
	--outFileNamePrefix ${base}.STAR. \
	--genomeDir Bcop_v2-chromosomes_Nmasked.STAR
done
```

6. Bin reads

Run SNPsplit on the RNAseq alignment files to bin the alignments into those containing the X and X' variants.

```
for file in $(ls *.out.bam)
do
	base=$(basename $file ".out.bam")
	SNPsplit \
	--paired \
	-o ${base}.out \
	--snp_file Xp_X_snp_file_X_only.txt \
	${base}.out.bam
done
```

Finally, use bamToFastq (bedtools) to convert these alignment files into reads:

```
for file in $(ls *.genome1.bam)
do
	base=$(basename ${file} ".STAR.Aligned.out.genome1.bam")
	bamToFastq -i ${base}.STAR.Aligned.out.genome1.bam -fq ${base}_X_reads.1.fq -fq2 ${base}_X_reads.2.fq
	gzip ${base}_X_reads.1.fq && gzip ${base}_X_reads.2.fq
done

for file in $(ls *.genome2.bam)
do
	base=$(basename ${file} ".STAR.Aligned.out.genome2.bam")
	bamToFastq -i ${base}.STAR.Aligned.out.genome2.bam -fq ${base}_Xp_reads.1.fq -fq2 ${base}_Xp_reads.2.fq
	gzip ${base}_Xp_reads.1.fq && gzip ${base}_Xp_reads.2.fq
done
```

Here's the proportion of reads assigned to each chromosome for each sample:

| Stage | Replicate | Prop reads assigned to X (%) | Prop reads assigned to X' (%) |
| - | - | - | - |
| Embryo | 1 | 5.96 | 4.05 |
| Embryo | 2 | 6.06 | 3.99 |
| Embryo | 3 | 6.25 | 4.04 |
| Larva | 1 | 7.49 | 2.33 |
| Larva | 2 | 7.31 | 2.36 |
| Pupa | 1 | 6.08 | 3.34 |
| Pupa | 2 | 6.25 | 3.42 |
| Adult | 1 | 7.22 | 2.31 |
| Adult | 2 | 7.62 | 2.44 |

This is more or less as expected. For larva, pupa and adult, there should be roughly 3 X copies for every 1 X' copy (because these are pooled XX and X'X females), and the proportion of reads assigned to the X is roughly in line with this (maybe not in pupa, can't think why). For embryos, these are 2h to 2-day-old samples, so bear in mind these are female embryos **laid by X'X mothers**. These samples may contain a decent proportion of maternal transcripts (we think ZGA probably occurs around 12h), which may explain why the proportion of reads assigned to the X' is more similar to the X. In other words, these embryo samples are more likly to be a mixture of pooled X'X females as well as pooled X'X and XX females, rather than just the latter (like in the other samples).

**Quantifying gene expression**

Now I have (i) X-X' single-copy homologs identified and (ii) RNAseq reads successfully binned as X or X' in origin, I can accurately quantify expression of these homologs.

```
# pull X-linked single-copy homologs to use as a reference
cat orthofinder_single_copy_assignments.tsv | tail -n 2321 | cut -f3 > X_single_copy_homologs.txt
seqtk subseq Xhom_longest_genes.fasta X_single_copy_homologs.txt > X_single_copy_homologs.fasta
# make kallisto index
kallisto index -i X_single_copy_homologs_kallisto.idx X_single_copy_homologs.fasta

# Run kallisto
for file in $(ls *_reads.1.fq.gz)
do
	base=$(basename ${file} "_reads.1.fq.gz")
	kallisto quant \
	--index=X_single_copy_homologs_kallisto.idx \
	-b 20 \
	--output-dir=${base}_out \
	--threads=16 \
	--plaintext ${base}_reads.1.fq.gz ${base}_reads.2.fq.gz \
	&& mv ${base}_out/abundance.tsv ${base}_vs_X_homs.tsv
done

# this loop modifies the .tsv files so I can read and merge them all into R at once painlessly:
for file in $(ls *.tsv)
do
	base=$(basename ${file} "_vs_X_homs.tsv")
	cat ${base}_vs_X_homs.tsv \
	| cut -f1,5 \
	| sed "s/tpm/${base}_tpm/" \
	> ${base}_vs_X_homs.modified.tsv
done
```

**Differential gene expression analysis**

Using EdgeR.

```
# Generate a counts table for EdgeR
setwd("/Users/robertbaird/Documents/analyses/functional_degradation/tsv_files_stage_data_X_vs_XX_no_bs")
all_files <- dir("/Users/robertbaird/Documents/analyses/functional_degradation/tsv_files_stage_data_X_vs_XX_no_bs")
file_names <- grep(all_files,pattern = "*.modified.tsv",value = TRUE)
Data_file <- map(file_names,read.delim, stringsAsFactors = FALSE, check.names = FALSE, row.names = NULL)
Merge_All_Samples <- Data_file %>% reduce(inner_join, by = "target_id")
# c1 = Xp, c2 = X... condition = chromosome the reads come from, because I care about comparing X vs X' expression
colnames(Merge_All_Samples) <- c('index', 'c1_r1', 'c2_r1', 'c1_r2', 'c2_r2', 'c1_r3', 'c2_r3', 'c1_r4', 'c2_r4', 'c1_r5', 'c2_r5', 'c1_r6', 'c2_r6', 'c1_r7', 'c2_r7', 'c1_r8', 'c2_r8', 'c1_r9', 'c2_r9')
write.table(Merge_All_Samples, file="counts.txt", append = FALSE, quote = FALSE, sep = "\t", row.names = F, col.names = T)

## READ IN COUNTS TABLE
counts <- read.delim('counts.txt', row.names='index')
counts <- counts[complete.cases(counts),]

## DEFINE GROUPS
group <- factor(c(rep(c("X", "Xp"), 9)),
                levels = c("Xp","X"))

## PREPARE DATA
counts_dgelist <- DGEList(counts=counts, group=group)

## FILTERING THE DATA - for now keep all counts because I'm interested in X' genes that are not expressed
dim(counts_dgelist)

## DATA VISUALISATION
# plot TPM distribution per sample
# if similar distributions, could justify using same arbitrary cut-off (e.g. TPM<1) for every sample...
counts_long <- gather(counts, sample, TPM, c1_r1:c2_r9)
ggplot(counts_long, aes(log10(TPM), colour=sample, group=sample)) + geom_density(adjust=2)

## NORMALISING THE DATA
counts_dgelist <- calcNormFactors(counts_dgelist)
head(counts_dgelist$counts)

# for each sample:
# 1. subset genes with zero TPM (as non-expressed)
# 2. subset genes with lowest 0.1% of non-zero TPMs (also as non-expressed, to allow for stochastic mismapping)
# 3. include only genes that are 'non-expressed' for all samples within a stage

# split by stage + chrom
embryo_X <- counts_dgelist$counts[,c(5,7,9)]
embryo_Xp <- counts_dgelist$counts[,c(6,8,10)]
larva_X <- counts_dgelist$counts[,c(11,13)]
larva_Xp <- counts_dgelist$counts[,c(12,14)]
pupa_X <- counts_dgelist$counts[,c(15,17)]
pupa_Xp <- counts_dgelist$counts[,c(16,18)]
adult_X <- counts_dgelist$counts[,c(1,3)]
adult_Xp <- counts_dgelist$counts[,c(2,4)]

# 1. count genes with TPM > 0 in each stage (change > to == to count genes with 0 counts across all samples for that stage)
nrow(embryo_X_zero_TPM <- embryo_X[rowSums(embryo_X[])>0,]) # 2046
nrow(embryo_Xp_zero_TPM <- embryo_Xp[rowSums(embryo_Xp[])>0,]) # 1961
nrow(larva_X_zero_TPM <- larva_X[rowSums(larva_X[])>0,]) # 2083
nrow(larva_Xp_zero_TPM <- larva_Xp[rowSums(larva_Xp[])>0,]) # 1927
nrow(pupa_X_zero_TPM <- pupa_X[rowSums(pupa_X[])>0,]) # 2098
nrow(pupa_Xp_zero_TPM <- pupa_Xp[rowSums(pupa_Xp[])>0,]) # 2031
nrow(adult_X_zero_TPM <- adult_X[rowSums(adult_X[])>0,]) # 2011
nrow(adult_Xp_zero_TPM <- adult_Xp[rowSums(adult_Xp[])>0,]) # 1849
```

Now I have counts for gene copies expressed at each stage. However, to account for stochastic mismapping of RNAseq reads, I also want to include genes that fall within the lowest 0.1% of non-zero TPMs as non-expressed too:

```
# subset remaining genes with <0.1% TPM within a sample and subtract these from non-zero TPMs
all_counts <- as.data.frame(counts_dgelist$counts)
all_counts$gene <- row.names(all_counts)
# This function will check if elements of each column are less than the 0.001 and quantile (0.1%) & keep rows where all columns lower than that quantile - apply it to all replicates
df <- all_counts[,c(5,19)] # subset sample
df <- df[rowSums(df[1])>0,] # exclude non-zero TPMs
keep <- Reduce(`&`, lapply(df[1], function(x) x <= quantile(x, 0.001)))
df[keep,] # 2
```

Here are the updated numbers:
| stage | chromosome | N gene expressed |
| - | - | - |
| Embryo | X | 2046 |
| Embryo | X' | 1961 |
| Larva | X | 2084 |
| Larva | X' | 1928 |
| Pupa | X | 2099 |
| Pupa | X' | 2031 |
| Adult | X | 2011 |
| Adult | X' | 1850 |


## IV. DOSAGE COMPENSATION

To look at potential dosage compensation of degraded X' genes, I will compare expression of X-linked genes between X'X and XX females. The only data I currently have for this is from early embryos (0-8h) for both sexes, where we are pretty certain ZGA has not begun and therefore these contain exclusively maternally-deposited transcripts (see supplementary materials for justification of this), i.e. we treat the male embryos as XX and the female embryos as X'X.

I can think of two ways to do this. The first way, which is more complicated, is the one I have included in the paper. The second approach is simpler but I only thought of it around the time I submitted the paper. It's likely I'll do it at some point as a complementary analysis, but I believe both approaches are reasonable.

The idea here is that X-linked genes that have degraded X' homologs, if dosage-compensated, should be upregulated in X'X females to match the expression of those X-linked genes in XX females. Here is how the two methods work:
- Method 1: Pull reads originating from the X from X'X and XX samples (using the same pipeline as in the silenced gene analysis). Quantify these expression of X-linked genes using these X reads. Conduct a differential gene expression analysis; DC'd genes should be upregulated in the X'X samples.
- Method 2: Map all RNAseq reads to Bcop_v2 (i.e. without the X' to force-map X' transcripts to X genes). We know genes where the X' copy is silenced, so look at these genes and compare expression between the X'X and XX samples. DC'd genes should have equal expression.

As I mentioned just now, I've only done method 1 for now. I won't repeat code below because it's essentially the same as above (RNAseq alignments, SNPsplit, Kallisto).

here are the proportions of reads for each file assigned to the X for each file. Here there are 12 files: 6 for each sex and 3 replicates within a sex for each time stage (0-4h and 4-8h):

| Sex | Stage | Replcatre | Prop reads assigned to X (%) |
| - | - | - | - |
| Female | 0-4h | 1 | 4.66 |
| Female | 0-4h | 2 | 4.71 |
| Female | 0-4h | 3 | 5.22 |
| Male | 0-4h | 1 | 9.36 |
| Male | 0-4h | 2 | 9.33 |
| Male | 0-4h | 3 | 9.24 |
| Female | 4-8h | 1 | 4.69 |
| Female | 4-8h | 2 | 6.24 |
| Female | 4-8h | 3 | 4.66 |
| Male | 4-8h | 1 | 8.83 |
| Male | 4-8h | 2 | 9.00 |
| Male | 4-8h | 3 | 8.87 |

This makes sense - XX females (i.e. male embryos) have twice as many X copies as X'X females, so we expect twice as many reads to be assigned to the X. The only exception is FL2 (i.e. female 4-8h), which appears to be an outlier, so I will exclude it from further analysis.

DGE code:

```
# Generate a counts table for EdgeR
setwd("/Users/robertbaird/Documents/analyses/dosage_compensation/tsv_files_embryo_data_X_vs_XX_no_bs")
all_files <- dir("/Users/robertbaird/Documents/analyses/dosage_compensation/tsv_files_embryo_data_X_vs_XX_no_bs")
file_names <- grep(all_files,pattern = "*.modified.tsv",value = TRUE)
Data_file <- map(file_names,read.delim, stringsAsFactors = FALSE, check.names = FALSE, row.names = NULL)
Merge_All_Samples <- Data_file %>% reduce(inner_join, by = "target_id")
head(Merge_All_Samples)
# Remove FL2 because it's an outlier
Merge_All_Samples <- Merge_All_Samples[,c(1:5,7:13)]
head(Merge_All_Samples)
# c1 = XX, c2 = XpX... here condition = genotype the sample comes from (number of X chromosomes)
colnames(Merge_All_Samples) <- c('index', 'c1_r1', 'c1_r2', 'c1_r3', 'c1_r4', 'c1_r5', 'c2_r1', 'c2_r2', 'c2_r3', 'c2_r4', 'c2_r5', 'c2_r6')
write.table(Merge_All_Samples, file="counts.txt", append = FALSE, quote = FALSE, sep = "\t", row.names = F, col.names = T)

## READ IN COUNTS TABLE
counts <- read.delim('counts.txt', row.names='index')
counts <- counts[complete.cases(counts),]

## DEFINE GROUPS
group <- factor(c('XpX', 'XpX', 'XpX', 'XpX', 'XpX', 'XX', 'XX', 'XX', 'XX', 'XX', 'XX'),
                levels = c("XpX","XX"))

## PREPARE DATA
counts_dgelist <- DGEList(counts=counts, group=group)

## FILTERING THE DATA - for now keep all counts because I'm interested in X' genes that are not expressed
dim(counts_dgelist)
head(cpm(counts_dgelist)) # calculate counts per million
apply(counts_dgelist$counts, 2, sum) # total gene counts per sample
keep <- rowSums(cpm(counts_dgelist)>100) >= 2 # genes must have at least 100 counts per million
counts_dgelist <- counts_dgelist[keep,]
dim(counts_dgelist) # 1128/2330 kept
counts_dgelist$samples$lib.size <- colSums(counts_dgelist$counts)
counts_dgelist$samples

## NORMALISING THE DATA
counts_dgelist <- calcNormFactors(counts_dgelist)
head(counts_dgelist)

## DATA EXPLORATION / ESTIMATING DISPERSION
plotMDS(counts_dgelist, method="bcv", col=as.numeric(counts_dgelist$samples$group))
legend("topright", as.character(unique(counts_dgelist$samples$group)), col=1:2, pch=20)
# 0-4h samples cluster together; F + M 4-8h samples cluster separately...  maybe there's some zygotic gene expression there then?

# estimating dispersion
counts_dgelist1 <- estimateCommonDisp(counts_dgelist, verbose=T)
names(counts_dgelist1)
# using empirical Bayes tagwise dispersions
counts_dgelist1 <- estimateTagwiseDisp(counts_dgelist1)
names(counts_dgelist1)
#  plots the tagwise biological coefficient of variation (square root of dispersions) against log2-CPM.
plotBCV(counts_dgelist1)

## DIFFERENTIAL EXPRESSION

# The function exactTest() conducts tagwise tests using the exact negative binomial test
et12 <- exactTest(counts_dgelist1, pair=c(1,2)) # compare groups 1 and 2 i.e. XpX vs XX
topTags(et12, n=10)
# Total number of differentially expressed genes at FDR < 0:05:
de1 <- decideTestsDGE(et12, adjust.method="BH", p.value=0.05)
summary(de1)
# Down: 4; NotSig: 1124, Up: 0
# In de1, 1 = significantly up-regulated in XX sample (0 genes), -1 = significantly upregulated in XpX sample (4 genes)

# the genes:
# jg4857.t1
# jg5988.t1
# jg6324.t1
# jg7315.t1
# Look at their functionality: grep jg5988 single_copy_homologs_functionality.tsv

# plotSmear generates a plot of the tagwise log-fold-changes against log-cpm
de1tags12 <- rownames(counts_dgelist1)[as.logical(de1)]
pdf(file="DC_smear.pdf", width=4, height=4)
plotSmear(et12, de.tags=de1tags12, cex=0.8, pch=1, col="black")
abline(h = c(-0.5, 0.5), col = c("dodgerblue"),  lty=c(2))
dev.off()
```

## IV. TE ACCUMULATION

Analysis of transposable elements was outsourced to Malte Grewoldt.


