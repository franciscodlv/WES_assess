#!/bin/bash
# Script to facilitate the clinical assesment of capture assays for NGS WES
# v1, 7/14/2020
# (c) 2020 by Francisco M. De La Vega, Stanford University
# All rights reserved.
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the  
# Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, 
# and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR 
# ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH 
# THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# Prerequisites:
# RTG core software
# Bedtools
#
SCRIPT_DIR=$(dirname "$0")  # The directory where this script resides, so we can look for local files here.
### Pointer to the executable and bedtools exectuables, required, this sets a default if not already defined in the environment
# change as needed
RTG=/Users/francisco.delavega/workdir/rtg_distributions/rtg-core-3.12.1/rtg
bedtools=/usr/local/bin/bedtools
###Locations of directory for fixed input files (change as needed)
#Human gene lists - exon BEDs - currently fixed and assumes properly prepared and with "Chr" interval prefix
exons_dir=/Users/delavefm/workdir/Human_Genes
#ClinVar location BEDs/usr/local/bin/bedtools
clinvar_dir=/Users/francisco.delavega/workdir/ClinVar
# Genome reference SDF for RTG analyzes 
template=/Users/francisco.delavega/workdir/data/references/human_g1k_v37_decoy

### Process Command line arguments
if [ "$#" -lt 1 ]; then
    cat <<EOF
USAGE: benchmark_wes_assays.sh /path_to_bam/sample.bam /path_to_assay_bed/assay_roi_hg19.bed 
e.g:
    benchmark_wes_assays.sh ../SEC/NA12878.bam ../BED/Assay_exome_targets_sorted_hg19.bed
	Note: BAM file should be indexed and BED file should be sorted and start with "assayname_"
	It is assumed assay BED file used "Chr" prefix for intevals.
	Execute this script on the directory where you want results to ouput.
EOF
    exit 1
fi
### Define variables for input files provided in the the command line
bam=$1
roi_bed=$2
THIS_DIR=$PWD
#Parse the assay name from ROI BED prefix (required with underscore separator)
assay=$(echo "$roi_bed" | awk 'BEGIN{FS="/"}{print $NF}' | awk -F'[_]' '{print $1}')

#create results subdirs under dirname
echo "Creating ouput directories"
mkdir ${THIS_DIR}/BED
mkdir ${THIS_DIR}/Clinvar

#Create a version of ROI BED files sorted & clean from garbage that does not have "Chr" prefix for use with RTG tools
cat ${roi_bed} | awk -F'\t' -v OFS="\t" '{ print $1, $2, $3 }' | ${bedtools} sort -i stdin | ${bedtools} merge | sed -e 's/chr//g' | sed -e '/7_gl000195_random/d' | sed -e '/Un_gl000228/d' | sed -e '/chr6_apd_hap1/d' | sed -e '/-AS/d' | tr -d $'\r' > ${THIS_DIR}/BED/${assay}_roi_bed_noChr.bed || exit 1
roi_bed_nochr=${THIS_DIR}/BED/${assay}_roi_bed_noChr.bed
#Create a version of ROI BED files sorted & clean of common garbage but with Chr to compare with medical genes BEDS with Chr
cat ${roi_bed} | awk -F'\t' -v OFS="\t" '{ print $1, $2, $3 }' | ${bedtools} sort -i stdin | ${bedtools} merge | sed -e '/7_gl000195_random/d' | sed -e '/Un_gl000228/d' | sed -e '/chr6_apd_hap1/d' | sed -e '/-AS/d' | tr -d $'\r' > ${THIS_DIR}/BED/${assay}_roi_bed_Chr.bed || exit 1
roi_bed_chr=${THIS_DIR}/BED/${assay}_roi_bed_Chr.bed 

### Intersects and counts of exons in gene lists
echo "Calculating intersections between assay ROI and exons from medical gene lists"
${bedtools} intersect -a ${exons_dir}/UCSC_Human_exons_medical_genes_hg19.bed -b ${roi_bed_chr} | sed -e 's/chr//g' > ${THIS_DIR}/BED/UCSC_Human_exons_medical_genes_in_${assay}_nochr.bed || exit 1
wc1="$(wc -l ${exons_dir}/UCSC_Human_exons_medical_genes_hg19.bed | awk '{print $1;}')"
wc2="$(wc -l ${THIS_DIR}/BED/UCSC_Human_exons_medical_genes_in_${assay}_nochr.bed | awk '{print $1;}')"
fraction=`echo "$wc2/$wc1" | bc -l`
echo " Ratio between Medical exons in ${assay}.bed (n=${wc2}) and total Medical exons list (n=${wc1}) = $fraction" > ${THIS_DIR}/BED/exon_overlap_counts.txt
echo " Ratio between Medical exons in ${assay}.bed (n=${wc2}) and total Medical exons list (n=${wc1}) = $fraction"

${bedtools}  intersect -a ${exons_dir}/UCSC_Human_exons_medical_genes_ngs_problem_ls_hg19.bed -b ${roi_bed_chr} | sed -e 's/chr//g' > ${THIS_DIR}/BED/UCSC_Human_exons_medical_genes_ngs_problem_ls_in_${assay}_nochr.bed || exit 1
wc1="$(wc -l ${exons_dir}/UCSC_Human_exons_medical_genes_ngs_problem_ls_hg19.bed | awk '{print $1;}')"
wc2="$(wc -l ${THIS_DIR}/BED/UCSC_Human_exons_medical_genes_ngs_problem_ls_in_${assay}_nochr.bed | awk '{print $1;}')"
fraction=`echo "$wc2/$wc1" | bc -l`
echo " Ratio between low stringency problem Medical exons in ${assay}.bed (n=${wc2}) and total exons in low stringency problem list (n=${wc1}) = $fraction" >> ${THIS_DIR}/BED/exon_overlap_counts.txt
echo " Ratio between low stringency problem Medical exons in ${assay}.bed (n=${wc2}) and total exons in low stringency problem list (n=${wc1}) = $fraction"

${bedtools}  intersect -a ${exons_dir}/UCSC_Human_exons_medical_genes_ngs_problem_hs_hg19.bed -b ${roi_bed_chr} | sed -e 's/chr//g' > ${THIS_DIR}/BED/UCSC_Human_exons_medical_genes_ngs_problem_hs_in_${assay}_nochr.bed || exit 1
wc1="$(wc -l ${exons_dir}/UCSC_Human_exons_medical_genes_ngs_problem_hs_hg19.bed | awk '{print $1;}')"
wc2="$(wc -l ${THIS_DIR}/BED/UCSC_Human_exons_medical_genes_ngs_problem_hs_in_${assay}_nochr.bed | awk '{print $1;}')"
fraction=`echo "$wc2/$wc1" | bc -l`
echo " Ratio between high stringency problem Medical exons in ${assay}.bed (n=${wc2}) and total exons in high stringency problem list (n=${wc1}) = $fraction" >> ${THIS_DIR}/BED/exon_overlap_counts.txt
echo " Ratio between high stringency problem Medical exons in ${assay}.bed (n=${wc2}) and total exons in high stringency problem list (n=${wc1}) = $fraction"

${bedtools}  intersect -a ${exons_dir}/UCSC_Human_exons_medical_genes_ngs_dead_zone_hg19.bed -b ${roi_bed_chr} | sed -e 's/chr//g' > ${THIS_DIR}/BED/UCSC_Human_exons_medical_genes_ngs_dead_zone_in_${assay}_nochr.bed || exit 1
wc1="$(wc -l ${exons_dir}/UCSC_Human_exons_medical_genes_ngs_dead_zone_hg19.bed | awk '{print $1;}')"
wc2="$(wc -l ${THIS_DIR}/BED/UCSC_Human_exons_medical_genes_ngs_dead_zone_in_${assay}_nochr.bed | awk '{print $1;}')"
fraction=`echo "$wc2/$wc1" | bc -l`
echo " Ratio between dead zone Medical exons in ${assay}.bed (n=${wc2}) and total dead zone exons in list (n=${wc1}) = $fraction" >> ${THIS_DIR}/BED/exon_overlap_counts.txt
echo " Ratio between dead zone Medical exons in ${assay}.bed (n=${wc2}) and total dead zone exons in list (n=${wc1}) = $fraction"

### Calculate overlaps between exons of Medical gene lists and assay equivalents for figures anmd histograms
# Results are saved as BED files that can be opened by Excel for analysis there
# Calculate overalps, incuding exons that are not covered in Assay ROI (these will have zero overlaps to be counted in the analysis)
echo "Now calculating bp overlaps between medical gene exons and assay ROIs for further analysis"
echo "overlaps for medical genes"
${bedtools}  intersect -a ${exons_dir}/UCSC_Human_exons_medical_genes_hg19.bed -b ${roi_bed_chr} -loj |  awk -f ${SCRIPT_DIR}/overlap.awk | sed -e 's/chr//g' | tr -d $'\r' > ${THIS_DIR}/BED/UCSC_Human_exons_medical_genes_in_${assay}_overlap.bed || exit 1
echo "overlaps for medical genes in ngs dead zone"
${bedtools}  intersect -a ${exons_dir}/UCSC_Human_exons_medical_genes_ngs_dead_zone_hg19.bed -b ${roi_bed_chr} -loj | awk -f ${SCRIPT_DIR}/overlap.awk | sed -e 's/chr//g' | tr -d $'\r' > ${THIS_DIR}/BED/UCSC_Human_exons_medical_genes_ngs_dead_zone_in_${assay}_overlap_nochr.bed || exit 1
echo "overlaps for medical genes in hs list"
${bedtools}  intersect -a ${exons_dir}/UCSC_Human_exons_medical_genes_ngs_problem_hs_hg19.bed -b ${roi_bed_chr} -loj | awk -f ${SCRIPT_DIR}/overlap.awk | sed -e 's/chr//g' | tr -d $'\r' > ${THIS_DIR}/BED/UCSC_Human_exons_medical_genes_ngs_problem_hs_in_${assay}_overlap_nochr.bed || exit 1
echo "overlaps for medical genes in ls list"
${bedtools}  intersect -a ${exons_dir}/UCSC_Human_exons_medical_genes_ngs_problem_ls_hg19.bed -b ${roi_bed_chr} -loj | awk -f ${SCRIPT_DIR}/overlap.awk | sed -e 's/chr//g' | tr -d $'\r' > ${THIS_DIR}/BED/UCSC_Human_exons_medical_genes_ngs_problem_ls_in_${assay}_overlap_nochr.bed || exit 1

# Gene lists column swaping the WES assay interval cooridnated of intersecting regions to estimate coverage with the Assay BED intervals
# colappsing record with identical cooridnates resulted from intersect to avpid redundancies, combining annotations into a single column to allow gropupby function
# and remove "Chr" prefix - these  files are only for use in RTG coverage
echo "Prep BED files for coverage analysis using Assay interval coordinates of overlapping exons"
echo "Medical gene overalp"
${bedtools} intersect -a ${exons_dir}/UCSC_Human_exons_medical_genes_hg19.bed -b ${roi_bed_chr} -wb | ${bedtools} sort -i stdin | awk -F'\t' '{ print $1"\t"$6"\t"$7"\t"$4","$5","$2","$3 }' - | bedtools groupby -g 1,2,3 -c 4 -o first | sed -e 's/chr//g' | tr -d $'\r' > ${THIS_DIR}/BED/${assay}_medical_genes_hg19_rev_nochr.bed || exit 1
echo "High stringency list overlaps"
${bedtools} intersect -a ${exons_dir}/UCSC_Human_exons_medical_genes_ngs_problem_hs_hg19.bed -b ${roi_bed_chr} -wb | ${bedtools} sort -i stdin | awk -F'\t' '{ print $1"\t"$6"\t"$7"\t"$4","$5","$2","$3 }' - | bedtools groupby -g 1,2,3 -c 4 -o first | sed -e 's/chr//g' | tr -d $'\r' > ${THIS_DIR}/BED/${assay}_medical_genes_hs_NGSproblem_hg19_rev_nochr.bed || exit 1
echo "Low stringency list overlaps"
${bedtools} intersect -a ${exons_dir}/UCSC_Human_exons_medical_genes_ngs_problem_ls_hg19.bed -b ${roi_bed_chr} -wb | ${bedtools} sort -i stdin | awk -F'\t' '{ print $1"\t"$6"\t"$7"\t"$4","$5","$2","$3 }' - | bedtools groupby -g 1,2,3 -c 4 -o first | sed -e 's/chr//g' | tr -d $'\r' > ${THIS_DIR}/BED/${assay}_medical_genes_ls_NGSproblem_hg19_rev_nochr.bed || exit 1
echo "NGS dead zone list overlaps"
${bedtools} intersect -a ${exons_dir}/UCSC_Human_exons_medical_genes_ngs_dead_zone_hg19.bed -b ${roi_bed_chr} -wb | ${bedtools} sort -i stdin | awk -F'\t' '{ print $1"\t"$6"\t"$7"\t"$4","$5","$2","$3 }' - | bedtools groupby -g 1,2,3 -c 4 -o first | sed -e 's/chr//g' | tr -d $'\r' > ${THIS_DIR}/BED/${assay}_medical_genes_dead_zone_NGSproblem_hg19_rev_nochr.bed || exit 1

####Analysis of coverage across medical genes by Assay
#Calculate coverage across medical lists exons from input BAM using rtg coverage
#all regions
echo "Calculating coverage aross all Assay ROI regions from alignment data"
${RTG} coverage --template=${template} --bed-regions=${roi_bed_nochr} --per-region -o ${THIS_DIR}/BED/${assay}_targets_noChr_hg19_coverage ${bam} || exit 1

echo "Analyze coverage across medical exons (original canonical transcript exon cooridnates)"
echo "Calculating coverage across Medical genes exons"
${RTG} coverage --template=${template} --bed-regions=${THIS_DIR}/BED/UCSC_Human_exons_medical_genes_in_${assay}_nochr.bed --per-region -o ${THIS_DIR}/BED/UCSC_Human_exons_medical_genes_in_${assay}_coverage ${bam} || exit 1

echo "Calculating coverage across medical Genes exons in NGS dead zone"
${RTG} coverage --template=${template} --bed-regions=${THIS_DIR}/BED/UCSC_Human_exons_medical_genes_ngs_dead_zone_in_${assay}_nochr.bed --per-region -o ${THIS_DIR}/BED/UCSC_Human_exons_medical_genes_ngs_dead_zone_in_${assay}_coverage ${bam} || exit 1

echo "Calculating coverage across medical Genes exons in high stringency NGS problem zone"
${RTG} coverage --template=${template} --bed-regions=${THIS_DIR}/BED/UCSC_Human_exons_medical_genes_ngs_problem_hs_in_${assay}_nochr.bed --per-region -o ${THIS_DIR}/BED/UCSC_Human_exons_medical_genes_ngs_problem_hs_in_${assay}_coverage ${bam} || exit 1

echo "Calculating coverage across medical Genes exons in low stringency NGS problem zone"
${RTG} coverage --template=${template} --bed-regions=${THIS_DIR}/BED/UCSC_Human_exons_medical_genes_ngs_problem_ls_in_${assay}_nochr.bed --per-region -o ${THIS_DIR}/BED/UCSC_Human_exons_medical_genes_ngs_problem_ls_in_${assay}_coverage ${bam} || exit 1

### Anayiss of coverage in WES assay target regions for Medical Genes (reversed targets vs above)
#Across Medical genes exons but only those alerady included in Assay - this is lenient for assay
#Could be over counting some segments because it is not unique? Not sure what RTG coverage does with repeted intrevals
echo "Analysis of coverage for medical gene exons but within overalping exon using Twist ROI coordinates"
echo "Calculating coverage across Medical genes exons"
${RTG} coverage --template=${template} --bed-regions=${THIS_DIR}/BED/${assay}_medical_genes_hg19_rev_nochr.bed --per-region -o ${THIS_DIR}/BED/${assay}_targets_medical_genes_coverage ${bam} || exit 1

echo "Calculating coverage across medical Genes exons in NGS dead zone, assay intervals"
${RTG} coverage --template=${template} --bed-regions=${THIS_DIR}/BED/${assay}_medical_genes_dead_zone_NGSproblem_hg19_rev_nochr.bed --per-region -o ${THIS_DIR}/BED/${assay}_targets_ngs_dead_zone_coverage ${bam} || exit 1

echo "Calculating coverage across medical Genes exons in high stringency NGS problem zone, assay intervals"
${RTG} coverage --template=${template} --bed-regions=${THIS_DIR}/BED/${assay}_medical_genes_hs_NGSproblem_hg19_rev_nochr.bed --per-region -o ${THIS_DIR}/BED/${assay}_targets_medical_genes_ngs_problem_hs_coverage ${bam} || exit 1

echo "Calculating coverage across medical Genes exons in low stringency NGS problem zone, assay intervals"
${RTG} coverage --template=${template} --bed-regions=${THIS_DIR}/BED/${assay}_medical_genes_ls_NGSproblem_hg19_rev_nochr.bed --per-region -o ${THIS_DIR}/BED/${assay}_targets_medical_genes_ngs_problem_ls_coverage ${bam} || exit 1

### Call variants for ClinVar loci analyses 
bamfile=$(basename "$bam")
bampath=$(dirname "$bam")
if [ -f ${bampath}/${bamfile}.calibration ]; then
   echo "Calibration file exists."
else
   echo "Calibrating BAM to call variants"
   ${RTG} calibrate --bed-regions=${roi_bed_nochr} -t ${template} ${bam}
fi

echo "Calling SNVs in all Clinvar positions emmiting all genotypes"
${RTG} snp --bed-regions=${clinvar_dir}/clinvar_20200520_vcf_merged_snv_indel.bed --avr-model=illumina-exome.avr -a -o ${THIS_DIR}/Clinvar/clinvar_all_calls -t ${template} ${bam} || exit 1

echo "Generate tables for histograms"
echo "Generating table for depth histogram from calls"
${RTG} bgzip -d -c ${THIS_DIR}/Clinvar/clinvar_all_calls/snps.vcf.gz | egrep -v "^#" - | cut -f 8 | sed 's/^.*;DP=\([0-9]*\);.*$/\1/' > ${THIS_DIR}/Clinvar/clinvar_all_calls/depth.txt || exit 1
echo "Generating table for Quality Histogram from calls"
${RTG} bgzip -d -c ${THIS_DIR}/Clinvar/clinvar_all_calls/snps.vcf.gz | egrep -v "^#" - | cut -f 6 |  sed 's/^\./0/' > ${THIS_DIR}/Clinvar/clinvar_all_calls/qual.txt || exit 1

### Metrics for ClinVar positions to ONLY include sites on Twist ROI
#Create overlap BED between all snv+indels in ClinVar and Twist ROI BED:
echo "Calculate coverage metrics only for ClinVar positions in Assay"
echo "intesrecting variants with Assay ROI intervals"
${bedtools}  intersect -a ${clinvar_dir}/clinvar_20200520_vcf_merged_snv_indel.bed -b ${roi_bed_nochr} > ${THIS_DIR}/BED/clinvar_20200520_vcf_snv_indel_merged_in_${assay}_nochr.bed || exit 1

echo "Calculating coverage wIth BED restricted to Assay ovelap (SNV and indels)"
${RTG} coverage --template=${template} --bed-regions=${THIS_DIR}/BED/clinvar_20200520_vcf_snv_indel_merged_in_${assay}_nochr.bed --per-region -o ${THIS_DIR}/Clinvar/Twist_Exome_RefSeq_ClinVar_snv_indels_coverage_in_${assay} ${bam} || exit 1

echo "Calculating coverage with variant calls filtered to intersect region"
${RTG} vcffilter -i ${THIS_DIR}/Clinvar/clinvar_all_calls/snps.vcf.gz -o ${THIS_DIR}/Clinvar/clinvar_all_calls/clinvar_all_in_${assay}_snps_indel.vcf.gz --bed-regions=${THIS_DIR}/BED/clinvar_20200520_vcf_snv_indel_merged_in_${assay}_nochr.bed|| exit 1

echo "Generate tables for histograms"
echo "Generating table for depth histogram from calls"
${RTG} bgzip -d -c  ${THIS_DIR}/Clinvar/clinvar_all_calls/clinvar_all_in_${assay}_snps_indel.vcf.gz | egrep -v "^#" - | cut -f 8 | sed 's/^DP=\([0-9]*\);.*$/\1/' > ${THIS_DIR}/Clinvar/clinvar_all_calls/${assay}_regions_depth.txt || exit 1

echo "Generating table for Quality Histogram from calls"
${RTG} bgzip -d -c  ${THIS_DIR}/Clinvar/clinvar_all_calls/clinvar_all_in_${assay}_snps_indel.vcf.gz | egrep -v "^#" - | cut -f 6 |  sed 's/^\./0/' > ${THIS_DIR}/Clinvar/clinvar_all_calls/${assay}_regions_qual.txt || exit 1

echo "DONE"

