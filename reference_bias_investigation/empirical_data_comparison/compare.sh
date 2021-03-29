#! /bin/bash 
# pipeline to examine basecall accuracy of RapUp, Snippy and Gon_phyling (de novo assembly) programs

outdir=${1}
ref_1=${2}
ref_2=${3}
rapup_path=${4}
loci_positions_1=${5}
loci_positions_2=${6}
program_path=${7}

mkdir ${outdir}/${ref_1}_${ref_2}_comparison_output
mkdir ${outdir}/${ref_1}_${ref_2}_comparison_output/locus_alignments
mkdir ${outdir}/${ref_1}_${ref_2}_comparison_output/${ref_1}_individ_taxa_loci
mkdir ${outdir}/${ref_1}_${ref_2}_comparison_output/${ref_2}_individ_taxa_loci

${rapup_path}/modules/locus_position_identifier.py \
--out_file_dir  ${outdir}/${ref_1}_${ref_2}_comparison_output/${ref_1}_locus_msa_files \
--position_csv_file ${loci_positions_1} \
--concatenated_fasta ${outdir}/${ref_1}_rapup_run/combine_and_infer/extended.aln

${rapup_path}/modules/locus_position_identifier.py \
--out_file_dir  ${outdir}/${ref_1}_${ref_2}_comparison_output/${ref_2}_locus_msa_files \
--position_csv_file ${loci_positions_2} \
--concatenated_fasta ${outdir}/${ref_2}_rapup_run/combine_and_infer/extended.aln


for i in $(ls ${outdir}/${ref_1}_${ref_2}_comparison_output/${ref_1}_locus_msa_files )
do

    cluster_name=$(basename ${i} | sed 's/.fasta//g')
    grep ">" ${outdir}/${ref_1}_${ref_2}_comparison_output/${ref_1}_locus_msa_files/${i} | sed 's/>//g' > ${outdir}/${ref_1}_${ref_2}_comparison_output/taxa_list.txt

    for j in $(cat ${outdir}/${ref_1}_${ref_2}_comparison_output/taxa_list.txt);
    do

        ${rapup_path}/modules/ref_producer.py -r --align_file ${outdir}/${ref_1}_${ref_2}_comparison_output/${ref_1}_locus_msa_files/${i} \
	    --out_file ${outdir}/${ref_1}_${ref_2}_comparison_output/${ref_1}_individ_taxa_loci/single_tax-${j}-${cluster_name}

    done

done

for i in $(ls ${outdir}/${ref_1}_${ref_2}_comparison_output/${ref_2}_locus_msa_files )
do

    cluster_name=$(basename ${i} | sed 's/.fasta//g')
    grep ">" ${outdir}/${ref_1}_${ref_2}_comparison_output/${ref_2}_locus_msa_files/${i} | sed 's/>//g' > ${outdir}/${ref_1}_${ref_2}_comparison_output/taxa_list.txt

    for j in $(cat ${outdir}/${ref_1}_${ref_2}_comparison_output/taxa_list.txt);
    do

        ${rapup_path}/modules/ref_producer.py -r --align_file ${outdir}/${ref_1}_${ref_2}_comparison_output/${ref_2}_locus_msa_files/${i} \
	    --out_file ${outdir}/${ref_1}_${ref_2}_comparison_output/${ref_2}_individ_taxa_loci/single_tax-${j}-${cluster_name}

    done

done

${program_path}/reference_bias_investigation/empirical_data_comparison/locus_match_and_combiner.py \
		--long_seqs_folder ${outdir}/${ref_1}_${ref_2}_comparison_output/${ref_1}_individ_taxa_loci \
		--manipulate_seqs_folder ${outdir}/${ref_1}_${ref_2}_comparison_output/${ref_1}_individ_taxa_loci \
		--output_dir ${outdir}/${ref_1}_${ref_2}_comparison_output/locus_alignments


# for j in $(ls -1 ${outdir}/${ref_1}_${ref_2}_comparison_output/${ref_1}_individ_taxa_loci/);
# 	do
# 		# locus=$(basename ${j} | sed -e 's/single_tax-//g')
# 		# taxon=$(basename ${j} | rev | cut -d "-" -f1 | rev)
# 		${program_path}/reference_bias_investigation/empirical_data_comparison/locus_match_and_combiner.py \
# 		--long_seqs_folder ${outdir}/${ref_1}_${ref_2}_comparison_output/${ref_1}_individ_taxa_loci \
# 		--manipulate_seqs_folder ${outdir}/${ref_1}_${ref_2}_comparison_output/${ref_1}_individ_taxa_loci \
# 		--output_dir ${outdir}/${ref_1}_${ref_2}_comparison_output/locus_alignments

# 	done
