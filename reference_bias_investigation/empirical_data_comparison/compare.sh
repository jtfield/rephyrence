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
mkdir ${outdir}/${ref_1}_${ref_2}_comparison_output/alignments_for_comparison
mkdir ${outdir}/${ref_1}_${ref_2}_comparison_output/alignment_comparison_results

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
		--output_dir ${outdir}/${ref_1}_${ref_2}_comparison_output/locus_alignments \
        --manip_ref_name ${ref_1} \
        --long_ref_name ${ref_2}

echo "Aligning loci with Mafft before performing comparison"

for i in $(ls ${outdir}/${ref_1}_${ref_2}_comparison_output/locus_alignments);
do
	mafft \
	--thread 2 \
	--op 5 \
	--lexp -0.5 \
	${outdir}/${ref_1}_${ref_2}_comparison_output/locus_alignments/${i} \
	> ${outdir}/${ref_1}_${ref_2}_comparison_output/alignments_for_comparison/aligned_${i}

	${program_path}/reference_bias_investigation/empirical_data_comparison/emp_snippy_gapped_align_compared.py \
	-t \
	--align_1 ${outdir}/${ref_1}_${ref_2}_comparison_output/alignments_for_comparison/aligned_${i} \
	--output_stub assessment_${i} \
	--output_dir ${outdir}/${ref_1}_${ref_2}_comparison_output/alignment_comparison_results
done