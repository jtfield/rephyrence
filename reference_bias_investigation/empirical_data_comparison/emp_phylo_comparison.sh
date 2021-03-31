#! /bin/bash 
# pipeline to examine basecall accuracy of RapUp, Snippy and Gon_phyling (de novo assembly) programs

source $1

# Function
# Accepts inputs from config file and assembles sequences using EP based on a given reference taxon
EP_run () {

	outdir=${1}
	snip_ref=${2}
	rapup_path=${3}
	r1_tail=${4}
	r2_tail=${5}
	loci_positions=${6}
	bootstrapping=${7}
	intermediate_files=${8}


	printf "\nBeginning small parsnp run for rapup reference.\n"
	ref_name=$(basename ${snip_ref} | sed 's/.fasta//g')

	echo "${ref_name}"

	mkdir ${outdir}/${ref_name}_short_ref

	mv ${outdir}/read_copies/${ref_name}${r1_tail} ${outdir}/${ref_name}_short_ref
	mv ${outdir}/read_copies/${ref_name}${r2_tail}	${outdir}/${ref_name}_short_ref

	for i in $(cat ${outdir}/taxa_list.txt | head -5);
	do
		ln -s $i ${outdir}/${ref_name}_short_ref/
		ln -s ${i%$r1_tail}$r2_tail ${outdir}/${ref_name}_short_ref/
	done

	

	#CALL GON_PHYLING ON ALL TAXA READS
	if [ $intermediate_files == "KEEP" ]; then
		${rapup_path}/gon_phyling.sh -d ${outdir}/${ref_name}_short_ref -1 $r1_tail -2 $r2_tail -o LOCUS -l $loci_positions -b $bootstrapping
	elif [ $intermediate_files == "CLEAN" ]; then
		${rapup_path}/gon_phyling.sh -d ${outdir}/${ref_name}_short_ref -1 $r1_tail -2 $r2_tail -o LOCUS -l $loci_positions -i $intermediate_files -b $bootstrapping
	fi

	cd ${outdir}/${ref_name}_short_ref

	${rapup_path}/modules/ref_producer.py --align_file ${outdir}/${ref_name}_short_ref/trimmed_reads/spades_output/genomes_for_parsnp/alignment_fixing/combo.fas \
	--ref_select ${ref_name} \
	--out_file ${outdir}/${ref_name}.fasta \
	-s

	${rapup_path}/multi_map.sh \
	-a ${outdir}/${ref_name}.fasta \
	-r ${ref_name} \
	-d ${read_copies_pwd} \
	-1 ${r1_tail} -2 ${r2_tail} \
	-o ${outdir}/${ref_name}_rapup_run \
	-g LOCUS \
	-f positional.csv \
	-i ${intermediate_files} \
	-b ${bootstrapping}

	cd ${outdir}/${ref_name}_rapup_run/combine_and_infer

	# sed -i 's/.fasta//g' extended.aln
	# sed -i 's/.ref//g' extended.aln

	#SET UP REFERENCE BY PULLING A SINGLE SEQUENCE FROM THE GON_PHYLING RUN ALIGNMENT
	cd $outdir


	cd $outdir


	# mv ${read_copies_pwd}/trimmed_reads $outdir/trimmed_reads


	###################################################################################################33
	# IMPORTANT!
	# Output all variables (including the ones generated and used by this script) into a file
	# This is to reduce the necessity of assembling sequences when basecall is the sectio under development
	( set -o posix ; set ) > $outdir/variables.txt
	${rapup_path}/modules/var_file_fixer.py \
	--var_file ${outdir}/variables.txt \
	--out_file ${outdir}/variables.txt

	if [[ $bootstrapping == "ON" ]];
	then
		# BOOTSTRAP RAPUP AND GONPHYLING FOR CONSENSUS TREE PRODUCTION
		# cd $outdir/trimmed_reads/spades_output/genomes_for_parsnp/alignment_fixing/ 
		# raxmlHPC -J MR -m GTRGAMMA -z RAxML_bootstrap.core_genome_run.out -n gon_phy_MR_CONS
		# cat RAxML_MajorityRuleConsensusTree.gon_phy_MR_CONS | sed -E 's/\[[0-9]+\]//g' > fixed_gon_phy_MR.tre
		# sed -i -e "s/;/:0.0;/g" fixed_gon_phy_MR.tre
		# sed -i -e 's/.ref//g' fixed_gon_phy_MR.tre
		cd ${outdir}/${ref_name}_rapup_run/combine_and_infer
		raxmlHPC -J MR -m GTRGAMMA -z RAxML_bootstrap.consensusFULL_bootstrap -n rapup_MR_CONS
		cat RAxML_MajorityRuleConsensusTree.rapup_MR_CONS | sed -E 's/\[[0-9]+\]//g' > fixed_rapup_MR.tre
		sed -i -e "s/;/:0.0;/g" fixed_rapup_MR.tre
		# printf "\n$program_path\n"

	elif [[ $bootstrapping == "OFF" ]];
	then
		echo "BOOTSTRAPPING TURNED OFF"
	fi

# End of EP_run function
}


if [[ $assembly == "ON" ]];
then 

	# update_reads=update_reads_dir
	# start_dir=start_tree_dir
	# update_dir=update_alignment_dir
	# loci_blast_indexes=loci_blast_indexes_dir
	
	# blast_output=blast_output_dir
	# gon_phy_basecall=gon_phy_basecall
	# rapup_basecall=rapup_basecall
	
	# short_ref=five_align
	# read_copies=read_copies

	mkdir -p $outdir

	printf "\n$outdir\n"

	cd $outdir

	# mkdir update_reads
	# mkdir start_dir
	# mkdir update_dir
	# mkdir loci_blast_indexes
	
	# blast_output
	# gon_phy_basecall
	# rapup_basecall
	
	# mkdir short_ref
	mkdir read_copies

	# mkdir -p $outdir

	# printf "\n$outdir\n"

	# cd $outdir

	# mkdir $snippy_reads
	# mkdir $read_copies

	cd ${outdir}/read_copies

	read_copies_pwd=$(pwd)

	cd ${outdir}

	# mkdir ${short_ref}

	# # Make folders holding reads for gon_phyling/rapup and Snippy
	# # This folder is for gon_phyling and rapup
	# ls ${master_reads}/*$r1_tail | sort -R > ${outdir}/taxa_list.txt

	# Make folders holding reads for gon_phyling/rapup and Snippy
	# This folder is for gon_phyling and rapup
	ls ${master_reads}/*$r1_tail | sort > ${outdir}/taxa_list.txt

	for i in $(cat ${outdir}/taxa_list.txt);
	do
		ln -s $i ${outdir}/read_copies/
	  	ln -s ${i%$r1_tail}$r2_tail ${outdir}/read_copies/
	done

	######################################
	# FUNCTION INSERT HERE ##

	# EP_run ${outdir} ${snip_ref} ${rapup_path} ${r1_tail} ${r2_tail} ${loci_positions} ${bootstrapping} ${intermediate_files}

	for ref in $(echo ${snip_ref} | sed 's/-/ /g');
	do

		EP_run ${outdir} ${ref} ${rapup_path} ${r1_tail} ${r2_tail} ${ref}_positions.csv ${bootstrapping} ${intermediate_files}

	done
	


elif [[ $assembly == "OFF" ]];
then
	echo "Skipping assembly and read alignment stage as per user request"
	echo "Reading variable file and beginning basecall stage"
fi




cd $outdir

# outdir=${1}
# ref_1=${2}
# ref_2=${3}
# rapup_path=${4}
# r1_tail=${5}
# r2_tail=${6}
# loci_positions_1=${7}
# loci_positions_2=${8}

# ${program_path}/reference_bias_investigation/empirical_data_comparison/compare.sh ${outdir}

##############################################################################################
# Seperating alignments into loci (if applicable) and individual taxa

if [[ $basecall == "ON" ]]; then

	source $outdir/variables.txt

	rm -r $outdir/$rapup_basecall
	rm -r $outdir/$snippy_basecall
	rm -r $outdir/$gon_phy_basecall
	rm -r $outdir/$snippy_rapup_basecall
	rm -r $outdir/gon_phy_to_rapup
	rm -r $outdir/gon_phy_to_snippy
	rm -r $outdir/rapup_to_snippy
	# rm -r $outdir/gon_phy_to_rapup/align_files
	# rm -r $outdir/gon_phy_to_snippy/align_files
	# rm -r $outdir/rapup_to_snippy/align_files
	# rm -r $outdir/gon_phy_to_rapup/assessment_output
	# rm -r $outdir/gon_phy_to_snippy/assessment_output
	# rm -r $outdir/rapup_to_snippy/assessment_output




	mkdir $outdir/$rapup_basecall
	mkdir $outdir/$snippy_basecall
	mkdir $outdir/$gon_phy_basecall
	mkdir $outdir/$snippy_rapup_basecall

	mkdir $outdir/gon_phy_to_rapup
	mkdir $outdir/gon_phy_to_snippy
	mkdir $outdir/rapup_to_snippy

	mkdir $outdir/gon_phy_to_rapup/align_files
	mkdir $outdir/gon_phy_to_snippy/align_files
	mkdir $outdir/rapup_to_snippy/align_files

	mkdir $outdir/gon_phy_to_rapup/assessment_output
	mkdir $outdir/gon_phy_to_snippy/assessment_output
	mkdir $outdir/rapup_to_snippy/assessment_output

	mkdir $outdir/gon_phy_to_snippy/blast_results
	mkdir $outdir/rapup_to_snippy/blast_results

	mkdir $outdir/$rapup_basecall/basecall_results
	mkdir $outdir/$snippy_basecall/basecall_results
	mkdir $outdir/$gon_phy_basecall/basecall_results
	mkdir $outdir/$snippy_rapup_basecall/basecall_results

	# # mkdir $outdir/$rapup_basecall/sep_loci
	# mkdir $outdir/$snippy_basecall/sep_loci
	# # mkdir $outdir/$gon_phy_basecall/sep_loci
	# mkdir $outdir/$snippy_rapup_basecall/sep_loci
	
	# mkdir $outdir/$snippy_rapup_basecall/sep_loci/individual_tax
	# mkdir $outdir/$rapup_basecall/sep_loci/individual_tax
	# mkdir $outdir/$snippy_basecall/sep_loci/individual_tax
	# mkdir $outdir/$gon_phy_basecall/sep_loci/individual_tax

	${rapup_path}/modules/locus_position_identifier.py --out_file_dir ${outdir}/${rapup_basecall}/sep_loci --position_csv_file ${outdir}/rapup_run/rapup_loci_positions.csv --concatenated_fasta ${outdir}/rapup_run/combine_and_infer/extended.aln
	# ${rapup_path}/modules/locus_position_identifier.py --out_file_dir ${outdir}/${snippy_basecall}/sep_loci --position_csv_file $loci_positions --concatenated_fasta ${outdir}/core.full.aln
	${rapup_path}/modules/locus_position_identifier.py --out_file_dir ${outdir}/${gon_phy_basecall}/sep_loci --position_csv_file $loci_positions --concatenated_fasta ${outdir}/trimmed_reads/spades_output/genomes_for_parsnp/alignment_fixing/combo.fas

	# mkdir $outdir/$rapup_basecall/sep_loci
	mkdir $outdir/$snippy_basecall/sep_loci
	# mkdir $outdir/$gon_phy_basecall/sep_loci
	mkdir $outdir/$snippy_rapup_basecall/sep_loci
	
	mkdir $outdir/$snippy_rapup_basecall/sep_loci/individual_tax
	mkdir $outdir/$rapup_basecall/sep_loci/individual_tax
	mkdir $outdir/$snippy_basecall/sep_loci/individual_tax
	mkdir $outdir/$gon_phy_basecall/sep_loci/individual_tax
	mkdir ${outdir}/${gon_phy_basecall}/chosen_loci
	mkdir ${outdir}/${gon_phy_basecall}/chosen_loci/individual_tax
	mkdir ${outdir}/${rapup_basecall}/chosen_loci
	mkdir ${outdir}/${rapup_basecall}/chosen_loci/individual_tax

	mkdir $outdir/$rapup_basecall/blast_results
	mkdir $outdir/$snippy_basecall/blast_results
	mkdir $outdir/$gon_phy_basecall/blast_results

	mkdir $outdir/$gon_phy_basecall/match_finding_loci
	mkdir $outdir/$gon_phy_basecall/loci_finding_results
	mkdir $outdir/$gon_phy_basecall/matched_loci_results

	ls ${outdir}/${gon_phy_basecall}/sep_loci/*.fasta | sort -R | head -${locus_num} > ${outdir}/${gon_phy_basecall}/loci_list.txt

	for i in $(cat ${outdir}/${gon_phy_basecall}/loci_list.txt);
	do
		cluster_name=$(basename ${i} | sed 's/.fasta//g')
		printf "\n$i\n"
		printf "\n$(basename $i)\n"
		${rapup_path}/modules/ref_producer.py -r --align_file $i \
		--out_file ${outdir}/${gon_phy_basecall}/match_finding_loci/single_tax_gon_phy-${cluster_name}
	done

	for i in $(ls $outdir/$rapup_basecall/sep_loci/*.fasta);
	# for i in $(ls $outdir/$rapup_basecall/sep_loci/individual_tax);
	do
		# makeblastdb -in $outdir/$rapup_basecall/sep_loci/individual_tax/$i -dbtype nucl -parse_seqids
		makeblastdb -in $i -dbtype nucl -parse_seqids
	done

	# BLAST THE SELECTED GON_PHYLING LOCI AGAINST THE RAPUP ALIGNMENTS
	# TO FIND THE MATCHING RAPUP LOCUS
	for i in $(ls ${outdir}/${gon_phy_basecall}/match_finding_loci/);
	do
		gon_phy_seq=$(basename $i | sed 's/.fasta//g')
		gon_phy_seq_path=$(realpath $i)
		for j in $(ls ${outdir}/${rapup_basecall}/sep_loci/*.fasta);
		do
			echo "${i}"
			echo "${j}"
			rapup_file=$(basename $j | sed 's/.fasta//g')
			rapup_file_path=$(realpath ${j})
			echo "blastn \
			-db ${rapup_file_path} \
			-query ${outdir}/${gon_phy_basecall}/match_finding_loci/$i \
			-out $outdir/$gon_phy_basecall/loci_finding_results/blast_output_${rapup_file}-${gon_phy_seq}.out \
			-max_hsps 1 \
			-outfmt 5"

			blastn \
			-db $j \
			-query ${outdir}/${gon_phy_basecall}/match_finding_loci/$i \
			-out ${outdir}/${gon_phy_basecall}/loci_finding_results/blast_output_${rapup_file}-${gon_phy_seq}.out \
			-max_hsps 1 \
			-outfmt 5
		done
	done

	#IDENTIFY EACH GON_PHY LOCUS AND ANALYZE ALL BLAST ALIGNMENTS TO FIND THE BEST ALIGNMENT
		${program_path}/empirical_data_comparison/emp_blast_matcher.py \
		--input_folder ${outdir}/${gon_phy_basecall}/loci_finding_results \
		--output_file ${outdir}/${gon_phy_basecall}/matched_loci_results/locus_file_matches_
	
	# SEPARATE EACH LOCUS ALIGNMENT INTO INDIVIDUAL SEQUENCES (SNIPPY SEQS STAY FULL SIZE)
	for j in $(ls ${outdir}/${gon_phy_basecall}/sep_loci/*.fasta);
	do
		for i in $(grep ">" $j | sed -e 's/>//g' | sed -e 's/.ref//g');
		do
			#printf "\n$j\n$i\n"
			msa_file=$(basename $j | sed -e 's/.fasta//g')
			#printf "\n$msa_file\n"
			${rapup_path}/modules/ref_producer.py \
			-s \
			--align_file $j \
			--out_file ${outdir}/${gon_phy_basecall}/sep_loci/individual_tax/single_tax_gon_phy-${msa_file}-${i} \
			--ref_select $i
		done
	done

	for j in $(ls ${outdir}/${rapup_basecall}/sep_loci/*.fasta);
	do
	        for i in $(grep ">" $j | sed -e 's/>//g' | sed -e 's/.ref//g');
	        do
	                #printf "\n$j\n$i\n"
	                msa_file=$(basename $j | sed -e 's/.fasta//g')
			#printf "\n$msa_file\n"
	                ${rapup_path}/modules/ref_producer.py -s --align_file $j --out_file ${outdir}/${rapup_basecall}/sep_loci/individual_tax/single_tax_rapup-${msa_file}-${i} --ref_select $i
	        done
	done

	#Create single taxa locus files for snippy output
	for j in $(ls ${outdir}/core.full.aln);
	#for j in $(ls $outdir/$snippy_basecall/sep_loci/*.fasta);
	do
	        for i in $(grep ">" $j | sed -e 's/>//g' | sed -e 's/.ref//g');
	        do
	                #printf "\n$j\n$i\n"
	                msa_file=$(basename $j | sed -e 's/.fasta//g')
			#printf "\n$msa_file\n"
	                ${rapup_path}/modules/ref_producer.py \
					-s --align_file $j \
					--out_file ${outdir}/${snippy_basecall}/sep_loci/individual_tax/single_tax_snippy-${msa_file}-${i} \
					--ref_select $i
	        done
	done

	snippy_seqs=$(ls $outdir/$snippy_basecall/sep_loci/individual_tax)
	rapup_seqs=$(ls $outdir/$rapup_basecall/sep_loci/individual_tax)
	gon_phyling_seqs=$(ls $outdir/$gon_phy_basecall/sep_loci/individual_tax)

	# NOW MAKE THE ACTUAL ALIGNMENTS AND COMPARISONS
	# BETWEEN RAPUP AND GON_PHYLING FIRST
	for j in $(cat ${outdir}/${gon_phy_basecall}/loci_list.txt);
	do
		# locus=$(basename ${j} | sed 's/.fasta//g')
		msa_file=$(basename ${j} | sed -e 's/.fasta//g')
		ln -s ${j} ${outdir}/${gon_phy_basecall}/chosen_loci/
		for i in $(grep ">" $j | sed -e 's/>//g' | sed -e 's/.ref//g');
		do
			#printf "\n$j\n$i\n"
			# msa_file=$(basename $j | sed -e 's/.fasta//g')
			#printf "\n$msa_file\n"
			${rapup_path}/modules/ref_producer.py \
			-s \
			--align_file $j \
			--out_file ${outdir}/${gon_phy_basecall}/chosen_loci/individual_tax/single_tax_gon_phy-${msa_file}-${i} \
			--ref_select $i
		done
	done

	for file in $(ls -1 ${outdir}/${gon_phy_basecall}/matched_loci_results/);
	do
		gon_phy_locus=$(head -1 ${outdir}/${gon_phy_basecall}/matched_loci_results/$file | \
		sed 's/single_tax_gon_phy-//g' | sed 's/-//g')
		
		rapup_locus=$(tail -1 ${outdir}/${gon_phy_basecall}/matched_loci_results/$file | \
		sed 's/.fasta//g' | sed 's/-//g')
		
		stripped_path_rapup_file=$(basename ${rapup_locus})
		rapup_locus=${stripped_path_rapup_file}
		echo $gon_phy_locus
		echo $rapup_locus

		# PROGRAM TAKES IN THE FOLDER OF SELECTED LOCI PRODUCED BY GON_PHYLING
		# FOLDER OF ALL RAPUP PRODUCED LOCI
		# THE NAME OF THE LOCUS FOR THE GON_PHYLING DATASET
		# THE NAME OF THE LOCUS FOR THE RAPUP DATASET
		# THE OUTPUT DIRECTORY FOR THE RESULTING COMBINED SEQUENCE FILES
		# THE OUTPUT DIRECTORY FOR THE MATCHED RAPUP SINGLE SEQUENCE FILES FOR FURTHER ALIGNMENT TO SNIPPY SEQS
		${program_path}/empirical_data_comparison/emp_seq_matcher.py \
		--folder_1 ${outdir}/${gon_phy_basecall}/chosen_loci/individual_tax \
		--folder_2 ${outdir}/${rapup_basecall}/sep_loci/individual_tax \
		--cluster_id_1 $gon_phy_locus \
		--cluster_id_2 $rapup_locus \
		--output_dir $outdir/gon_phy_to_rapup/align_files/ \
		--matched_seq_output_dir ${outdir}/${rapup_basecall}/chosen_loci/individual_tax

		for i in $(cat $outdir/gon_phy_to_rapup/align_files/taxa_name_list.txt);
		do
			mafft --op 5 --lexp -0.5 --thread $align_threads $outdir/gon_phy_to_rapup/align_files/combined_original_${gon_phy_locus}_${rapup_locus}--${i}-- > $outdir/gon_phy_to_rapup/align_files/aligned_combine-${gon_phy_locus}-${rapup_locus}--${i}-original.fasta
			mafft --op 5 --lexp -0.5 --thread $align_threads $outdir/gon_phy_to_rapup/align_files/combined_reverse_complement_${gon_phy_locus}_${rapup_locus}--${i}-- > $outdir/gon_phy_to_rapup/align_files/aligned_combine-${gon_phy_locus}-${rapup_locus}--${i}-reverse_complement.fasta
			mafft --op 5 --lexp -0.5 --thread $align_threads $outdir/gon_phy_to_rapup/align_files/combined_reverse_${gon_phy_locus}_${rapup_locus}--${i}-- > $outdir/gon_phy_to_rapup/align_files/aligned_combine-${gon_phy_locus}-${rapup_locus}--${i}-reverse.fasta
			mafft --op 5 --lexp -0.5 --thread $align_threads $outdir/gon_phy_to_rapup/align_files/combined_complement_${gon_phy_locus}_${rapup_locus}--${i}-- > $outdir/gon_phy_to_rapup/align_files/aligned_combine-${gon_phy_locus}-${rapup_locus}--${i}-complement.fasta

			${program_path}/empirical_data_comparison/emp_align_compare.py \
			--align_1 $outdir/gon_phy_to_rapup/align_files/aligned_combine-${gon_phy_locus}-${rapup_locus}--${i}-reverse_complement.fasta \
			--align_2 $outdir/gon_phy_to_rapup/align_files/aligned_combine-${gon_phy_locus}-${rapup_locus}--${i}-complement.fasta \
			--align_3 $outdir/gon_phy_to_rapup/align_files/aligned_combine-${gon_phy_locus}-${rapup_locus}--${i}-reverse.fasta \
			--align_4 $outdir/gon_phy_to_rapup/align_files/aligned_combine-${gon_phy_locus}-${rapup_locus}--${i}-original.fasta \
			--output_stub $outdir/gon_phy_to_rapup/assessment_output/gon_rap_results-${gon_phy_locus}-${rapup_locus}--${i}

		done
		wait

	done

	

	for j in $(cat ${outdir}/${gon_phy_basecall}/loci_list.txt);
	do
		locus=$(basename ${j} | sed -e 's/.fasta//g')
		${program_path}/empirical_data_comparison/blast_location_finder.py \
		--long_seqs_folder ${outdir}/snippy_basecall/sep_loci/individual_tax \
		--manipulate_seqs_folder ${outdir}/gon_phy_basecall/chosen_loci/individual_tax \
		--output_dir ${outdir}/gon_phy_to_snippy/align_files

	done

	for i in $(ls ${outdir}/gon_phy_to_snippy/align_files);
	do
		mafft \
		--thread $align_threads \
		--op 5 \
		--lexp -0.5 \
		${outdir}/gon_phy_to_snippy/align_files/${i} \
		> \
		${outdir}/gon_phy_to_snippy/align_files/aligned_${i}


		${program_path}/empirical_data_comparison/emp_snippy_gapped_align_compared.py \
		-t \
		--align_1 ${outdir}/gon_phy_to_snippy/align_files/aligned_${i} \
		--output_stub snip_gon_assessment_${i} \
		--output_dir ${outdir}/gon_phy_to_snippy/assessment_output

	done

#####################################################################################################
# RAPUP TO SNIPPY SECTION

	for j in $(ls -1 ${outdir}/${rapup_basecall}/chosen_loci/individual_tax/);
	do
		locus=$(basename ${j} | sed -e 's/single_tax-//g')
		taxon=$(basename ${j} | rev | cut -d "-" -f1 | rev)
		${program_path}/empirical_data_comparison/blast_location_finder.py \
		--long_seqs_folder ${outdir}/snippy_basecall/sep_loci/individual_tax \
		--manipulate_seqs_folder ${outdir}/rapup_basecall/chosen_loci/individual_tax \
		--output_dir ${outdir}/rapup_to_snippy/align_files

	done

	for i in $(ls ${outdir}/rapup_to_snippy/align_files);
	do
		mafft \
		--thread $align_threads \
		--op 5 \
		--lexp -0.5 \
		${outdir}/rapup_to_snippy/align_files/${i} \
		> \
		${outdir}/rapup_to_snippy/align_files/aligned_${i}


		${program_path}/empirical_data_comparison/emp_snippy_gapped_align_compared.py \
		-t \
		--align_1 ${outdir}/rapup_to_snippy/align_files/aligned_${i} \
		--output_stub snip_rap_assessment_${i} \
		--output_dir ${outdir}/rapup_to_snippy/assessment_output
	done

	


elif [[ $basecall == "OFF" ]]; then
printf "\nSKIPPING BASECALL SECTION OF PIPELINE PER USER SETTING\n"

fi
