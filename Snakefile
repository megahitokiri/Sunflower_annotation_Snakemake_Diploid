#--------------------------------------------------------------------------------
# MAIN CONFIGURATION PARAMETERS (INTERNAL CONFIG FILE)
#Cancel all jobs if fail: squeue -u $USER -h | awk '{print $1}' | xargs scancel
#--------------------------------------------------------------------------------

CHRS = '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35'.split()
PROJECT = "Harg2202"
REFERENCE = "MAIN_FASTAs/Harg2202r1.0-20210824.genome.fasta"
NEW_NAMES_REFERENCE = "Harg2202r1.0-20210824.genome.new_names.fasta"
GFF3_FILE = "MAIN_GFF3s/Harg2202r1.0-20210824.gff3"

#--------------------------------------------------------------------------------
# TargetRule FINAL_GFF3
#--------------------------------------------------------------------------------

rule FINAL_GFF3:
	input:
		expand("{Main_Reference}",Main_Reference=REFERENCE),
		expand("{Project}/Ref/scaffold_{Chrs}.fasta",Project=PROJECT,Chrs = CHRS),
		expand("{Project}/EDTA_Files/scaffold_{Chrs}.fasta.mod.EDTA.TEanno.gff3",Project=PROJECT,Chrs = CHRS),
		expand("{Project}/EDTA_Files/{Project}_chr{Chrs}.masked.fasta",Project=PROJECT,Chrs = CHRS),
		expand("{Gff3_file}",Gff3_file=GFF3_FILE),
		expand("{Project}/Annotation_steps/{Project}_step1_chr{Chrs}.gff3",Project=PROJECT,Chrs = CHRS),
		expand("{Project}/Annotation_steps/{Project}_step2_chr{Chrs}.gff3",Project=PROJECT,Chrs = CHRS),
		expand("{Project}/FINAL_ANNOTATION/FINAL_{Project}_v1_1.sorted.gff3",Project=PROJECT),
		expand("{Project}/Summary_data/{Project}.protein.fasta",Project=PROJECT),
		expand("{Project}/Summary_data/JML.{Project}_GFF3_summary.txt",Project=PROJECT),
		expand("{Project}/Summary_data/Original.{Project}_GFF3_summary.txt",Project=PROJECT),
		expand("{Project}/Summary_data/{Project}_busco/short_summary.specific.eudicots_odb10.{Project}_busco.txt",Project=PROJECT),
		expand("{Project}/Summary_data/{Project}_busco/run_eudicots_odb10/missing_busco_list.tsv",Project=PROJECT),
		expand("{Project}/Summary_data/{Project}.original.protein.fasta",Project=PROJECT),
		expand("{Project}/Summary_data/{Project}_busco_original/short_summary.specific.eudicots_odb10.{Project}_busco_original.txt",Project=PROJECT),
		expand("{Project}/Summary_data/{Project}_busco_original/run_eudicots_odb10/missing_busco_list.tsv",Project=PROJECT),
		expand("{Project}/Summary_data/Missing_buscos/{Project}_Missing_buscos.txt",Project=PROJECT),
		expand("{Project}/FINAL_ANNOTATION/FINAL_{Project}_v1_2.sorted.flagged.gff3",Project=PROJECT),
	params:
		project=PROJECT,
	shell:
		"""
		echo Pipeline Finished correctly for {params.project}..... CONGRATS!!! > {params.project}.Pipeline_complete.txt

		mv *.err {params.project}/logs
		mv *.out {params.project}/logs

		date "+DATE: %D%nTIME: %T" >> {params.project}.Pipeline_complete.txt
		
		cp -v {params.project}/FINAL_ANNOTATION/FINAL_{params.project}_v1_2.sorted.flagged.gff3 {params.project}/{params.project}.Rieseberg.v1_2.gff3
		echo final annotation file located at: {params.project}/{params.project}.Rieseberg.v1_2.gff3 >> {params.project}.Pipeline_complete.txt

		echo Pipeline Finished correctly for {params.project}..... CONGRATS!!!
		"""	
#--------------------------------------------------------------------------------
# Init: Initializing files and folder
#--------------------------------------------------------------------------------
rule Init:
	input:
		reference={REFERENCE},
		Assembly_spliter="Assembly_Chr_splitter.R",
	output:
		expand("{Project}/Ref/{New_names_ref}",Project=PROJECT,New_names_ref=NEW_NAMES_REFERENCE),
	params:
		project=PROJECT,
	shell:
		"""
		snakemake --dag | dot -Tsvg > dag.svg
		mkdir {params.project}
		cd {params.project}
		
			mkdir logs
			mkdir EDTA_Files
			mkdir Ref
			mkdir Annotation_steps
			mkdir FINAL_ANNOTATION
			mkdir Summary_data
		cd ..
		cp  {input.reference} {output}
		ml samtools
		#Index FASTA file
		samtools faidx {output} 
		
		cd {params.project}
		
		ml unload samtools
		"""

#--------------------------------------------------------------------------------
# Chr_splitting: Split the Main Assembly in Chromosomes for easy handling
#--------------------------------------------------------------------------------

rule Chr_splitting:
	input:
		rules.Init.output,
	output:
		expand("{Project}/Ref/scaffold_{Chrs}.fasta",Project=PROJECT,Chrs=CHRS),
		
	params:
		project=PROJECT,
	shell:
		"""
		echo Assembly_Chr_splitter.R --args -f {input} 
		ml r/4.1.0
		echo Assembly split into Chromosomes
		R --vanilla < Assembly_Chr_splitter.R --args -f {input} &&
		mv *scaffold*.fasta {params.project}/Ref/
		ml unload r/4.1.0
		ml unload samtools
		"""

#--------------------------------------------------------------------------------
# EDTA_individual: Look for TE elements on individual Fasta Chr
#--------------------------------------------------------------------------------
rule EDTA_individual:
	input:
		rules.Chr_splitting.output,
	output:
		gff3_file="{Project}/EDTA_Files/scaffold_{Chrs}.fasta.mod.EDTA.TEanno.gff3",
	params:
		project=PROJECT,

	shell:
		"""
		cp {params.project}/Ref/scaffold_{wildcards.Chrs}.fasta {params.project}/EDTA_Files
		cp -v Name_checker_pre.sh {params.project}/EDTA_Files
		cp -v Name_checker_post.sh {params.project}/EDTA_Files
		
		cd {params.project}/EDTA_Files
		
		echo "Name_cheker pre and post correct EDTA bigger than 15 characters name error on FASTA."
		bash Name_checker_pre.sh scaffold_{wildcards.Chrs}.fasta {params.project}
		
		eval "$(conda shell.bash hook)"
		conda activate EDTA
			echo starting EDTA process on: scaffold_{wildcards.Chrs}.fasta
			EDTA.pl --overwrite 0 --genome scaffold_{wildcards.Chrs}.fasta --sensitive 0 --anno 1 --evaluate 0 --threads 16 --force 1
		conda deactivate
		
		bash Name_checker_post.sh scaffold_{wildcards.Chrs}.fasta {params.project} scaffold_{wildcards.Chrs}.fasta.mod.EDTA.TEanno.gff3
		"""
		
#--------------------------------------------------------------------------------
# Masked_FASTA: Create masked fasta for further analysis from EDTA results.
#--------------------------------------------------------------------------------

rule Masked_FASTA:
	input:
		EDTA_repeats_file=rules.EDTA_individual.output.gff3_file,
		reference=rules.Chr_splitting.output,
	output:
		masked_fasta_file="{Project}/EDTA_Files/{Project}_chr{Chrs}.masked.fasta",
	params:
		project=PROJECT,
	shell:
		"""
		ml bedtools
		cd {params.project}/EDTA_Files
		echo "Creating Masked Reference Genome for scaffold_{wildcards.Chrs}.masked.fasta"
		bedtools maskfasta -fi scaffold_{wildcards.Chrs}.fasta -bed scaffold_{wildcards.Chrs}.fasta.mod.EDTA.TEanno.gff3 -fo {params.project}_chr{wildcards.Chrs}.masked.fasta
		ml unload bedtools	
		"""

#--------------------------------------------------------------------------------
# STEP1_annotation: Transform Analysis into parsed GFF3 files .
#--------------------------------------------------------------------------------

rule STEP1_annotation:
	input:
		Masked_FASTA_file=rules.Masked_FASTA.output.masked_fasta_file,
		Gff3_file={GFF3_FILE},
	output:
		gff3_step1="{Project}/Annotation_steps/{Project}_step1_chr{Chrs}.gff3",
	params:
		project=PROJECT,
	shell:
		"""
		ml r/4.1.0
		cp {input.Gff3_file} {params.project}/Annotation_steps/{params.project}.gff3
		echo "Processing Step1 of annotation"
		R --vanilla < Step1_Annotation.R --args -a {params.project} -c {wildcards.Chrs}
		ml unload r/4.1.0
		"""
#------------------------------------------------------------------------------------
# STEP2_annotation: Uses EDTA maked fasta file and ORF analysis to determine viable genes
#------------------------------------------------------------------------------------

rule STEP2_annotation:
	input:
		GFF3_File=rules.STEP1_annotation.output,
		Ref_File=rules.Chr_splitting.output,
		Masked_FASTA_File=rules.Masked_FASTA.output,
		Original_Gff3_file={GFF3_FILE},
	output:
		gff3_step2="{Project}/Annotation_steps/{Project}_step2_chr{Chrs}.gff3",
	params:
		project=PROJECT,
	shell:
		"""
		BASEDIR=$PWD 
		pwd
		cp -v $BASEDIR/Step2_Filtering.R {params.project}/Annotation_steps/
		cd {params.project}/Annotation_steps/
		ml r/4.1.0
		R --vanilla < Step2_Filtering.R --args -g $BASEDIR/{params.project}/Annotation_steps/{params.project}_step1_chr{wildcards.Chrs}.gff3 -a $BASEDIR/{params.project}/Ref/scaffold_{wildcards.Chrs}.fasta -m $BASEDIR/{params.project}/EDTA_Files/{params.project}_chr{wildcards.Chrs}.masked.fasta -o {params.project}_step2_chr{wildcards.Chrs} -s $BASEDIR/{input.Original_Gff3_file}
		ml unload r/4.1.0
		cd ..
		"""
		
#------------------------------------------------------------------------------------
# Chr_merge: Fuse all gff3 individual chromosomes into complete assembly again
#------------------------------------------------------------------------------------

rule Chr_merge:
	input:
		expand("{Project}/Annotation_steps/{Project}_step2_chr{Chrs}.gff3",Project=PROJECT,Chrs = CHRS),
	output:
		"{Project}/FINAL_ANNOTATION/FINAL_{Project}_v1_1.sorted.gff3",
	params:
		project=PROJECT,
		Chrs=CHRS,
	shell:
		"""
		cat {params.project}/Annotation_steps/{params.project}_step2_chr1.gff3 > {params.project}/FINAL_ANNOTATION/FINAL_{params.project}_v1.gff3
		
		for i in {{2..35}}
			do
			tail -n +4 {params.project}/Annotation_steps/{params.project}_step2_chr$i.gff3 >> {params.project}/FINAL_ANNOTATION/FINAL_{params.project}_v1.gff3
			done
		
		/home/jmlazaro/github/gff3sort/gff3sort.pl --chr_order original {params.project}/FINAL_ANNOTATION/FINAL_{params.project}_v1.gff3 > {params.project}/FINAL_ANNOTATION/FINAL_{params.project}_v1.sorted.gff3
	
		echo awk '{{gsub("character\\\(0\\\)", '0');print}}' {params.project}/FINAL_ANNOTATION/FINAL_{params.project}_v1.sorted.gff3 
		awk '{{gsub("character\\\(0\\\)", "0");print}}' {params.project}/FINAL_ANNOTATION/FINAL_{params.project}_v1.sorted.gff3 > {params.project}/FINAL_ANNOTATION/FINAL_{params.project}_v1_1.sorted.gff3
		"""		

#------------------------------------------------------------------------------------
# Summary_statistics: Get summary of the new gff3 file
#------------------------------------------------------------------------------------

rule Summary_statistics:
	input:
		GFF3_file=rules.Chr_merge.output,
		Ref_file=rules.Init.output,
	output:
		Protein_FASTA="{Project}/Summary_data/{Project}.protein.fasta",
	params:
		project=PROJECT,
	shell:
		"""
		BASEDIR=$PWD
		ml StdEnv/2020  gcc/9.3.0
		ml nixpkgs/16.09  gcc/5.4.0
		ml nixpkgs/16.09  gcc/7.3.0
		ml transdecoder/5.5.0
		
		gff3_file_to_proteins.pl --gff3 {input.GFF3_file} --fasta $BASEDIR/{input.Ref_file} --seqType prot > $BASEDIR/{params.project}/Summary_data/{params.project}.protein.fasta
		ml unload transdecoder/5.5.0
		
		echo Protein File Generated correctly CORRECTLY ....
		ml unload perl
		"""
#------------------------------------------------------------------------------------
# GFF3_statistics: Calculate Statistics for Original and filtered GFF3 files
#------------------------------------------------------------------------------------

rule GFF3_statistics:
	input:
		GFF3_file=rules.Chr_merge.output,
		Original_Gff3_file={GFF3_FILE},
	output:
		New_GFF3_summary="{Project}/Summary_data/JML.{Project}_GFF3_summary.txt",
		Original_GFF3_summary="{Project}/Summary_data/Original.{Project}_GFF3_summary.txt",
	params:
		project=PROJECT,
	shell:
		"""
		BASEDIR=$PWD
		cp -v GFF3_Summary_Statistics.R {params.project}/Summary_data/
		cd {params.project}/Summary_data/
		
		ml r/4.1.0
		echo Calculating Statistics on New processed file
		R --vanilla < GFF3_Summary_Statistics.R --args --gff $BASEDIR/{input.GFF3_file} -o JML.{params.project}
		
		echo Calculating Statistics on Original file
		R --vanilla < GFF3_Summary_Statistics.R --args --gff $BASEDIR/{input.Original_Gff3_file} -o Original.{params.project}
		cd ../..
		"""
		
#------------------------------------------------------------------------------------
# BUSCO: Evaluate the Cognate results into BUSCO protein mode
#------------------------------------------------------------------------------------

rule BUSCO:
	input:
		Protein_fasta=rules.Summary_statistics.output.Protein_FASTA,
	output:
		Busco_New_results="{Project}/Summary_data/{Project}_busco/short_summary.specific.eudicots_odb10.{Project}_busco.txt",
		Busco_New_Missing="{Project}/Summary_data/{Project}_busco/run_eudicots_odb10/missing_busco_list.tsv",
	params:
		project=PROJECT,
	shell:
		"""
		BASEDIR=$PWD
		ml StdEnv/2020
		ml gcc/9.3.0
		ml openmpi/4.0.3
		ml busco/5.2.2
		mkdir {params.project}/Summary_data/busco_downloads
		mkdir {params.project}/Summary_data/busco_downloads/lineages
		
		cp /home/jmlazaro/BUSCO/eudicots_odb10.tar.gz {params.project}/Summary_data/busco_downloads/lineages
		cd {params.project}/Summary_data/busco_downloads/lineages
		tar -xvf eudicots_odb10.tar.gz
		cd $BASEDIR/{params.project}/Summary_data/
		
		busco -f -c 4 -m protein -i $BASEDIR/{params.project}/Summary_data/{params.project}.protein.fasta -o {params.project}_busco -l eudicots_odb10 --offline --download_path $BASEDIR/{params.project}/Summary_data/busco_downloads
		echo done BUSCO analysis
		"""

#------------------------------------------------------------------------------------
# Summary_statistics_Original: Get summary of the Original gff3 file
#------------------------------------------------------------------------------------

rule Summary_statistics_Original:
	input:
		Original_Gff3_file={GFF3_FILE},
		Ref_file=rules.Init.output,
	output:
		Protein_FASTA="{Project}/Summary_data/{Project}.original.protein.fasta",
	params:
		project=PROJECT,
	shell:
		"""
		BASEDIR=$PWD
		ml StdEnv/2020  gcc/9.3.0
		ml nixpkgs/16.09  gcc/5.4.0
		ml nixpkgs/16.09  gcc/7.3.0
		ml transdecoder/5.5.0
		
		gff3_file_to_proteins.pl --gff3 {input.Original_Gff3_file} --fasta $BASEDIR/{input.Ref_file} --seqType prot > $BASEDIR/{params.project}/Summary_data/{params.project}.original.protein.fasta
		ml unload transdecoder/5.5.0
		
		echo Original Protein File Generated correctly CORRECTLY ....
		ml unload perl
		"""

#---------------------------------------------------------------------------------------
# BUSCO_Original: Evaluate the Cognate results into BUSCO protein mode from Orignal GFF3
#---------------------------------------------------------------------------------------

rule BUSCO_Original:
	input:
		Protein_fasta=rules.Summary_statistics_Original.output.Protein_FASTA,
	output:
		Busco_Original_results="{Project}/Summary_data/{Project}_busco_original/short_summary.specific.eudicots_odb10.{Project}_busco_original.txt",
		Busco_Original_Missing="{Project}/Summary_data/{Project}_busco_original/run_eudicots_odb10/missing_busco_list.tsv",
	params:
		project=PROJECT,
	shell:
		"""
		BASEDIR=$PWD
		ml StdEnv/2020
		ml gcc/9.3.0
		ml openmpi/4.0.3
		ml busco/5.2.2
		mkdir {params.project}/Summary_data/busco_downloads_original
		mkdir {params.project}/Summary_data/busco_downloads_original/lineages
		
		cp /home/jmlazaro/BUSCO/eudicots_odb10.tar.gz {params.project}/Summary_data/busco_downloads_original/lineages
		cd {params.project}/Summary_data/busco_downloads_original/lineages
		tar -xvf eudicots_odb10.tar.gz
		cd $BASEDIR/{params.project}/Summary_data/
		
		busco -f -c 4 -m protein -i $BASEDIR/{params.project}/Summary_data/{params.project}.original.protein.fasta -o {params.project}_busco_original -l eudicots_odb10 --offline --download_path $BASEDIR/{params.project}/Summary_data/busco_downloads_original
		echo done BUSCO analysis
		rm -rf $BASEDIR/{params.project}/Summary_data/busco_downloads_original
		echo done removing termporal BUSCO_downloads_Original Folder 
		"""		

#---------------------------------------------------------------------------------------
# Missing_BUSCO: Evaluate Missing Buscos
#---------------------------------------------------------------------------------------

rule Missing_BUSCO:
	input:
		Busco_Original=rules.BUSCO_Original.output.Busco_Original_Missing,
		Busco_New=rules.BUSCO.output.Busco_New_Missing,
		Original_GFF3_file={GFF3_FILE},
		GFF3_file=rules.Chr_merge.output,
	output:
		Missing_Busco_List="{Project}/Summary_data/Missing_buscos/{Project}_Missing_buscos.txt",
		New_BUSCO_GFF3_file="{Project}/FINAL_ANNOTATION/FINAL_{Project}_v1_2.sorted.flagged.gff3",
	params:
		project=PROJECT,
	shell:
		"""
		BASEDIR=$PWD
		mkdir {params.project}/Summary_data/Missing_buscos
		cp -v {input.Busco_Original} {params.project}/Summary_data/Missing_buscos/{params.project}_busco_original.txt
		cp -v {input.Busco_New} {params.project}/Summary_data/Missing_buscos/{params.project}_busco_new.txt
		cp -v {params.project}/Summary_data/{params.project}_busco_original/run_eudicots_odb10/full_table.tsv {params.project}/Summary_data/Missing_buscos
		
		cd {params.project}/Summary_data/Missing_buscos
		echo "processing diff file"
		diff -y {params.project}_busco_new.txt {params.project}_busco_original.txt > $BASEDIR/{params.project}/Summary_data/Missing_buscos/{params.project}_Missing_preprocessed.txt
		grep "<" $BASEDIR/{params.project}/Summary_data/Missing_buscos/{params.project}_Missing_preprocessed.txt | awk '{{print $1}}' > $BASEDIR/{params.project}/Summary_data/Missing_buscos/{params.project}_Missing_buscos.txt
		echo "diff file correctly processed"
		grep -F -f {params.project}_Missing_buscos.txt full_table.tsv > $BASEDIR/{params.project}/Summary_data/Missing_buscos/{params.project}_Missing_Busco_genes_list.txt
		echo Missing Busco genes list ready
		awk '{{gsub("mRNA:", "");print}}' $BASEDIR/{params.project}/Summary_data/Missing_buscos/{params.project}_Missing_Busco_genes_list.txt | awk '{{print $3}}' > $BASEDIR/{params.project}/Summary_data/Missing_buscos/{params.project}_Missing_Busco_genes_list.tsv

		echo extracting missing buscos from original GFF3
#		head -n 1 $BASEDIR/{input.Original_GFF3_file} > $BASEDIR/{params.project}/Summary_data/Missing_buscos/{params.project}_Missing_Busco_genes.gff3
		grep -F -f {params.project}_Missing_Busco_genes_list.tsv $BASEDIR/{input.Original_GFF3_file} > $BASEDIR/{params.project}/Summary_data/Missing_buscos/{params.project}_Missing_Busco_genes.tsv
		awk '{{gsub("locus_tag=", "locus_tag=BUSCO_");print}}' {params.project}_Missing_Busco_genes.tsv > $BASEDIR/{params.project}/Summary_data/Missing_buscos/{params.project}_Missing_Busco_genes.gff3
		echo missing busco GFF3 created correctly.

		cat $BASEDIR/{params.project}/FINAL_ANNOTATION/FINAL_{params.project}_v1_1.sorted.gff3 {params.project}_Missing_Busco_genes.gff3 > $BASEDIR/{params.project}/Summary_data/Missing_buscos/FINAL_{params.project}_v1_2.flagged.gff3
		echo Merging GFF3 files to reincorporate missing BUSCOs gene
		
		/home/jmlazaro/github/gff3sort/gff3sort.pl --chr_order original $BASEDIR/{params.project}/Summary_data/Missing_buscos/FINAL_{params.project}_v1_2.flagged.gff3 > $BASEDIR/{params.project}/FINAL_ANNOTATION/FINAL_{params.project}_v1_2.sorted.flagged.gff3
		echo Sorted completed correcty.
		"""			
