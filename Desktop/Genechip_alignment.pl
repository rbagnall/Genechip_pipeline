#!/usr/bin/perl

#This script aligns paired end NGS data to the illumina TruSeq or Nextera Exome or sureselect_EZ_Exome_v4 + UTR (+50 nucleotides of padding either side of exome intervals)
#written by richard bagnall (r.bagnall@centenary.org.au)
#Usage: Exome_align.pl -fastq [path/2/user] -cohort [cohort name] -exome [truseq or sureselectv4]
#save raw reads as [name1].1.fastq.gz, [name1].2.fastq.gz , [name2].1.fastq.gz , [name2].2.fastq.gz, etc...
#save the raw reads in a directory called Rawdata e.g. /home/shared/NGS/human/[USER]/Rawdata

use strict; use warnings;
use File::Copy;
use Getopt::Long;
use Parallel::ForkManager;
use List::Util qw( min max );

####################################################
# Setting up command line options and usage errors #
####################################################

# set a command line option variable (-cohort); is required and used to assign a foldername
my $cohort = ''; 
# set a command line option variable (-fastq); is required and is the path to the rawdata
my $input ='';
# set a command line option variable (-L); is required and is the intervals list: genechip101, genechip69, genechip98 or trusight_cardiomyopathy
my $intervals ='';

GetOptions ("cohort=s" => \$cohort,
			"fastq=s" => \$input,
            "L=s" => \$intervals);

my $offset = length($input); # use this extensively in subroutines for getting names of samples, files and folders

# -cohort is required, else print usage and die
if ($cohort eq '') {
	print "\n\n\t*** ERROR: You need to define a name for this cohort with the -cohort option\n";
	usage();
    exit;
}
# -fastq is required, else print usage and die
if ($input eq '') {
	print "\n\n\t*** ERROR: You need to define the path to the raw fastq files with the -fastq option\n";
	usage();
    exit;
}
# pre-empt common error
elsif ($input =~ m/\/$/) {
	print "\n\n\t*** ERROR: Please remove the last / from the -fastq option\n";
	usage();
    exit;
}

if ($intervals !~ m/^(genechip101|genechip69|genechip98|trusight_cardiomyopathy)$/) {
	print "\n\n\t*** ERROR: You need to define the intervals used with the -L option (genechip101|genechip69|genechip98|trusight_cardiomyopathy)\n";
	usage();
    exit;
}

print "\n\n\n\t\t---------GENE PANEL ALIGNMENT-----------\n";
print "\t\tAnalysing $cohort cohort exome data \n";
print "\t\tAligning Paired-end sequence reads (Illumina 1.8+)\n";
print "\t\tUsing the $intervals panel\n";
print "\t\tRealigning and recalibrating over the $intervals intervals\n";
print "\t\tUsing BWA MEM for read alignment\n";
print "\t\tUsing Picard to remove duplicate reads\n";
print "\t\tUsing GATK v3 for realigning\n";
print "\t\tUsing GATK v3 for recalibration\n";
print "\t\t-------------------------------------------------------------\n\n";
print "\n";

# create temporary folder
mkdir ("$input/tempSAMs") or die "Unable to create tempSAMs directory: <$!>\n";
mkdir ("$input/tempBAMs") or die "Unable to create tempBAMs directory: <$!>\n";
mkdir ("$input/tempLOGs") or die "Unable to create tempLOGs directory: <$!>\n";

# paths
my $path2rawdata = "$input/Rawdata";
my $path2bam = "$input/tempBAMs";
my $path2sam = "$input/tempSAMs";
my $path2log = "$input/tempLOGs";
my $path2gatk = '/home/groups/cardio/Applications/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar';
my $path2indel1kg = '/home/groups/cardio/References/INDELS/1000G_phase1.indels.b37.vcf';
my $path2indelmills = '/home/groups/cardio/References/INDELS/Mills_and_1000G_gold_standard.indels.b37.sites.vcf';
my $path2mergedexome50 = '/home/groups/cardio/References/Exome/Genechip_seq_pipeline/merged_plus50.intervals';
my $path2mergedexome10 = '/home/groups/cardio/References/Exome/Genechip_seq_pipeline/merged_plus10.intervals';
my $path2targetregions = "/home/groups/cardio/References/Exome/Genechip_seq_pipeline/$intervals.targetregions.bed";
my $path2SNP =   '/home/groups/cardio/References/SNP/dbsnp135.b37.vcf';
#my $path2omni = '/home/groups/cardio/References/SNP/1000G_omni2.5.b37.sites.vcf';
#my $path2hapmap = '/home/groups/cardio/References/SNP/hapmap_3.3.b37.sites.vcf';
my $path2ref = '/home/groups/cardio/References/Bwa_b37/hs37d5.fa';

#####################
#   Start BWA mem   #
#####################

my @fq = glob("$path2rawdata/*.fastq.gz"); # make array of fastq.gz

# check that there are the same number of forward and reverse files, and they have correct extensions
my @f_fq = grep(/1.fastq.gz$/i, @fq);
my @r_fq = grep( /2.fastq.gz$/i, @fq );
if (scalar(@f_fq) != scalar(@r_fq)) {
    rmdir "$path2bam";
    rmdir "$path2log";
    rmdir "$path2sam";
    die "*** error: There must be a 1.fastq.gz and 2.fastq.gz file for each sample\n\n";
}

my $eightfork_manager = Parallel::ForkManager->new(8);

for (my $i = 0; $i < @fq; $i = $i ++) {
    my @fq_pair = splice(@fq, $i, 2);
    $eightfork_manager->start and next;
    bwa_mem(join(" ", @fq_pair));
    $eightfork_manager->finish;
}

$eightfork_manager-> wait_all_children;

print "\n\n*** BWA mapping complete ***\n";
clock();

rmdir "$path2sam"; # remove tempSAMs directory

#####################
#   Sort Bamfile    #
#####################

my $tenfork_manager = Parallel::ForkManager->new(10);

my @bamfiles = glob("$path2bam/*unsorted.bam"); # make array of bam files
for (my $i = 0; $i < @bamfiles; $i++) {  # loop through them and pass to sortbam subroutine
    $tenfork_manager->start and next;
    sortbam($bamfiles[$i]);
    $tenfork_manager->finish;
}

$tenfork_manager-> wait_all_children;

print "\n\n*** Sorting initial bamfiles complete ***\n";
clock();

############################
#   Remove Duplicates      #
############################

my $fourfork_manager = Parallel::ForkManager->new(4);

my @sortedbamfiles = glob("$path2bam/*sorted.bam"); #make array of sorted.bam files
for (my $i = 0; $i < @sortedbamfiles; $i++) {
    $fourfork_manager->start and next;
    mark_dups($sortedbamfiles[$i]); # loop throught and pass to picard mark_dups subroutine
    $fourfork_manager->finish;
}

$fourfork_manager-> wait_all_children;

print "\n\n*** Removing duplicates complete ***\n";
clock();

#####################
#   Sort Bamfile    #
#####################

my $fsttenfork_manager = Parallel::ForkManager->new(10);

my @ddbamfiles = glob("$path2bam/*dedupped.bam"); #make array of removed_duplicates (dedupped) bam files
for (my $i = 0; $i < @ddbamfiles; $i++) { # loop through array of bamfiles and pass to sortbam subroutine
    $fsttenfork_manager->start and next;
    sortbam($ddbamfiles[$i]);
    $fsttenfork_manager->finish;
}

$fsttenfork_manager-> wait_all_children;

print "\n\n*** Sorting dedupped bamfiles complete ***\n";
clock();

#######################
#  Realign Bamfile    #
#######################

my $sixfork_manager = Parallel::ForkManager->new(6);

my @sort_ddbamfiles = glob("$path2bam/*.sorted.bam"); #make array of dedupped.sorted.bam files
for (my $i = 0; $i < @sort_ddbamfiles; $i++) {
    $sixfork_manager->start and next;
	indexbam($sort_ddbamfiles[$i]); # loop through array of bamfiles and pass to indexbam subroutine
	realigner_target_creator($sort_ddbamfiles[$i]); # pass to realigner_target_creator subroutine
	indel_realigner($sort_ddbamfiles[$i]); # pass to indel_realigner subroutine
    #~#    unlink ($sort_ddbamfiles[$i]); # cleanup: delete bam files
    #~#    unlink ("$sort_ddbamfiles[$i].bai"); # cleanup: delete bam index files
    $sixfork_manager->finish;
}

$sixfork_manager-> wait_all_children;

my @indel_intervals = glob("$path2bam/*.indel.intervals"); # loop through array of .intervals files and delete them
for (my $i = 0; $i < @indel_intervals; $i++) { unlink ("$indel_intervals[$i]") }


print "\n\n*** GATK indel realigner complete ***\n";
clock();

#####################################
# base quality score recalibration  #
#####################################

my $fstsixfork_manager = Parallel::ForkManager->new(6);

my @realigned = glob("$path2bam/*realigned.bam"); #make array of .realigned.bam files
for (my $i = 0; $i < @realigned; $i++) { # loop through array of bamfiles
    $fstsixfork_manager->start and next;
	base_recalibration($realigned[$i]); # pass to base_recalibration subroutine
    $fstsixfork_manager->finish;
}

$fstsixfork_manager-> wait_all_children;

my @grp = glob("$path2bam/*.grp"); # loop through array of .grp files and delete them
for (my $i = 0; $i < @grp; $i++) { unlink ("$grp[$i]") }

print "\n\n*** GATK BQSR complete ***\n";
clock();

#################################
# Depth of coverage calculation #
#################################

my $sndsixfork_manager = Parallel::ForkManager->new(6);

my @recalibrated = glob("$path2bam/*recalibrated.bam"); #make array of .recalibrated.bam files
for (my $i = 0; $i < @recalibrated; $i++) { # loop through array of bamfiles
    $sndsixfork_manager->start and next;
	coverage($recalibrated[$i]); # pass to coverage sub
    $sndsixfork_manager->finish;
}

$sndsixfork_manager-> wait_all_children;

print "\n\n*** Depth of coverage calculation complete ***\n";
clock();

cleanup();

################## SUBROUTINES ###################

sub clock {
	my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
	my @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
	my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
	my $year = 1900 + $yearOffset;
	my $theTime = "$hour:$minute:$second, $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year";
	print "*** $theTime ***\n\n";
}

sub usage {
	print "\n\n\n\t---------EXOME SEQUENCING ALIGNMENT PIPELINE------------------\n\n";
	print "\tUsage:\tExome_alignment_v1.pl [-fastq -cohort -exome]\n";
	print "\t-fastq - is required. The path to the directory containing the fastq files. NB: do not add trailing / \n";
	print "\t-cohort - is required. Define a name for the final VCF file / \n";
    print "\t-exome - is required. Define an exome used from truseq, nextera or sureselectV4 / \n";
	print "\tUse this script to align, realign and recalibrate a cohort of exome sequencing data. Requires paired end fastq.gz compressed reads\n";
	print "\tReads must be stored in your /path/2/Rawdata folder and labelled:\n\n";
	print "\t\tsamplename1.f.fastq.gz\n";
	print "\t\tsamplename1.r.fastq.gz\n";
	print "\t\tsamplename2.f.fastq.gz\n";
	print "\t\tsamplename2.r.fastq.gz\n";
	print "\t\tetc..  \n\n";
	print "\tReplace samplename with a unique ID (e.g. blood code)\n";
	print "\tContact r.bagnall\@centenary.org.au\n";
	print "\t-------------------------------------------------------------\n\n";
	exit;
}

sub bwa_mem {
	my ($current_read_pair) = shift(@_);
	my @read_pair = split(" ", $current_read_pair);
	my $current_samplename = substr $read_pair[0], ($offset + 9), -11; # i.e. get IO2 from $input/Rawdata/IO2.f.fastq.gz
	my $RG ='@RG'; #need this to be able to print out @RG
	print "\n\n*** Aligning $current_samplename reads ***\n";
	clock();

    my $bwa_mem = system("bwa mem -PM -t 3 -R '$RG\tID:$current_samplename\tSM:$current_samplename\tPL:ILLUMINA' $path2ref $current_read_pair > $path2sam/$current_samplename.sam");
	print "\n\n*** Converting $current_samplename sam to bam ***\n";
	clock();
	my $sam2bam = system("samtools view -bS -o $path2bam/$current_samplename.unsorted.bam $path2sam/$current_samplename.sam");
    unlink ("$path2sam/$current_samplename.sam") or die "Could not delete $current_samplename.sam file: <$!>\n"; # cleanup: delete samfile
}

sub sortbam {
	my ($current_bam) = shift(@_);
	my $current_name = substr $current_bam, ($offset + 10), -13; # eg $input/tempBAMs/IO2.unsorted.bam or $input/tempBAMs/IO2.dedupped.bam
	print "\n\n*** Sorting $current_name bam file ***\n";
	clock();
	my $sortbam = system("samtools sort -o $path2bam/$current_name.sorted.bam -T $path2bam/$current_name -O bam $current_bam ");
	unlink($current_bam) or die "Could not delete $current_bam\n"; # delete the unsorted bamfile
}

sub indexbam {
	my ($current_bam) = shift(@_);
	my $current_name = substr $current_bam, ($offset + 10), -4; # eg $input/tempBAMs/IO2.sorted.bam
	print "\n\n*** Indexing $current_name file ***\n";
	clock();
	my $indexbam = system("samtools index $current_bam");
}

sub mark_dups {
	my ($current_sortedbam) = shift(@_);
	my $current_name = substr $current_sortedbam, ($offset + 10), -11; # eg $input/tempBAMs/IO2.sorted.bam
	print "\n\n*** Removing duplicate reads in $current_name bam file ***\n";
	clock();
	my $mark_dups = system("java -jar /home/groups/cardio/Applications/picard-tools-1.106/picard-tools-1.106/MarkDuplicates.jar INPUT=$current_sortedbam OUTPUT=$path2bam/$current_name.dedupped.bam METRICS_FILE=$path2log/$current_name.dedup.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT");
	
    unlink($current_sortedbam) or die "Could not delete $current_sortedbam\n"; # delete the sorted bamfile, which has duplicates
}

sub realigner_target_creator {
	my ($current_bam) = shift(@_);# $input/tempBAMs/sample1.sorted.bam
    my $current_name = substr $current_bam, ($offset + 10), -11;
	print "\n\n*** Creating $current_name target indel sites ***\n";
	clock();
	my $realigner_target_creator = system("java -jar $path2gatk -T RealignerTargetCreator -R $path2ref -o $path2bam/$current_name.indel.intervals -I $current_bam -L $path2mergedexome50 -known $path2indel1kg -known $path2indelmills");
}

sub indel_realigner {
	my ($current_bam) = shift(@_);# $input/tempBAMs/sample1.sorted.bam
    my $current_name = substr $current_bam, ($offset + 10), -11;
	print "\n\n*** Realigning $current_name bam file ***\n";
	clock();
	my $local_realignment = system("java -jar $path2gatk -T IndelRealigner -R $path2ref -o $path2bam/$current_name.realigned.bam -I $current_bam -L $path2mergedexome50 -known $path2indel1kg -known $path2indelmills -targetIntervals $path2bam/$current_name.indel.intervals");
}

sub base_recalibration {
	my ($current_bam) = shift(@_);# $input/tempBAMs/sample1.realigned.bam
	my $current_name = substr $current_bam, ($offset + 10), -14;
	print "\n\n*** Performing pre base quality score recalibration for $current_name ***\n";
	clock();
	my $bqsr_pre = system("java -jar $path2gatk -T BaseRecalibrator -R $path2ref -I $current_bam -knownSites $path2SNP -knownSites $path2indel1kg -knownSites $path2indelmills -L $path2mergedexome50 -o $path2bam/$current_name.prerecal_data.grp");
	print "\n\n*** Performing post base quality score recalibration for $current_name ***\n";
	clock();
    my $bqsr_post = system("java -jar $path2gatk -T BaseRecalibrator -R $path2ref -I $current_bam -knownSites $path2SNP -knownSites $path2indel1kg -knownSites $path2indelmills -L $path2mergedexome50 -BQSR $path2bam/$current_name.prerecal_data.grp -o $path2bam/$current_name.postrecal_data.grp");
    print "\n\n*** Printing recalibration plots for $current_name ***\n";
	clock();
    my $plots = system("java -jar $path2gatk -T AnalyzeCovariates -R $path2ref -before $path2bam/$current_name.prerecal_data.grp -after $path2bam/$current_name.postrecal_data.grp -plots $path2bam/$current_name.postrecal.pdf");
	print "\n\n*** Printing recalibrated reads for $current_name ***\n";
	clock();
	my $bqsr = system("java -jar $path2gatk -T PrintReads -R $path2ref -I $current_bam -BQSR $path2bam/$current_name.prerecal_data.grp -L $path2mergedexome50 -o $path2bam/$current_name.recalibrated.bam");
    
	unlink ($current_bam); # cleanup: delete bam files
	unlink ("$input/tempBAMs/$current_name.realigned.bai"); # cleanup: delete bam index files
}

sub coverage {
	my ($current_bam) = shift(@_);# $input/tempBAMs/sample1.recalibrated.bam
	my $current_name = substr $current_bam, ($offset + 10), -17;
	print "\n\n*** Calculating coverage for $current_name ***\n";
	clock();
	my $coverage = system("samtools view -b $current_bam | coverageBed -abam stdin -b $path2targetregions -d | awk -F\"\t\" '{print\$5}' > $input/tempBAMs/$current_name.coverage");
    my $coverage1 = system("Rscript /home/groups/cardio/Applications/Rscripts/genechip_coverage.R $input/tempBAMs/$current_name.coverage");
    unlink ("$input/tempBAMs/$current_name.coverage");
    print "\n\n*** Finished $current_name ***\n\n";
}

sub cleanup{
    
    #create folder using $cohort, move contents of tempBAMs and tempLOGs
    
    mkdir ("$input/$cohort") or die "Unable to create new $cohort directory: <$!>\n";
    mkdir ("$input/$cohort/Logs") or die "Unable to create new $cohort Logs directory: <$!>\n";
    mkdir ("$input/$cohort/Bamfiles") or die "Unable to create new $cohort Bamfiles directory: <$!>\n";
    
    system("mv $path2bam/* $input/$cohort/Bamfiles/");
    system("mv $path2log/* $input/$cohort/Logs/");

    
    rmdir "$path2bam" or die "Unable to delete $path2bam folder: <$!>\n";
    rmdir "$path2log" or die "Unable to delete $path2log folder: <$!>\n";
}