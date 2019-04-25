#!usr/bin/perl -w
use strict;
use Getopt::Long qw(GetOptions);

sub show_version{
    warn <<END;
    Version 1.0.0 (copyright R. van Boxtel 2017)
    Version 1.0.1 (modified to avoid permission problems. R.Oka 2018)
    Version 1.2.0 (skips lines without info in ref/sample columns. R.Oka 2018)
    Version 1.3.0 (fixed paths to binaries removed. A. Huber 2018)
    Version 1.4.0 (X chromosome is no longer removed)
    Version 1.5.0 (GQ Control sample fixed at 60. R. van Boxtel 2019)
END
    exit;
}

sub usage{
    warn <<END;

    Usage:
    Run by typing: perl INDELFI.pl -i [input VCF] -s [column test sample] -c [column control sample]

    Dependencies:
    The following packages need to be installed under GUIX.

          guixr package -i bedtools
          guixr package -i vcftools

    Your GUIX profile should be loaded.

    Required parameters:
    -i | input VCF (including the PATH)
    -s | column number, which contains the genotype field of the TEST sample in the VCF (0-based)
    -c | column number, which contains the genotype field of the CONTROL sample in the VCF (0-based)

    Optional parameters:
    -q | minimal GATK quality score for filtering (Default: 250)
    -RD | minimal read depth at a variant position in both TEST and CONTROL sample to pass the filter (Default: 20)
    -VAF | minimal variant allele frequency in the TEST sample for a variant to pass the filter (Default: 0.1)
    -GQ | sample-specific genotype quality as determined by GATK (Default: 99)
    -f | flanks of control sample to include in the filtering process (Default: 100 bp)
    -mq | minimal mapping quality (MQ) score (Default: 60)
    -bl | Option to add a blacklist for filtering (in VCF format)
    -call | Option to add a GATK callableLoci file for filtering (in BED format)
    -o | output DIR PATH (if not specified, will write in the input folder)

END
    exit;
}

die usage() if @ARGV == 0;

my $vcf;
my $sample;
my $control;
my $qual = 250; #optional parameter
my $min_rd = 20;
my $min_vaf = 0.1; #optional parameter
my $GQ = 99; #optional parameter
my $flank = 100;
my $min_mq = 60;
my $blacklist = 0;
my $callable = 0;
my $out_dir; #optional parameter
my $help;
my $version;
GetOptions(
    'i=s' => \$vcf,
    'o=s' => \$out_dir,
    's=s' => \$sample,
    'c=s' => \$control,
    'q=s' => \$qual,
    'RD=s' => \$min_rd,
    'VAF=s' => \$min_vaf,
    'GQ=s' => \$GQ,
    'f=s' => \$flank,
    'mq=s' => \$min_mq,
    'bl=s' => \$blacklist,
    'call=s' => \$callable,
    'help' => \$help,
    'version' => \$version,
);

if ($help){
    usage();
}

if ($version){
    show_version();
}


#--get names of the samples
my $sample_name = "Error";
my $control_name = "Error";

my %filter = ();

open (IN, $vcf);
while (my $line = <IN>){
    chomp $line;
    next if ($line =~ m/\#\#/);
    my @data = split("\t", $line);
    if ($line =~ m/\#/){
	$sample_name = $data[$sample];
	$control_name = $data[$control];
    }
    else{
	my $start = $data[1] - $flank;

	# if start coordinate is negative, set it to 0
	if( $start < 0 ){
		$start = 0 ;
	}

	next if ($data[$control] eq "./.");
        next if ($data[$sample] eq "./.");
		# skipping line without ref and sample information (ROka, 20/06/2018)

	my $end = $data[1] + $flank;
	my @info_control = split(":", $data[$control]); #get specs for control sample

	next if ($info_control[1] eq "\.");
		# skipping line if no allele count is available (ROka 20/06/2018)

	if ($info_control[0] ne "0/0"){
	    $filter{$line} = $data[0]."\t".$start."\t".$end;
	}
	else{
	    my @alleles_control = split(",", $info_control[1]);
	    my $RD_control = 0;
	    foreach my $allele (@alleles_control){
		if ($RD_control == 0){
		    $RD_control = $allele;
		}
		else{
		    $RD_control = $RD_control + $allele;
		}
	    }
	    my $AA_control = $RD_control - $alleles_control[0];
	    if ($AA_control > 0){
		$filter{$line} = $data[0]."\t".$start."\t".$end;
	    }
	}
    }
}
close IN;

#--generate ASC-specific output file
my @path = split("/", $vcf);
pop @path;

	# if out_dir not specified in the parameter, output in the input directory (R. Oka 2018)
if ( !$out_dir ){
	$out_dir = join("/", @path)."/".$sample_name."/";
}

mkdir $out_dir;

my $out = $out_dir.$sample_name."_".$control_name."_QUAL".$qual."_RD".$min_rd."_VAF".$min_vaf."_GQsample".$GQ."_GQcontrol60_flank".$flank."_MQ".$min_mq."_INDELs_autosomal_noEvidenceControl_temp.vcf";
my $bed = $out_dir.$sample_name."_flank".$flank.".bed";

#--generate filter
open (BED, '>'.$bed);
foreach my $line (sort(keys%filter)){
    print BED $filter{$line}, "\n";
}
close BED;

my $bed_sorted = $bed;
$bed_sorted =~ s/.bed/_sorted.bed/;
system ("sortBed -i $bed > $bed_sorted");
unlink $bed;

my $bed_merged = $bed_sorted;
$bed_merged =~ s/_sorted.bed/_merged.bed/;
system ("mergeBed -i $bed_sorted > $bed_merged");
unlink $bed_sorted;

#--filtering
open (IN, $vcf);
open (OUT, '>'.$out);

while (my $line = <IN>){
    chomp $line;
    my @data = split("\t", $line);
    if ($line =~ m/\#/){
	print OUT $line, "\n";
    }
    else{
	next if ($data[0] eq "Y" or $data[0] eq "MT"); #get autosomal genome;
	next if ($data[2] =~ m/rs/ and $data[2] !~ m/COSM/); #remove SNP_ids w/o COSMIC id
	next if ($data[4] =~ m/\,/);
	next if ($data[5] < $qual); #minimal quality requirement
	next if ($data[6] ne "PASS"); #only consider passed variants

# -- Filter MQ score
	my @info = split(";", $data[7]);
	my $mapq = 0;

	foreach my $element (@info){
	    if ($element =~ m/MQ\=/){
		my @mapqs = split("=", $element);
		$mapq = $mapqs[1];
	    }
	}
	next if ($mapq < $min_mq);
	next if ($data[$control] eq "./.");

	my @info_control = split(":", $data[$control]); #get specs for control sample
	next if ($info_control[0] ne "0/0"); #remove calls in ref sample
	next if ($info_control[3] < 60); #remove calls with low GQ-score in control sample
	my @alleles_control = split(",", $info_control[1]);
	my $RD_control = 0;
	foreach my $allele (@alleles_control){
	    if ($RD_control == 0){
	    $RD_control = $allele;
	    }
	    else{
		$RD_control = $RD_control + $allele;
	    }
	}
	my $AA_control = $RD_control - $alleles_control[0];
	next if ($RD_control < $min_rd); #remove positions with less than 20 informative reads in ref sample
	next if ($AA_control > 0); #remove positions with any evidence in ref sample

	my @info_sample = split(":", $data[$sample]); #get specs for test sample
	next if ($info_sample[0] eq "0/0" or $info_sample[0] eq "./."); #remove lines w/o call in test sample
	next if ($info_sample[3] < $GQ); #remove calls with low GQ-score in test sample
	my @alleles_sample = split(",", $info_sample[1]);
	my $RD_sample = 0;
	foreach my $allele (@alleles_sample){
	    if ($RD_sample == 0){
	    $RD_sample = $allele;
	    }
	    else{
	    $RD_sample = $RD_sample + $allele;
	    }
	}
	next if ($RD_sample < $min_rd); #remove positions with less than 20 informative reads in test sample
	my $AA = $RD_sample - $alleles_sample[0];
	my $VAF = $AA/$RD_sample;

	next if ($AA == 0); #remove calls with no alternative reads, but a call in sample
	next if ($VAF < $min_vaf); #remove samples with low VAF

	print OUT $line, "\n";
    }
}
close IN;
close OUT;

my $final_out = $out;
$final_out =~ s/_temp.vcf/.vcf/;
system ("intersectBed -header -v -a $out -b $bed_merged > $final_out");

unlink $out;
unlink $bed_merged;

if ($blacklist ne 0){
    print "Processing blacklist: ", $blacklist, "\n";
    my $out_blacklist = $final_out;
    $out_blacklist =~ s/.vcf/_noBlacklist/;
    system ("vcftools --gzvcf $final_out --exclude-positions $blacklist --recode --recode-INFO-all --out $out_blacklist");

    my $out_blacklist_new = $out_blacklist.".vcf";
    $final_out = $out_blacklist_new;
    my $log = $out_blacklist.".log";

    $out_blacklist = $out_blacklist.".recode.vcf";
    rename $out_blacklist, $out_blacklist_new;

    unlink $log;
}

if ($callable ne 0){
    print "Intersecting callableLoci file: ", $callable, "\n";
    my $out_call = $final_out;
    $out_call =~ s/.vcf/_CALLABLE.vcf/;
    system ("intersectBed -header -u -a $final_out -b $callable > $out_call");

}
