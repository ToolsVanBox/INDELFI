# INDELFI
This script filters low quiality INDEL calls from VCF files.  



# Dependency

The followig packages need to be installed under GUIX.

bedtools
vcftools



# Before running (script requirement)

This is a stand-alone perl script.  Except for the dependencies, no installation required.  


# Usage:

	perl INDELFI.pl -i [input VCF] -s [column test sample] -c [column control sample]


##Required parameters:
    -i | input VCF (including the PATH)
    -s | column number, which contains the genotype field of the TEST sample in the VCF (0-based)
    -c | column number, which contains the genotype field of the CONTROL sample in the VCF (0-based)

##Optional parameters:
    -q | minimal GATK quality score for filtering (Default: 250)
    -RD | minimal read depth at a variant position in both TEST and CONTROL sample to pass the filter (Default: 20)
    -VAF | minimal variant allele frequency in the TEST sample for a variant to pass the filter (Default: 0.1)
    -GQ | sample-specific genotype quality as determined by GATK (Default: 99)
    -f | flanks of control sample to include in the filtering process (Default: 100 bp)
    -mq | minimal mapping quality (MQ) score (Default: 60)
    -bl | Option to add a blacklist for filtering (in VCF format)
    -call | Option to add a GATK callableLoci file for filtering (in BED format)
    -o | output DIR PATH (if not specified, will write in the input folder)



# Output

VCF file containing INDELS that passed the specified quality filters.  
