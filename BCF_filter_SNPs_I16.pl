#!/usr/perl
use strict;
use warnings;
use Getopt::Long;
use List::Util qw(shuffle);

#usage:
#cat file.vcf | Clean_VCF.pl
#cat VCF_From_SamtoolsPileup.vcf |perl Clean_VCF.pl 
#to get the alternative bases from two bases in codon 12 of KRAS
#samtools mpileup -uf $REFERENCE $BAM -r 12:25289550-25289551 | bcftools view - | perl BCF_filter_SNPs_I16.pl | grep -v '^##' > $BAM.KRAScodon12.vcf

#my $allele_freq_cutoff = 0.05;
#my $allele_freq_cutoff = 0.01;
#my $allele_freq_cutoff = $ARGV[0];
my $allele_freq_cutoff = -1;
#my $rawdepth_cutoff = 30;
my $rawdepth_cutoff = -1;


while(<STDIN>){
    if ($_ =~ /^#CHROM/){
	chomp;
        print $_ . "\trawdepth\tref_reads\talt_reads\tallele_freq\n";
        next;
    }


    if ($_ =~ /^#/){
	print $_;
	next;
    }
    chomp;
    my @values = split (/\t/,$_);
    #alt_base is pos 5 == values[4]
    $values[4] =~ s/,X//g;
    $values[4] =~ s/X/\./g;
    my $info = $values[7];
    $info =~  /DP=([0-9,]+);/;
    my $rawdepth = $1;
    $info =~  /I16=([0-9,]+)/;
    #print $1 . "\n";
    my ($ref_forward,$ref_reverse,$alt_forward,$alt_reverse,$REST) = split (/,/,$1);
    #print $ref_forward . "\n";;
    my $ref_count = $ref_forward + $ref_reverse;
    my $alt_count = $alt_forward + $alt_reverse;
    my $allele_freq = 0; #maybe should have 'NA' but this is only used to separate
    if ($ref_count + $alt_count > 0){
        $allele_freq = $alt_count / ($alt_count + $ref_count);
    }
    if ($allele_freq > $allele_freq_cutoff && $rawdepth > $rawdepth_cutoff){
    #if(1){
	push(@values,$rawdepth);
	push(@values,$ref_count);
	push(@values,$alt_count);
	push(@values,$allele_freq);
	
	print join("\t", @values) . "\n";
    }

}

print STDERR "finished\n";


