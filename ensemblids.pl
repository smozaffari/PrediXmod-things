#!usr/bin/perl
use strict;
use warnings;

my %gene;
my ($file1, $file2) =  @ARGV;
#input file first ,output file second
#fed in my pred_obs.R script.

open (LIST, "/group/im-lab/nas40t2/egamazon/annotation/gencode.v12.summary.gene") || die "nope: $!";
#open file that has ensembl ids and genes listed
while (my $line = <LIST>) {
    chomp $line;
    my @line = split("\t", $line);
    my @dec = split(/\./, $line[4]);
    my $ens = $dec[0];
    print @dec;
    $gene{$ens} = $line[5];
    #save the name of each gene to be called by ensembl id
}
close (LIST);

open(NEW, ">$file2") || die "nope: $!";

#my @file = `zcat ../haky/main/Transcriptome/GEUVADIS/E-GEUV-3/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz`;
#foreach my $line (@file) {
open (OBS, $file1) || die "nope: $!";
my $firstline = <OBS>;
print NEW ("Gene ",$firstline);
#read and write the first line without any changes since it indludes the individual IDs

while (my $line = <OBS>) {
    chomp $line;
    my @line = split(" ", $line);
    my @id = split(/\./, $line[0]);
    my $ens = $id[0];
    #if genename for ensembl id exists, output that gene name with the rest of gene expression for that gene
    if ($gene{$ens}) {
        print NEW (join (" ", $gene{$ens}, @line[1..$#line]), "\n");
    }
}
close (OBS);
close (NEW);
