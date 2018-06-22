#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;

my $AUGUSTUS = $ARGV[0];

if($ARGV[0] eq "-h" || $ARGV[0] eq "--help" || !$ARGV[0]){

print "usage: [Augustus output] [Protein-FASTA output] [GFF output]\n" and exit;
}

open (AUGUSTUS, "< $AUGUSTUS")  or die "Could not open file $AUGUSTUS $!";


my $count;
my %proteins;
my %gff;
my @gff;
my $sequence;
my @sequence;
my $ID;
while(my $line =<AUGUSTUS>){
        chop $line;
	if ($line !~ m/^\#/g){
		my $sequence;
		$line =~ s/CDS/exon/g;
		my @split = split("\t", $line);
 	 	push (@gff,$line);
	}
	if($line =~ m/protein sequence = \[([A-Z]+)/g){
	       push (@sequence,$1);
	}
	if($line =~ m/\# ([A-Z]{98})/g){
		push (@sequence,$1);
	}
	if($line =~ m/\# ([A-Z]+)\]/g){
		push (@sequence,$1);
	}


	if($line =~ m/^\# end/){
		$count++;
		foreach my $member (@gff){
			$member =~ s/ID=g[0-9]+/ID=g$count/g;
	                $member =~ s/Parent=g[0-9]+/Parent=g$count/g;
                }		
		
                if( $gff[0] =~ m/gene/g){
                        my @ID_1 = split("\t", $gff[0]);
                        $ID_1[8] =~ s/ID=//g;
                        $ID = join("_",$ID_1[0],$ID_1[8]);
                }
		$sequence = join("",@sequence);
		$sequence =~ s/ //g;
		
		if(exists $gff{$ID}){
			print "already exists ! and exit\n";
		}else{
			@{$gff{$ID}} = @gff;

		}
		$proteins{$ID}=$sequence;	
		@gff =();
		$sequence =();
		@sequence=();

	}
}

my $PROTEIN = $ARGV[1];
my $NEW_GFF = $ARGV[2];

open (PROTEIN, "> $PROTEIN")  or die "Could not open file $PROTEIN $!";
open (NEW_GFF, "> $NEW_GFF")  or die "Could not open file $NEW_GFF $!";

if (keys %proteins ne keys %gff){
	print "error: number of elements in Gff and Protein file do not match!\n" and exit;
}

$count = 0;
foreach my $keys(sort natural_p keys %proteins){
        $count++;
	my $new_key = $keys;
	$new_key =~ s/_g[0-9]+/_g$count/;
	my $tmp =$proteins{$keys};
	$tmp  =~ s/(.{1,80})/$1\n/g;
	print PROTEIN "\>$new_key\n$tmp\n";
}

$count = 0;
#use Data::Dumper; die Dumper \%gff;
foreach my $keys(sort natural_gff keys %gff){
	$count++;
	foreach(@{$gff{$keys}}){
		s/ID=g[0-9]+/ID=g$count/g;
		s/Parent=g[0-9]+/Parent=g$count/g;
	}	
	my $result = join("\n", @{$gff{$keys}});
	print NEW_GFF "$result\n";
}

########

sub natural_p{
	my @a = split("\_g",$a);
	my @b = split("\_g",$b);
	my $contig_a = $a[0];
        my $contig_b = $b[0];
	my $gene_a = $a[1];
	my $gene_b = $b[1];
        $contig_a =~ s/[^0-9*.]//g;
        $contig_b =~ s/[^0-9*.]//g;
	if ($contig_a != $contig_b ){
	       return $contig_a <=> $contig_b;
	}elsif($contig_a == $contig_b){
		return $gene_a <=> $gene_b;
	}



}

sub natural_gff{

        my @a = split("\_g",$a);
        my @b = split("\_g",$b);
        my $contig_a = $a[0];
        my $contig_b = $b[0];
        my $gene_a = $a[1];
        my $gene_b = $b[1];
        $contig_a =~ s/[^0-9*.]//g;
        $contig_b =~ s/[^0-9*.]//g;
        if ($contig_a != $contig_b ){
               return $contig_a <=> $contig_b;
        }elsif($contig_a == $contig_b){
                return $gene_a <=> $gene_b;
        }


}

