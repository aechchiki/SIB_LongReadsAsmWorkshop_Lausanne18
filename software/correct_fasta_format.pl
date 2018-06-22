#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
use Data::Dumper;
my @blocks;
my $FASTA;
my $renum;
my $UP;
my $One;
my $divide;
my $prefix;
my $msa;
my $clean;
my $explode;
my $pacbio;
my $quiver;
my $verbose;
my $ncbi;
my $leading;
my $trailing;
#added quiver polishing cutoff and dump of alternative sequences
my $help = "Help:
		 -i [FASTA file] [options]

		 -r is optional and will rename all header by 1,2,3....
		 -p is optional and provides a prefix for the renaming
		 -u is optional as well and will replace lowercase DNA-letters by uppercase
		 -d is optional and will verify that all sequences are a multiple of 3 (removing bases)
		 -o is optional and will present the entire sequence in one line
		 -m is optional and will extract sequences from a MSA alignment
		 -n NCBI depositing conform (no seq <200bp, gap presence at start/end not checked automaticall- see other options)
		 -t trim trailing Ns
		 -l trim leading Ns
		 -b make header PacBio conform for blasr (>header/0/0_X, X being the length of the sequence)
		 -c is optional and will remove any FASTA with * (STOP CODON)
		 -q quiver downstream, remove contigs with >50% not polished (small letters)
		 -v verbose, provides additional information in STDERR
		 -x is optional and will explode it --> each entry a file with its name $!\n";

my %arguments;
getopts('i:hxrcnuodp:mbqvlt', \%arguments);
if ($arguments{h}){
	print $help and die;
}
if($arguments{n}){
	$ncbi = 1;
	print STDERR "ncbi cutoff 200bp chosen\n";

}
if ($arguments{r}){
    $renum  = "1";
	print STDERR "renumbering chosen\n";
}
if ($arguments{u}){
    $UP  = "1";
	print STDERR "upscaling chosen \n";
}
if ($arguments{v}){
    $verbose  = "1";
	print STDERR "verbose reporting chosen \n";
}
if($arguments{o}){
	$One = "1";
	print STDERR "sequence output in one line \n";
}else{
	$One = "0";
}
if($arguments{c}){
	$clean = "1";
	print STDERR "removing * containing sequences \n";
}else{
	$clean = "0";
}
if($arguments{d}){
	$divide = "1";
	print STDERR "sequence will be made a multiple of 3\n";
}else{
	$divide = "0";
}
if($arguments{b}){
	$pacbio = "1";
	print STDERR "header will be made PacBio conform\n";
}
if($arguments{q}){
	if($arguments{u}){
		print STDERR "ERROR:options q and u should/cant be combined !\n" and die;
	}
	$quiver = "1";
	print STDERR "sequences with <50% polished will be removed\n";
}
if($arguments{x}){
	$explode = "1";
	print STDERR "each FASTA will give rise to a new file\n";
}else{
	$explode = "0";
}
if(!$arguments{i}){
	print $help and die;
}else{
	$FASTA = $arguments{i};
}
if(!$arguments{p}){
	$prefix = "sequence";
}else{
	$prefix = $arguments{p};
	print STDERR "manually prefix $prefix chosen\n";
}
if($arguments{m}){
    $msa = "1";
    print STDERR "assuming file is MSA and extracting sequences\n";
}else{
    $msa = "0";
}
if($arguments{t}){
  	$trailing = 1;
    print STDERR "trailing gaps are removed\n";
}
if($arguments{l}){
  	$leading = 1;
    print STDERR "leading gaps are removed\n";
}
if($explode && $pacbio){
	print STDERR "ERROR: PacBio format and explode cant be chosen together\n" and die;
}


open (SEQUENCE, "< $FASTA")  or die "Could not open file 'FASTA' $!";


my @seq;
my $seq;
my $old_ID;
my $new_ID;
my %parsed_fasta;
my $counter =0;
my $header;
my @fasta;
my $fasta;

if($verbose && $quiver){
	print STDERR "Dumping non-polished contigs into dump.fasta \n";
        open (DUMP, "> dump.fasta") or print "ERROR: could no write to file dump.fasta \n" and die;
}



while (my $line = <SEQUENCE>){
	chomp $line;
	if($line =~ /^\>(.+)/ && $old_ID){
		$new_ID = $1;
		$header = $old_ID;
		$fasta = join("",@fasta);
		if($trailing){
			$fasta =~ s/[N]+$//;
		}
		if($leading){
			$fasta =~ s/^[N]+//;
		}
		if($fasta !~ m/\w+/){
			print STDERR "ERROR: $header had empty sequence and was removed!\n";
			$old_ID = $new_ID;
			@fasta=();
			next;
		}
		$counter++;
		verify_fasta($header,$fasta,"0");
		if(!$parsed_fasta{$header}){
			$parsed_fasta{$header} = 1;
		}else{
			print STDERR "WARNING: you might have duplicated FASTA entries, $header already encountered\n";
		}
		#######finished , just cleaning up
		$old_ID = $new_ID;
		@fasta=();
	}elsif($line =~ /^\>\s+$/){
		print STDERR "ERROR: empty header encountered and removed!\n";
		@fasta=();
		next;
	}elsif($line =~ /^\>(.+)/ && !$old_ID){
		print STDERR "This should be first FASTA $line\n";
		$old_ID = $1;
	}else{
		push(@fasta,$line);
	}
}


#cleaning up the last element
$counter++;
$header = $old_ID;
$fasta = join("\n",@fasta);
if($trailing){
	$fasta =~ s/[N]+$//;
}
if($leading){
	$fasta =~ s/^[N]+//;
}
if($fasta !~ m/\w+/){
	 print STDERR "ERROR: $header had empty sequence and was removed!\n";
         $old_ID = $new_ID;
         @fasta=();
         next;
 }
verify_fasta($header,$fasta,"0");

sub verify_fasta {
	##here we expect a proper fasta containing some sort of header and a fasta sequence
	my ($header,$fasta,$return) = @_;
	my $seqLength = length($fasta);
	if(($ncbi && $seqLength <200 )|| $seqLength==0){return;}
	if(defined $verbose){
		print STDERR "$header\t";
	}
	if(defined $quiver){
		#now, if >50% of a sequence is in lower-case = not enough evidence for polishing, we remove it
		my $fasta2 = $fasta;
		$fasta2 =~ s/[A-Z]+//g;
		if(defined $verbose){
			print STDERR "fasta length: $seqLength bp\t";
		}
		my $seqLength_low = length($fasta2)-2;
		if(defined $verbose){
			print STDERR "lower-case: $seqLength_low bp\t";
		}
		if($seqLength_low != 0){
			my $ratio = int(($seqLength_low/$seqLength*100)+0.5);
			if(defined $verbose){
				print STDERR "polished: $ratio %\n";
			}
			if($ratio>=50){
				$fasta =~ s/(.{1,80})/$1\n/g;
				print DUMP ">$header\n$fasta\n";
				return;
			}else{
				if(defined $verbose){
					print STDERR "polished: $ratio %\n";
				}
			}
		}
	}
	if($clean ==1 && $fasta=~m/\*/){
		next;
	}
	if( $msa == 1){
		$fasta =~ s/-//g;
	}
	if( defined $renum ){
		$header = join("",$prefix,$counter);
	}
	if(defined $pacbio){
		my $tmp = join("/",$header,"0","0_$seqLength");
		$header = $tmp;
	}
	if($UP){
		$fasta =~ tr/a-z/A-Z/;
	}
	$header =~ s/^ //;
	$header =~ s/\s+$//;
	$header =~ s/\s+/_/g;
	$fasta =~ s/\s+//g;
	if($divide eq "1"){
		if(length($fasta) % 3 != 0){
			until(length($fasta) % 3 == 0){
				chop $fasta;
			}
		}
	}
	if($One ne "1"){
		$fasta =~ s/(.{1,80})/$1\n/g;
	}
	if($return == 1){
		my @result = join("\n",$header,$fasta);
		return @result;
		next;
	}
	if($explode==0){
		print STDOUT  ">$header\n$fasta";
	}else{
		my $out=join("",$header,".fasta");
		$out =~ s/^\>//;
		open (FILE, "> $out") or print "ERROR: could no write to file $out \n" and die;
		print FILE  ">$header\n$fasta";
		close FILE;
	}
}
