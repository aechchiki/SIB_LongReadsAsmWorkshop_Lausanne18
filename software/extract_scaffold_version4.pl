#!/software/bin/perl
use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;


#changelog:
#09.01.2018 - fix perl bug with scalar
#28.03.2014 - fixing bug in print buffer which led to loss of long sequences - eschmid
#	    - adding DAZZLER option
#	    - adding STRICT  option
#01.05.2015 - add extracting of regions
#24.07.2015 - adding reverse mode for reads
#14.06.2016 - adding optional 3rd column for regions to define the extracted name
#	      adding optional 4th column for the directions "+ or -"
#	      changed from 0 based to 1 based system
my @blocks;


my $ID;
my $FASTA;
my $DAZZ;
my $STRICT;
my $COORDS;
my $verb;
my $reverse;
my $help=  "usage: -i [query] -f [FASTA] [option] 
                query \= FASTA ID's to extract (wihtout leading \>)
                FASTA \= FASTA sequence from which to extract
                option: -d \= specific for reads from the DAZZLER dextract tool
                        -s \= strict treatment of FASTA header (no empty spaces) 
			-r \= inverse extraction, will extract everything BUT the provided reads
                        -c \= extract regions from the sequences (either 3 columns or 6)
                                  scaffold <\t> start <\t> stop <\t> name(optional) <\t> direction(optional, +/-) <\t> introns 
                                  1st  base =  1
                                  last base = -1
			          introns e.g.: 163123,163560..163561,163710..163713

				IMPORTANT: NO ! allowed in the sequence name

			-v \= activate verbose mode (needs module Data::Dumper)\n";
		
&argus;
sub argus{
        my %arguments;
        getopts('rdvsci:f:', \%arguments);
        if ($arguments{h}){
                die print $help;
	}
	if(!$arguments{i}){
		die print $help;
	}else{
	        $ID = $arguments{i};
	}
	if(!$arguments{f}){
		die print $help;
	}else{
		$FASTA = $arguments{f};
	}
        if ($arguments{d}){
                $DAZZ = 1;
        }
        else{
                $DAZZ = 0;
        }
        if ($arguments{s}){
                $STRICT = 1;
        }else{
                $STRICT = 0;
        }
	if($arguments{v}){
		$verb = "1";	
		use Data::Dumper;
	}else{
		$verb = "0";
	}
	if($arguments{r}){
		$reverse = "1";
	}else{
		$reverse = "0";
	}
	if($arguments{c}){
		$COORDS = 1;
	}else{
		$COORDS = 0;
	}
}

print STDERR "  query: $ID
	        FASTA: $FASTA
		-d: $DAZZ
		-s: $STRICT
		-c: $COORDS
		-r: $reverse
		-v: $verb
		\n";

open (QUERY, "< $ID")  or die "Could not open file 'query' $!";
open (SEQUENCE, "< $FASTA")  or die "Could not open file 'FASTA' $!";

#first lets quickly grep all FASTA ID's before we start selecting anything
my $quick_ID = `grep "^>" $FASTA`;
chomp $quick_ID;
my @all_IDs = split("\n",$quick_ID);
my %target_IDs;
my $counter_target = 0;
my $counter_keep = 0;
my $counter_FASTA = 0;

foreach my $tmp_ID(@all_IDs){
	$tmp_ID =~ s/^\>//;
	
	 if($STRICT eq 1){
                        $tmp_ID =~ s/(^\S+).*/$1/;
         }elsif($DAZZ eq 1){
                        $tmp_ID =~ s/ .*$//;
         }
	$target_IDs{$tmp_ID} =$tmp_ID;
	$counter_FASTA++;
}

print STDERR "Number of FASTA sequences: $counter_FASTA\n";

if($verb eq "1"){
	print STDERR Dumper \%target_IDs;
}


#now lets generate a hash with all the ones we are actually interested in
my %query;
my @region;

while(<QUERY>){
	 chomp;
	
	 my @line = split("\t",$_);
	 my $ID = $line[0];
	 if($COORDS eq 1){
		my $start   = $line[1];
		my $stop    = $line[2];
		my $region;
		if($line[3] && $line[4] && $line[5] && $line[5] ne 0){
			 $region  = join("!", $start,$stop,$line[3],$line[4],$line[5]);
		#		print STDERR "Name & direction & introns encountered and incorporated\n";
		}elsif($line[3] && $line[4]){
		#		print STDERR "Name & direction encountered and incorporated\n";
			$region  = join("!", $start,$stop,$line[3],$line[4]);
		}else{
			$region  = join("!", $start,$stop);
				print STDERR "Only 3 valid columns detected, no Name & direction or region encountered and incorporated\n";
		}
	 	push(@{$query{$ID}{$region}},$region);
	 }else{
		$query{$ID}="1";
	}
	$counter_target++;
}

print STDERR "Number of target sequences: $counter_target\n";

if($verb eq "1"){
	print STDERR Dumper \%query;
}

#now we check whether they match and generate the ones we want to keep
my %keep_IDS;
foreach my $entries (keys %query){
	if($target_IDs{$entries}){
		$keep_IDS{$entries}=$query{$entries};
	}
	$counter_keep++;
}


print STDERR "Number of sequences matching: $counter_keep\n";

if($verb eq "1"){
	print STDERR Dumper \%keep_IDS;
}

#this reads now the entire FASTA file but only the ones which match with our previous selection
#thereby the RAM consumption is way lower as if we read everything into a hash and then treat it
my @seq;
my $seq;
my $old_ID;
my $new_ID;
my %parsed_fasta;
while (my $line = <SEQUENCE>){
	chomp $line;
	if($line =~ /^\>(.+)/){
		$new_ID = $1;
		$seq = join("\n",@seq);
		if($STRICT eq 1){
			$new_ID =~ s/(^\S+).*/$1/;
			#print STDERR "$new_ID\n";
		}elsif($DAZZ eq 1){
			$new_ID =~ s/ .*$//;
			#print STDERR "$new_ID\n";
		}

		#push the new sequence into hash if matching
		if($reverse == "0"){
			if($old_ID && $keep_IDS{$old_ID}){
				$parsed_fasta{$old_ID}=$seq;
				#print STDERR "keeping: $old_ID\n";
			}
		}elsif($reverse == "1"){
			if($old_ID && !$keep_IDS{$old_ID}){
				$parsed_fasta{$old_ID}=$seq;
				#print STDERR "keeping: $old_ID\n";
			}
		}else{
			print "reverse option problem encountered\n" and die;
		}
		$old_ID = $new_ID;
		@seq=();
	}else{
		push(@seq,$line);
	}

}

$seq = join("\n",@seq);
if($STRICT eq 1){ 
     $new_ID =~ s/(^\S+).*/$1/;
}elsif($DAZZ eq 1){
    $new_ID =~ s/ .*$//;
}
if($reverse == "0"){
	if($old_ID && $keep_IDS{$old_ID}){
		$parsed_fasta{$old_ID}=$seq;
		#print STDERR "keeping: $old_ID\n";
	}
}elsif($reverse == "1"){
	if($old_ID && !$keep_IDS{$old_ID}){
		$parsed_fasta{$old_ID}=$seq;
		#print STDERR "keeping: $old_ID\n";
	}
}else{
	print "reverse option problem encountered\n" and die;
}       




                                                                       



my $counter_extracted = 0;
foreach my $key (sort keys %parsed_fasta){
	my $fasta = $parsed_fasta{$key};
	$fasta =~ s/\n//g;
	if($COORDS eq 1){
		foreach my $reg(sort keys %{$keep_IDS{$key}}){
			my @region= split("!",$reg);
			my $start = $region[0]-1;
			my $stop = -1;
			if ($region[1] ne -1){
				$stop = $region[1] ;
			}
			#sanity check:
			if( $stop ne -1 && $stop < $start){
				die "not valid bed file, start !< stop : @region \n" ;
			}
			my $fasta_region;
			if ($stop eq -1){
				$fasta_region = substr($fasta,$start);
			}else{
				my $end = $stop - $start;
				$fasta_region = substr($fasta,$start,$end);
			}
			my $new_name;	
			my $frag_length = length($fasta_region);
			#adding a name if present
			if($region[2]){
				$new_name = $region[2];
			}else{
				$new_name = ("$key\/$start\_$stop");
			}
                        #adding introns
                        #000018F MAN_ANNOT_PH    CDS     162750  165309  .       +       0       ID=ClasWithCRJE1c18 
                        #000180F MAN_ANNOT_PH    CDS     354     3484    .       -       0       ID=ClasWithCRJE1c180 
                        #
                        #
                        my @introns;  
                        if($region[4]){
                                my $new_fasta;
				my $exon_start = 0;
				my $exon_stop;
				my $extract;
                                @introns = split(",",$region[4]);
                                foreach my $deletions(@introns){
#					print STDERR "detected intron $deletions\n";
                                        my @coords = split(/\.\./, $deletions);
					my $intron_start = $coords[0] -1 - $start;
					my $intron_stop;
					#problem: I might have single nucleotide ones too...
					if($coords[1]){
						$intron_stop = $coords[1] - $start ;
					}else{
						$intron_stop = $coords[0] - $start;
					}
					$exon_stop  = $intron_start; 
#					print STDERR "extracting exon from $exon_start to $exon_stop\n";
					my $extract = substr($fasta_region,$exon_start,$exon_stop - $exon_start);
#					print STDERR "region extracted: $extract\n";
					#next we add it to the previous CDS extract
					$new_fasta .= $extract;
					#so we pass on the start of the exon to the next one
					$exon_start = $intron_stop ;
                                }
				#now lets get the last intron_stop to end of the read
				$extract = substr($fasta_region,$exon_start) or die "Could not extract $exon_start from @region from $new_name\n";
				$new_fasta .= $extract;
				$fasta_region = $new_fasta;
                        }
			
			#adding strang direction
			if($region[3] && $region[3] eq "-"){
			        $fasta_region = scalar reverse $fasta_region;
		                $fasta_region =~tr/ATGC/TACG/;
                		$fasta_region =~tr/atgc/tacg/;

			}elsif($region[3] && $region[3] ne  "+"){
				print STDERR "not valid bed file, directon not provided\n" and die;
			}
			$fasta_region =~ s/(.{1,80})/$1\n/g;
			print STDOUT "\>$new_name\n$fasta_region\n";
		}
	}else{
		$fasta =~ s/(.{1,80})/$1\n/g;
	        print STDOUT "\>$key\n$fasta";
	}
	$counter_extracted++;
}	

print STDERR "Number of extracted entries: $counter_extracted\n";

