#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Std;
use Data::Dumper;

############################################################
##############							####################
#############   FASTA scaffold analysis ####################
############							####################
###########   Emanuel Schmid  			####################
##########								####################
#########   02.10.2013  				####################
########	update 12.12.2013           ####################
#######changed from array to hash: 14.11.2014 ##############
######	bugfixed 08.01.14 				####################
#####   completely new organized Jan2015####################
####	contig measure altered	Nov16  #####################
###		07.07.17 add plotting method . #####################
###										####################
##										####################	
#										####################		
############################################################


my(@blocks, @c_size_list, @c_i_size_list,@split,$fasta_out);
my ($size,$thres) = (0,100);
my $print_s;
my $plot;
my $file;
my $contigize = 0;
my $min_gap = 10;

&argus;
sub argus{
	my %arguments;
	getopts('g:s:t:lpchi:', \%arguments);
	if ($arguments{h} || !%arguments){
		die "
		usage: assesV(*).pl -i test.fasta

		Note: this slurps everything in memory, be aware!

		-i DNA seq in *.fasta (mandatory)
		-s estimated genome size(MB)
		-t size threshold (default:100bp) 
		-l to print a list with ordered sequence sizes (default: no)
		-p plot the size distribution (needs R in ENV)
		-g define minimal considered gap-length (default:10bp)$!\n";
        }
	$file = $arguments{i};
    if (!$arguments{s}){
    	print STDERR " WARNING: no genome size provided 1000000 (1Mbp) set$!";
    	$size = 1000000;
	}else{ 
    	$size = ($arguments{s}*1000000);
    }
	if ($arguments{l}){
		$print_s = 1;
	}else{
		$print_s = 0;
	}
	if($arguments{g}){
		$min_gap = $arguments{g};
	}
	if($arguments{p}){
		$plot = 1;
		#we need to generate the file, can purge it later if output not desired
		$print_s = 1;
		#verify that R is present
		system "Rscript --version" and print "ERROR: R not found in environment\n" and die;
	}else{
		$plot = 0;
	}
	if($arguments{c}){
		$contigize = 1;
		print STDERR "scaffolds are broken into contigs for contig stats section only$! "
	}
	if (!$arguments{t}){
		print STDERR "No length cut-off for minimal sequence length chosen - default 100bp used$!";
		$thres = 100;
	}else{
		$thres = $arguments{t};
	}
}


#if proper name remove extension and make a proper file, else just concatenate
my $hist_file;
my $name = $file;
$name =~ s/\.\w+$//;
print "$name\n";
if($name){
	$hist_file = join(".",$name,"hist");
}else{
	$hist_file = join(".",$file,"hist");
}

#open (INPUT,"< $file") or print "no input defined\n" and die;
if($print_s==1){
	open (OUTPUT,"> $hist_file") or print "no input defined\n" and die;
}

my $data = slurp($file);


#all needed variables
my %contig_db;
my %sequences_db;
my $contig_counter=0;
my $tot_filled;
my $totn_filled;
my @Ns;
my @filled;
my $GC;
my $tot_ATGC;
my $scaf_ATGC;
my $cont_ATGC;
my @FASTA;
my $FASTA;
#now lets get already the first line, first ID
my $counter;


my @collection = split(">",$data);

shift @collection;
foreach my $FASTA (@collection){
	$FASTA =~ s/.*\n//;
	$FASTA =~ s/\R//g;
	#print "$FASTA\n";
	$counter++;
	my $ID = $counter;
	my $MyLength = length $FASTA;	
	if($MyLength < $thres){
		next;
	}else{
		#1st gaps
		$sequences_db{$ID}=$MyLength;
		while($FASTA =~ /(N+)/g){
			if(length($1)>=$min_gap){
       			push (@Ns, length($1));
       		}
        }
		while($FASTA =~ /([g|c|a|t]+)/g){
				push (@filled, length($1));
        }
		#next GC content
		my $DNA = $FASTA;
		$GC += $DNA =~ tr/[G|C]/[G|C]/;
		$GC += $DNA =~ tr/[g|c]/[g|c]/;
		$tot_filled += $DNA =~ tr/[g|c|a|t]/[g|c|a|t]/;
		$tot_ATGC += $DNA =~ tr/[G|C|A|T]/[G|C|A|T]/;		
		$tot_ATGC += $DNA =~ tr/[g|c|a|t]/[g|c|a|t]/;
		#now fragmenting in case we have scaffolds
		if($contigize == 1){
			if( $FASTA =~ m/[Nn]/) {
				#now we want to know for the scaffolded ones as well what the contigs looks like
				#only define gaps for a length>10
				my $gap_string = $min_gap x "N";
				my @split_contig = split($gap_string,$FASTA);
				foreach my $subcontigs (@split_contig){
					$contig_counter++;
					my $tmp_ID = join("_","\>Contig",$contig_counter);
					my $tmp_length= length $subcontigs;
					$contig_db{$tmp_ID}= $tmp_length;
				}
			}else{
				$contig_db{$ID}= $MyLength;
			}
		}
	}
}


####now lets produce a list of the sizes (eg. necessary for NG50)
my @t_SIZE_CUT_scaf;
my @t_SIZE_CUT_cont;


foreach (keys %sequences_db){
	push (@t_SIZE_CUT_scaf, $sequences_db{$_});
}

foreach (keys %contig_db){
    push (@t_SIZE_CUT_cont, $contig_db{$_});
}


#######sort them by size  (descending)
my @SIZE_CUT_scaf = sort { $b <=> $a } @t_SIZE_CUT_scaf;
my @SIZE_CUT_cont = sort { $b <=> $a } @t_SIZE_CUT_cont;

if($print_s==1){
	foreach (sort keys %sequences_db){
		print OUTPUT "$sequences_db{$_}\n";
	}
	close OUTPUT;
}
if($plot==1 && @SIZE_CUT_scaf>1){
	system "Rscript /home/eschmid/tools/R_scripts/plotAssessHist.R" ;
}

##how many sequences?
my $scaf_sequences  = @SIZE_CUT_scaf;
my $cont_sequences  = @SIZE_CUT_cont;

##SOmebody said N50?
my $scaf_N50  = &N50(@SIZE_CUT_scaf);
my $cont_N50  = &N50(@SIZE_CUT_cont);

##bah better NG50
my $scaf_NG50  = &NG50(@SIZE_CUT_scaf);
my $cont_NG50  = &NG50(@SIZE_CUT_cont);

#calculate the GC
my $GC_ratio = ($GC /$tot_ATGC)*100;
my $tot_pfilled  = ($tot_filled/&sum(@SIZE_CUT_scaf)) * 100;

## percentage of counts in regard to total counts
my @result_distr = &distr(@SIZE_CUT_scaf);
my @presult_distr = @result_distr;
foreach my $member(@presult_distr){
        $member=($member/$scaf_sequences*100);
}


# percentage of bp-sum in regard to total bp size
my @csum_distr = &collect_distr(@SIZE_CUT_scaf);
my @ptcsum_distr= @csum_distr;
foreach my $member(@ptcsum_distr){
	$member=($member/&sum(@SIZE_CUT_scaf)*100);
}

# percentage of bp-sum in regard to genome size
my @pgcsum_distr = @csum_distr;
foreach my $member(@pgcsum_distr){
        $member=($member/$size*100);
}



##get all regions which have been filled gaps
my ($gap_sum,$gap_average,$gap_bpsum) = (0,0,0);
my $filled_sum = @filled;
if (@Ns > 0){
    @Ns = sort {$a <=> $b} @Ns;
    #sum the bp of gaps
    foreach (@Ns){
        $gap_bpsum += $_;
    }
    $gap_sum = @Ns;
	$gap_average = $gap_bpsum / $gap_sum;
;
}else{
        print "    no gaps detected\n";
}


###Print all results
print "\n\n\n";
print "\n\n\n";
print "    INFORMATION:
    A threshold of $thres bp was applied for ALL measures and genome size estimated $size Bp\n";
print "    Gaps are considered from a size of  $min_gap bp on, below ignored\n";
print "    lower case bases are often either low confidence bases (PB) or filled gaps (PBHoney,SSPACE)\n";
print "\n\n\n";
print "    overall stats\n\n";
printf  "%-40s %12.2f\n", "    GC ratio (in \%):", $GC_ratio ;
printf "%-40s %12.0f\n", "    Shortest sequence length:", &min(@SIZE_CUT_scaf);
printf "%-40s %12.0f\n", "    Longest sequence length:",&max(@SIZE_CUT_scaf);
printf "%-40s %12.0f\n", "    bp-sum of all sequences :", &sum(@SIZE_CUT_scaf);
printf "%-40s %12.0f\n", "    Number of sequences >" . $thres . ":", $scaf_sequences;
printf "%-40s %12.0f\n", "    N50 sequence length :", $scaf_N50;
printf "%-40s %12.0f\n", "    NG50 sequence length :", $scaf_NG50;
printf "%-30s %14.0s %7s %10s %7.2f\n", "    lower Case (\#\/bps/\%):","", $filled_sum, $tot_filled , $tot_pfilled ;
printf "%-30s %13.0s  %7s %10s %7.2f\n", "   gaps (\#\/bps/\%):","",$gap_sum, $gap_bpsum , ($gap_bpsum/&sum(@SIZE_CUT_scaf))*100  ;
print "\n";
print "-----" x 20 . "\n";
print "\n";
if(@SIZE_CUT_cont && $contigize ==1){
	print " Breaking scaffolds into contigs stats\n";
	print "\n";
	printf "%-40s %12.0f\n", "    Shortest sequence length:", &min(@SIZE_CUT_cont);
	printf "%-40s %12.0f\n", "    Longest sequence length:",&max(@SIZE_CUT_cont);
	printf "%-40s %12.0f\n", "    bp-sum of all sequences :", &sum(@SIZE_CUT_cont);
	printf "%-40s %12.0f\n", "    Number of sequences >" . $thres . ":", $cont_sequences;
	printf "%-40s %12.0f\n", "    N50 sequence length :", $cont_N50;
	printf "%-40s %12.0f\n", "    NG50 sequence length :", $cont_NG50;
	print "\n";
	print "-----" x 20 . "\n";
}
print "\n";
print "    overall sequence size distribution (\# of contigs \& scaffolds)\n\n";
printf "%-25s %10s %18s \n","",  "    absolute    \|", "% total sequ.";
printf "%-25s %10s %20.2f\n", "    100-500:",$result_distr[0], $presult_distr[0];
printf "%-25s %10s %20.2f\n", "    500-1000:",$result_distr[1], $presult_distr[1];
printf "%-25s %10s %20.2f\n", "    1000-10000:",$result_distr[2], $presult_distr[2];
printf "%-25s %10s %20.2f\n", "    10000-100000:",$result_distr[3], $presult_distr[3];
printf "%-25s %10s %20.2f\n", "    100000-1000000:",$result_distr[4], $presult_distr[4];
printf "%-25s %10s %20.2f\n", "    >1000000:",$result_distr[5], $presult_distr[5];
print  "\n\n\n";
print "    overall sequence size distribution (bp SUM of contigs \& scaffolds)\n\n";
printf "%-25s %10s %20s %15s\n", "", "    absolute Bps \|", " % total Bps \|", "% genome size";
printf "%-25s %10s %22.2f %17.2f\n", "    100-500:", $csum_distr[0],$ptcsum_distr[0], $pgcsum_distr[0]; 
printf "%-25s %10s %22.2f %17.2f\n", "    500-1000:", $csum_distr[1],$ptcsum_distr[1],$pgcsum_distr[1];
printf "%-25s %10s %22.2f %17.2f\n", "    1000-10000:", $csum_distr[2],$ptcsum_distr[2],$pgcsum_distr[2];
printf "%-25s %10s %22.2f %17.2f\n", "    10000-100000:", $csum_distr[3],$ptcsum_distr[3],$pgcsum_distr[3];
printf "%-25s %10s %22.2f %17.2f\n", "    100000-1000000:", $csum_distr[4],$ptcsum_distr[4],$pgcsum_distr[4];
printf "%-25s %10s %22.2f %17.2f\n", "    >1000000:", $csum_distr[5],$ptcsum_distr[5],$pgcsum_distr[5];
print  "\n";
print  "\n";



############################sub-routines############################################################

#slurp in thx to https://perlmaven.com/slurp
sub slurp {
    my $file = shift;
    open my $fh, '<', $file or die;
    local $/ = undef;
    my $cont = <$fh>;
    close $fh;
    return $cont;
}


#distribution of size - how many sequences in the size bin
sub distr{
	my ($cut100,$cut500,$cut1000,$cut10000,$cut100000,$cut1000000) = (0,0,0,0,0,0);;
	foreach (@_){
	    if ($_ >= 100 and $_ <500 ){
	        $cut100++;
	    }
	    if ($_>= 500 and $_ <1000){
	        $cut500++;
	    }
	    if ($_ >= 1000 and $_<10000){
	        $cut1000++;
	    } 
	    if ($_ >= 10000 and $_ <100000){
	        $cut10000++;
	    }
	    if ($_ >= 100000 and $_ < 1000000){
	        $cut100000++;
	    }
	    if ($_ >= 1000000){
	        $cut1000000++;
	    }
	}
	return $cut100,$cut500,$cut1000,$cut10000,$cut100000,$cut1000000;
	     
}

#distribution of size - what total BP-sum in each bin
sub collect_distr{
	my ($cut100,$cut500,$cut1000,$cut10000,$cut100000,$cut1000000) = (0,0,0,0,0,0);;
	foreach (@_){
		if ($_ >= 100 and $_ <500 ){
			$cut100 += $_;
		}
		if ($_>= 500 and $_ <1000){
			$cut500 += $_;
		}
		if ($_ >= 1000 and $_<10000){
			$cut1000 += $_;
		}
		if ($_ >= 10000 and $_ <100000){
			$cut10000 += $_;
		}
		if ($_ >= 100000 and $_ < 1000000){
			$cut100000 += $_;
		}
		if ($_ >= 1000000){
			$cut1000000 += $_;
		}
	}

	return $cut100,$cut500,$cut1000,$cut10000,$cut100000,$cut1000000;

}


#calculate Median of obtained sequences:
sub median{
    my $count = @_;
    if ($count % 2){
        return $_ [int($count / 2)];
    }else {
        return ($_ [int($count /2)] + $_ [int($count /2)-1])/2;
    }
}


#caluclate minimum of sequences
sub min{
    my $min = $_[$#_];
    return $min;
}

#calculate maximum of sequences
sub max{
    my $max = $_[0];
    return $max;
}

#calculate BP-sum of sequences
sub sum{
    my $sub_sum = 0;
	foreach my $seq(@_){
		$sub_sum += $seq;
	}
    return $sub_sum;
}	


sub mean{
        my $count = @_;
        my $mean = &sum(@_) /$count;
        return $mean;
}

 #N50
sub N50{
    my ($counts,$sum) =(0,0);
	my $tot_sum = &sum(@_);
    if (@_ > 1){
		until ($sum >= $tot_sum/2 ){
            $sum += $_[$counts];
           	$counts++;
        }
		if($_[$counts] > 0){
	    	return ($_[$counts-1]);
		}else{
			return "0";
		}
	}else{
		return "0";
	}
}


#NG50
sub NG50{
    my ($counts,$sum) =(0,0);
    if (&sum(@_) > $size/2 ){
       	until ($sum >= $size/2 ){
           	$sum += $_[$counts];
           	$counts++;
        	}
       	return $_[$counts-1];
    }else { 
		return 0; 
	}
}
