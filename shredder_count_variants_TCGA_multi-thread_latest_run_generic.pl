#!/usr/bin/perl -w

use strict;
use warnings;
use threads;

###########################################
# Function Prototypes
###########################################

sub parse_file($$@);
sub find_repeat_block(\%\@\@$$$$);
sub count_repeats(\%$$$);

###########################################
# Initialisations
###########################################

my $max_num_of_threads = 8;
my %thread_samples= ();
my $thread;
my @thread_list;
my $portion = 0;
my $CRITERIA = "4TTAGGGu";

###########################################
# Main
###########################################


my $OUT = "telomere_variants_count_trimmed_${CRITERIA}.txt";


#print "$#ARGV\n";

my %hash = ();

if($#ARGV ne 0){
	print "File containing repeat units require\n";
	exit;
}

my $pattern = "";

open(PATTERN, "<", "$ARGV[0]") or die "cannot open < $ARGV[0]: $!";

while (my $line = <PATTERN>){
	chomp $line;
	$line =~ s/^\s*(.*)\s*$/$1/;
	$line =~ s/\s+/\|/g;
	if($pattern =~ /.+/){
		$pattern .= "|";
	}
	$pattern .= $line;
}

close PATTERN;

my $pattern_inv = join("|", reverse(split(/\|/,reverse($pattern))));

$pattern_inv =~ tr/ACGT/TGCA/;

my $pattern_combined = "$pattern|$pattern_inv";

#my $files = `ls *.bam`;
my $files = `ls *4TTAGGGu.ini.*.bam`;

my $dir = `pwd`;

chomp $files;

my @files = split("\n", $files);

open(OUT, ">", "$OUT") or die "cannot open < $OUT: $!";

my @keys = split(/\|/, $pattern);
my @keys_rc = split(/\|/, $pattern_inv);

print OUT "Sample_ID,Sample_Type,Repeat,Strand,Count\n";

foreach my $file (@files){
	my @input = split("\n", `samtools view $file | cut -f10,11`);
	my $name;
	my $type;

#	if($file =~ /.*(.RR\d+)/){
#		$name = $1;
#	}
#	elsif($file =~ /_([ACGT]{6})_/){
#		$name = $1;
#	}
#	els
	#if($file =~ /^[^_]*_([^_]*)_/){
	#	$name = $1;
	#	if($file =~ /^.*Tumour.*/){
	#		$type = 2;
	#	} elsif($file =~ /^.*Normal.*/){
	#		$type = 1;
	#	}else{
	#		$type = 0;
	#	}

	#}
	#else{
		$name = $file;
	#}

	#$hash{$name}{$type}=parse_file($file, $pattern_combined, @input);

	$thread = threads->create(\&parse_file, $file, $pattern_combined, @input);
	$thread_samples{$thread->tid()}{"name"}=$name;
	$thread_samples{$thread->tid()}{"type"}=$type;

	print "$name\n";

	while(scalar(threads->list(threads::running)) >= $max_num_of_threads){
		foreach(threads->list(threads::joinable)){
			$hash{$thread_samples{$_->tid()}{"name"}}{$thread_samples{$_->tid()}{"type"}} = $_->join();
		}
#		sleep(1);
	}

	if(scalar(threads->list(threads::joinable)) >= 10)
	{
		foreach(threads->list(threads::joinable)){
			$hash{$thread_samples{$_->tid()}{"name"}}{$thread_samples{$_->tid()}{"type"}} = $_->join();
		}
	}


}

foreach(threads->list()){
	$hash{$thread_samples{$_->tid()}{"name"}}{$thread_samples{$_->tid()}{"type"}} = $_->join();
}

foreach my $ID (sort {$a cmp $b} keys %hash)
{
        if(defined $hash{$ID}{'1'} && defined $hash{$ID}{'2'}){
	        #print OUT "$ID";

        	foreach my $type ('1', '2')
        	{
		        my $index = 0;

		        foreach my $key (@keys)
		        {
					if(defined $hash{$ID}{$type}{$key}){
						print OUT "$ID,$type,$key,G,$hash{$ID}{$type}{$key}\n";
					}else{
						print OUT "$ID,$type,$key,G,0\n";
					}
					if(defined $hash{$ID}{$type}{$keys_rc[$index]}){
						print OUT "$ID,$type,$key,C,$hash{$ID}{$type}{$keys_rc[$index]}\n";
					}else{
						print OUT "$ID,$type,$key,C,0\n";
					}
       		        $index++;
        		}
        	}
        }
        #print OUT "\n";
}


close OUT;


# TEST

#parse_file("test", "TTAGGG|GTAGGG", "AGAGTTAGGGTTAGGGATTAGGGCGTAGGGCCA\tABCDABCDABCDABCDABCDABCDABCDABCDA");

###########################################
# Fuction definitions
###########################################

sub parse_file($$@){
	my ($file_name, $repeat_pattern ,@input) = @_;
	my %hash = ();
	my $seq = "";
	my $qual = "";

	my $other = "";
	my $repeat = "";
	my @other_s = ();
	my @other_q = ();


	#open(DUMP, ">", "${file_name}_other.txt") or die "cannot open > ${file_name}_other_short.txt: $!";

	foreach my $line (@input){
		chomp $line;
		($seq, $qual) = split("\t", $line);
		chomp $seq;
		chomp $qual;
#		if($portion > 0){
#			$seq = substr $seq, 0, $portion;
#			$qual = substr $qual, 0, $portion;
#		}

		my $count = 0;
		my $phred_score = 0;
		my $i = 0;
		my @qual_array = split("", $qual);

		while($i < 6){
			$phred_score += (ord($qual_array[$i])-33);
			$i++;
		}

		while($i < length $qual && ($phred_score/6) >= 30){
			#print "$i\n";
			$phred_score -= (ord($qual_array[$i-6])-33);
			$phred_score += (ord($qual_array[$i])-33);
			$i++;
		}

		if($i < length $qual){
			$seq = substr $seq, 0, $i;
			$qual = substr $qual, 0, $i;
		}

		if($i >= 24){
			find_repeat_block(%hash, @other_s, @other_q, $seq, $qual, $repeat_pattern, "3");
		}
	}


	#for(my $i=0; $i < scalar @other_s; $i++){
	#	my $q_score = 0;
	#	my $j = 0;
	#	my @o_seq = ();
	#	my @o_qual = ();
	#	#my @window = ();
	#	#my $window_qual = 0;
	#	my @window_qual = ();
	#
	#	print DUMP $other_s[$i]."\t".$other_q[$i]."\n";
	#
	#	if(length($other_s[$i]) < 6){
	#		foreach (split(//, $other_q[$i])){
	#			$q_score += (ord($_)-33);
	#			$j++;
	#		}
	#		$q_score /= $j;
	#
	#		#print STDERR "$other_s[$i]\t$j\n";
	#
	#		if($q_score >= 20){
	#			print DUMP $other_s[$i] . "\t" . $other_q[$i] . "\t2\n";
	#			$hash{"other"}{2}+= $j;
	#		}else{
	#			print DUMP $other_s[$i] . "\t" . $other_q[$i] . "\t1\n";
	#			$hash{"other"}{1}+= $j;
	#		}
	#	}else{
	#		print DUMP $other_s[$i] . "\t" . $other_q[$i] . "\t3\n";
	#
	#		@o_seq = split(//, $other_s[$i]);
	#		@o_qual =split(//, $other_q[$i]);
	#
	#		#print STDERR length($other_s[$i])."\t$#o_qual\n";
	#
	#		for(my $k=0; $k<6; $k++){
	#			$q_score += (ord($o_qual[$k])-33);
	#		}
	#		$q_score /= 6;
	#		$window_qual[0] = $q_score;
	#
	#		my $wqc = 0;
	#
	#		for(my $l=6; $l<scalar @o_qual; $l++){
	#			$q_score -= ((ord($o_qual[$l-6])-33)/6);
	#			$q_score += ((ord($o_qual[$l])-33)/6);
	#			$window_qual[$wqc] = $q_score;
	#			$wqc++;
	#		}
	#
	#		for(my $m=0; $m<scalar @o_qual; $m++){
	#			my $q = 0;
	#			my $count = 0;
	#			my $start = $m-5;
	#			my $end = $m;
	#
	#			if($start < 0){$start = 0;}
	#			if($end > $#window_qual){$end = $#window_qual;}
	#
	#			for(my $h=$start; $h <= $end; $h++){
	#
	#				if((ord($window_qual[$h])-33)>=20){
	#					$q++;
	#				} else {
	#					$q--;
	#				}
	#			}
	#
	#			if($q >= 0){
	#				$hash{"other"}{4}++;
	#			}else{
	#				$hash{"other"}{3}++;
	#			}
	#
	#		}
	#	}
	#}

	#close(DUMP);

	return \%hash;
}


sub find_repeat_block(\%\@\@$$$$){
	my ($hash_ref, $other_seq_ref, $other_qual_ref, $seq, $qual, $pattern, $flag) = @_;
	#my %hash = %$hash_ref;
	#my @other_seq = @$other_seq_ref;
	#my @other_qual = @$other_qual_ref;

# $flag: 0=middle of read, 1=start of read, 2=end of read, 3=whole read #

	if($seq =~ /$pattern/){
		$seq =~ /^(.*?)((?:$pattern)+)(.*?)$/;

		my $o1 = $1;
		my $p = $2;
		my $o2 = $3;

		if($o1 =~ /.+/){
			my $o1_q = substr($qual, 0, length($o1));
			if($flag =~ /[13]/){
				find_repeat_block(%$hash_ref, @$other_seq_ref, @$other_qual_ref, $o1, $o1_q, $pattern, "1");
			} else {
				find_repeat_block(%$hash_ref, @$other_seq_ref, @$other_qual_ref, $o1, $o1_q, $pattern, "0");
			}
		}

		if($o2 =~ /.+/){
			my $o2_q = substr($qual, -length($o2));

			if($flag =~ /[23]/){
				find_repeat_block(%$hash_ref, @$other_seq_ref, @$other_qual_ref, $o2, $o2_q, $pattern, "2");
			}else{
				find_repeat_block(%$hash_ref, @$other_seq_ref, @$other_qual_ref, $o2, $o2_q, $pattern, "0");
			}
		}

		my $p_q = substr($qual, length($o1), length($p));

		count_repeats(%$hash_ref, $p, $p_q, $pattern);

		return;

	}
	#else {
	#	if($flag =~ /[12]/ && length($seq) < 6){
	#		$hash_ref->{"other"}->{0}+= length($seq);
	#		return;
	#	}else{
	#		push(@$other_seq_ref, $seq);
	#		push(@$other_qual_ref, $qual);
	#		return;
	#	}
	#}

}

sub count_repeats(\%$$$){
	my ($hash_ref, $seq, $qual, $pattern) = @_;
	my $repeat;
	my $i;
	my $phred_score;

	while($seq =~ /^($pattern)(?=(?:$pattern|$))(.*)/){
		$phred_score = 100;
		#$phred_score = 0;

		$repeat = $1;
		$seq = $2;

		foreach my $x (split("", substr($qual, 0, length($repeat)))){
			if($phred_score > (ord($x)-33)){
				$phred_score = (ord($x)-33);
			}

			#$phred_score += (ord($x)-33);
			#$hash_ref->{$repeat}->{"qual"}{$i}{int((ord($x)-33)/int(10))}++;
			#$qual_temp[$i]=(ord($x)-33);
			#$i++;
		}

		#$phred_score /= $i;

		$hash_ref->{$repeat}++;

	}
}
