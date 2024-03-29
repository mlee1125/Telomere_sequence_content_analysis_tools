#!/usr/bin/perl -w
use strict;
use warnings;

my $MAX_RL = 150;
my $BIN_SIZE = 50;
my $WINDOW_SIZE = 6;

my $num_args = $#ARGV + 1;
if($num_args != 1){
    print STDERR "\nUsage: sort_telo_ctrls.pl <fastq>\n";
    exit;
}

my $input_file = $ARGV[0];
unless($input_file =~ /^(.*)\.(?:fastq|fq)$/){
    print STDERR "\nUsage: sort_telo_ctrls.pl <fastq>\n";
    exit;
}

my $file_prefix = $1;
my $out_file = "$file_prefix.variants_analysis_6.csv";

open(INPUT, "<", $input_file) or die "Cannot open < $input_file: $!\n";
open(OUT, ">", $out_file) or die "Cannot open > $out_file: $!\n";

my %hash = ();
my $phase = "";
my $end = "";
my $id = "";
my $seq = "";
my $qual = "";
my @seq = "";
my $seq_len = "";
my $strand_type = "";

while(my $line = <INPUT>){
	unless($line =~ /^@.*/){
		print STDERR "Error: File not in fastq format, header incorrect!\n";
		exit;
	}
	$id = $line;
	$strand_type = "";

	$seq = <INPUT>;
	unless($seq =~ /^[ACGTN]+$/){
		print STDERR "Error: File not in fastq format, seq incorrect!\n";
		exit;
	}
	
	$line = <INPUT>;
	unless($line =~ /^\+$/){
		print STDERR "Error: File not in fastq format, missing +!\n";
		exit;
	}
	$qual = <INPUT>;
	unless($qual =~ /^[!-N]+$/){
		print STDERR "Error: File not in fastq format, qual incorrect!\n";
		exit;
	}
	
	if($seq =~ /^([ACGT]{0,5})(?:TCAGGG|ACAGGG|CCAGGG|GCAGGG|TAAGGG|TGAGGG|TTAGGG|TCCGGG|TCGGGG|TCTGGG|TCAAGG|TCACGG|TCATGG|TCATGAG|TCAGCG|TCAGTG|TCAGGA|TCAGGC|TCAGGT)+([ACGT]{0,5})$/){
		if(defined $1){
			$phase = length($1);
		}else{
			$phase = 0;
		}
		if(defined $2){
			$end = length($2);
		}else{
			$end = 0;
		}
		$strand_type = "G";
	}elsif($seq =~ /^([ACGT]{0,5})(?:CCCTGA|CCCTGT|CCCTGG|CCCTGC|CCCTTA|CCCTCA|CCCTAA|CCCGGA|CCCCGA|CCCAGA|CCTTGA|CCGTGA|CCATGA|CTCTGA|CGCTGA|CACTGA|TCCTGA|GCCTGA|ACCTGA)+([ACGT]{0,5})$/){
		if(defined $1){
			$phase = length($1);
		}else{
			$phase = 0;
		}
		if(defined $2){
			$end = length($2);
		}else{
			$end = 0;
		}
		$strand_type = "C";
	}elsif($seq =~ /^([ACGT]{0,11}?)(?:(?:[ACGT](?:GGGTCA|GGGACA|GGGCCA|GGGGCA|GGGTAA|GGGTGA|GGGTTA|GGGTCC|GGGTCG|GGGTCT|AGGTCA|CGGTCA|TGGTCA|GAGTCA|GCGTCA|GTGTCA|GGATCA|GGCTCA|GGTTCA)(?:GGGTCA|GGGACA|GGGCCA|GGGGCA|GGGTAA|GGGTGA|GGGTTA|GGGTCC|GGGTCG|GGGTCT|AGGTCA|CGGTCA|TGGTCA|GAGTCA|GCGTCA|GTGTCA|GGATCA|GGCTCA|GGTTCA))|(?:(?:GGGTCA|GGGACA|GGGCCA|GGGGCA|GGGTAA|GGGTGA|GGGTTA|GGGTCC|GGGTCG|GGGTCT|AGGTCA|CGGTCA|TGGTCA|GAGTCA|GCGTCA|GTGTCA|GGATCA|GGCTCA|GGTTCA)(?:GGGTC|GGGAC|GGGCC|GGGGC|GGGTA|GGGTG|GGGTT|AGGTC|CGGTC|TGGTC|GAGTC|GCGTC|GTGTC|GGATC|GGCTC|GGTTC))|(?:(?:GACCC|GTCCC|GGCCC|GCCCC|TACCC|CACCC|AACCC|GACCT|GACCG|GACCA|GACTC|GACGC|GACAC|GATCC|GAGCC|GAACC)(?:TGACCC|TGTCCC|TGGCCC|TGCCCC|TTACCC|TCACCC|TAACCC|GGACCC|CGACCC|AGACCC|TGACCT|TGACCG|TGACCA|TGACTC|TGACGC|TGACAC|TGATCC|TGAGCC|TGAACC))|(?:(?:TGACCC|TGTCCC|TGGCCC|TGCCCC|TTACCC|TCACCC|TAACCC|GGACCC|CGACCC|AGACCC|TGACCT|TGACCG|TGACCA|TGACTC|TGACGC|TGACAC|TGATCC|TGAGCC|TGAACC)(?:TGACCC|TGTCCC|TGGCCC|TGCCCC|TTACCC|TCACCC|TAACCC|GGACCC|CGACCC|AGACCC|TGACCT|TGACCG|TGACCA|TGACTC|TGACGC|TGACAC|TGATCC|TGAGCC|TGAACC)[ACGT]))+([ACGT]{0,11}?)$/){
		if(defined $1){
			$phase = length($1);
			$seq = substr($seq, $phase, length($seq));
			$qual = substr($qual, $phase, length($qual));
		}else{
			$phase = 0;
		}
		if(defined $2){
			#print STDERR "$seq\n$2\n\n";
			$end = length($2);
			$seq = substr($seq, 0, length($seq) - ($end+1));
			$qual = substr($qual, 0, length($qual) - ($end+1));
		}else{
			$end = 0;
		}
		$strand_type = "mixed";
	}else{
		print STDERR "Problem!\n";
		next;
	}
	
	$seq_len = length($seq);
	
	if($strand_type eq "G"){
		chomp $seq;
		@seq = split("", $seq);
		for(my $i = $phase; $i < (scalar @seq - $end); $i++){
			$hash{"G"}{int($seq_len / $BIN_SIZE)}{int($i / $WINDOW_SIZE)}{($i - $phase)%6}{$seq[$i]}++;
		}
	}elsif($strand_type eq "C"){
		chomp $seq;
		@seq = split("", $seq);
		for(my $i = $phase; $i < (scalar @seq - $end); $i++){
			$hash{"C"}{int($seq_len / $BIN_SIZE)}{int($i / $WINDOW_SIZE)}{($i - $phase)%6}{$seq[$i]}++;
		}
	}elsif($strand_type eq "mixed"){
		my $sub_seq = "";
		my $global_pos = $phase;
		while($seq =~ /^((?:[ACGT](?:GGGTCA|GGGACA|GGGCCA|GGGGCA|GGGTAA|GGGTGA|GGGTTA|GGGTCC|GGGTCG|GGGTCT|AGGTCA|CGGTCA|TGGTCA|GAGTCA|GCGTCA|GTGTCA|GGATCA|GGCTCA|GGTTCA)(?:GGGTCA|GGGACA|GGGCCA|GGGGCA|GGGTAA|GGGTGA|GGGTTA|GGGTCC|GGGTCG|GGGTCT|AGGTCA|CGGTCA|TGGTCA|GAGTCA|GCGTCA|GTGTCA|GGATCA|GGCTCA|GGTTCA))|(?:(?:GGGTCA|GGGACA|GGGCCA|GGGGCA|GGGTAA|GGGTGA|GGGTTA|GGGTCC|GGGTCG|GGGTCT|AGGTCA|CGGTCA|TGGTCA|GAGTCA|GCGTCA|GTGTCA|GGATCA|GGCTCA|GGTTCA)(?:GGGTC|GGGAC|GGGCC|GGGGC|GGGTA|GGGTG|GGGTT|AGGTC|CGGTC|TGGTC|GAGTC|GCGTC|GTGTC|GGATC|GGCTC|GGTTC))|(?:(?:GACCC|GTCCC|GGCCC|GCCCC|TACCC|CACCC|AACCC|GACCT|GACCG|GACCA|GACTC|GACGC|GACAC|GATCC|GAGCC|GAACC)(?:TGACCC|TGTCCC|TGGCCC|TGCCCC|TTACCC|TCACCC|TAACCC|GGACCC|CGACCC|AGACCC|TGACCT|TGACCG|TGACCA|TGACTC|TGACGC|TGACAC|TGATCC|TGAGCC|TGAACC))|(?:(?:TGACCC|TGTCCC|TGGCCC|TGCCCC|TTACCC|TCACCC|TAACCC|GGACCC|CGACCC|AGACCC|TGACCT|TGACCG|TGACCA|TGACTC|TGACGC|TGACAC|TGATCC|TGAGCC|TGAACC)(?:TGACCC|TGTCCC|TGGCCC|TGCCCC|TTACCC|TCACCC|TAACCC|GGACCC|CGACCC|AGACCC|TGACCT|TGACCG|TGACCA|TGACTC|TGACGC|TGACAC|TGATCC|TGAGCC|TGAACC)[ACGT]))((?:(?:[ACGT](?:GGGTCA|GGGACA|GGGCCA|GGGGCA|GGGTAA|GGGTGA|GGGTTA|GGGTCC|GGGTCG|GGGTCT|AGGTCA|CGGTCA|TGGTCA|GAGTCA|GCGTCA|GTGTCA|GGATCA|GGCTCA|GGTTCA)(?:GGGTCA|GGGACA|GGGCCA|GGGGCA|GGGTAA|GGGTGA|GGGTTA|GGGTCC|GGGTCG|GGGTCT|AGGTCA|CGGTCA|TGGTCA|GAGTCA|GCGTCA|GTGTCA|GGATCA|GGCTCA|GGTTCA))|(?:(?:GGGTCA|GGGACA|GGGCCA|GGGGCA|GGGTAA|GGGTGA|GGGTTA|GGGTCC|GGGTCG|GGGTCT|AGGTCA|CGGTCA|TGGTCA|GAGTCA|GCGTCA|GTGTCA|GGATCA|GGCTCA|GGTTCA)(?:GGGTC|GGGAC|GGGCC|GGGGC|GGGTA|GGGTG|GGGTT|AGGTC|CGGTC|TGGTC|GAGTC|GCGTC|GTGTC|GGATC|GGCTC|GGTTC))|(?:(?:GACCC|GTCCC|GGCCC|GCCCC|TACCC|CACCC|AACCC|GACCT|GACCG|GACCA|GACTC|GACGC|GACAC|GATCC|GAGCC|GAACC)(?:TGACCC|TGTCCC|TGGCCC|TGCCCC|TTACCC|TCACCC|TAACCC|GGACCC|CGACCC|AGACCC|TGACCT|TGACCG|TGACCA|TGACTC|TGACGC|TGACAC|TGATCC|TGAGCC|TGAACC))|(?:(?:TGACCC|TGTCCC|TGGCCC|TGCCCC|TTACCC|TCACCC|TAACCC|GGACCC|CGACCC|AGACCC|TGACCT|TGACCG|TGACCA|TGACTC|TGACGC|TGACAC|TGATCC|TGAGCC|TGAACC)(?:TGACCC|TGTCCC|TGGCCC|TGCCCC|TTACCC|TCACCC|TAACCC|GGACCC|CGACCC|AGACCC|TGACCT|TGACCG|TGACCA|TGACTC|TGACGC|TGACAC|TGATCC|TGAGCC|TGAACC)[ACGT]))*)$/){
			$sub_seq = $1;
			if(defined $2){
				$seq = $2;
			}else{
				$seq = "";
			}

			chomp $sub_seq;
			@seq = split("", $sub_seq);
				
			if($sub_seq =~ /^[ACGT](?:GGGTCA|GGGACA|GGGCCA|GGGGCA|GGGTAA|GGGTGA|GGGTTA|GGGTCC|GGGTCG|GGGTCT|AGGTCA|CGGTCA|TGGTCA|GAGTCA|GCGTCA|GTGTCA|GGATCA|GGCTCA|GGTTCA)(?:GGGTCA|GGGACA|GGGCCA|GGGGCA|GGGTAA|GGGTGA|GGGTTA|GGGTCC|GGGTCG|GGGTCT|AGGTCA|CGGTCA|TGGTCA|GAGTCA|GCGTCA|GTGTCA|GGATCA|GGCTCA|GGTTCA)$/){
				$phase = 4;
				for(my $i = 0; $i < (scalar @seq); $i++){
					$hash{"G"}{int($seq_len / $BIN_SIZE)}{int(($i + $global_pos) / $WINDOW_SIZE)}{($i - $phase)%6}{$seq[$i]}++;
				}
			}elsif($sub_seq =~ /^(?:GGGTCA|GGGACA|GGGCCA|GGGGCA|GGGTAA|GGGTGA|GGGTTA|GGGTCC|GGGTCG|GGGTCT|AGGTCA|CGGTCA|TGGTCA|GAGTCA|GCGTCA|GTGTCA|GGATCA|GGCTCA|GGTTCA)(?:GGGTC|GGGAC|GGGCC|GGGGC|GGGTA|GGGTG|GGGTT|AGGTC|CGGTC|TGGTC|GAGTC|GCGTC|GTGTC|GGATC|GGCTC|GGTTC)$/){
				$phase = 3;
				for(my $i = 0; $i < (scalar @seq); $i++){
					$hash{"G"}{int($seq_len / $BIN_SIZE)}{int(($i + $global_pos) / $WINDOW_SIZE)}{($i - $phase)%6}{$seq[$i]}++;
				}			
			}elsif($sub_seq =~ /^(?:GACCC|GTCCC|GGCCC|GCCCC|TACCC|CACCC|AACCC|GACCT|GACCG|GACCA|GACTC|GACGC|GACAC|GATCC|GAGCC|GAACC)(?:TGACCC|TGTCCC|TGGCCC|TGCCCC|TTACCC|TCACCC|TAACCC|GGACCC|CGACCC|AGACCC|TGACCT|TGACCG|TGACCA|TGACTC|TGACGC|TGACAC|TGATCC|TGAGCC|TGAACC)$/){
				$phase = 2;
				for(my $i = 0; $i < (scalar @seq); $i++){
					$hash{"C"}{int($seq_len / $BIN_SIZE)}{int(($i + $global_pos) / $WINDOW_SIZE)}{($i - $phase)%6}{$seq[$i]}++;
				}			
			}elsif($sub_seq =~ /^(?:TGACCC|TGTCCC|TGGCCC|TGCCCC|TTACCC|TCACCC|TAACCC|GGACCC|CGACCC|AGACCC|TGACCT|TGACCG|TGACCA|TGACTC|TGACGC|TGACAC|TGATCC|TGAGCC|TGAACC)(?:TGACCC|TGTCCC|TGGCCC|TGCCCC|TTACCC|TCACCC|TAACCC|GGACCC|CGACCC|AGACCC|TGACCT|TGACCG|TGACCA|TGACTC|TGACGC|TGACAC|TGATCC|TGAGCC|TGAACC)[ACGT]$/){
				$phase = 3;
				for(my $i = 0; $i < (scalar @seq); $i++){
					$hash{"C"}{int($seq_len / $BIN_SIZE)}{int(($i + $global_pos) / $WINDOW_SIZE)}{($i - $phase)%6}{$seq[$i]}++;
				}			
			}
			$global_pos += length($sub_seq);
		}
		if(length($seq) > 0){
			print STDERR "Problem! Mixed strand contained leftover sequence: $seq\n"
		}
	}
}

print OUT "strand,seq_len,seq_pos,base_pos,base_type,count\n";

foreach my $key_strand (sort {$b cmp $a} keys %hash){
	foreach my $key_sl (sort {$a <=> $b} keys $hash{$key_strand}){
		for(my $key_sp = 0; $key_sp < $MAX_RL/$WINDOW_SIZE; $key_sp++){
			for(my $key_bp = 0; $key_bp < 6; $key_bp++){
				foreach my $key_bt ("A", "C", "G", "T"){
					if(defined $hash{$key_strand}{$key_sl}{$key_sp}{$key_bp}{$key_bt}){
						print OUT $key_strand.",".$key_sl*$BIN_SIZE.",".$key_sp*$WINDOW_SIZE.",$key_bp,$key_bt,$hash{$key_strand}{$key_sl}{$key_sp}{$key_bp}{$key_bt}\n";
					}else{
						print OUT $key_strand.",".$key_sl*$BIN_SIZE.",".$key_sp*$WINDOW_SIZE.",$key_bp,$key_bt,0\n";
					}
				}
			}
		}
	}
}

close INPUT;
close OUT;
