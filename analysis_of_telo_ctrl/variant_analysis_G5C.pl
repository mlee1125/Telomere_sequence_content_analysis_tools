#!/usr/bin/perl -w
use strict;
use warnings;

my $MAX_RL = 150;
my $BIN_SIZE = 50;
my $WINDOW_SIZE = 6;
my $PERFECT_SEQ = "TTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGG";
my $PERFECT_SEQ2 = "CCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA";
my @PERFECT_SEQ = split("", $PERFECT_SEQ);

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
	
	if($seq =~ /^([ACGT]{0,5})(?:TTAGCG|ATAGCG|CTAGCG|GTAGCG|TAAGCG|TCAGCG|TGAGCG|TTCGCG|TTGGCG|TTTGCG|TTAACG|TTACCG|TTATCG|TTAGAG|TTAGGG|TTAGTG|TTAGCA|TTAGCC|TTAGCT)+([ACGT]{0,5})$/){
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
	}elsif($seq =~ /^([ACGT]{0,5})(?:CGCTAA|CGCTAT|CGCTAG|CGCTAC|CGCTTA|CGCTGA|CGCTCA|CGCGAA|CGCCAA|CGCAAA|CGTTAA|CGGTAA|CGATAA|CTCTAA|CCCTAA|CACTAA|TGCTAA|GGCTAA|AGCTAA)+([ACGT]{0,5})$/){
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
	}elsif($seq =~ /^([ACGT]{0,11}?)(?:(?:[ACGT](?:GCGTTA|GCGATA|GCGCTA|GCGGTA|GCGTAA|GCGTCA|GCGTGA|GCGTTC|GCGTTG|GCGTTT|ACGTTA|CCGTTA|TCGTTA|GAGTTA|GGGTTA|GTGTTA|GCATTA|GCCTTA|GCTTTA)(?:GCGTTA|GCGATA|GCGCTA|GCGGTA|GCGTAA|GCGTCA|GCGTGA|GCGTTC|GCGTTG|GCGTTT|ACGTTA|CCGTTA|TCGTTA|GAGTTA|GGGTTA|GTGTTA|GCATTA|GCCTTA|GCTTTA))|(?:(?:GCGTTA|GCGATA|GCGCTA|GCGGTA|GCGTAA|GCGTCA|GCGTGA|GCGTTC|GCGTTG|GCGTTT|ACGTTA|CCGTTA|TCGTTA|GAGTTA|GGGTTA|GTGTTA|GCATTA|GCCTTA|GCTTTA)(?:GCGTT|GCGAT|GCGCT|GCGGT|GCGTA|GCGTC|GCGTG|ACGTT|CCGTT|TCGTT|GAGTT|GGGTT|GTGTT|GCATT|GCCTT|GCTTT))|(?:(?:AACGC|TACGC|GACGC|CACGC|AACGC|AACGC|AACGC|AACGT|AACGG|AACGA|AACTC|AACCC|AACAC|AATGC|AAGGC|AAAGC)(?:TAACGC|TATCGC|TAGCGC|TACCGC|TTACGC|TGACGC|TCACGC|GAACGC|CAACGC|AAACGC|TAACGT|TAACGG|TAACGA|TAACTC|TAACCC|TAACAC|TAATGC|TAAGGC|TAAAGC))|(?:(?:TAACGC|TATCGC|TAGCGC|TACCGC|TTACGC|TGACGC|TCACGC|GAACGC|CAACGC|AAACGC|TAACGT|TAACGG|TAACGA|TAACTC|TAACCC|TAACAC|TAATGC|TAAGGC|TAAAGC)(?:TAACGC|TATCGC|TAGCGC|TACCGC|TTACGC|TGACGC|TCACGC|GAACGC|CAACGC|AAACGC|TAACGT|TAACGG|TAACGA|TAACTC|TAACCC|TAACAC|TAATGC|TAAGGC|TAAAGC)[ACGT]))+([ACGT]{0,11}?)$/){
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
	
	
	#print "$phase\n";
		
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
		while($seq =~ /^((?:[ACGT](?:GCGTTA|GCGATA|GCGCTA|GCGGTA|GCGTAA|GCGTCA|GCGTGA|GCGTTC|GCGTTG|GCGTTT|ACGTTA|CCGTTA|TCGTTA|GAGTTA|GGGTTA|GTGTTA|GCATTA|GCCTTA|GCTTTA)(?:GCGTTA|GCGATA|GCGCTA|GCGGTA|GCGTAA|GCGTCA|GCGTGA|GCGTTC|GCGTTG|GCGTTT|ACGTTA|CCGTTA|TCGTTA|GAGTTA|GGGTTA|GTGTTA|GCATTA|GCCTTA|GCTTTA))|(?:(?:GCGTTA|GCGATA|GCGCTA|GCGGTA|GCGTAA|GCGTCA|GCGTGA|GCGTTC|GCGTTG|GCGTTT|ACGTTA|CCGTTA|TCGTTA|GAGTTA|GGGTTA|GTGTTA|GCATTA|GCCTTA|GCTTTA)(?:GCGTT|GCGAT|GCGCT|GCGGT|GCGTA|GCGTC|GCGTG|ACGTT|CCGTT|TCGTT|GAGTT|GGGTT|GTGTT|GCATT|GCCTT|GCTTT))|(?:(?:AACGC|TACGC|GACGC|CACGC|AACGC|AACGC|AACGC|AACGT|AACGG|AACGA|AACTC|AACCC|AACAC|AATGC|AAGGC|AAAGC)(?:TAACGC|TATCGC|TAGCGC|TACCGC|TTACGC|TGACGC|TCACGC|GAACGC|CAACGC|AAACGC|TAACGT|TAACGG|TAACGA|TAACTC|TAACCC|TAACAC|TAATGC|TAAGGC|TAAAGC))|(?:(?:TAACGC|TATCGC|TAGCGC|TACCGC|TTACGC|TGACGC|TCACGC|GAACGC|CAACGC|AAACGC|TAACGT|TAACGG|TAACGA|TAACTC|TAACCC|TAACAC|TAATGC|TAAGGC|TAAAGC)(?:TAACGC|TATCGC|TAGCGC|TACCGC|TTACGC|TGACGC|TCACGC|GAACGC|CAACGC|AAACGC|TAACGT|TAACGG|TAACGA|TAACTC|TAACCC|TAACAC|TAATGC|TAAGGC|TAAAGC)[ACGT]))((?:(?:[ACGT](?:GCGTTA|GCGATA|GCGCTA|GCGGTA|GCGTAA|GCGTCA|GCGTGA|GCGTTC|GCGTTG|GCGTTT|ACGTTA|CCGTTA|TCGTTA|GAGTTA|GGGTTA|GTGTTA|GCATTA|GCCTTA|GCTTTA)(?:GCGTTA|GCGATA|GCGCTA|GCGGTA|GCGTAA|GCGTCA|GCGTGA|GCGTTC|GCGTTG|GCGTTT|ACGTTA|CCGTTA|TCGTTA|GAGTTA|GGGTTA|GTGTTA|GCATTA|GCCTTA|GCTTTA))|(?:(?:GCGTTA|GCGATA|GCGCTA|GCGGTA|GCGTAA|GCGTCA|GCGTGA|GCGTTC|GCGTTG|GCGTTT|ACGTTA|CCGTTA|TCGTTA|GAGTTA|GGGTTA|GTGTTA|GCATTA|GCCTTA|GCTTTA)(?:GCGTT|GCGAT|GCGCT|GCGGT|GCGTA|GCGTC|GCGTG|ACGTT|CCGTT|TCGTT|GAGTT|GGGTT|GTGTT|GCATT|GCCTT|GCTTT))|(?:(?:AACGC|TACGC|GACGC|CACGC|AACGC|AACGC|AACGC|AACGT|AACGG|AACGA|AACTC|AACCC|AACAC|AATGC|AAGGC|AAAGC)(?:TAACGC|TATCGC|TAGCGC|TACCGC|TTACGC|TGACGC|TCACGC|GAACGC|CAACGC|AAACGC|TAACGT|TAACGG|TAACGA|TAACTC|TAACCC|TAACAC|TAATGC|TAAGGC|TAAAGC))|(?:(?:TAACGC|TATCGC|TAGCGC|TACCGC|TTACGC|TGACGC|TCACGC|GAACGC|CAACGC|AAACGC|TAACGT|TAACGG|TAACGA|TAACTC|TAACCC|TAACAC|TAATGC|TAAGGC|TAAAGC)(?:TAACGC|TATCGC|TAGCGC|TACCGC|TTACGC|TGACGC|TCACGC|GAACGC|CAACGC|AAACGC|TAACGT|TAACGG|TAACGA|TAACTC|TAACCC|TAACAC|TAATGC|TAAGGC|TAAAGC)[ACGT]))*)$/){
			$sub_seq = $1;
			if(defined $2){
				$seq = $2;
			}else{
				$seq = "";
			}

			chomp $sub_seq;
			@seq = split("", $sub_seq);
				
			if($sub_seq =~ /^[ACGT](?:GCGTTA|GCGATA|GCGCTA|GCGGTA|GCGTAA|GCGTCA|GCGTGA|GCGTTC|GCGTTG|GCGTTT|ACGTTA|CCGTTA|TCGTTA|GAGTTA|GGGTTA|GTGTTA|GCATTA|GCCTTA|GCTTTA)(?:GCGTTA|GCGATA|GCGCTA|GCGGTA|GCGTAA|GCGTCA|GCGTGA|GCGTTC|GCGTTG|GCGTTT|ACGTTA|CCGTTA|TCGTTA|GAGTTA|GGGTTA|GTGTTA|GCATTA|GCCTTA|GCTTTA)$/){
				$phase = 4;
				for(my $i = 0; $i < (scalar @seq); $i++){
					$hash{"G"}{int($seq_len / $BIN_SIZE)}{int(($i + $global_pos) / $WINDOW_SIZE)}{($i - $phase)%6}{$seq[$i]}++;
				}
			}elsif($sub_seq =~ /^(?:GCGTTA|GCGATA|GCGCTA|GCGGTA|GCGTAA|GCGTCA|GCGTGA|GCGTTC|GCGTTG|GCGTTT|ACGTTA|CCGTTA|TCGTTA|GAGTTA|GGGTTA|GTGTTA|GCATTA|GCCTTA|GCTTTA)(?:GCGTT|GCGAT|GCGCT|GCGGT|GCGTA|GCGTC|GCGTG|ACGTT|CCGTT|TCGTT|GAGTT|GGGTT|GTGTT|GCATT|GCCTT|GCTTT)$/){
				$phase = 3;
				for(my $i = 0; $i < (scalar @seq); $i++){
					$hash{"G"}{int($seq_len / $BIN_SIZE)}{int(($i + $global_pos) / $WINDOW_SIZE)}{($i - $phase)%6}{$seq[$i]}++;
				}			
			}elsif($sub_seq =~ /^(?:AACGC|TACGC|GACGC|CACGC|AACGC|AACGC|AACGC|AACGT|AACGG|AACGA|AACTC|AACCC|AACAC|AATGC|AAGGC|AAAGC)(?:TAACGC|TATCGC|TAGCGC|TACCGC|TTACGC|TGACGC|TCACGC|GAACGC|CAACGC|AAACGC|TAACGT|TAACGG|TAACGA|TAACTC|TAACCC|TAACAC|TAATGC|TAAGGC|TAAAGC)$/){
				$phase = 2;
				for(my $i = 0; $i < (scalar @seq); $i++){
					$hash{"C"}{int($seq_len / $BIN_SIZE)}{int(($i + $global_pos) / $WINDOW_SIZE)}{($i - $phase)%6}{$seq[$i]}++;
				}			
			}elsif($sub_seq =~ /^(?:TAACGC|TATCGC|TAGCGC|TACCGC|TTACGC|TGACGC|TCACGC|GAACGC|CAACGC|AAACGC|TAACGT|TAACGG|TAACGA|TAACTC|TAACCC|TAACAC|TAATGC|TAAGGC|TAAAGC)(?:TAACGC|TATCGC|TAGCGC|TACCGC|TTACGC|TGACGC|TCACGC|GAACGC|CAACGC|AAACGC|TAACGT|TAACGG|TAACGA|TAACTC|TAACCC|TAACAC|TAATGC|TAAGGC|TAAAGC)[ACGT]$/){
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
