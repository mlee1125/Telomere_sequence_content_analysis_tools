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
	
	if($seq =~ /^([ACGT]{0,5})(?:TTAGGG|ATAGGG|CTAGGG|GTAGGG|TAAGGG|TCAGGG|TGAGGG|TTCGGG|TTGGGG|TTTGGG|TTAAGG|TTACGG|TTATGG|TTAGAG|TTAGCG|TTAGTG|TTAGGA|TTAGGC|TTAGGT)+([ACGT]{0,5})$/){
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
	}elsif($seq =~ /^([ACGT]{0,5})(?:CCCTAA|CCCTAT|CCCTAG|CCCTAC|CCCTTA|CCCTGA|CCCTCA|CCCGAA|CCCCAA|CCCAAA|CCTTAA|CCGTAA|CCATAA|CTCTAA|CGCTAA|CACTAA|TCCTAA|GCCTAA|ACCTAA)+([ACGT]{0,5})$/){
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
	}elsif($seq =~ /^([ACGT]{0,11}?)(?:(?:[ACGT](?:GGGTTA|GGGATA|GGGCTA|GGGGTA|GGGTAA|GGGTCA|GGGTGA|GGGTTC|GGGTTG|GGGTTT|AGGTTA|CGGTTA|TGGTTA|GAGTTA|GCGTTA|GTGTTA|GGATTA|GGCTTA|GGTTTA)(?:GGGTTA|GGGATA|GGGCTA|GGGGTA|GGGTAA|GGGTCA|GGGTGA|GGGTTC|GGGTTG|GGGTTT|AGGTTA|CGGTTA|TGGTTA|GAGTTA|GCGTTA|GTGTTA|GGATTA|GGCTTA|GGTTTA))|(?:(?:GGGTTA|GGGATA|GGGCTA|GGGGTA|GGGTAA|GGGTCA|GGGTGA|GGGTTC|GGGTTG|GGGTTT|AGGTTA|CGGTTA|TGGTTA|GAGTTA|GCGTTA|GTGTTA|GGATTA|GGCTTA|GGTTTA)(?:GGGTT|GGGAT|GGGCT|GGGGT|GGGTA|GGGTC|GGGTG|AGGTT|CGGTT|TGGTT|GAGTT|GCGTT|GTGTT|GGATT|GGCTT|GGTTT))|(?:(?:AACCC|TACCC|GACCC|CACCC|AACCC|AACCC|AACCC|AACCT|AACCG|AACCA|AACTC|AACGC|AACAC|AATCC|AAGCC|AAACC)(?:TAACCC|TATCCC|TAGCCC|TACCCC|TTACCC|TGACCC|TCACCC|GAACCC|CAACCC|AAACCC|TAACCT|TAACCG|TAACCA|TAACTC|TAACGC|TAACAC|TAATCC|TAAGCC|TAAACC))|(?:(?:TAACCC|TATCCC|TAGCCC|TACCCC|TTACCC|TGACCC|TCACCC|GAACCC|CAACCC|AAACCC|TAACCT|TAACCG|TAACCA|TAACTC|TAACGC|TAACAC|TAATCC|TAAGCC|TAAACC)(?:TAACCC|TATCCC|TAGCCC|TACCCC|TTACCC|TGACCC|TCACCC|GAACCC|CAACCC|AAACCC|TAACCT|TAACCG|TAACCA|TAACTC|TAACGC|TAACAC|TAATCC|TAAGCC|TAAACC)[ACGT]))+([ACGT]{0,11}?)$/){
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
		print STDERR "Problem!\n$seq\n";
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
		while($seq =~ /^((?:[ACGT](?:GGGTTA|GGGATA|GGGCTA|GGGGTA|GGGTAA|GGGTCA|GGGTGA|GGGTTC|GGGTTG|GGGTTT|AGGTTA|CGGTTA|TGGTTA|GAGTTA|GCGTTA|GTGTTA|GGATTA|GGCTTA|GGTTTA)(?:GGGTTA|GGGATA|GGGCTA|GGGGTA|GGGTAA|GGGTCA|GGGTGA|GGGTTC|GGGTTG|GGGTTT|AGGTTA|CGGTTA|TGGTTA|GAGTTA|GCGTTA|GTGTTA|GGATTA|GGCTTA|GGTTTA))|(?:(?:GGGTTA|GGGATA|GGGCTA|GGGGTA|GGGTAA|GGGTCA|GGGTGA|GGGTTC|GGGTTG|GGGTTT|AGGTTA|CGGTTA|TGGTTA|GAGTTA|GCGTTA|GTGTTA|GGATTA|GGCTTA|GGTTTA)(?:GGGTT|GGGAT|GGGCT|GGGGT|GGGTA|GGGTC|GGGTG|AGGTT|CGGTT|TGGTT|GAGTT|GCGTT|GTGTT|GGATT|GGCTT|GGTTT))|(?:(?:AACCC|TACCC|GACCC|CACCC|AACCC|AACCC|AACCC|AACCT|AACCG|AACCA|AACTC|AACGC|AACAC|AATCC|AAGCC|AAACC)(?:TAACCC|TATCCC|TAGCCC|TACCCC|TTACCC|TGACCC|TCACCC|GAACCC|CAACCC|AAACCC|TAACCT|TAACCG|TAACCA|TAACTC|TAACGC|TAACAC|TAATCC|TAAGCC|TAAACC))|(?:(?:TAACCC|TATCCC|TAGCCC|TACCCC|TTACCC|TGACCC|TCACCC|GAACCC|CAACCC|AAACCC|TAACCT|TAACCG|TAACCA|TAACTC|TAACGC|TAACAC|TAATCC|TAAGCC|TAAACC)(?:TAACCC|TATCCC|TAGCCC|TACCCC|TTACCC|TGACCC|TCACCC|GAACCC|CAACCC|AAACCC|TAACCT|TAACCG|TAACCA|TAACTC|TAACGC|TAACAC|TAATCC|TAAGCC|TAAACC)[ACGT]))((?:(?:[ACGT](?:GGGTTA|GGGATA|GGGCTA|GGGGTA|GGGTAA|GGGTCA|GGGTGA|GGGTTC|GGGTTG|GGGTTT|AGGTTA|CGGTTA|TGGTTA|GAGTTA|GCGTTA|GTGTTA|GGATTA|GGCTTA|GGTTTA)(?:GGGTTA|GGGATA|GGGCTA|GGGGTA|GGGTAA|GGGTCA|GGGTGA|GGGTTC|GGGTTG|GGGTTT|AGGTTA|CGGTTA|TGGTTA|GAGTTA|GCGTTA|GTGTTA|GGATTA|GGCTTA|GGTTTA))|(?:(?:GGGTTA|GGGATA|GGGCTA|GGGGTA|GGGTAA|GGGTCA|GGGTGA|GGGTTC|GGGTTG|GGGTTT|AGGTTA|CGGTTA|TGGTTA|GAGTTA|GCGTTA|GTGTTA|GGATTA|GGCTTA|GGTTTA)(?:GGGTT|GGGAT|GGGCT|GGGGT|GGGTA|GGGTC|GGGTG|AGGTT|CGGTT|TGGTT|GAGTT|GCGTT|GTGTT|GGATT|GGCTT|GGTTT))|(?:(?:AACCC|TACCC|GACCC|CACCC|AACCC|AACCC|AACCC|AACCT|AACCG|AACCA|AACTC|AACGC|AACAC|AATCC|AAGCC|AAACC)(?:TAACCC|TATCCC|TAGCCC|TACCCC|TTACCC|TGACCC|TCACCC|GAACCC|CAACCC|AAACCC|TAACCT|TAACCG|TAACCA|TAACTC|TAACGC|TAACAC|TAATCC|TAAGCC|TAAACC))|(?:(?:TAACCC|TATCCC|TAGCCC|TACCCC|TTACCC|TGACCC|TCACCC|GAACCC|CAACCC|AAACCC|TAACCT|TAACCG|TAACCA|TAACTC|TAACGC|TAACAC|TAATCC|TAAGCC|TAAACC)(?:TAACCC|TATCCC|TAGCCC|TACCCC|TTACCC|TGACCC|TCACCC|GAACCC|CAACCC|AAACCC|TAACCT|TAACCG|TAACCA|TAACTC|TAACGC|TAACAC|TAATCC|TAAGCC|TAAACC)[ACGT]))*)$/){
			$sub_seq = $1;
			if(defined $2){
				$seq = $2;
			}else{
				$seq = "";
			}

			chomp $sub_seq;
			@seq = split("", $sub_seq);
				
			if($sub_seq =~ /^[ACGT](?:GGGTTA|GGGATA|GGGCTA|GGGGTA|GGGTAA|GGGTCA|GGGTGA|GGGTTC|GGGTTG|GGGTTT|AGGTTA|CGGTTA|TGGTTA|GAGTTA|GCGTTA|GTGTTA|GGATTA|GGCTTA|GGTTTA)(?:GGGTTA|GGGATA|GGGCTA|GGGGTA|GGGTAA|GGGTCA|GGGTGA|GGGTTC|GGGTTG|GGGTTT|AGGTTA|CGGTTA|TGGTTA|GAGTTA|GCGTTA|GTGTTA|GGATTA|GGCTTA|GGTTTA)$/){
				$phase = 4;
				for(my $i = 0; $i < (scalar @seq); $i++){
					$hash{"G"}{int($seq_len / $BIN_SIZE)}{int(($i + $global_pos) / $WINDOW_SIZE)}{($i - $phase)%6}{$seq[$i]}++;
				}
			}elsif($sub_seq =~ /^(?:GGGTTA|GGGATA|GGGCTA|GGGGTA|GGGTAA|GGGTCA|GGGTGA|GGGTTC|GGGTTG|GGGTTT|AGGTTA|CGGTTA|TGGTTA|GAGTTA|GCGTTA|GTGTTA|GGATTA|GGCTTA|GGTTTA)(?:GGGTT|GGGAT|GGGCT|GGGGT|GGGTA|GGGTC|GGGTG|AGGTT|CGGTT|TGGTT|GAGTT|GCGTT|GTGTT|GGATT|GGCTT|GGTTT)$/){
				$phase = 3;
				for(my $i = 0; $i < (scalar @seq); $i++){
					$hash{"G"}{int($seq_len / $BIN_SIZE)}{int(($i + $global_pos) / $WINDOW_SIZE)}{($i - $phase)%6}{$seq[$i]}++;
				}			
			}elsif($sub_seq =~ /^(?:AACCC|TACCC|GACCC|CACCC|AACCC|AACCC|AACCC|AACCT|AACCG|AACCA|AACTC|AACGC|AACAC|AATCC|AAGCC|AAACC)(?:TAACCC|TATCCC|TAGCCC|TACCCC|TTACCC|TGACCC|TCACCC|GAACCC|CAACCC|AAACCC|TAACCT|TAACCG|TAACCA|TAACTC|TAACGC|TAACAC|TAATCC|TAAGCC|TAAACC)$/){
				$phase = 2;
				for(my $i = 0; $i < (scalar @seq); $i++){
					$hash{"C"}{int($seq_len / $BIN_SIZE)}{int(($i + $global_pos) / $WINDOW_SIZE)}{($i - $phase)%6}{$seq[$i]}++;
				}			
			}elsif($sub_seq =~ /^(?:TAACCC|TATCCC|TAGCCC|TACCCC|TTACCC|TGACCC|TCACCC|GAACCC|CAACCC|AAACCC|TAACCT|TAACCG|TAACCA|TAACTC|TAACGC|TAACAC|TAATCC|TAAGCC|TAAACC)(?:TAACCC|TATCCC|TAGCCC|TACCCC|TTACCC|TGACCC|TCACCC|GAACCC|CAACCC|AAACCC|TAACCT|TAACCG|TAACCA|TAACTC|TAACGC|TAACAC|TAATCC|TAAGCC|TAAACC)[ACGT]$/){
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
