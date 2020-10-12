#!/usr/bin/perl -w
use strict;
use warnings;

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

my $out_file_pure_G = "$file_prefix.G-strands.fastq";
my $out_file_pure_C = "$file_prefix.C-strands.fastq";
my $out_file_mixed = "$file_prefix.mixed-strands.fastq";
my $out_file_other = "$file_prefix.others.fastq";

open(INPUT, "<", $input_file) or die "Cannot open < $input_file: $!\n";

open(OUT_G, ">", $out_file_pure_G) or die "Cannot open > $out_file_pure_G: $!\n";
open(OUT_C, ">", $out_file_pure_C) or die "Cannot open > $out_file_pure_C: $!\n";
open(OUT_mix, ">", $out_file_mixed) or die "Cannot open > $out_file_mixed: $!\n";
open(OUT_other, ">", $out_file_other) or die "Cannot open > $out_file_other: $!\n";

my $id = "";
my $seq = "";
my $qual = "";

while(my $line = <INPUT>){
	unless($line =~ /^@.*/){
		print STDERR "Error: File not in fastq format, header incorrect!\n";
		exit;
	}
	$id = $line;
	
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
	
	if($seq =~ /^[ACGT]{0,5}(?:TTAGGG|ATAGGG|CTAGGG|GTAGGG|TAAGGG|TCAGGG|TGAGGG|TTCGGG|TTGGGG|TTTGGG|TTAAGG|TTACGG|TTATGG|TTAGAG|TTAGCG|TTAGTG|TTAGGA|TTAGGC|TTAGGT)+[ACGT]{0,5}$/){
		print OUT_G $id;
		print OUT_G $seq;
		print OUT_G "+\n";
		print OUT_G $qual;
	}elsif($seq =~ /^[ACGT]{0,5}(?:CCCTAA|CCCTAT|CCCTAG|CCCTAC|CCCTTA|CCCTGA|CCCTCA|CCCGAA|CCCCAA|CCCAAA|CCTTAA|CCGTAA|CCATAA|CTCTAA|CGCTAA|CACTAA|TCCTAA|GCCTAA|ACCTAA)+[ACGT]{0,5}$/){
		print OUT_C $id;
		print OUT_C $seq;
		print OUT_C "+\n";
		print OUT_C $qual;
	}elsif($seq =~ /^[ACGT]{0,11}(?:(?:[ACGT](?:GGGTTA|GGGATA|GGGCTA|GGGGTA|GGGTAA|GGGTCA|GGGTGA|GGGTTC|GGGTTG|GGGTTT|AGGTTA|CGGTTA|TGGTTA|GAGTTA|GCGTTA|GTGTTA|GGATTA|GGCTTA|GGTTTA)(?:GGGTTA|GGGATA|GGGCTA|GGGGTA|GGGTAA|GGGTCA|GGGTGA|GGGTTC|GGGTTG|GGGTTT|AGGTTA|CGGTTA|TGGTTA|GAGTTA|GCGTTA|GTGTTA|GGATTA|GGCTTA|GGTTTA))|(?:(?:GGGTTA|GGGATA|GGGCTA|GGGGTA|GGGTAA|GGGTCA|GGGTGA|GGGTTC|GGGTTG|GGGTTT|AGGTTA|CGGTTA|TGGTTA|GAGTTA|GCGTTA|GTGTTA|GGATTA|GGCTTA|GGTTTA)(?:GGGTT|GGGAT|GGGCT|GGGGT|GGGTA|GGGTC|GGGTG|AGGTT|CGGTT|TGGTT|GAGTT|GCGTT|GTGTT|GGATT|GGCTT|GGTTT))|(?:(?:AACCC|TACCC|GACCC|CACCC|AACCC|AACCC|AACCC|AACCT|AACCG|AACCA|AACTC|AACGC|AACAC|AATCC|AAGCC|AAACC)(?:TAACCC|TATCCC|TAGCCC|TACCCC|TTACCC|TGACCC|TCACCC|GAACCC|CAACCC|AAACCC|TAACCT|TAACCG|TAACCA|TAACTC|TAACGC|TAACAC|TAATCC|TAAGCC|TAAACC))|(?:(?:TAACCC|TATCCC|TAGCCC|TACCCC|TTACCC|TGACCC|TCACCC|GAACCC|CAACCC|AAACCC|TAACCT|TAACCG|TAACCA|TAACTC|TAACGC|TAACAC|TAATCC|TAAGCC|TAAACC)(?:TAACCC|TATCCC|TAGCCC|TACCCC|TTACCC|TGACCC|TCACCC|GAACCC|CAACCC|AAACCC|TAACCT|TAACCG|TAACCA|TAACTC|TAACGC|TAACAC|TAATCC|TAAGCC|TAAACC)[ACGT]))+[ACGT]{0,11}$/){
		print OUT_mix $id;
		print OUT_mix $seq;
		print OUT_mix "+\n";
		print OUT_mix $qual;
	}else{
		print OUT_other $id;
		print OUT_other $seq;
		print OUT_other "+\n";
		print OUT_other $qual;
	}
}

#(?:TTAGGG|ATAGGG|CTAGGG|GTAGGG|TAAGGG|TCAGGG|TGAGGG|TTCGGG|TTGGGG|TTTGGG|TTAAGG|TTACGG|TTATGG|TTAGAG|TTAGCG|TTAGTG|TTAGGA|TTAGGC|TTAGGT)
#(?:CCCTAA|CCCTAT|CCCTAG|CCCTAC|CCCTTA|CCCTGA|CCCTCA|CCCGAA|CCCCAA|CCCAAA|CCTTAA|CCGTAA|CCATAA|CTCTAA|CGCTAA|CACTAA|TCCTAA|GCCTAA|ACCTAA)

#(A(GGGTTA)(GGGTTA))|((GGGTTA)(GGGTT))|((AACCC)(TAACCC))|((TAACCC)(TAACCC)T)
#(?:GGGTTA|GGGATA|GGGCTA|GGGGTA|GGGTAA|GGGTCA|GGGTGA|GGGTTC|GGGTTG|GGGTTT|AGGTTA|CGGTTA|TGGTTA|GAGTTA|GCGTTA|GTGTTA|GGATTA|GGCTTA|GGTTTA)
#(?:GGGTT|GGGAT|GGGCT|GGGGT|GGGTA|GGGTC|GGGTG|AGGTT|CGGTT|TGGTT|GAGTT|GCGTT|GTGTT|GGATT|GGCTT|GGTTT)
#(?:TAACCC|TATCCC|TAGCCC|TACCCC|TTACCC|TGACCC|TCACCC|GAACCC|CAACCC|AAACCC|TAACCT|TAACCG|TAACCA|TAACTC|TAACGC|TAACAC|TAATCC|TAAGCC|TAAACC)
#(?:AACCC|ATCCC|AGCCC|ACCCC|TACCC|GACCC|CACCC|AACCT|AACCG|AACCA|AACTC|AACGC|AACAC|AATCC|AAGCC|AAACC)
#(?:(?:[ACGT](?:GGGTTA|GGGATA|GGGCTA|GGGGTA|GGGTAA|GGGTCA|GGGTGA|GGGTTC|GGGTTG|GGGTTT|AGGTTA|CGGTTA|TGGTTA|GAGTTA|GCGTTA|GTGTTA|GGATTA|GGCTTA|GGTTTA)(?:GGGTTA|GGGATA|GGGCTA|GGGGTA|GGGTAA|GGGTCA|GGGTGA|GGGTTC|GGGTTG|GGGTTT|AGGTTA|CGGTTA|TGGTTA|GAGTTA|GCGTTA|GTGTTA|GGATTA|GGCTTA|GGTTTA))|(?:(?:GGGTTA|GGGATA|GGGCTA|GGGGTA|GGGTAA|GGGTCA|GGGTGA|GGGTTC|GGGTTG|GGGTTT|AGGTTA|CGGTTA|TGGTTA|GAGTTA|GCGTTA|GTGTTA|GGATTA|GGCTTA|GGTTTA)(?:GGGTT|GGGAT|GGGCT|GGGGT|GGGTA|GGGTC|GGGTG|AGGTT|CGGTT|TGGTT|GAGTT|GCGTT|GTGTT|GGATT|GGCTT|GGTTT))|(?:(?:AACCC|ATCCC|AGCCC|ACCCC|TACCC|GACCC|CACCC|AACCT|AACCG|AACCA|AACTC|AACGC|AACAC|AATCC|AAGCC|AAACC)(?:TAACCC|TATCCC|TAGCCC|TACCCC|TTACCC|TGACCC|TCACCC|GAACCC|CAACCC|AAACCC|TAACCT|TAACCG|TAACCA|TAACTC|TAACGC|TAACAC|TAATCC|TAAGCC|TAAACC))|(?:(?:TAACCC|TATCCC|TAGCCC|TACCCC|TTACCC|TGACCC|TCACCC|GAACCC|CAACCC|AAACCC|TAACCT|TAACCG|TAACCA|TAACTC|TAACGC|TAACAC|TAATCC|TAAGCC|TAAACC)(?:TAACCC|TATCCC|TAGCCC|TACCCC|TTACCC|TGACCC|TCACCC|GAACCC|CAACCC|AAACCC|TAACCT|TAACCG|TAACCA|TAACTC|TAACGC|TAACAC|TAATCC|TAAGCC|TAAACC)[ACGT]))


close INPUT;

close OUT_G;
close OUT_C;
close OUT_mix;
close OUT_other;