#!/usr/bin/perl -w

# Copyright (C) 2015 Adam Ameur

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


#########################################################################################
##
## File: cava.pl
##
## CAVA: Count-based analysis of variants in long amplicons
##
## Author: Adam Ameur
##
## Input: Two parameters are required
##
## 1) FASTA file containg full-length amplicon sequence data,
##    produced for example by PacBio sequencing
##
## 2) A pattern string containing the DNA sequence, with mutation indicated by '[]'
##    Example: TGGAACGCACGGACATCACC[A/G]TGAAGCACAAGCTGGGCGGG
##
## The number of allowed mismatches (sequence errors) in the pattern string is passed as
## an optional parameter.
##
## Output: All matching hits of wild-type and mutation sequences and their orientation
## within the reads are reported in the output. Optionally, detailed information about
## for each individual read is printed.
##
##########################################################################################

################
##
## Main program
##
################

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;

# Global hash variables. These are used for speeding up the matching code.
my %seen_substrings;  # Hash with number of mismatches in substrings of same length as pattern
my %unmatched_reads;  # Hash with all reads that do not contain any match to pattern

# Read command line arguments
my ($in_file, $mutation, $max_mismatch, $output_read_ids, $hash_off, $silent) = get_args_and_error_check();

# Execute mutation matching
cava();


################
##
## Functions
##
################

# Function that performs the mutation screening
sub cava{

	unless($silent){
		print STDOUT "\n######### Screening of $mutation in file: $in_file  #########\n\n";
	}

	my %reads=read_fasta_file($in_file); # Read sequence data from FASTA file
	my %updated_reads=%reads; # Hash to keep track of reads that remain to be screened

	my $seq_up = "undefined";
	my $seq_down = "undefined";
	my $mut = "undefined";
	my $wt = "undefined";
	my $seq_up_rc = "undefined";
	my $seq_down_rc = "undefined";
	my $mut_rc = "undefined";
	my $wt_rc = "undefined";
	my $len_up=-1;
	my $len_down=-1;
	my $ref_count_fwd=0;
	my $ref_count_rev=0;
	my $alt_count_fwd=0;
	my $alt_count_rev=0;
	my $other_count_fwd=0;
	my $other_count_rev=0;

	if($mutation =~ /^\s*([ACGT]+)\[([ACGT]*)\/([ACGT]*)\]([ACGT]+)\s*$/){
		$seq_up = $1;
		$wt = $2;
		$mut = $3;
		$seq_down = $4;
		$seq_up_rc = reverse_complement_IUPAC($seq_up);
		$seq_down_rc = reverse_complement_IUPAC($seq_down);
		$mut_rc = reverse_complement_IUPAC($mut);
		$wt_rc = reverse_complement_IUPAC($wt);

		$len_up=length $seq_up;
		$len_down=length $seq_down;
	}


	my $pat = "$seq_up$wt$seq_down";

	## First count all identical matches to wt/mut...
	foreach my $key (keys %reads) {
		my $mutation_result="undefined";
		my $match_sequence="undefined";
		my $nr_mismatches=-1;

		my $sequence = $reads{$key};

		if($sequence =~/($seq_up[^$mut]$seq_down)/){
			$match_sequence=$1;
			if($sequence =~/$seq_up$wt$seq_down/){
				$ref_count_fwd++;
				$nr_mismatches = 0;
				$mutation_result = "ref_fwd";
			}
			else{
				$other_count_fwd++;
				$nr_mismatches = 1;
				$mutation_result = "other_fwd";
			}
		}
		if($sequence =~/($seq_down_rc[^$mut_rc]$seq_up_rc)/){
			$match_sequence=$1;
			if($sequence =~/$seq_down_rc$wt_rc$seq_up_rc/){
				$ref_count_rev++;
				$nr_mismatches = 0;
				$mutation_result = "ref_rev";
			}
			else{
				$other_count_rev++;
				$nr_mismatches = 1;
				$mutation_result = "other_rev";
			}
		}
		if($sequence =~/($seq_up$mut$seq_down)/){
			$match_sequence=$1;
			$alt_count_fwd++;
			$mutation_result = "alt_fwd";
			$nr_mismatches = 0;
		}
		if($sequence =~/($seq_down_rc$mut_rc$seq_up_rc)/){
			$match_sequence=$1;
			$alt_count_rev++;
			$mutation_result = "alt_rev";
			$nr_mismatches = 0;
		}

		## Print results and remove from hash if there is a match
		if(!($mutation_result eq "undefined")){
			if($output_read_ids){
				print "$key\t$match_sequence\t$mutation_result\t$nr_mismatches\n";
			}
			delete $updated_reads{$key};
		}
	}


	## Go through all remaing reads and count occurences of ref/alt by fuzzy matching...
	if($max_mismatch > 0){
		while (my ($key) = each(%updated_reads)) {

			my $mutation_result="undefined";
			my $match_sequence="undefined";
			my $orientation="undefined";
			my $nr_mismatches=-1;
			my $sequence = $updated_reads{$key};

			## Check if read is in a global hash of 'unmatching' reads. Execute the commands
			## only if this is not true.
			if( !(defined $unmatched_reads{$sequence}) ){
                ## make a fuzzy matching on fwd strand
				my %best_alignment_fwd = get_best_alignment($sequence, $pat, $max_mismatch, $hash_off);

				my @mismatch_vec_fwd = (keys %best_alignment_fwd);
				my $min_mismatch = $mismatch_vec_fwd[0];

				if($min_mismatch <= $max_mismatch){
					$match_sequence = $best_alignment_fwd{$min_mismatch};
					$orientation = "fwd";
					$nr_mismatches = $min_mismatch;
				}

				## make a fuzzy matching on rev strand
				if($min_mismatch > $max_mismatch){
					my $sequence_rev = reverse_complement_IUPAC($sequence);
					my %best_alignment_rev = get_best_alignment($sequence_rev, $pat, $max_mismatch, $hash_off);
					my @mismatch_vec_rev = (keys %best_alignment_rev);
					my $min_mismatch_rev = $mismatch_vec_rev[0];

					if($min_mismatch_rev <= $max_mismatch){
						$match_sequence = $best_alignment_rev{$min_mismatch_rev};
						$orientation = "rev";
						$min_mismatch = $min_mismatch_rev;
						$nr_mismatches = $min_mismatch;
					}
				}

				## If a match is found either on fwd or rev, go through hash and report all other reads
				## having the exact same matching sequence
				if(($nr_mismatches >= 0) && ($nr_mismatches <= $max_mismatch)){

					my $tmp_trimmed = substr $match_sequence, $len_up; # trim bases upstream
					my $trimmed_len = length($tmp_trimmed);
					my $variable_seq = substr $tmp_trimmed, 0, $trimmed_len-$len_down; # trim bases downstream

					if($orientation eq "fwd"){
					if($variable_seq eq $wt){
						$ref_count_fwd++;
						$mutation_result = "ref_fwd";
					}
					if($variable_seq eq $mut){
						$alt_count_fwd++;
						$mutation_result = "alt_fwd";
					}
					if($mutation_result eq "undefined"){
						$ref_count_fwd++;
						$mutation_result = "other_fwd";
					}
					}
					if($orientation eq "rev"){
						if($variable_seq eq $wt){
							$ref_count_rev++;
							$mutation_result = "ref_rev";
						}
						if($variable_seq eq $mut){
							$alt_count_rev++;
							$mutation_result = "alt_rev";
						}
						if($mutation_result eq "undefined"){
							$other_count_rev++;
							$mutation_result = "other_rev";
						}
					}

					if($orientation eq "rev"){
						$match_sequence = reverse_complement_IUPAC($match_sequence);
					}

					## Loop through, report all matching reads and remove from hash
					foreach my $k (keys %updated_reads) {
						my $s = $reads{$k};
						if($s =~/$match_sequence/){

							if($mutation_result eq "alt_fwd"){
								$alt_count_fwd++;
							}
							if($mutation_result eq "ref_fwd"){
								$ref_count_fwd++;
							}
							if($mutation_result eq "other_fwd"){
								$other_count_fwd++;
							}
							if($mutation_result eq "alt_rev"){
								$alt_count_rev++;
							}
							if($mutation_result eq "ref_rev"){
								$ref_count_rev++;
							}
							if($mutation_result eq "other_rev"){
								$other_count_rev++;
							}

							if($output_read_ids){
								print "$k\t$match_sequence\t$mutation_result\t$nr_mismatches\n";
							}

							delete $updated_reads{$k};
						}
					}
				}
				else{
					if(!(defined $hash_off)){
						$unmatched_reads{$sequence}=$nr_mismatches;
					}
					delete $updated_reads{$key};
				}
			}

			if(defined $updated_reads{$key}){
				delete $updated_reads{$key};
			}
		}
	}

	if(!$output_read_ids){
		my $freq_fwd = 0;
		if(($ref_count_fwd+$alt_count_fwd)>0){
			$freq_fwd = $alt_count_fwd/($ref_count_fwd+$alt_count_fwd+$other_count_fwd);
		}
		my $freq_rev = 0;
		if(($ref_count_rev+$alt_count_rev)>0){
			$freq_rev = $alt_count_rev/($ref_count_rev+$alt_count_rev+$other_count_rev);
		}
		my $ref_count_total=$ref_count_fwd+$ref_count_rev;
		my $alt_count_total=$alt_count_fwd+$alt_count_rev;
		my $other_count_total=$other_count_fwd+$other_count_rev;
		my $freq = 0;
		if(($ref_count_total+$alt_count_total+$other_count_total)>0){
			$freq = $alt_count_total/($ref_count_total+$alt_count_total+$other_count_total);
		}
		print "$mutation\t$ref_count_fwd\t$alt_count_fwd\t$other_count_fwd\t$freq_fwd\t$ref_count_rev\t$alt_count_rev\t$other_count_rev\t$freq_rev\t$ref_count_total\t$alt_count_total\t$other_count_total\t$freq\n";
	}
}

# Function for reading a FASTA sequence file
sub read_fasta_file{
	my $fasta_file=shift;

	my %id2seq = ();
	my $id = '';

	open(IN_FILE, "<", $fasta_file) or die "cannot open: $!";

	while(<IN_FILE>){
		chomp;
		if($_ =~ /^>(.+)/){
			$id = $1;
		}
		else{
			$id2seq{$id} .= $_;
		}
	}

	close(IN_FILE);

	return(%id2seq);
}


# Function that reports the first position in a sequence that contains a match to a given
# pattern, allowing for a given number of mismatches.
sub get_best_alignment{
	my $seq = shift;
	my $pat = shift;
	my $max_mismatches = shift;
	my $hash_off = shift;

	my $wt_str="undefined";
	my $mut_str="undefided";

	my $len_in  = length $seq;
	my $len_pat = length $pat;

	my %best_alignment_result;
	my $best_alignment = "undefined";
	my $best_mm = $len_pat;

	for my $i ( 0 .. $len_in - $len_pat ) {

		#Extract k-mer of length = $len_pat from position $i
		my $substr = substr $seq, $i, $len_pat;

		my $current_mm = $len_pat;

		## Only make a new matching if substring has not been seen before!
		## Store all seen substrings in a global hash
		if(!(defined $seen_substrings{$substr})){
			$current_mm = ($pat ^ $substr) =~ tr/\0//c;
			if(!(defined $hash_off)){
				$seen_substrings{$substr} = $current_mm;
			}
		}
		else{
			$current_mm = $seen_substrings{$substr};
		}

		if($current_mm < $best_mm){
			$best_alignment = $substr;
			$best_mm = $current_mm;
			if($best_mm <= $max_mismatches){
				$best_alignment_result{$best_mm}=$best_alignment;
				return(%best_alignment_result);
			}
		}

	}

	$best_alignment_result{$best_mm}=$best_alignment;

	return(%best_alignment_result);
}


# Function to perform a reverse complement of an IUPAC sequence string
sub reverse_complement_IUPAC{

	my $dna = shift;

	# reverse the DNA sequence
	my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
	$revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;

	return $revcomp;
}


# Function for reading arguments and error handling
sub get_args_and_error_check{

	if (@ARGV == 0) {pod2usage(-exitval => 2, -verbose => 0);}

	my ($in_file, $mutation, $max_mismatch, $output_read_ids, $hash_off, $silent);

	my $result = GetOptions("--help"           => sub{local *_=\$_[1];
							                      pod2usage(-exitval =>2, -verbose => 1)},
		                    "-f=s"             =>\$in_file,
							"-m=s"               =>\$mutation,
							"-max_mm=i"          =>\$max_mismatch,
							"-d!"                =>\$output_read_ids,
							"-h_off!"            =>\$hash_off,
		                    "-silent!"           =>\$silent)|| pod2usage(-exitval => 2, -verbose => 1);

	my $error_to_print;

	unless(defined($in_file)) {
		$error_to_print .= "\tNo input FASTA file specified.\n";
	}

    unless(defined($mutation)) {
		$error_to_print .= "\tNo mutation specified.\n";
	}
	else{
		if(!($mutation =~ /^\s*([ACGT]+)\[([ACGT]*)\/([ACGT]*)\]([ACGT]+)\s*$/)){
			$error_to_print .= "\tMutation string (-m) not correctly formatted.\n";
		}
	}

	unless(defined($max_mismatch)) {
		$max_mismatch=5;
	}

	if(defined $error_to_print) {
		my $error_msg="ERROR(s):\n$error_to_print\n";
		pod2usage(-message => $error_msg, -exitval => 2, -verbose => 0);
	}

	else{
		return ($in_file, $mutation, $max_mismatch, $output_read_ids, $hash_off, $silent);
	}
}


__END__

=head1 NAME


cava.pl

=head1 SYNOPSIS

./cava.pl [options] B<--help> B<-f> B<-m> B<-max_mm> B<-d> B<-h_off> B<--silent>

=head1 OPTIONS

=over 8

=item [REQUIRED]

=item B<-f>

FASTA formatted file containing reads.

=item B<-m>

Mutation string to use for screening FASTA file.

=item [OPTIONAL]

=item B<-max_mm>

Maximum number of mismatches allowed in mutation string (default = 5).

=item B<-d>

Print detailed matching info for each read (default = no).

=item B<-h_off>

If set, use of hash tables is reduced in the matching. This slows down the program but reduces memory consumption (default = hash on).

=item B<--silent>

Do not print status to stdout.

=back

=head1 DESCRIPTION

B<This program> performs a matching of mutations in amplicon sequence data.

=cut


