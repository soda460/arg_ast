#!/usr/bin/perl
use strict;
use List::MoreUtils qw(firstidx);


# Ouverture des fichiers passes en arguments dans la ligne de commande
unless (open(CLE, "$ARGV[0]")) {
	die ("Cannot open file with the corresponing names (key values)!\n");
}
unless (open(ACTUTAL_NAMES, "$ARGV[1]")) {
	die ("Cannot open  litteral card amr names!\n");
}

my @clef = <CLE>;
chomp @clef;

my @litteral = <ACTUTAL_NAMES>;
chomp @litteral;

my %dict;

for my $line (@clef) {
	my @key_value = split("\t", $line);
	$dict{$key_value[1]} = $key_value[0];
}

# tests
print($dict{'tet(W/N/W)'}, "\n");
print($dict{'MIR-12'}, "\n");

my $my_output_string = "";

# Step 2 , iterate on litteral ARG
for my $actual_name (@litteral){

	if (exists $dict{ $actual_name }) {
		print ($actual_name, " will be replaced by ", $dict{ $actual_name }, "\n");
		$my_output_string .= $dict{ $actual_name } . "\t";
	} else {
		print("not found ", $actual_name, "\n");
		$my_output_string .= '!' . $actual_name . "\t";
	}

}

print "----\n";
print $my_output_string;


open (OUT, ">edited_lane");

# Header
print OUT ($my_output_string);



