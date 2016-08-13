#!/usr/bin/perl -w
##########################################################
# Author  :  Aurelie Kapusta
# version :  see below / see change log
# email   :  4urelie.k@gmail.com
# PURPOSE :  Parser for the Taxonomer Web version output: convert json to a 2 columns table
##########################################################
use warnings;
use strict;
use Carp;
use Getopt::Long;
use JSON::Parse 'json_file_to_perl';
use Data::Dumper;

my $version = "2.0";
my $scriptname = "TaxoWeb_json-to-tab.pl";
my $changelog = "
#	- v1.0 = 20-24 June 2016
#	- v2.0 = 12 Aug 2016
#            Bug fix - output was incomplete (was not going back up the levels while looping)
 

\n";

my $usage = "\nUsage [v$version]: 
    perl $scriptname 
            -i <taxonomer_output.json> [-d] [-p] 
           [-o <X>] [-s <X>] [-c <X>] [-v] [-h] [-l]
    			
    PURPOSE:
    Convert the taxonomer .json output from the www.taxonomer.com web interface
    to a 2 columns table containing the number of reads classified at each level (=lineages)
    as well as a lighter output (shorter lineage strings; everything before the phylum, when any, is kicked out, and all \"no_rank:\" are removed
		
    MANDATORY ARGUMENTS:	
     -i,--in      => (STRING) input files = output file in json format of taxonomer.com
                              If several files to parse, you can separate file names with commas: \",\"
                              or if -d is set, -i can designate a directory that contains the files to parse 
                              (only .json files in that directory will be parsed)
 
    OPTIONAL ARGUMENTS
     -d,--dir     => (BOOL)   if -in is a directory and not a file or a list of file
     -p,--pool    => (BOOL)   Relevant only if several input files
                              To concatenate all outputs (individual files are still created)
     -o,--out     => (STRING) Use this to define the core name of output files for --pool [irrelevant otherwise]
                              By default, first file of the list provided as input is used, and there will be CAT in the name.
     
     -v,--v       => (BOOL)   verbose, makes the script talk to you
     -v,--v       => (BOOL)   print version if only option
     -h,--help    => (BOOL)   print this usage
     -l,--log     => (BOOL)   print the change log (updates)
\n";

################################################################################
### Get arguments/options
### check some of them, print details of the run if verbose chosen
################################################################################
my ($tied,$score,$fconf,$outname) = ("na","na","na","na");
my $tieout = 1;
my ($in,$d,$cat,$cluster,$fa,$key,$lin_are_ref,$help,$chlog,$v);
GetOptions ('in=s'      => \$in,  
            'out=s'     => \$outname,
            'dir'       => \$d, 
            'pool'      => \$cat, 
            'score=s'   => \$score, 
            'conf=s'    => \$fconf, 
            'log'       => \$chlog, 
            'help'      => \$help, 
            'v'         => \$v);

#check step for options
die "\n $scriptname version: $version\n\n" if ((! $in) && (! $help) && (! $chlog) && ($v));
die $changelog if ($chlog);
die $usage if ($help);
die "\n Please provide input file (-i, use -h to see the usage)\n\n" if (! $in);
die "\n -i $in does not exist?\n\n" if (($in !~ /,/) && (! -e $in));
$in = $1 if ($in =~ /^(.*)\/$/); #remove the / at the end

################################################################################
### MAIN
################################################################################
########################################
# VERBOSE STUFF
########################################
print STDERR "\n--------------------------------------------------\n" if ($v);
print STDERR " --- $scriptname started (v$version), with:\n" if ($v);
print STDERR "       Input files are located in $in\n" if (($v) && ($d));
print STDERR "       --in  => $in\n" if (($v) && (! $d));
print STDERR "             -> if there were any previous output files in $in, they will be deleted\n" if (($v) && ($d));
print STDERR "             -> outputs will be concatenated\n" if (($v) && ($cat));
print STDERR "                in files with core name = $outname\n" if (($v) && ($cat) && ($outname ne "na"));
print STDERR "     Filtering options:\n" if ($v);
print STDERR "          Threshold for maximum score is set to $score; if the score < $score, classification is ignored\n" if (($v) && ($score ne "na"));
print STDERR "          Threshold for confidence score is set to $fconf; if the confidence < $fconf, classification is ignored\n" if (($v) && ($fconf ne "na"));
print STDERR "          no filter on maximum score or confidence score\n" if (($v) && ($score eq "na") && ($fconf eq "na"));

########################################
# ACTUAL PARSING
########################################
print STDERR "--------------------------------------------------\n" if ($v);
print STDERR " --- Get list of input file(s) and delete previous output files\n" if ($v);
($d)?($d="y"):($d="n");
#get list of input files + delete previous outputs
my $listin = ();
($listin,$outname) = input_files($in,$outname,$d,$v);

# now loop through input file(s) to convert them
################################################################################
print STDERR " --- Input file(s) set with -i will be converted from json to tab\n" if ($v);
INPUT: foreach my $in (@{$listin}) {
	chomp $in;
	print STDERR " --- Parsing the classifier output $in \n" if ($v);
	unless (-e $in) {
		print STDERR "     WARN: (main): $in does not exist? skipping\n";
		next INPUT;
	}
	print STDERR "     Loading the data\n" if ($v);
	my $p = json_file_to_perl($in);		
	my $data = load_json($p,$v);
	print STDERR "     Printing the data in $outname.tab\n" if ($v);
	print STDERR "     Printing lighter data in $outname.light.tab\n" if ($v);
	print STDERR "        it is \"light\" because to avoid super long lineage strings,\n" if ($v);
	print STDERR "        everything before the phylum (when any) is kicked out, and all \"no_rank:\" are removed\n" if ($v);
	print_json($data,$outname,$v);
}	
print STDERR " --- Done\n" if ($v);
print STDERR "--------------------------------------------------\n\n" if ($v);
exit;



################################################################################
### SUBROUTINES
################################################################################
#----------------------------------------------------------------------------
# from a filename or a directory keep only its path - independently of the current dir
# my $path = path($filename);
#----------------------------------------------------------------------------
sub get_path {
	my($file) = shift;
	($file =~ /\//)?($file =~ s/(.*)\/.*$/$1/):($file = ".");
	return $file;
}

#----------------------------------------------------------------------------
# get list of input files + delete previous outputs
# ($listin,$outname) = input_files($in,$outname,$d,$v);
#----------------------------------------------------------------------------
sub input_files {
	my($in,$outname,$d,$v) = @_;
	my $path;
	my @listin = ();
	if ($d eq "y") {
		$path = $in;
		my @namelist = `ls $path`;
		NAME: foreach my $name (@namelist) {
			chomp($name);
			$name = "$path/$name";
			if (substr($name,-4,4) eq ".tab") {
				`rm -Rf $name`;
				print STDERR "        rm -Rf $name\n" if ($v);
				next NAME; 
			} elsif (substr($name,-5,5) eq ".json") {
				push(@listin,$name);
			}	
		}
		$outname = $path if ($outname eq "na");
	} elsif ($in =~ /,/) {		
		@listin = split(",",$in);
		$path = get_path($listin[0]);
		foreach my $f (@listin) {
			`rm -Rf $f.*.tab`;
			print STDERR "        rm -Rf $f.*.tab\n" if ($v);
		}
		$outname = $listin[0] if ($outname eq "na");
	} else {
		if ($outname eq "na") {
			$outname = $in;
			$outname =~ s/\.json//;
		}
		`rm -Rf $outname.tab`;
		`rm -Rf $outname.light.tab`;
		print STDERR "        rm -Rf $outname.tab\n" if ($v);
		print STDERR "        rm -Rf $outname.light.tab\n" if ($v);
		push(@listin,$in);
		
	}
	`rm -Rf $outname.CAT.*.parsed.tab` if ($cat);
	print STDERR "        rm -Rf $outname.CAT.*.parsed.tab\n" if (($cat) && ($v));
	$outname =~ s/\s/_/g;
	return (\@listin,$outname);
}

#----------------------------------------------------------------------------
# load the json
# my $data = load_json($p,$v);$data = load_json($p,0);
#----------------------------------------------------------------------------
sub load_json {
	my($p,$v) = @_;
	my $data = ();
	my $lineage = ();
	my $level = 0;
	foreach my $bin (keys $p->{'binnerResult'}) {
		my $lin = $level."__".$bin.";";
		$lineage->{$level}{$bin}{$bin}=$lin;
		$data->{$level}{$bin}{$lin}=$p->{'binnerResult'}{$bin}; #get counts for that lineage = the binning level
	}
 	foreach my $bin (keys $p->{'classifierResult'}) {
 		($data,$lineage) = dump_hash($p->{'classifierResult'}{$bin},$bin,$data,$lineage,++$level);	
 		--$level;
 	}	
	return ($data);
} 

sub dump_hash {
	my ($ref,$bin,$data,$lineage,$level) = @_;
	my $name = "undef";
	my $count = 0;
	foreach my $key ( keys %{$ref} ) {
		if ( ref $ref->{$key} eq 'ARRAY' ) { #that means the children => enter
			my $i=0;
			until (! $ref->{'children'}[$i]) { #loop on ALL the children
				#get the full lineage of the child beforehand
				my $lin;
				my $childname = $ref->{'children'}[$i]{'name'};
				$childname =~ s/\s/_/g;
				$name =~ s/\s/_/g;
				my $lvl = sprintf("%02d", $level);
				my $lpar = sprintf("%02d", $level-1);
				($lineage->{$bin}{$lpar}{$name})?($lin = $lineage->{$bin}{$lpar}{$name}.$lvl."__".$name.";"):($lin=$bin."\t".$lvl."__".$name.";");
				$lineage->{$bin}{$lvl}{$childname}=$lin; #store that lineage => will be the parental
				$data->{$lin}=$count;
				#Now go deeper
				($data,$lineage) = dump_hash($ref->{'children'}[$i],$bin,$data,$lineage,++$level);
				--$level;
				$i++;
			}	 
		} else {
			#not children => the rest
			$name = $ref->{$key} if ($key eq "name");
			$count = $ref->{$key} if ($key eq "count");
		}
	}
	return ($data,$lineage);
  
# 
# #NB, subroutine above based on this:
# #https://www.daniweb.com/programming/software-development/threads/461091/how-to-find-the-depth-of-perl-hash-which-has-no-consistent-structure
# sub dump_data {
#     my ( $ref, $level ) = @_;
#     die "You can only pass an HASH reference data"
#       unless ref $ref eq 'HASH';
#     for ( keys %{$ref} ) {
#         --$level;
#         if ( ref $ref->{$_} eq 'HASH' ) {
#             dump_data( $ref->{$_}, ++$level );    ## sub called itself
#             --$level;
#         }
#         else { print $ref->{$_}, $/ }
#     }
# }

}

#----------------------------------------------------------------------------
# print the json file in a tab format
# print_json($data,$outname,$v);
#----------------------------------------------------------------------------
sub print_json {
	my($data,$outname,$v) = @_;
	open(my $fh, ">", "$outname.tab") or confess "     \nERROR (sub print_json): could not open to write $outname.tab $!\n";
	open(my $fh2, ">", "$outname.light.tab") or confess "     \nERROR (sub print_json): could not open to write $outname.light.tab $!\n";
	foreach my $lin (sort keys %{$data}){
		print $fh "$lin\t$data->{$lin}\n";	
		my $lightlin = $lin;
		$lightlin =~ s/^([A-Za-z]+?\t).+?([0-9][0-9]__phylum:)/$1$2/;
		$lightlin =~ s/no_rank://g;
		print $fh2 "$lightlin\t$data->{$lin}\n";
	}
	close $fh;
	close $fh2;
	return 1;
}
