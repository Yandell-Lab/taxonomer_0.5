#!/usr/bin/perl -w
##########################################################
# Author  :  Aurelie Kapusta
# version :  see below + see changelog
# email   :  4urelie.k@gmail.com
# PURPOSE :  Part of Taxonomer make Classifier db util.
#            Add lineages using a list of sequence IDs associated with the lineage
##########################################################
use warnings;
use strict;
use Carp;
use Getopt::Long;

#-----------------------------------------------------------------------------
#------------------------------- DESCRIPTION ---------------------------------
#-----------------------------------------------------------------------------
my $version = "1.0";
my $scriptname = "fasta_add_lineage_using_key.pl";
my $changelog = "
# Change log for $scriptname
#	- v1.0 = 23 May 2016
\n";

my $usage = "\nUsage [v$version]: 

    perl $scriptname -f input.fasta[,input2.fasta...] -t taxonomy.txt [-o <output>] [-d] [-v] [-c] [-h]
    
    DESCRIPTION:
     Add lineages to sequences, with a file containing the correspondance
     If any issue, warnings are issued with the tag WARN 
    
    MANDATORY ARGUMENTS:	    	
     -f,--fasta   => (STRING) input fasta file(s)
                              if several input files, they can be comma separated, 
                              or placed in a directory if -d is used
     -t,--taxo    => (STRING) text file with at least 2 columns, tab separated
                              sequence_id \\t lineage

    OPTIONAL ARGUMENTS:
     -o,--out     => (STRING) output name of the folder with the processed files
                              default = input file name (without extension if any)
     -d,--dir     => (BOOL)   if -i corresponds to a directory containing the fasta files
                              (any file without a proper extension will be ignored, and
                              proper extensions are .fa .faa .fasta .fas)
     -m,--merge   => (BOOL)   Concatenate all files
     -v,--version => (BOOL)   print version
     -c,--chlog   => (BOOL)   print change log (updates)
     -h,--help    => (BOOL)   print this usage\n\n";    

#-----------------------------------------------------------------------------
#------------------------------ LOAD AND CHECK -------------------------------
#-----------------------------------------------------------------------------
my ($fasta,$taxofile,$out,$dir,$merge,$v,$chlog,$help);
my $opt_success = GetOptions(
			 	  'fasta=s'   => \$fasta,
			 	  'taxo=s'    => \$taxofile,
			 	  'out=s'     => \$out,
			 	  'dir'       => \$dir,
			 	  'merge'     => \$merge,
			 	  'chlog'     => \$chlog,
			 	  'version'   => \$v,
			 	  'help'      => \$help,);

#Check options, if files exist, etc
die "\n --- $scriptname version $version\n\n" if ($v);
die $usage if ($help);
die $changelog if ($chlog);
die "\n SOME MANDATORY ARGUMENTS MISSING, CHECK USAGE with -h or --help\n\n" if (! $fasta || ! $taxofile);
die "\n -l $fasta does not exist?\n\n" if (! -e $fasta);
die "\n -l $taxofile does not exist?\n\n" if (! -e $taxofile);
die "\n -f $fasta is not a fasta file? Maybe you forgot to use --dir?\n\n" unless (($dir) || ($fasta =~ /\.(fa$|fasta$|faa$|fas$)/));

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
print STDERR "\n --- $scriptname (v$version)\n";
print STDERR "     input = $fasta\n";
print STDERR "        (is a directory)\n" if ($dir);
print STDERR "     taxonomy info = $taxofile\n";
$fasta = $1 if (($dir) && ($fasta =~ /^(.*)\/$/)); #remove / at the end if any
if (! $out) {
	$out = $fasta;
	$out =~ s/(.*?).(fa$|fasta$|faa$|fas$)/$1/;
}
if ($dir) {
	print STDERR " --- Making output directory $out\n";
	my $outdir = $out.".lineages";
	`rm -Rf $outdir` if (-e $outdir);
	`mkdir $outdir` if ($dir);
}
($dir)?($dir="y"):($dir="n");
#get input files
my @files = ();
($dir)?(@files = `ls $fasta`):(push(@files,$fasta));

#Loading superfamilies from -s
print STDERR " --- Loading taxonomy from $taxofile\n";
my $taxo = get_taxo($taxofile);

#Add lineages
print STDERR " --- Adding lineage information\n";
my $check = add_lineage(\@files,$fasta,$dir,$out,$taxo);

#check if nd stuff
my ($totalfa,$totalnd) = (0,0);
foreach my $check_nd (keys %{$check}) {
	$totalfa++ unless ($check->{$check_nd}==0);
	$totalnd+=$check->{$check_nd};
}
if ($totalfa > 0) {
	print STDERR "     WARN: there are $totalfa files with sequences without lineage, listed below:\n";
	foreach my $check_nd (keys %{$check}) {
		print STDERR "           $check_nd\n" unless ($check->{$check_nd}==0);
	}
	print STDERR "     => for a total of $totalnd sequences without lineage\n";
}

#Concatenate if relevant
if (($dir eq "y") && ($merge)) {
	print STDERR " --- Concatenating fasta files\n";
	`cat $out.lineages/*.fa > $out.lineages.fa`;
	print STDERR "     => $out.lineages.fa\n";
}
print STDERR " --- Script done\n\n";
exit;


#-----------------------------------------------------------------------------
#-------------------------------- SUBROUTINES --------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# get taxonomy, 1 column file
# my $taxo = get_taxo($taxofile);
#-----------------------------------------------------------------------------
sub get_taxo {
	my($in) = shift;
	my %taxo = ();	
	open(my $fh, "<", $in) or confess "ERROR (sub get_taxo): could not open taxonomy file $in $!\n";
	LIN: while(<$fh>) {
		chomp (my $lin = $_);
		next LIN if (($lin !~ /\w/) || (substr($lin,0,1) eq "#")); #skip blank lines and comments
		my @lin = split(/\t/,$lin);
		$taxo{$lin[0]} = $lin[1];
	}
	close($fh);
	return \%taxo;
}

#-----------------------------------------------------------------------------
# Add lineages for all files
# my $check = add_lineage(\@files,$fasta,$dir,$out,$taxo);
#-----------------------------------------------------------------------------
sub add_lineage {
	my ($files,$fasta,$dir,$out,$taxo) = @_;		
	my $check = ();
	FA: foreach my $fa (@{$files}) {
		chomp($fa);
		#check fasta format
		unless ($fa =~ /(.*?).(fa$|fasta$|faa$|fas$)/) { #I could be less stringent, but better for lineage finding that they have proper extension
#			my $head = `$head -n 1 $fa`;
#			if (substr($head,0,1) ne ">") {
				print STDERR "     WARN: $fa not a fasta file? Skipped.\n";
				next FA;
#			}
		}
#		print STDERR "     -> $fa\n";
		$check->{$fa}=0;
		#Deal with path
		my $fullfa = $fa;
		$fullfa = $fasta."/".$fa if ($dir eq "y");
		#output file
		my $outfa;
		$outfa = $out.".lineages.fa" if ($dir eq "n");
		$outfa = $out.".lineages/".$fa if ($dir eq "y");
		$outfa =~ s/(.*?).(fa$|fasta$|faa$|fas$)/$1.lineages.fa/ if ($dir eq "y");
		#now open and lopp
#		print STDERR "     -> $outfa\n";
		open(my $fh, "<", $fullfa) or confess "ERROR (sub add_lineage): could not open to read $fullfa $!\n";
		open(my $fhout, ">", $outfa) or confess "ERROR (sub add_lineage): could not open to write $outfa $!\n";
		LINE: while(<$fh>) {
			chomp (my $line = $_);
			next LINE if ($line !~ /\w/); #skip blank lines
			if (substr($line,0,1) eq ">") {
				$line =~ s/^>//;
				$line =~ s/\s+/\t/ if ($line !~ /\t/); #make the first space a tab if no tabs detected
				my ($seqid,$desc) = split(/\t/,$line);
 				$desc = "" unless ($desc);
				$seqid =~ s/_hg38//; #since not in the taxonomy file
				my $lineage = "nd";
				$lineage = $taxo->{$seqid} if ($taxo->{$seqid});
				$check->{$fa}++ if ($lineage eq "nd");
				my @seqid = split(/\s+/,$line);
 				my $fullseqid = join("_",@seqid);
				print $fhout ">$fullseqid\t$lineage\n";
				print STDERR "No lineage for $seqid ($fullseqid)\n" if ($lineage eq "nd");
			} else {
				print $fhout $line."\n";
			}	
		}
		close $fh;
		close $fhout;
	}
	return ($check);
}



