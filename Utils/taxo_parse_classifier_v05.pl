#!/usr/bin/perl -w
##########################################################
# Author  :  Aurelie Kapusta
# email   :  4urelie.k@gmail.com
# PURPOSE :  Parser for taxonomer classifier with custom database
##########################################################
use warnings;
use strict;
use Carp;
use Getopt::Long;
#use Data::Dumper;

my $version = "1.5";
my $scriptname = "taxo_parse_classifier_v05.pl";
my $changelog = "
#	- v1.0 = 19-21 Jun 2017
#            Version of the parser for taxonomer_0.5 classifier & custom databases 
#	- v1.1 = 22 Jun 2017
#            taxid output was missing the taxid and only had the lineage
#	- v1.2 = 23 Jun 2017
#            --agg option
#	- v1.3 = 07 Aug 2017
#            fix --agg option... matrix was not initialized with 0s
#	- v1.4 = 09 Aug 2017
#            Add option --ties in mandatory; rename the filter that was named -t to -u
#              => fix count of classified sequences for --ties 1
#                 if --ties 0 then ties have to be skipped for refid output
#            Deal with full tax IDs => lineage as taxid in the output
#	- v1.5 = 13 Jan 2018
#            Small bug fix and adapt with the new -a option from the prep script
\n";

my $usage = "\nUsage [v$version]: 
    perl $scriptname 
             -i <classifier_output> [-d] [-f <file.fa>] [-k <db_key.txt>] 
             [-n] [-a] [-p <X>] [-b] [-o <X>] [-t <X>] [-s <X>] [-v] [-h] [-l]
			
    PURPOSE:
    Parser for taxonomer outputs (custom databases only), GitHub version of the tool.
    It outputs read counts for each taxid, and each reference if the database fasta file 
    is provided. Note that for this output, when reads are tied to several reference sequences
    and taxonomer was run with --ties 1, the counts are divided by the number of ties 
    for each reference. However, if --ties 0 was used then only the taxid output is reliable.
    See the usage for filtering and normalization.

    To concatenate existing parsed outputs, use the -b option. Besides -i, other 
    mandatory arguments become non mandatory with this option. The -n option can 
    still be used, but then -f needs to be set, and only the pooled output will 
    be normalized. Typically:
    perl $scriptname -i <classifier_output_directory> -d -p 5 -b -v [other filtering options]
		
    MANDATORY ARGUMENTS:	
     -i,--in      => (STRING) input files = output file(s) of the taxonomer classification (classify_reads)
                              If several files to parse, you can separate file names with commas: \",\"
                              or if -d is set, -in can designate a directory that contains the files to parse 
                              (do not put other files than the ones to parse in it)
                              To get a matrix, use -p           
     -t,--ties    => (INT)    Same as when classify_reads.py was run; if --ties 1, then there is one line
                              per tied reference; with --ties 0 there is one line per read.
     -m,--map     => (STRING) Provide the root of the database .tri and .sti:
                              if name.tri and name.sti, then use -m name
                              This allows to get the parent taxid of each refname 
                              (also needed to get the lineage in the refname and aggrname outputs)
                                                            
    OPTIONAL ARGUMENTS:
     BEHAVIOR
     -f,--fa      => (STRING) To get counts at the reference ID level, set -f with the fasta file
                              used to build the database (it has the correspondance between seqID 
                              and reference name in its headers). Also required for -a and -n.
     -k,--key     => (STRING) To get the lineages and not the taxids in the taxid output,
                              load the key file (_key.txt) from the taxo_prep_classifier.pl script             
     -n,--norm    => (BOOL)   To normalize counts for each reference by ref length (requires -f)
                              and by total amount of classified reads of each sample for all outputs
                              Also x 10e9 for readability -> formula of RPKM (Mortazavi 2008)
                              Note that the taxid output will be unchanged since reference length is needed
     -a,--agg     => (BOOL)   To aggregate the refname counts to the NAME part of the header, 
                              if it is formmated as: >NAME__details    [description]
                              This is done once the numbers for the refname output are obtained,
                              in order to have proper normalization (since all copies of a 'NAME' 
                              are of different length). It will work with -b, but not if -n is set as well.
                              The fasta file (set with -f) is also required.
   
     IF SEVERAL INPUT FILES        
     -d,--dir     => (BOOL)   Set if -i is a directory and not a file or a list of file
     -p,--pool    => (FLOAT   To concatenate all outputs (individual files are still created)
                              Relevant only if several input files: the read counts
                              for each file will be in different columns.
                              Set the min. % of samples that should have > 0 read count
                              to a reference (for it to be kept in the concatenated output);
                              use -p 0 to keep all reference sequences.
     -b,--before  => (BOOL)   To load existing parsed outputs. Use with -d only.
                              Make sure to set a different -o if you wish to keep 
                              a previous concatenated file
     -o,--out     => (STRING) Use this to define the core name of output files for -p 
                              [irrelevant otherwise]

     FOR FILTERING
     -u,--uniq    => (INT)    To skip sequences that have more than (>) X ties
                              Typically: -u 10
     -s,--score   => (INT)    To skip sequences that have a max score lower than (<) the value X
                              Typically: -s 40                          
                                   
     OTHERS
     -v,--v       => (BOOL)   verbose, makes the script talk to you
     -v,--v       => (BOOL)   print version if only option
     -h,--help    => (BOOL)   print this usage
     -l,--log     => (BOOL)   print changelog
\n";

#-----------------------------------------------------------------------------
#--------------------------- CHECK & VERBOSE ---------------------------------
#-----------------------------------------------------------------------------
my ($in,$ties);
my ($d,$key,$name,$fa,$norm,$agg,$cat,$before,$outname);
my ($tied,$score);
my ($help,$chlog,$v);
GetOptions ('in=s'      => \$in, 
            'ties=s'    => \$ties, 
            'key=s'     => \$key,
            'map=s'     => \$name,
            'fa=s'      => \$fa,
            'norm'      => \$norm,
            'agg'       => \$agg,
            'dir'       => \$d, 
            'pool=s'    => \$cat,
            'before'    => \$before, 
            'out=s'     => \$outname, 
            'uniq=s'    => \$tied, 
            'score=s'   => \$score, 
            'log'       => \$chlog, 
            'help'      => \$help, 
            'v'         => \$v);

#check step for options
check_options();

#Print log
print_log() if ($v);

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
#get list of input files + delete previous outputs unless $cat
my $summaryfile;
my @listin = ();
input_files();

#load the key file if any
my %map = ();
load_key() if ($key && ! $before);

#load the sti and tri files as one hash to get the parent of a reference
my %sti = ();
my $FILE;
load_sti_tri();

#Load the headers
my %rfinfo = ();
load_fa() unless (! $fa);

#now parse
my %totals = ();
init_totals() unless ($before && ! $norm);
my %cat_tx = ();
my %cat_ref = ();
my %cat_agg = ();
my %allkeys = ();
my ($txheader,$rfheader,$rnheader);
set_headers();
print STDERR " --- Parsing input file(s):\n" if (! $before && $v);
print STDERR " --- Loading previous output files:\n" if ($before && $v);
for my $file (@listin) {
	if ($before) {
		load_parsed($file,"taxid") if ($file =~ /taxid/);
		load_parsed($file,"refname") if ($file =~ /refname/ && $fa);
		load_parsed($file,"aggrname") if ($file =~ /aggrname/ && $fa); #this is why having the 3 options -d -n -a is not okay
	} else {
		my ($counts,$refcounts) = parse_taxo($file); #keys are the integer IDs
		#Now print & load concat stuff
		print_parsed($file,"taxid",$counts,$txheader);
		my $rncounts = print_parsed($file,"refname",$refcounts,$rfheader) if ($fa);
		print_parsed($file,"aggrname",$rncounts,$rnheader) if ($fa && $agg);
	}	
}

#Now print cat if relevant
if (defined $cat) {
	fill_zeros_cat();
	print STDERR " --- Printing concatenated output for all samples\n" if ($v);
	print_cat("taxid",\%cat_tx,$txheader);
	print STDERR "     -> $outname.CAT.taxid.tab\n" if ($v);
	print_cat("refname",\%cat_ref,$rfheader) if ($fa);
	print STDERR "     -> $outname.CAT.refname.tab\n" if ($fa && $v);
	print_cat("aggrname",\%cat_agg,$rnheader) if ($fa && $agg);
	print STDERR "     -> $outname.CAT.aggrname.tab\n" if ($fa && $agg && $v);
}

#Now print summary
print_summary() if (! $before);
 
print STDERR " --- DONE\n\n" if ($v);
exit;

#-----------------------------------------------------------------------------
#------------------------------- SUBROUTINES ---------------------------------
#-----------------------------------------------------------------------------

#----------------------------------------------------------------------------
sub input_files {
	print STDERR " --- Get list of input file(s) and delete previous output files\n" if (! $before && $v);
	print STDERR " --- Get list of previous output files\n" if ($before && $v);
	my $path;
	if ($d) {
		$path = $in;
		my @namelist = `ls $path`;
		NAME: foreach my $name (@namelist) {
			chomp($name);
			$name = "$path/$name";
			if ($name =~ /(taxid|refname|aggrname)\.tab$/) {
				`rm -Rf $name` if (! $before);
				push(@listin,$name) if ($before && $name !~ /summary\.tab$/);
				next NAME; 
			}
			push(@listin,$name) if (! $before);
		}
		$outname = $path unless ($outname);
		$summaryfile = $outname.".summary.tab";
	} elsif ($in =~ /,/) { #no need to check with $before here - die check upstream
		@listin = split(",",$in);
		$path = get_path($listin[0]);
		foreach my $f (@listin) {
			`rm -Rf $f.*.taxid.tab`;
			`rm -Rf $f.*.refname.tab`;
		}
		$outname = $listin[0] unless ($outname);
		$summaryfile = $outname.".summary.tab";
	} else {
		`rm -Rf $in.*.tab`;
		push(@listin,$in);
		$path = get_path($in);
		$summaryfile = $in;
		$summaryfile =~ s/\.out$//;
		$summaryfile = $summaryfile.".summary.tab";
		$outname = $in unless ($outname);
	}
	`rm -Rf $summaryfile` unless ($before);
	`rm -Rf $outname.CAT.*.tab` if (defined $cat);
	return;
}

#----------------------------------------------------------------------------
sub load_key {
	print STDERR " --- Loading correspondance: lineages <=> TaxIDs from $key\n" if ($v);
	open(my $fh, "<", $key) or confess "     \nERROR (sub load_key): could not open to read $key $!\n";
	LINE: while(<$fh>) {
		chomp (my $l = $_);
		next LINE if ($l !~ /\w/ || $l =~ /^#/);
		my ($tx, $parent, $rank, $lineage) = split(/\s+/,$l);
		$map{$tx} = $lineage;
	}
	close ($fh);
	return;
}

#----------------------------------------------------------------------------
sub load_sti_tri {
	my $child;	
	#check first if sti contains any info or if it has id under itself
	$FILE = $name.".sti";	
	my @lines = `head -n 2 $FILE`;
	my $ch = ($lines[0]);
	my $pa = ($lines[1]);
	$ch =~ s/^>//;
	chomp($ch);
	chomp($pa);
	$FILE = $name.".tri" if ($ch == $pa);
	print STDERR " --- Loading correspondance: seqid <=> TaxIDs from $FILE\n" if ($v);

	#Now load - tri has all info if -a was used to prep the db files	
	open(my $fh, "<", $FILE) or confess "     \nERROR (sub load_sti_tri): could not open to read $FILE $!\n";
	LINE: while(defined(my $l = <$fh>)) {
		chomp $l;
		next LINE if ($l !~ /\w/);		
		if (substr($l,0,1) eq ">") {
			$child = $l;
			$child =~ s/^>//;
		} else {		
			$sti{$child}=$l;
		}
	}
	close ($fh);
	return;
}

#----------------------------------------------------------------------------
sub load_fa {
	print STDERR " --- Getting ref_names & length from $fa\n" if ($v);
	my ($id,$name);
	my $length = 0;
	my $c = 0;
	open (my $fa_fh, "<", $fa) or confess "     \nERROR (sub load_fa): could not open to read $fa $!\n";
	while (<$fa_fh>) {
		chomp (my $l = $_);
		if (substr($l,0,1) eq ">") {	
			#First thing, store length of previous sequence
			if ($c != 0) {
				#save with the $name so that -n can be done with the before option
				$rfinfo{$name}{'ln'}=$length; 
				$length = 0;
			}
			$c=1;	
			#get the id
			my @id = split(/\s+/,$l);
			$id = shift @id;
			$id =~ s/>//;
			$name = shift @id;			
			#load refname			
			$rfinfo{$id}{'id'}=$name;
			#rebuild description = what's left
			my $desc = join(" ",@id);			
			$rfinfo{$id}{'dc'}=$desc;
		} else {
			$length+=length($l);
		}
	}
	close ($fa_fh);
	#store length for the last sequence
	$rfinfo{$name}{'ln'}=$length;
	return;
}

#----------------------------------------------------------------------------
sub init_totals {
	foreach my $file (@listin) {
		$totals{$file}{'c'} = 0;
		$totals{$file}{'u'} = 0;
		$totals{$file}{'p'} = 0;
	}
	return;
}

#----------------------------------------------------------------------------
sub set_headers {
	$txheader = "classified_taxid";
	$rfheader = "ref_id\tref_name\tref_len\tref_desc\ttaxid";
	$rnheader = "ref_name_aggregated";
	if ($key) {
		if ($FILE =~ /\.tri$/) {
			$txheader = "lineage"; #not one taxid per lineage in that case
		} else {
			$txheader = $txheader."\tlineage";
		}
		$rfheader = $rfheader."\tlineage";
		$rnheader = $rnheader."\tlineage";
	}
	return;
}

#----------------------------------------------------------------------------
sub load_parsed {
	my $file = shift;
	my $type = shift;
	print STDERR "      - $file...\n" if ($v);
	open(my $fh, "<", $file) or confess "     \nERROR (sub load_taxo): could not open to read $file $!\n";	
	LINE: while(<$fh>) {
		chomp (my $line = $_);
		next LINE if ($line !~ /\w/ || substr($line,0,1) eq "#");			
		my @l = split(/\t/,$line);
		my $c = pop @l;
		my $k = join("\t",@l);
		$cat_tx{$k}{$file} = $c if ($type eq "taxid");
		$cat_ref{$k}{$file} = $c if ($type eq "refname");
		$cat_agg{$k}{$file} = $c if ($type eq "aggrname");
		#Save all keys
		$allkeys{$type}{$k}++;
		#Count totals if needed
		$totals{$file}{'c'}+=$c if ($norm);
	}
	close ($fh);
	return;
}

#----------------------------------------------------------------------------
sub parse_taxo {
	my $file = shift;
	my %counts = ();
	my %refcounts = ();
	print STDERR "      - $file\n" if ($v);
	open(my $in_fh, "<", $file) or confess "     \nERROR (sub load_taxo): could not open to read $file $!\n";	
##CLASSIFIER OUTPUT (COMMAND LINE):
#	U	query_name
#	C	query_name	taxid	refid	rank	ties	max_score	avg_score	query_len
	LINE: while(<$in_fh>) {
		chomp (my $line = $_);
		next LINE if ($line !~ /\w/);			
		#load values
		my ($class,$query,$taxid,$ref,$depth,$nt,$maxsc,$avgsc,$rlen,$kc) = split(/\t/,$line);
		if ($class eq "U") { 
			$totals{$file}{'u'}++;
			next LINE;
		}
		#OK, now classified stuff:
		if ($ties) {
			$totals{$file}{'c'}+= 1/$nt;
		} else {
			$totals{$file}{'c'}++;
		}
		#filter if relevant
		next LINE if ($score && $maxsc < $score);
		next LINE if ($tied && $nt > $tied);
		#went through filters => count reads + by taxid and refid
		my $k = $taxid;
		#Now get the lineage corresponding to that taxID if -k set		
		if ($key) {
			if ($map{$k}) {
				$k = $map{$k};
			} elsif ($map{$sti{$k}}) {
				$k = $map{$sti{$k}};
			}					
		}
		if ($ties) {
			#normalize if tied since one line per tie
			$totals{$file}{'p'}+= 1/$nt;
			$counts{$k}+=1/$nt; #counts per taxids
			$refcounts{$ref}+=1/$nt; #counts per refids
		} else {
			$totals{$file}{'p'}++;
			$counts{$k}++;
			if ($nt == 1) {
				#then the refid is OK
				$refcounts{$ref}+=1/$nt;
			}
		}
	}
	close ($in_fh);
	$totals{$file}{'tot'} = $totals{$file}{'c'} + $totals{$file}{'u'};
	return(\%counts,\%refcounts);
}

#----------------------------------------------------------------------------
sub print_parsed {
	my $file = shift;
	my $type = shift;
	my $counts = shift;
	my $header = shift;
	my %rncounts = ();
	my $out = $file;
	$out =~ s/\.out$//;
	$out = $out.".".$type.".tab";
	print STDERR "        -> in $out\n" if ($v);
	open(my $fh, ">", $out) or confess "     \nERROR (sub print_parsed): could not open to write $out $!\n";
	print $fh "#$header\tread_count";
	print $fh "_normalized" if ($norm);
	print $fh "\n";
	foreach my $k (keys %{$counts}) {
		my $c = $counts->{$k};
		$c=$c/$totals{$file}{'c'} if ($norm); #normalize if asked
		my $fullkey;
		($fullkey,$c) = get_fullkey_counts($k,$type,$c);	
		#for cat - save all + initialize all keys
		$allkeys{$type}{$fullkey}++ if (defined $cat);
		$cat_tx{$fullkey}{$file} = $c if (defined $cat && $type eq "taxid");
		$cat_ref{$fullkey}{$file} = $c if (defined $cat && $type eq "refname");
		$cat_agg{$fullkey}{$file} = $c if (defined $cat && $type eq "aggrname");
		print $fh "$fullkey\t$c\n";
		#for --agg, if parsing the refname right now:
		if ($agg && $type eq "refname") {
			my ($name,$details) = split("__",$rfinfo{$k}{'id'});
			$rncounts{$name}+=$c;
			#add Rnames to the key hash => store lineages
			my $tx = $sti{$k};
			my $lin = $map{$tx};
			$map{$name} = $lin if (! $map{$name});
		}
	}	
	close $fh;
	return \%rncounts;
}

#----------------------------------------------------------------------------
sub get_fullkey_counts {
	my $k = shift;
	my $type = shift;
	my $c = shift;
	my $fullkey = "";
	if ($type eq "taxid") {
		$fullkey = $k;
		$fullkey = $fullkey."\t".$map{$k} if ($key && $map{$k});
	} elsif ($type eq "refname") {
		my $name = $rfinfo{$k}{'id'};
		my $len = "nd";
		$len = $rfinfo{$name}{'ln'} if ($rfinfo{$name}{'ln'});
		$c = $c/$len*10e9 if ($norm);
		my $tx = $sti{$k};
		my $desc = $rfinfo{$k}{'dc'};
		#clean up the description with the lineage if it is in it
		if ($key) {
			my $lineage = $map{$tx};
			$desc =~ s/$lineage\s//;
		}
		$fullkey = "$k\t$name\t$len\t$desc\t$tx";
		$fullkey = $fullkey."\t".$map{$tx} if ($key);
	} elsif ($type eq "aggrname") {
		$fullkey = $k;
		$fullkey = $fullkey."\t".$map{$k} if $map{$k};
	}
	return($fullkey,$c);
}	

#----------------------------------------------------------------------------
sub fill_zeros_cat {
	foreach my $file (@listin) {
		TYPE: foreach my $type (keys %allkeys) {
			next TYPE if ($before && $file !~ /$type/);
			foreach my $fullkey (keys %{$allkeys{$type}}) {
				if ($type eq "taxid") {				
					$cat_tx{$fullkey}{$file}=0 unless ($cat_tx{$fullkey}{$file});
				} elsif ($type eq "refname") {
					$cat_ref{$fullkey}{$file}=0 unless ($cat_ref{$fullkey}{$file});
				} else {	
					$cat_agg{$fullkey}{$file}=0 unless ($cat_agg{$fullkey}{$file});
				}
			}
		}
	}
	return;
}

#----------------------------------------------------------------------------
sub print_cat {
	my $type = shift;
	my $all = shift;
	my $header = shift;
		
	#get name of outputs
	my $out = $outname.".CAT.".$type.".tab";

	open(my $fh, ">", $out) or confess "     \nERROR (sub print_cat): could not open to write $out $!\n";
	print $fh "$header";
	KEYS: foreach my $k (keys %{$all}) {
		foreach my $file (sort keys %{$all->{$k}}) {
			my $filename = filename($file);
			$filename =~ s/\.$type\.tab$// if ($before);
			print $fh "\t$filename";
		}
		print $fh "\n";
		last KEYS;
	}	
	#now print values
	my $tot = scalar(@listin);
	$tot = grep {$_ =~ /$type/} @listin if ($before);
	REFS: foreach my $fullkey (keys %{$all}) {
		if (defined $cat && $cat != 0) {
			my $per = $allkeys{$type}{$fullkey} / $tot * 100;
			next REFS if ($per < $cat);
		}
		print $fh "$fullkey";
		foreach my $file (sort keys %{$all->{$fullkey}}) {
			my $c = 0;
			$c = $all->{$fullkey}{$file} if ($all->{$fullkey}{$file});	
			if ($before && $norm) {	#divide the counts if needed, for refs only
				$c = $c/$totals{$file}{'c'};
				my ($k,$stuff)=split(/\t/,$fullkey) unless ($type eq "taxid");
				$c = $c/$rfinfo{$k}{'id'}*10e9 unless ($type eq "taxid");
			}
			print $fh "\t$c";
		}
		print $fh "\n";
	}
	close $fh;
	return;
}

#----------------------------------------------------------------------------
sub print_summary {
	print STDERR " --- Read count summary in $summaryfile\n" if ($v);
	open(my $fh, ">", $summaryfile) or confess "     \nERROR (main): could not open to write $summaryfile $!\n";
	print $fh "#File\ttotal_nb\tclassified\tunclassified\tclassified_and_passed_filters\n";
	foreach my $file (@listin) {
		print $fh "$file\t$totals{$file}{'tot'}\t$totals{$file}{'c'}\t$totals{$file}{'u'}\t$totals{$file}{'p'}\n";
	}
	close $fh;
	return;
}

#----------------------------------------------------------------------------
sub get_path {
	my($file) = shift;
	($file =~ /\//)?($file =~ s/(.*)\/.*$/$1/):($file = ".");
	return $file;
}

#----------------------------------------------------------------------------
sub filename {
	my($name) = shift;
	$name =~ s/.*\/(.*)$/$1/;
	return $name;
}

#----------------------------------------------------------------------------
sub check_options {
	die "\n $scriptname version: $version\n\n" if (! $in && ! $help && ! $chlog && $v);
	die $changelog if ($chlog);
	die $usage if ($help);
	die "\n Please provide input file (-i, use -h to see the usage)\n\n" if (! $in);
	die "\n -i $in does not exist?\n\n" if ($in !~ /,/ && ! -e $in);
	die "\n Please set --ties to 0 or 1 (use -h to see the usage)\n\n" if (! defined $ties || $ties !~ /^(0|1)$/);
	die "\n -p is set but is not a float number?\n\n" if (defined $cat && $cat !~ /[\.0-9]+?/);
	die "\n Please provide core name of sti and tri file (-m, use -h to see the usage)\n\n" if (! $name);
	die "\n -m $name.sti or $name.tri do not exist?\n\n" if (! -e $name.".sti" || ! -e $name.".tri");
	die "\n -k $key does not exist?\n\n" if (! $before && $key && ! -e $key);
	die "\n -f $fa does not exist?\n\n" if (! $before && $fa && ! -e $fa);
	die "\n -b is set but -d is not - please check usage\n\n" if ($before && ! $d);
	die "\n -b is set but -p is not - please check usage\n\n" if ($before && ! defined $cat);
	die "\n -n is set but -f is not - please check usage\n\n" if ($norm && ! $fa);
	die "\n -a is set but -f is not - please check usage\n\n" if ($agg && ! $fa);
	$in = $1 if ($in =~ /^(.*)\/$/);
	return;
}

#----------------------------------------------------------------------------
sub print_log {
	print STDERR "\n --- $scriptname started (v$version), with:\n";
	print STDERR "       Input files are located in $in\n" if ($d);
	print STDERR "       --in  => $in\n" if (! $d);
	print STDERR "             -> classsify_reads.py was run with --ties "; 
	print STDERR "1\n" if ($ties);
	print STDERR "0\n" if (! $ties);
	print STDERR "             -> previous output files in $in will be loaded\n" if ($d && $before);
	print STDERR "             -> if there were any previous output files in $in, they will be deleted\n" if ($d && ! $before);
	print STDERR "             -> outputs will be concatenated\n" if (defined $cat);
	print STDERR "                in files with core name = $outname\n" if (defined $cat && $outname);
	print STDERR "             -> references with less than $cat % samples with classified reads will be skipped\n" if (defined $cat);
	print STDERR "       --map [sti or tri file] => $name.sti or $name.tri\n";
	print STDERR "       --key => $key\n" if ($key && ! $before);
	print STDERR "       --fa  => $fa\n" if (! $before && $fa);
	print STDERR "             -> additional output for counts per reference\n" if (! $before && $fa);
	print STDERR "       --norm => counts will be divided by reference length and total number of classified reads\n" if ($norm);
	if ($agg && $before && $norm) {
		print STDERR "       --agg is set, but so are -b and -n => ignored\n";
		undef $agg;
	}
	print STDERR "       --agg  => counts will be aggreated to the NAME part of a refid formatted as NAME__DETAILS\n" if ($agg);
	print STDERR "       Filtering options:\n";
	print STDERR "         Threshold for number of ties is set to $tied, so if a sequence has a number\n" if ($tied && ! $before);
	print STDERR "         of references tied for max read weight > $tied, classification is ignored\n" if ($tied && ! $before);
	print STDERR "         Threshold for score is set to $score; if the score < $score, classification is ignored\n" if (! $before && $score);
	print STDERR "         No filter on score, confidence or ties number\n" if (! $before && ! $tied && ! $score);
	return;
}
