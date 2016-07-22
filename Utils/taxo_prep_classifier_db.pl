#!/usr/bin/perl -w
##########################################################
# Author  :  Aurelie Kapusta; GitHub 4ureliek
# version :  see below + see changelog
# email   :  4urelie.k@gmail.com
# PURPOSE :  Part of Taxonomer make Classifier db util
#            From a fasta file with lineages in their headers, create all files required to make a classifier db
##########################################################
use warnings;
use strict;
use Carp;
use Getopt::Long;

#-----------------------------------------------------------------------------
#------------------------------- DESCRIPTION ---------------------------------
#-----------------------------------------------------------------------------
#flush buffer
$| = 1;

my $version = "1.4";
my $scriptname = "wd_prep_classifier_db.pl";
my $changelog = "
# Change log for $scriptname
#	- v1.0 = 22 March 2016
#	- v1.1 = 07 April 2016
#            print sequences in uppercase if -u
#	- v1.2 = 12 April 2016
#            Fix bug to get the whole description
#	- v1.3 = 21 July 2017
#            Allow lineages to not have the same length in the input file
#            => they will be put at the same length when populated down
#            Be more flexible on the lineage formatting
#	- v1.4 = 22 July 2017
#            Bug fix in writing .fa and .sti (was expecting __ in lineages)
\n";

my $usage = "\nUsage [v$version]: 
    perl $scriptname -f input.fasta [-o coredbname] [-s ;] [-t] [-w] [-u] [-v] [-c] [-h]
    
    PURPOSE:
     Create the files required to build a classification database with Taxonomer 

    DESCRIPTION:
     Input = fasta file containing lineage information
     (one way to add these lineages is to use the script fasta_add_lineage_using_key.pl)
     The fasta file will be parsed to obtain the files required, 
     so the formatting of the headers and descriptions is important:
       - NO SPACES ALLOWED in the seqID, or in the lineage
       - lineage should be a serie of AT LEAST ONE node, with a semicolon (;) as the default separator (see -s)
       - sequences may have lineages of different length - they will be populated down
         if there is a rank (example 1) but not a rank value, the rank value will be populated down
         if only the first node has a rank value, a rank value will be added to all nodes
         it is OK if for some sequences nodes are empty, and others only have rank values
       - Examples that will be processed properly:
           example 1: >Seq_ID1    00__rankValue0;01__rankValue1;02__rankValue2;03__rankValue3
                      >Seq_ID2    00__rankValue0;01__rankValue1;02__;03__
                      >Seq_ID3    00__rankValue0;01__rankValue1;02__
                      >Seq_ID4    00__rankValue0;01__rankValue1
                      >Seq_ID5    00__rankValue0;rankValue1;rankValue2;rankValue3
           example 2: >Seq_ID1    Prokaryotes;Bacteria;Proteobacteria;Ecoli
                      >Seq_ID2    Prokaryotes;Bacteria;Proteobacteria;Senterica

       - Examples that will NOT be processed properly:
           example 1: >Seq_ID1    rankValue0;01__rankValue1;02__rankValue2
                      >Seq_ID2    rankValue0;01__rankValue1

     The fasta file will be parsed and create files required to build a classification database:
      -> get the taxid map = .tri
           >1
           0
           >child_ID
           parent_ID
      -> get the seq_id map = .sti
           >seq_id
           tax_id 
      -> a modified fasta file = _taxonomer.fa
         Will have a shifted header, to have a numerical ID matching the .sti 
         Use grep \"^>\" <coredb>_taxonomer.fa to obtain the correspondance between ID / header 
         Note that lineages will be populated down, so that all sequences have a lineage of the same length
    
     Two files that can be useful to parse the outputs will be created:
      -> a key file _key.txt file, with 4 columns seperated by a tab: 
         taxid \t parent_taxid \t rank \t lineage_string
      -> an index file _index.txt, with 3 columns seperated by a tab: 
         taxid \t rank \t label=last_value_of_lineage 
    
    MANDATORY ARGUMENTS:	    	
     -f,--fasta   => (STRING) input fasta file with lineage information, see above for its formatting
     
    OPTIONAL ARGUMENTS:
     -o,--out     => (STRING) output core name of the db (corename.tri, corename.sti...)
                              default = input file without fasta extension
     -s,--sep     => (STRING) to set the separator for the nodes in the lineage
                              default = ;
                              special characters will be annoying to use, better to avoid
     -t,--test    => (BOOL)   will exit if missing lineages or duplicated sequences
     -w,--warn    => (BOOL)   print warnings (details of missing lineages or duplicated sequences)
     -u,--uc      => (BOOL)   print sequences in uc for the file _taxonomer.fa
                              (kmers with bases in lower cases are ignored when the db is built)
     -v,--version => (BOOL)   print version
     -c,--chlog   => (BOOL)   print change log (updates)
     -h,--help    => (BOOL)   print this usage\n\n";    

#-----------------------------------------------------------------------------
#------------------------------ LOAD AND CHECK -------------------------------
#-----------------------------------------------------------------------------
my ($fasta,$out,$test,$warn,$uc,$v,$chlog,$help);
my $sep = ";";
my $opt_success = GetOptions(
			 	  'fasta=s'   => \$fasta,
			 	  'out=s'     => \$out,
			 	  'sep=s'     => \$sep,
			 	  'test'      => \$test,
			 	  'warn'      => \$warn,
			 	  'uc'        => \$uc,
			 	  'chlog'     => \$chlog,
			 	  'version'   => \$v,
			 	  'help'      => \$help,);

#Check options, if files exist, etc
die "\n --- $scriptname version $version\n\n" if $v;
die $changelog if ($chlog);
die "\n SOME MANDATORY ARGUMENTS MISSING, CHECK USAGE:\n$usage" if ($help || ! $fasta);
die "\n -f $fasta is not a fasta file?\n\n" unless ($fasta =~ /\.fa|faa|fasta|fas$/);
die "\n -l $fasta does not exist?\n\n" if (! -e $fasta);


#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
print STDERR "\n --- $scriptname (v$version)\n";
print STDERR "     input = $fasta\n";
if (! $out) {
	$out = $fasta;
	$out =~ s/(.*?)\.fa$|\.faa$|\.fasta$|\.fas$/$1/;
}
print STDERR "     output core name = $out\n";
print STDERR " --- Getting fasta IDs\n"; #could be all done at the same time, but that way taxIDs will be smaller numbers than the seqIDs
my @headers = `grep "^>" $fasta`;
my $nb_seqs = scalar(@headers);
print STDERR "     => total of $nb_seqs sequence headers to process\n";
print STDERR " --- Checking lineages for each seqIDs to get the max. number of nodes\n";
my $maxlvl = get_maxlvl(\@headers,$sep);
print STDERR "     => $maxlvl\n";

print STDERR " --- Parsing headers to load lineages for each seqIDs with filled empty nodes (up to $maxlvl)\n";
my ($taxonomy,$nb_duplicated_seqids,$nb_no_lineage,$indexmap) = load_headers(\@headers,$sep,$maxlvl,$warn);
if (($nb_duplicated_seqids > 0) || ($nb_no_lineage > 0)) {
	print STDERR " --- Parsing done, but:\n";
	print STDERR "     there were $nb_duplicated_seqids duplicated sequences (same id)\n" if ($nb_duplicated_seqids > 0);
	print STDERR "     there were $nb_no_lineage sequences without lineage information\n" if ($nb_no_lineage > 0);
	die "     => EXITING (--test chosen)\n\n" if ($test);
}
print STDERR " --- Obtaining tax id map (tri map) + taxids of each sequence\n";
my ($taxIDmap,$seqIDmap,$counts) = get_tri_sti($taxonomy);
print STDERR " --- Print .tri, _key.txt and _index.txt files\n";
print_tri($taxIDmap,$indexmap,$out);
print STDERR " --- Obtaining the sti map + rewriting the fasta file with populated down lineage and numerical IDs\n";
fasta_and_sti($fasta,$seqIDmap,$taxonomy,$out,$counts,$uc);
print STDERR " --- DONE\n\n";
exit;


#-----------------------------------------------------------------------------
#-------------------------------- SUBROUTINES --------------------------------
#-----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# get number of level (depth)
# my $maxlvl = get_maxlvl(\@headers,$sep);
#----------------------------------------------------------------------------
sub get_maxlvl {
	my ($headers,$sep) = @_;
	my $maxlvl = 0;
	HEADER: foreach my $header (@{$headers}) {
		chomp($header);
		$header =~ s/^>//;
		my ($seqID,$lin) = split(/\s+/,$header);
		next HEADER if ((! $lin) || ($lin eq "nd")); #missing ones will be counted later
		my @nodes;
		($lin =~ /;/)?(@nodes = split($sep,$lin)):($nodes[0]=$lin);
		$maxlvl = $#nodes+1 unless ($#nodes+1 < $maxlvl);
	}
	return ($maxlvl);
}

#----------------------------------------------------------------------------
# Load headers
# my ($taxonomy,$nb_duplicated_seqids,$nb_no_lineage,$indexmap) = load_headers(\@headers,$sep,$maxlvl,$warn);
#----------------------------------------------------------------------------
sub load_headers{
	my ($headers,$sep,$maxlvl,$warn) = @_;
	my $nb_duplicated_seqids = 0;
	my $nb_no_lineage = 0;
	my $taxonomy = ();
	my $indexmap = ();
	HEADER: foreach my $header (@{$headers}) {
		chomp($header);
		$header =~ s/^>//;
		my ($seqid,$lin) = split(/\s+/,$header); #description useless here, but could be shifted to be as the lineage, so check that as well:
		if ((! $lin) || ($lin eq "nd")) {
			$nb_no_lineage++;
			print STDERR "        WARN: no lineage found for $seqid -> skipping sequence\n" if ($warn);
			next HEADER;			
		} 
		my $poplin;
		($poplin,$indexmap) = pop_down($lin,$sep,$maxlvl,$indexmap); #$poplin = array, $indexmap = hash		
		if ($taxonomy->{$seqid}) {
			$nb_duplicated_seqids++;
			print STDERR "        WARN: duplicated entry for $seqid -> skipping sequence\n" if ($warn);
			next HEADER;
		} else {
			$taxonomy->{$seqid}=$poplin;#$poplin = array
		}
	}	
	return ($taxonomy,$nb_duplicated_seqids,$nb_no_lineage,$indexmap);
}

#----------------------------------------------------------------------------
# Populate down the lineage + get the index map, with rank for each name
# my ($poplin,$indexmap) = pop_down($lin,$sep,$maxlvl,$indexmap); #$poplin = array, $indexmap = hash
#----------------------------------------------------------------------------
sub pop_down {
	my ($lin,$sep,$maxlvl,$indexmap) = @_;
	my @nodes = ();	
	if ($lin =~ /$sep/) {
		@nodes = split($sep,$lin);
	} else {
		$nodes[0]=$lin; #in case it's a one node taxonomy => no need to carry on; still need to enter this sub to get the indexmap		
		$indexmap->{$lin}{'rank'} = 0;
		$indexmap->{$lin}{'label'} = $nodes[0];
		$indexmap->{$lin}{'label'} = $1 if ($indexmap->{$lin}{'label'} =~ /^[0-9]+?__(.+?)$/);
		return (\@nodes,$indexmap);
	}
	my $lineage_string = "";
	for (my $i = 0; $i < $maxlvl; $i++){
		#account for empty nodes values but that have rank values => $nodes[$i] would be defined
		if (($nodes[$i]) && ($nodes[$i] =~ /__/)) {
			my ($rank,$label) = split("__",$nodes[$i]);
			my ($prevrank,$prevlabel) = split("__",$nodes[$i-1]) unless ($label);
			$nodes[$i]=$rank."__".$prevlabel unless ($label);
		} else {
			#need to check if the previous node has a rank value => increment that
			if ($nodes[$i-1] =~ /__/) {
				my ($prevrank,$prevlabel) = split("__",$nodes[$i-1]);
				my $newrank;
				(substr($prevrank,0,1) eq 0)?($newrank="0".$i):($newrank = $prevrank++);
				($nodes[$i])?($nodes[$i]=$newrank."__".$nodes[$i]):($nodes[$i]=$newrank."__".$prevlabel);
			} else {
				$nodes[$i]=$nodes[$i-1] unless ($nodes[$i]); #basically assign the parent one to that nodes, but only if empty
			}	
		}
		($lineage_string)?($lineage_string = $lineage_string.";".$nodes[$i]):($lineage_string = $nodes[$i]); #current lineage - may miss levels, but for the indexmap I need all of them
		$indexmap->{$lineage_string}{'rank'} = $i;
		$indexmap->{$lineage_string}{'label'} = $nodes[$i];
		$indexmap->{$lineage_string}{'label'} = $1 if ($indexmap->{$lineage_string}{'label'} =~ /^[0-9]+?__(.+?)$/);
		}
	return (\@nodes,$indexmap); #the @nodes array has the populated down lineages
}

#----------------------------------------------------------------------------
# Get taxID map from the lineages that are stored
# my ($taxIDmap,$seqIDmap,$counts) = get_tri_sti($taxonomy);
#----------------------------------------------------------------------------
sub get_tri_sti {
	my $taxonomy = shift;	
	my $taxid = 1;
	my $seq_id = 1;
	my %taxIDmap = ();
	my %seqIDmap = ();
	foreach my $seqID (keys %{$taxonomy}) {
		my $lineage = "";
		#For each level, check if exists, otherwise create
		foreach my $rank (@{$taxonomy->{$seqID}}) {
			($lineage eq "")?($lineage = $rank):($lineage = $lineage.";".$rank);
			if (! $taxIDmap{$lineage}) { #loop through all nodes, basically
				$taxIDmap{$lineage}{'own'}=$taxid;
				$taxIDmap{$lineage}{'parent'}=$taxid-1; #the parent will necessarily be just one value above
				$taxid++;
			}
			$seqIDmap{$seqID}=$taxIDmap{$lineage}{'own'}; #each seqID, the last taxid (before the ++) is the one to remember; no possible duplicated seqid at this point           
		}	
	}
	return(\%taxIDmap,\%seqIDmap,$taxid); #lasttaxid
}

#----------------------------------------------------------------------------
# Print files now
# print_tri($taxIDmap,$indexmap,$out);
#----------------------------------------------------------------------------
sub print_tri {
	my ($taxIDmap,$indexmap,$out) = @_;	
	my $key = $out."_key.txt";
	my $idx = $out."_index.txt";	
	my %trimap = ();
	open(my $keyfh, ">", $key) or confess "ERROR (sub print_tri): could not open to write $key $!\n";
	print $keyfh "#taxID\tparent_taxID\trank\tlineage_string\n";
	open(my $idxfh, ">", $idx) or confess "ERROR (sub print_tri): could not open to write $idx $!\n";
	foreach my $lineage (sort keys %{$taxIDmap}) {
		#NB: tri does not need to be sorted, but that's 'prettier'
		$trimap{$taxIDmap->{$lineage}{'own'}}=$taxIDmap->{$lineage}{'parent'};
		my ($rank,$label) = ($indexmap->{$lineage}{'rank'},$indexmap->{$lineage}{'label'});
		print $keyfh "$taxIDmap->{$lineage}{'own'}\t$taxIDmap->{$lineage}{'parent'}\t$rank\t$lineage\n"; #taxid \t parent_taxid \t lineage_string		
		print $idxfh "$taxIDmap->{$lineage}{'own'}\t$rank\t$label\n"; #taxid \t rank \t label
	}
	close ($keyfh);
	close ($idxfh);
	
	#print a sorted tri
	my $tri = $out.".tri";
	open(my $trifh, ">", $tri) or confess "ERROR (sub print_tri): could not open to write $tri $!\n";
	foreach my $child (sort {$a<=>$b} keys %trimap) {
		print $trifh ">$child\n$trimap{$child}\n";
	}
	close ($trifh);
	return 1;
}

#----------------------------------------------------------------------------
# sti stuff now
# fasta_and_sti($fasta,$seqIDmap,$taxonomy,$out,$counts,$uc);
#----------------------------------------------------------------------------
sub fasta_and_sti {
	my ($fa,$seqIDmap,$taxonomy,$out,$taxid,$uc) = @_;
	#$taxid+=1;
	my $sti = $out.".sti";
	$out = $out."_taxonomer.fa";
	my $skip;
	my %check = ();
	open(my $fh, "<", $fa) or confess "ERROR (sub fasta_and_sti): could not open to read $fa $!\n";
	open(my $stifh, ">", $sti) or confess "ERROR (sub fasta_and_sti): could not open to write $sti $!\n";
	open(my $fafh, ">", $out) or confess "ERROR (sub fasta_and_sti): could not open to write $out $!\n";
	LINE: while(<$fh>) {
		chomp (my $line = $_);
		if (substr($line,0,1) eq ">") {
			$line =~ s/^>//;
			my @header = split(/\s+/,$line);
			my ($seqid,$lin) = ($header[0],$header[1]);
			my @desc = @header;
			splice(@desc,0,2); #remove the first 2 elements; id_name & lineage
			my $desc = ""; #just in case			
			$desc = join(" ",@desc) if ($desc[0]); #re-build the description if any
			#checks
			if ((! $lin) || ($lin eq "nd") || ($check{$seqid})){
				$skip = 1;
			} else {	
				my $lineage = join(";",@{$taxonomy->{$seqid}});
				print $fafh ">$taxid\t$seqid\t$lineage\t$desc\n";
				print $stifh ">$taxid\n$seqIDmap->{$seqid}\n";
				$check{$seqid}=1;
				$skip = 0;
				$taxid++;
			}
		} else {
			if ($uc) {
				print $fafh uc($line)."\n" unless ($skip == 1);
			} else {
				print $fafh $line."\n" unless ($skip == 1);
			}	
		}	
	}
	close $fh;
	close $stifh;
	close $fafh;
	return 1;
}

