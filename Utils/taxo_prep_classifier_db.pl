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

my $version = "1.9";
my $scriptname = "taxo_prep_classifier_db.pl";
my $changelog = "
# Change log for $scriptname
#	- v1.0 = 22 March 2016
#	- v1.1 = 07 April 2016
#            print sequences in uppercase if -u
#	- v1.2 = 12 April 2016
#            Fix bug to get the whole description
#	- v1.3 = 21 July 2016
#            Allow lineages to not have the same length in the input file
#            => they will be put at the same length when populated down
#            Be more flexible on the lineage formatting
#	- v1.4 = 22 July 2016
#            Bug fix in writing .fa and .sti (was expecting __ in lineages)
#	- v1.5 = 27 Oct 2016
#            new output = tgs file
#	- v1.6 = 31 Jan 2017
#            Option to get sti and tri with a tab instead
#            Bug fix in the tri map...!!
#	- v1.7 = 16 Mar 2017
#            Bug fix for lineages without 00__ etc
#	- v1.8 = 22 Mar 2017
#            Populating down lineages is now an option
#            Bug fix in the tri, again
#	- v1.9 = 21 Jun 2017
#            Bug fix in the tri, again; was adding a layer

# TO DO = detect if same seqIDs have different lineages, output warnings, and allow correction / chose 1 by deleting lines in a file, then rerun script with an option loading that file to correct lineages
\n";

my $usage = "\nUsage [v$version]: 
    perl $scriptname 
         -f input.fasta [-o coredbname] [-n] [-m <X>] [-s ;]
         [-g map.tab] [-p] [-t] [-w] [-u] [-v] [-c] [-h]
    
    PURPOSE:
     From a fasta file with lineages in their headers,
     creates the files required to build a Taxonomer classification database
     (.tri, .sti, .tgs) 

    INPUT FILE:
     A fasta file containing lineage information
     (one way to add these lineages is to use the script fasta_add_lineage_using_key.pl)
     The fasta file will be parsed to obtain the files required, so the formatting 
     of the headers and descriptions is important and should be:
        >seqID    lineage    description if any
     Also, note that:
       - lineage is anything that is after the first space or tab,
         so no spaces are allowed in seqID or lineages
       - Lineage levels (or nodes) are identified by a separator (set with -s), for ex.:
           Prokaryotes;Bacteria;Proteobacteria;Ecoli [-s ;]
           T,HLA,A*01 [-s ,]
           gene01:transcript0001:exon1 [-s :]
       - sequences may have lineages of different length, use -p to populate them down.
         Note that one type of lineage format will be specifically recognized:
            00__rankValue0;01__rankValue1;02__rankValue2;03__rankValue3
            (the ranks will be incremented when populated down)

    OUTPUTS:
     The fasta file will be parsed and create files required to build a classification database:
      -> get the taxid map = .tri
           >1
           0
           >child_ID
           parent_ID
      -> get the seq_id map = .sti
           >seq_id
           tax_id
      -> if -g is set, get the flexible taxonomy map = .tgs
           taxID \\t geneID \\t seqID \\t info[string_no_spaces]
      -> a modified fasta file = _taxonomer.fa
         Will have a shifted header, to have a numerical ID matching the .sti 
         Use grep \"^>\" <coredb>_taxonomer.fa to obtain the correspondance between ID / header 
    
     Two files that can be useful to parse the outputs will be created:
      -> a key file _key.txt file, with 4 columns seperated by a tab: 
         taxid \t parent_taxid \t rank \t lineage_string
      -> an index file _index.txt, with 3 columns seperated by a tab: 
         taxid \t rank \t label=last_value_of_lineage 
         (required for iobio visualization)
    
    MANDATORY ARGUMENTS:	    	
     -f,--fasta   => (STRING) input fasta file with lineage information, see above for its formatting
     
    OPTIONAL ARGUMENTS:
     -o,--out     => (STRING) output core name of the db (corename.tri, corename.sti...)
                              default = input file without fasta extension
                              (should be the database name)
     -n,--new     => (BOOL)   Get tabulated tri and sti files
                              sti: <id> \\t <taxid>
                              tri: <child_id> \\t <parent_id> 
                                   (with a mandatory line = 1 \\t 0)
     -m,--mirr    => (INT)    To set the first seqID number
                              0 is reserved, so if 0 is set, it will be 1
                              Default: last taxID+1
     -s,--sep     => (STRING) to set the separator for the nodes in the lineage
                              (special characters will be annoying to use, better to avoid)
                              default = ;
     -g,--genes   => (STRING) to also get the .tgs file.
                              getting the geneID requires a file with the geneID for each seqID:
                                 seqID \\t geneID
                              The file may contain only some seqIDs; 
                              the others will get a geneID = seqID
                              if no file at all, use -g na                             
     -p,--pop     => (BOOL)   To populate down the lineages, so that all sequences have a lineage of the same length
                              if there is a rank but not a rank value, the rank value will be populated down:
                                 00__rankValue0;01__rankValue1;02__;03__
                              will be:
                                 00__rankValue0;01__rankValue1;02__rankValue1;03__rankValue1
                              The lineage format is very flexible except for one thing:
                              if any node has a __ in it, the ones before need to as well
                              (because it expects that __ is the separator for rank and value)                        
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
my ($fasta,$out,$tab,$genes,$pop,$setid,$test,$warn,$uc,$v,$chlog,$help);
my $sep = ";";
my $opt_success = GetOptions(
			 	  'fasta=s'   => \$fasta,
			 	  'out=s'     => \$out,
			 	  'new'       => \$tab,
			 	  'mirr=s'    => \$setid,
			 	  'genes=s'   => \$genes,
			 	  'sep=s'     => \$sep,
			 	  'pop'       => \$pop,
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
die "\n -f $fasta does not exist?\n\n" if (! -e $fasta);
die "\n -m $setid is not an integer!\n\n" if (($setid) && ($setid =~ /\D/));
die "\n -g $genes does not exist?\n\n" if (($genes) && ($genes ne "na") && (! -e $genes));
($tab)?($tab = "y"):($tab = "n");
($pop)?($pop = "y"):($pop = "n");
($uc)?($uc = "y"):($uc = "n");

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
print STDERR "\n --- $scriptname (v$version)\n";
print STDERR "     input = $fasta\n";
if (! $out) {
	$out = $fasta;
	$out =~ s/\.fa$|\.faa$|\.fasta$|\.fas$//;
}
print STDERR "     output core name = $out\n";

print STDERR " --- Getting fasta IDs\n"; #could be all done at the same time, but that way taxIDs will be smaller numbers than the seqIDs
my @headers = `grep "^>" $fasta`;
my $nb_seqs = scalar(@headers);
print STDERR "     => total of $nb_seqs sequence headers to process\n";
my $maxlvl = 1;
if ($pop eq "y") {
print STDERR " --- Checking lineages for each seqIDs to check their depth\n";
	$maxlvl = get_maxlvl(\@headers,$sep);
	print STDERR "     => largest number of nodes in any lineage = $maxlvl\n";
	print STDERR "        (separator for nodes = $sep)\n";
}
print STDERR " --- Parsing headers to load lineages for each seqIDs \n";
print STDERR "     (separator for nodes = $sep)\n" if ($pop eq "n");
print STDERR "     (with filled empty nodes up to $maxlvl = populating down)\n" if ($pop eq "y");
my ($taxonomy,$nb_duplicated_seqids,$nb_no_lineage,$indexmap) = load_headers(\@headers,$sep,$maxlvl,$pop,$warn);
if (($nb_duplicated_seqids > 0) || ($nb_no_lineage > 0)) {
	print STDERR " --- Parsing done, but:\n";
	print STDERR "     there were $nb_duplicated_seqids duplicated sequences (same id)\n" if ($nb_duplicated_seqids > 0);
	print STDERR "     there were $nb_no_lineage sequences without lineage information\n" if ($nb_no_lineage > 0);
	die "     => EXITING (--test chosen)\n\n" if ($test);
}
print STDERR " --- Obtaining tax id map (tri map) + taxids of each sequence\n";
my ($taxIDmap,$seqIDmap,$counts) = get_tri_sti($taxonomy,$indexmap,$sep);
$counts = $setid if ($setid);
$counts++ if $counts == 0;

print STDERR " --- Print .tri, _key.txt and _index.txt files\n";
print STDERR "     (tri in tabulated format)\n" if ($tab eq "y");
print_tri_key_idx($taxIDmap,$indexmap,$out,$tab);

print STDERR " --- Obtaining the sti map and writing it + rewriting the fasta file with integer IDs\n";
print STDERR "     (with lineages populated down)\n" if ($pop eq "y");
print STDERR "     (sti in tabulated format)\n" if ($tab eq "y");
print_fasta_sti_tgs($fasta,$seqIDmap,$taxonomy,$out,$counts,$uc,$tab,$genes);
print STDERR " --- DONE\n\n";
exit;


#-----------------------------------------------------------------------------
#-------------------------------- SUBROUTINES --------------------------------
#-----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# Load the geneIDs
# my $genID = load_geneIDs($genes) if ($genes) && ($genes ne "na");
#----------------------------------------------------------------------------
sub load_geneIDs {
	$genes = shift;
	my %ids = ();
	open(my $fh, "<", $genes) or confess "ERROR (sub load_geneIDs): could not open to read $genes $!\n";
	while (<$fh>) {
		chomp(my $l = $_);
		next LINE if ((substr($l,0,1) eq "#") || ($l =~ /\W/));
		my @l = split('\s+',$l);
		$ids{$l[0]}=$l[1];
	}
	close ($fh);
	return(\%ids);	
}

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
		($lin =~ /$sep/)?(@nodes = split($sep,$lin)):($nodes[0]=$lin);
		$maxlvl = $#nodes+1 unless ($#nodes+1 < $maxlvl);
	}	
	return ($maxlvl);
}

#----------------------------------------------------------------------------
# Load headers
# my ($taxonomy,$nb_duplicated_seqids,$nb_no_lineage,$indexmap) = load_headers(\@headers,$sep,$maxlvl,$pop,$warn);
#----------------------------------------------------------------------------
sub load_headers{
	my ($headers,$sep,$maxlvl,$pop,$warn) = @_;
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
		($poplin,$indexmap) = get_map($lin,$sep,$maxlvl,$indexmap,$pop); #$poplin = array, $indexmap = hash		
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
# my ($poplin,$indexmap) = get_map($lin,$sep,$maxlvl,$indexmap); #$poplin = array, $indexmap = hash
#----------------------------------------------------------------------------
sub get_map {
	my ($lin,$sep,$maxlvl,$indexmap,$pop) = @_;
	my @nodes = ();	
	($lin =~ /$sep/)?(@nodes = split($sep,$lin)):($nodes[0]=$lin);
	$maxlvl = $#nodes +1 if ($pop eq "n"); #Won't populate down
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
		#current lineage - may miss levels, but for the indexmap I need all of them
		($lineage_string)?($lineage_string = $lineage_string.$sep.$nodes[$i]):($lineage_string = $nodes[$i]); 
		$indexmap->{$lineage_string}{'rank'} = $i;
		$indexmap->{$lineage_string}{'label'} = $nodes[$i];
		$indexmap->{$lineage_string}{'label'} = $1 if ($indexmap->{$lineage_string}{'label'} =~ /^[0-9]+?__(.+?)$/);
	}
	return (\@nodes,$indexmap); #the @nodes array has the populated down lineages
}

#----------------------------------------------------------------------------
# Get taxID map from the lineages that are stored
# my ($taxIDmap,$seqIDmap,$counts) = get_tri_sti($taxonomy,$indexmap,$sep);
#----------------------------------------------------------------------------
sub get_tri_sti {
	my ($taxonomy,$indexmap,$sep) = @_;	
	my $taxid = 1; #0 is the root
	my $seq_id = 1;
	my %taxIDmap = ();
	my %seqIDmap = ();
	my $prev_lineage = "";
	foreach my $seqID (keys %{$taxonomy}) {
		my $lineage = "";				
		foreach my $rank (@{$taxonomy->{$seqID}}) { #loop through all nodes
			($lineage eq "")?($lineage = $rank):($lineage = $lineage.$sep.$rank);
			if (! $taxIDmap{$lineage}{'own'}) { #no taxid yet
				$taxIDmap{$lineage}{'own'}=$taxid; #set as the previously incremented taxid
				if ((! $taxIDmap{$prev_lineage}{'own'}) || ($indexmap->{$lineage}{'rank'} == 0)) {
					$taxIDmap{$lineage}{'parent'}=0; #when under the root
				} else {
					$taxIDmap{$lineage}{'parent'}=$taxIDmap{$prev_lineage}{'own'};
				}
				$taxid++;
			}
			$prev_lineage = $lineage;	
			#each seqID, the last taxid (before the ++) is the one to remember; no possible duplicated seqid at this point 	
			$seqIDmap{$seqID}=$taxIDmap{$lineage}{'own'};   
		}		
	}
	return(\%taxIDmap,\%seqIDmap,$taxid); #lasttaxid+1 returned => will be the first seqID
}

#----------------------------------------------------------------------------
# Print files now
# print_tri_key_idx($taxIDmap,$indexmap,$out);
#----------------------------------------------------------------------------
sub print_tri_key_idx {
	my ($taxIDmap,$indexmap,$out) = @_;	
	my $key = $out."_key.txt";
	my $idx = $out."_index.txt";	
	my %trimap = ();
	$trimap{1}=0;
	open(my $keyfh, ">", $key) or confess "ERROR (sub print_tri_key_idx): could not open to write $key $!\n";
	print $keyfh "#taxID\tparent_taxID\trank\tlineage_string\n";
	open(my $idxfh, ">", $idx) or confess "ERROR (sub print_tri_key_idx): could not open to write $idx $!\n";
	LIN: foreach my $lineage (sort keys %{$taxIDmap}) {
		next LIN unless ($lineage);
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
	open(my $trifh, ">", $tri) or confess "ERROR (sub print_tri_key_idx): could not open to write $tri $!\n";
	foreach my $child (sort {$a<=>$b} keys %trimap) {
		print $trifh ">$child\n$trimap{$child}\n" if ($tab eq "n");
		print $trifh "$child\t$trimap{$child}\n" if ($tab eq "y");
	}
	close ($trifh);
	return 1;
}

#----------------------------------------------------------------------------
# sti stuff now
# print_fasta_sti_tgs($fasta,$seqIDmap,$taxonomy,$out,$counts,$uc,$tab,$genes);
#----------------------------------------------------------------------------
sub print_fasta_sti_tgs {
	my ($fa,$seqIDmap,$taxonomy,$out,$seqid_int,$uc,$tab,$genes) = @_;
	
	#Deal with the tgs
	print STDERR "     Getting seqID <-> geneID correspondance to get the tgs\n" if (($genes) && ($genes ne "na"));
	my $gID = load_geneIDs($genes) if (($genes) && ($genes ne "na"));
	my $sti = $out.".sti";
	my $tgs = $out.".tgs";
	$out = $out."_taxonomer.fa";	
	my $skip;
	my %check = ();
	open(my $fh, "<", $fa) or confess "ERROR (sub print_fasta_sti_tgs): could not open to read $fa $!\n";
	open(my $stifh, ">", $sti) or confess "ERROR (sub print_fasta_sti_tgs): could not open to write $sti $!\n";
	open(my $fafh, ">", $out) or confess "ERROR (sub print_fasta_sti_tgs): could not open to write $out $!\n";
	open(my $tgsfh, ">", $tgs) or confess "ERROR (sub print_fasta_sti_tgs): could not open to write $tgs $!\n" if ($genes);
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
				my $lineage = join("$sep",@{$taxonomy->{$seqid}});
				print $fafh ">$seqid_int\t$seqid\t$lineage\t$desc\n";
				print $stifh ">$seqid_int\n$seqIDmap->{$seqid}\n" if ($tab eq "n");
				print $stifh "$seqid_int\t$seqIDmap->{$seqid}\n" if ($tab eq "y");
				if ($genes) {
					my $g;
					($gID->{$seqid})?($g = $gID->{$seqid}):($g = $seqid_int);
					print $tgsfh "$seqIDmap->{$seqid}\t$g\t$seqid_int\t$seqid\n";
				}
				$check{$seqid}=1;
				$skip = 0;
				$seqid_int++;
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
	close $tgsfh if ($genes);
	return 1;
}

