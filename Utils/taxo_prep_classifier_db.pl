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

#flush buffer
$| = 1;
#-------------------------------------------------------------------------------
#-------------------------------- DESCRIPTION ----------------------------------
#-------------------------------------------------------------------------------
my $VERSION = "1.10";
my $SCRITPNAME = "taxo_prep_classifier_db.pl";
my $CHANGELOG = "
# Change log for $SCRITPNAME
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
#	- v1.10 = 3 Jan 2018
#            Add option -a to make each sequence its own node
#            RW with convention upercase, and make them global (no need to pass to sub)
#            bug fix uc, it was always doing it

# TO DO = detect if same seqIDs have different lineages, output warnings, and allow correction / chose 1 by deleting lines in a file, then rerun script with an option loading that file to correct lineages
\n";

my $usage = "\nUsage [v$VERSION]: 
    perl $SCRITPNAME 
         -f input.fasta [-o coredbname] [-n] [-m <X>] [-s ;]
         [-g map.tab] [-p] [-t] [-w] [-u] [-a] [-v] [-c] [-h]
    
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
     -a,--all     => (BOOL)   each sequence will also be a node (to call counts to the reference in the public version)
                              this only affects the tri file.
     -m,--mirr    => (INT)    To set the first seqID number
                              0 is reserved, so if 0 is set, it will be 1
                              Default: last taxID+1
                              Will be ignored if -a is set
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

#-------------------------------------------------------------------------------
#------------------------------- LOAD AND CHECK --------------------------------
#-------------------------------------------------------------------------------
my ($FASTA,$OUT,$TAB,$GENES,$POP,$SETID,$TEST,$WARN,$UC,$ALL,$V,$CHLOG,$HELP);
my $SEP = ";";
my $opt_success = GetOptions(
			 	  'fasta=s'   => \$FASTA,
			 	  'out=s'     => \$OUT,
			 	  'new'       => \$TAB,
			 	  'mirr=s'    => \$SETID,
			 	  'genes=s'   => \$GENES,
			 	  'sep=s'     => \$SEP,
			 	  'pop'       => \$POP,
			 	  'test'      => \$TEST,
			 	  'warn'      => \$WARN,
			 	  'uc'        => \$UC,
			 	  'all'       => \$ALL,
			 	  'chlog'     => \$CHLOG,
			 	  'version'   => \$V,
			 	  'help'      => \$HELP,);

#Check options, if files exist, etc
die "\n --- $SCRITPNAME version $VERSION\n\n" if $V;
die $CHANGELOG if ($CHLOG);
die "\n SOME MANDATORY ARGUMENTS MISSING, CHECK USAGE:\n$usage" if ($HELP || ! $FASTA);
die "\n -f $FASTA is not a fasta file?\n\n" unless ($FASTA =~ /\.fa|faa|fasta|fas$/);
die "\n -f $FASTA does not exist?\n\n" if (! -e $FASTA);
die "\n -m $SETID is not an integer!\n\n" if ($SETID && $SETID =~ /\D/);
die "\n -g $GENES does not exist?\n\n" if ($GENES && $GENES ne "na" && ! -e $GENES);
undef($SETID) if ($ALL);
($TAB)?($TAB="\t"):($TAB="\n");

#-------------------------------------------------------------------------------
#------------------------------------ MAIN -------------------------------------
#-------------------------------------------------------------------------------
print STDERR "\n --- $SCRITPNAME (v$VERSION)\n";
print STDERR "     input = $FASTA\n";
if (! $OUT) {
	$OUT = $FASTA;
	$OUT =~ s/\.fa$|\.faa$|\.fasta$|\.fas$//;
}
print STDERR "     output core name = $OUT\n";

print STDERR " --- Getting fasta IDs\n"; #could be all done at the same time, but that way taxIDs will be smaller numbers than the seqIDs
my @HEADERS = `grep "^>" $FASTA`;
my $nb_seqs = scalar(@HEADERS);
print STDERR "     => total of $nb_seqs sequence headers to process\n";
my $MAXLVL = 1;
if ($POP) {
print STDERR " --- Checking lineages for each seqIDs to check their depth\n";
	$MAXLVL = get_maxlvl();
	print STDERR "     => largest number of nodes in any lineage = $MAXLVL\n";
	print STDERR "        (separator for nodes = $SEP)\n";
}
print STDERR " --- Parsing headers to load lineages for each seqIDs \n";
print STDERR "     (separator for nodes = $SEP)\n" if (! $POP);
print STDERR "     (with filled empty nodes up to $MAXLVL = populating down)\n" if ($POP);
my $TAXONOMY = ();
my $INDEXMAP = ();
my $NB_DUP_ID = 0;
my $NB_NO_LIN = 0;
load_headers();
if (($NB_DUP_ID > 0) || ($NB_NO_LIN > 0)) {
	print STDERR " --- Parsing done, but:\n";
	print STDERR "     there were $NB_DUP_ID duplicated sequences (same id)\n" if ($NB_DUP_ID > 0);
	print STDERR "     there were $NB_NO_LIN sequences without lineage information\n" if ($NB_NO_LIN > 0);
	die "     => EXITING (--test chosen)\n\n" if ($TEST);
}
print STDERR " --- Obtaining tax id map (tri map) + taxids of each sequence\n";
my $TAXIDMAP = ();
my $SEQIDMAP = ();
my $TAXID = 1; #0 is the root
get_tri_sti();
$TAXID = $SETID if ($SETID);
$TAXID++ if $TAXID == 0;

print STDERR " --- Print .tri, _key.txt and _index.txt files\n";
print STDERR "     (tri in tabulated format)\n" if ($TAB eq "\t");
print STDERR "     (-a set, so tri also has each sequence as its own node in the tri)\n" if ($ALL);
print STDERR "     (WARN: -m was set, with is incompatible with -a, so was ignored)\n" if ($ALL && $SETID);
print_tri_key_idx();

print STDERR " --- Obtaining the sti map and writing it + rewriting the fasta file with integer IDs\n";
print STDERR "     (with lineages populated down)\n" if ($POP);
print STDERR "     (sti in tabulated format)\n" if ($TAB eq "\t");
print_fasta_sti_tgs();
print STDERR " --- DONE\n\n";
exit;


#-------------------------------------------------------------------------------
#--------------------------------- SUBROUTINES ---------------------------------
#-------------------------------------------------------------------------------
sub get_maxlvl {
	HEADER: foreach my $header (@HEADERS) {
		chomp($header);
		$header =~ s/^>//;
		my ($id,$lin) = split(/\s+/,$header);
		next HEADER if ((! $lin) || ($lin eq "nd")); #missing ones will be counted later
		my @nodes;
		($lin =~ /$SEP/)?(@nodes = split($SEP,$lin)):($nodes[0]=$lin);
		$MAXLVL = $#nodes+1 unless ($#nodes+1 < $MAXLVL);
	}	
	return 1;
}

#-------------------------------------------------------------------------------
sub load_headers{
	HEADER: foreach my $header (@HEADERS) {
		chomp($header);
		$header =~ s/^>//;
		my ($id,$lin) = split(/\s+/,$header); #description useless here, but could be shifted to be as the lineage, so check that as well:
		if ((! $lin) || ($lin eq "nd")) {
			$NB_NO_LIN++;
			print STDERR "        WARN: no lineage found for $id -> skipping sequence\n" if ($WARN);
			next HEADER;			
		} 
		my $poplin = ();
		$poplin = get_map($lin);
		if ($TAXONOMY->{$id}) {
			$NB_DUP_ID++;
			print STDERR "        WARN: duplicated entry for $id -> skipping sequence\n" if ($WARN);
			next HEADER;
		} else {
			$TAXONOMY->{$id}=$poplin; #$POPlin = array
		}
	}	
	return 1;
}

#-------------------------------------------------------------------------------
sub get_map {
	my $lin = shift;
	my @nodes = ();	
	($lin =~ /$SEP/)?(@nodes = split($SEP,$lin)):($nodes[0]=$lin);
	$MAXLVL = $#nodes +1 if (! $POP); #Won't populate down
	my $lineage_string = "";
	for (my $i = 0; $i < $MAXLVL; $i++){
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
		($lineage_string)?($lineage_string = $lineage_string.$SEP.$nodes[$i]):($lineage_string = $nodes[$i]); 
		$INDEXMAP->{$lineage_string}{'rank'} = $i;
		$INDEXMAP->{$lineage_string}{'label'} = $nodes[$i];
		$INDEXMAP->{$lineage_string}{'label'} = $1 if ($INDEXMAP->{$lineage_string}{'label'} =~ /^[0-9]+?__(.+?)$/);
	}
	return (\@nodes); #the @nodes array has the populated down lineages
}

#-------------------------------------------------------------------------------
sub get_tri_sti {
	my $prev_lineage = "";
	foreach my $seqid (keys %{$TAXONOMY}) {
		my $lineage = "";				
		foreach my $rank (@{$TAXONOMY->{$seqid}}) { #loop through all nodes
			($lineage eq "")?($lineage = $rank):($lineage = $lineage.$SEP.$rank);
			if (! $TAXIDMAP->{$lineage}{'own'}) { #no taxid yet
				$TAXIDMAP->{$lineage}{'own'}=$TAXID; #set as the previously incremented taxid
				if ((! $TAXIDMAP->{$prev_lineage}{'own'}) || ($INDEXMAP->{$lineage}{'rank'} == 0)) {
					$TAXIDMAP->{$lineage}{'parent'}=0; #when under the root
				} else {
					$TAXIDMAP->{$lineage}{'parent'}=$TAXIDMAP->{$prev_lineage}{'own'};
				}
				$TAXID++;
			}
			$prev_lineage = $lineage;	
			#each seqID, the last taxid (before the ++) is the one to remember; no possible duplicated seqid at this point 	
			$SEQIDMAP->{$seqid}=$TAXIDMAP->{$lineage}{'own'};   
		}		
	}
	return 1;
}

#-------------------------------------------------------------------------------
sub print_tri_key_idx {
	my $key = $OUT."_key.txt";
	my $idx = $OUT."_index.txt";	
	my %trimap = ();
	$trimap{1}=0;
	open(my $keyfh, ">", $key) or confess "ERROR (sub print_tri_key_idx): could not open to write $key $!\n";
	print $keyfh "#taxID\tparent_taxID\trank\tlineage_string\n";
	open(my $idxfh, ">", $idx) or confess "ERROR (sub print_tri_key_idx): could not open to write $idx $!\n";
	LIN: foreach my $lineage (sort keys %{$TAXIDMAP}) {
		next LIN unless ($lineage);
		#NB: tri does not need to be sorted, but that's 'prettier'
		$trimap{$TAXIDMAP->{$lineage}{'own'}}=$TAXIDMAP->{$lineage}{'parent'};
		my ($rank,$label) = ($INDEXMAP->{$lineage}{'rank'},$INDEXMAP->{$lineage}{'label'});
		print $keyfh "$TAXIDMAP->{$lineage}{'own'}\t$TAXIDMAP->{$lineage}{'parent'}\t$rank\t$lineage\n"; #taxid \t parent_taxid \t lineage_string		
		print $idxfh "$TAXIDMAP->{$lineage}{'own'}\t$rank\t$label\n"; #taxid \t rank \t label
	}
	close ($keyfh);
	close ($idxfh);
	
	#print a sorted tri
	my $tri = $OUT.".tri";
	open(my $trifh, ">", $tri) or confess "ERROR (sub print_tri_key_idx): could not open to write $tri $!\n";
	foreach my $child (sort {$a<=>$b} keys %trimap) {
		print $trifh ">".$child.$TAB.$trimap{$child}."\n";
	}
	close ($trifh);
	return 1;
}

#-------------------------------------------------------------------------------
sub print_fasta_sti_tgs {
	#Deal with the tgs
	print STDERR "     Getting seqID <-> geneID correspondance to get the tgs\n" if ($GENES && $GENES ne "na");
	my $gID = load_geneIDs() if ($GENES && $GENES ne "na");
	my $sti = $OUT.".sti";
	my $tri = $OUT.".tri";
	my $tgs = $OUT.".tgs";
	$OUT = $OUT."_taxonomer.fa";	
	my $skip;
	my %check = ();
	open(my $fh, "<", $FASTA) or confess "ERROR (sub print_fasta_sti_tgs): could not open to read $FASTA $!\n";
	open(my $stifh, ">", $sti) or confess "ERROR (sub print_fasta_sti_tgs): could not open to write $sti $!\n";
	open(my $fafh, ">", $OUT) or confess "ERROR (sub print_fasta_sti_tgs): could not open to write $OUT $!\n";
	open(my $tgsfh, ">", $tgs) or confess "ERROR (sub print_fasta_sti_tgs): could not open to write $tgs $!\n" if ($GENES);
	open(my $trifh, ">>", $tri) or confess "ERROR (sub print_tri_key_idx): could not open to write $tri $!\n" if ($ALL);
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
			if (! $lin || $lin eq "nd" || $check{$seqid}){
				$skip = 1;
			} else {	
				my $lineage = join("$SEP",@{$TAXONOMY->{$seqid}});
				#fasta file
				print $fafh ">$TAXID\t$seqid\t$lineage\t$desc\n";
				#sti file
				if ($ALL) {
					#then the sti file will be in tri and sti just has seqs as children of themselves
					print $stifh ">".$TAXID.$TAB.$TAXID."\n";
				} else {
					print $stifh ">".$TAXID.$TAB.$SEQIDMAP->{$seqid}."\n";
				}
				#Deal with tgs file if relevant
				if ($GENES) {
					my $g;
					($gID->{$seqid})?($g = $gID->{$seqid}):($g = $TAXID);
					print $tgsfh "$SEQIDMAP->{$seqid}\t$g\t$TAXID\t$seqid\n";
				}
				#now if -a option; add the sti in the tri 
				if ($ALL) {
					print $trifh ">".$TAXID.$TAB.$SEQIDMAP->{$seqid}."\n"; #what would be in sti otherwise
				}	
				$check{$seqid}=1;
				$skip = 0;
				$TAXID++;
			}
		} else {
			if ($UC) {
				print $fafh uc($line)."\n" unless ($skip == 1);
			} else {
				print $fafh $line."\n" unless ($skip == 1);
			}	
		}	
	}
	close $fh;
	close $stifh;
	close $fafh;
	close $tgsfh if ($GENES);
	close ($trifh) if ($ALL);
	return 1;
}

#-------------------------------------------------------------------------------
sub load_geneIDs {
	my %ids = ();
	open(my $fh, "<", $GENES) or confess "ERROR (sub load_geneIDs): could not open to read $GENES $!\n";
	while (<$fh>) {
		chomp(my $l = $_);
		next LINE if ((substr($l,0,1) eq "#") || ($l =~ /\W/));
		my @l = split('\s+',$l);
		$ids{$l[0]}=$l[1];
	}
	close ($fh);
	return(\%ids);	
}

