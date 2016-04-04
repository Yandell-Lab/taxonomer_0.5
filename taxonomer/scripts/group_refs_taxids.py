from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description="Group reference sequences by taxonomic info and output reorganized fasta file")
parser.add_argument("reference", type=str, help="reference fasta sequences")
parser.add_argument("sti", type=str, help=".sti file")
parser.add_argument("output", type=str, help="output fasta file")
parser.add_argument("-protein", help="flag to mark if protein sequences", action="store_true", default=False)
args = parser.parse_args()

cur_seqid = None
seqid_taxid = {} #{seqid} -> taxid
taxid_seqid = {} #{taxid} -> [seqids]
f = open(args.sti,'r')
for line in f:
    if line[0] == ">":
        cur_seqid = line[1:].strip()
    else:
        taxid = int(line.strip())
        seqid_taxid[cur_seqid] = taxid
        if taxid not in taxid_seqid:
            taxid_seqid[taxid] = []
        taxid_seqid[taxid].append(cur_seqid)
f.close()


ref_count = 0
reference_seqs = {} #taxid->sequence
handle = open(args.reference,'r')
for record in SeqIO.parse(handle, "fasta"):
    if ref_count % 10000 == 0:
        print ref_count
    ref_count += 1
    if record.id not in seqid_taxid:
        print "seqid not found in sti: %s"%(record.id)
        continue
    taxid = seqid_taxid[record.id]
    if taxid not in reference_seqs:
        reference_seqs[taxid] = []
    reference_seqs[taxid].append(str(record.seq))
handle.close()

print "output references"
out = open(args.output,'w')
for taxid,seq in reference_seqs.items():
    if args.protein == False:
        out.write(">%s\t%d\n%s\n" % (taxid_seqid[taxid][0],len(seq),"N".join(seq)) )
    else:
        out.write(">%s\t%d\n%s\n" % (taxid_seqid[taxid][0],len(seq),"$".join(seq)) )
out.close()
