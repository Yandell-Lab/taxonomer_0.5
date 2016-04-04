import argparse

parser = argparse.ArgumentParser(description="Generate data for random forest based on taxonomer input from simulated reads")
parser.add_argument("taxonomer_input", type=str, help="taxonomer input file name")
parser.add_argument("sti_file", type=str, help="sti file")
parser.add_argument("tri_file", type=str, help="tri file")
parser.add_argument("output", type=str, help="output file")
args = parser.parse_args()

f = open(args.sti_file,'r')
cur_seqid = None
seqid_taxid = {} #seqid->taxid
for line in f:
  if line[0] == ">":
    cur_seqid = line[1:].strip()
  else:
    seqid_taxid[cur_seqid] = line.strip()
f.close()

f = open(args.tri_file,'r')
taxid_rel = {} #child->parent
cur_child = None
for line in f:
  if line[0] == ">":
    cur_child = line[1:].strip()
  else:
    taxid_rel[cur_child] = line.strip()
f.close()

def get_tax_path(taxid):
  path = {}
  cur_taxid = taxid
  while cur_taxid != "0":
    path[cur_taxid] = 0
    cur_taxid = taxid_rel[cur_taxid]
  return path

out = open(args.output,'w')
f = open(args.taxonomer_input,'r')
for line in f:
  if line[0] == "U":
    continue
  data = line.strip().split("\t")
  correct_taxid = seqid_taxid[data[1].split("_")[0]]
  correct_path = get_tax_path(correct_taxid)
  if data[2] in correct_path:
    out.write("1\t%s\n"%("\t".join(data[4:])))
  else:
    out.write("0\t%s\n"%("\t".join(data[4:])))
f.close()
out.close()
