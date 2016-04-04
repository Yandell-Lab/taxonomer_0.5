
def merge_counts(char * kmer_counts_1, unsigned long int id_1, char * kmer_counts_2, unsigned long int id_2, char * outfile):
  out = open(outfile,'w')

  f1 = open(kmer_counts_1,'r')
  f2 = open(kmer_counts_2,'r')
  line1 = f1.readline()
  line2 = f2.readline()

  cdef unsigned long long int kmer1 = 0
  cdef unsigned long long int kmer2 = 0
  kmer1_hex = ""
  kmer2_hex = ""
  cdef unsigned long int kmer1_id = 0
  cdef unsigned long int kmer2_id = 0

  cdef int flag = 1
  if line1 and line2:
    line1_data = line1.strip().split("\t")
    line2_data = line2.strip().split("\t")
    kmer1_hex = line1_data[0]
    kmer2_hex = line2_data[0]
    kmer1 = long(kmer1_hex,16)
    kmer2 = long(kmer2_hex,16)
    kmer1_id = id_1
    kmer2_id = id_2
    if id_1 == 0:
      kmer1_id = int(line1_data[1])
    kmer2_id = id_2
    if id_2 == 0:
      kmer2_id = int(line2_data[1])
  else:
    flag = 0

  while flag == 1:
    if kmer1 == kmer2:
      out.write("%s\t%d\n"%(kmer1_hex, kmer1_id | kmer2_id))
      line1 = f1.readline()
      line2 = f2.readline()
      if line1 and line2:
        line1_data = line1.strip().split("\t")
        line2_data = line2.strip().split("\t")
        kmer1_hex = line1_data[0]
        kmer1 = long(kmer1_hex,16)
        kmer2_hex = line2_data[0]
        kmer2 = long(kmer2_hex,16)
        kmer1_id = id_1
        kmer2_id = id_2
        if id_1 == 0:
          kmer1_id = int(line1_data[1])
        if id_2 == 0:
          kmer2_id = int(line2_data[1])
      else:
        flag = 0
        break
    elif kmer1 > kmer2:
      out.write("%s\t%d\n"%(kmer2_hex, kmer2_id))
      line2 = f2.readline()
      if line2:
        line2_data = line2.strip().split("\t")
        kmer2_hex = line2_data[0]
        kmer2 = long(kmer2_hex,16)
        kmer2_id = id_2
        if id_2 == 0:
          kmer2_id = int(line2_data[1])
      else:
        flag = 0
        break
    else:
      if id_1 != 0:
        out.write("%s\t%d\n"%(kmer1_hex, id_1))
      else:
        out.write("%s\t%d\n"%(kmer1_hex, kmer1_id))
      line1 = f1.readline()
      if line1:
        line1_data = line1.strip().split("\t")
        kmer1_hex = line1_data[0]
        kmer1 = long(kmer1_hex,16)
        kmer1_id = id_1
        if id_1 == 0:
          kmer1_id = int(line1_data[1])
      else:
        flag = 0
        break

  if line1:
    while line1:
      if id_1 != 0:
        out.write("%s\t%d\n"%(line1.strip().split("\t")[0],id_1))
      else:
        out.write("%s"%(line1))
      line1 = f1.readline()
  elif line2:
    while line2:
      if id_2 != 0:
        out.write("%s\t%d\n"%(line2.strip().split("\t")[0],id_2))
      else:
        out.write("%s"%(line2))
      line2 = f2.readline()
  out.close()
  f1.close()
  f2.close()
