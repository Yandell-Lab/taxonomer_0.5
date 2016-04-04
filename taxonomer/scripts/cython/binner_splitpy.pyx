
cdef long int * db_ids = [0x1,0x2,0x4,0x8,0x10,0x20,0x40,0x80,0x100,0x200,0x400,0x800,0x1000,0x2000,0x4000,0x8000,0x10000,0x20000,0x40000,0x80000,0x100000,0x200000,0x400000,0x800000,0x1000000,0x2000000,0x4000000,0x8000000,0x10000000,0x20000000,0x40000000,0x80000000]
cdef int num_dbs = 32

def parse_line(char * line,int kc, long int db_id, long int collisions):
  cdef int dbs[32]
  for i in xrange(num_dbs):
    dbs[i] = 0
  cdef long int tmp_db = 0
  data = line.strip().split("\t")
  for dbinfo in data[1].split(";"):
    dbi = dbinfo.split("-")
    tmp_db = long(dbi[0])
    tmp_count = int(dbi[1])
    for i in xrange(num_dbs):
      if db_ids[i] & tmp_db != 0:
        dbs[i] += tmp_count

  cdef long int max_dbs = 0 #number representing bitflags for databases
  cdef int max_count = 0

  for i in xrange(num_dbs):
    if dbs[i] > max_count:
      max_dbs = 0
      max_dbs = max_dbs | db_ids[i]
      max_count = dbs[i]
    elif dbs[i] == max_count:
      max_dbs = max_dbs | db_ids[i]

  if (db_id & max_dbs > 0) and (max_count >= kc) and (max_dbs | collisions == collisions):
    return ">%s\t%s\t%d\n%s\n"%(data[0],data[1],max_count,data[2])
  elif db_id == 0 and max_count < kc:
    return ">%s\t%s\t%d\n%s\n"%(data[0],data[1],max_count,data[2])
  else:
    return None
