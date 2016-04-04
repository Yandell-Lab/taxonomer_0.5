import os
path = os.path.dirname(os.path.abspath(__file__))
import sys
sys.path.append("%s/cython/"%(path))
import argparse
import socket
import classify_READpy
import binner_PY

parser = argparse.ArgumentParser(description="Create query server on a particular socket for a single database, for INTERNAL use")
parser.add_argument("port", type=int, help="port number to listen at")
parser.add_argument("db_prefix", type=str, help="full path for database prefix")
args = parser.parse_args()

host=""
backlog = 5
s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
s.bind((host,args.port))
s.listen(backlog)

query_object = classify_READpy.classify_READpy()
query_object.load_binner_db(1, args.db_prefix+".bmi", args.db_prefix+".btbi")
client, address = s.accept()
data_size = 1024

while 1:
    #cr.classify_reads_stream(client.makefile())
    read_data = client.recv(data_size)
    data = read_data.strip().split("\t")
    bin = query_object.bin_read_stream(data[0], len(data[0]), data[1])
    client.send("%d"%(bin))



