import os
from Bio import SeqIO
from Bio import Entrez


Entrez.email = "cmzmasek@yahoo.com"
#ids = ['MW738412', 'MW738362', 'MW696434', 'MW594041', 'MW644278', 'MW725704', 'MW596208', 'MW709028', 'MT834682',
#       'MT846481', 'MT846486', 'MW521724', 'MW338785', 'MW528761', 'MW694473'
#       ]


filename = 'out.txt'

ids = []

with open('ext_node_ids_reorg.txt') as f:
    for line in f:
        s = line.split(',')
        for ss in s:
            sss = ss.strip()
            if len(sss) > 0:
                ids.append(sss)

# ids = ['MT334544'] <=== multiple
print(ids)
print(str(len(ids)))

# Downloading...
out_handle = open(filename, "w")
c = 0
for myid in ids:

    print(str(c))
    print(myid)
    c += 1
    net_handle = Entrez.efetch(
        db="nucleotide", id=myid, rettype="record", retmode="text"
    )
    result = net_handle.read()
    #print(result)
    for line in result.splitlines():
        # print(line)
        if 'title "surface glycoprotein' in line:
            print(line)

    # records = Entrez.read(net_handle)
    # print(str(records))
    # xx = records['Bioseq-set_seq-set']

    # seq desc title "surface glycoprotein [Severe acute respiratory syndrome coronavirus 2]"
    # inst
    # print(str(xx))
# print(records['Seq-entry'])
# for record in xx:
#    print(str(record))
#    xxx = xx['Seq-entry_set']
# print(str(record))
# y = json.loads(records)
# each record is a Python dictionary or list.


# print(net_handle.read())
# out_handle.write(net_handle.read())

out_handle.close()
net_handle.close()
print("Saved")

# print("Parsing...")
# record = SeqIO.read(filename, "record")
# print(record)
