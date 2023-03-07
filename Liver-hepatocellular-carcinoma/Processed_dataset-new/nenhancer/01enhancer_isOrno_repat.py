from Bio import SeqIO
import pysam

enhancers = []
for record in SeqIO.parse("/data/lyli/Liver/NoEnhancer/test_nenseq634_300.fa", "fasta"):
    name = record.id
    enhancers.append(name)

enhancers1484 = []
for record in SeqIO.parse("/data/lyli/Liver/NoEnhancer/NoEnhancer1484_300.fa", "fasta"):
    name = record.id
    enhancers1484.append(name)

enhancer_del = []
for enhancer in enhancers:
    for i in enhancers1484:
        if enhancer == i:
            # enhancer_del.append(enhancer)
            print(enhancer)

# for i in enhancer_del:
#     enhancers.remove(i)
# 
# newfile = open('/data/lyli/Liver/NoEnhancer/test_nenseq299_300-.fa', 'w')
# for j in enhancers:
#     chromosome = j.split(":")[0]
#     start_position = j.split(":")[1].split("-")[0]
#     end_position = j.split("-")[1].split("\n")[0]
#     start_position = int(start_position)
#     end_position = int(end_position)
#     ref = '/data/lyli/Hg19/tolhg19.fa'
#     with pysam.FastaFile(ref) as fa_obj:
#         seq = fa_obj.fetch(chromosome, start_position, end_position)
# 
#     newfile.write(">" + chromosome + ":" + str(start_position) + "-" + str(end_position) + "\n" + seq + "\n")
# 
# newfile.close()