#!/usr/bin/python3

from collections import defaultdict

from Bio import SeqIO
import sys

if len(sys.argv) > 4:
    print("Incorrect number of arguments\nUsage: ./pa.py reads_file transcriptome_file k\nUsage: ./pa.py".format(len(sys.argv) - 1))
    exit()
elif len(sys.argv) == 4:
    reads_file = sys.argv[1]
    transcriptome_file = sys.argv[2]
    k = int(sys.argv[3])
    if k < 21 or k > 36:
        print("Please input a k-mer length between 21-36")
        exit()
elif len(sys.argv) == 2 or len(sys.argv) == 3:
    print("Incorrect number of arguments\nUsage: ./pa.py reads_file transcriptome_file k\nUsage: ./pa.py")
    exit()
else:
    reads_file = "reads.fasta"
    transcriptome_file = "chr11_transcriptome.fasta"
    k = 21

reads = []
# for reads that have base N, skip the base N
for record in SeqIO.parse(reads_file, "fasta"):
    seq = record.seq.replace('N', '')
    reads.append(seq)

# given gene annotation, create index
isoforms = {}
for record in SeqIO.parse(transcriptome_file, "fasta"):
    isoforms[record.id] = record.seq

kmer_map = defaultdict(lambda: set())
for iso in isoforms:
    seq = isoforms[iso]
    for i in range(len(seq)-k+1):
        kmer = seq[i:i+k]
        kmer_map[kmer].add(iso)

def gen_eq_classes(read_seq):
    eq_classes = None
    for i in range(len(read_seq)-k+1):
        kmer = read_seq[i:i+k]
        if kmer not in kmer_map:
            return None
        elif eq_classes is None:
            eq_classes = kmer_map[kmer]
        else:
            eq_classes = kmer_map[kmer] & eq_classes
        if not eq_classes: # if eq_classes is empty, then return (will be empty set)
            break
    return eq_classes

def get_rev_comp(read_seq):
    rev_comp = ''
    for base in read_seq[::-1]:
        if base == 'A':
            c = 'T'
        elif base == 'T':
            c = 'A'
        elif base == 'C':
            c = 'G'
        elif base == 'G':
            c = 'C'
        else:
            c = 'N' # there are reads that have base N, we just treat these reads as bad reads backwards - their eq class will be empty
        rev_comp += c
    return rev_comp

reads_eq_classes = []
for read_seq in reads:
    read_seq_rc = get_rev_comp(read_seq)
    eq_class = gen_eq_classes(read_seq)
    eq_class_rev = gen_eq_classes(read_seq_rc)

    if eq_class is None and eq_class_rev is None:
        reads_eq_classes.append(None)
    elif eq_class is None:
        reads_eq_classes.append(eq_class_rev)
    elif eq_class_rev is None:
        reads_eq_classes.append(eq_class)
    else:
        reads_eq_classes.append(eq_class | eq_class_rev)

eq_class_labels_dict = defaultdict(lambda: 0)
for read_eq in reads_eq_classes:
    if read_eq is None:
        eq_class_labels_dict[tuple('')] += 1
        continue
    eq_classes = sorted(read_eq)
    eq_class_labels_dict[tuple(eq_classes)] += 1

sorted_labels = sorted(eq_class_labels_dict.items(), key=lambda kv: kv[0])
sorted_labels = sorted(sorted_labels, key=lambda x: len(x[0]))

f = open("output.tsv", "w")
f.write("counts\tnumber of items in equivalence class\tisoforms in equivalence class\n")
for classes, counts in sorted_labels:
    label = 'NA' if len(classes) == 0 else ','.join(classes)
    f.write("{}\t{}\t{}\n".format(counts, len(classes), label))
f.close()
