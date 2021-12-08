import os,sys,gzip

dic = {}
for line in gzip.open('data/uniprot_sprot.fasta.gz', 'rt'):
    if line[0] == '>':
        acc = line.split('|')[1]
        if 'GN=' in line:
            gene = line.split('GN=')[1].split(' ')[0]
            dic[acc] = gene

l = ''
for line in open('data/GPCRDB_fasta_SCOP_7TM_348_hmmsearch_GPCRDB_numberings.txt', 'r'):
    acc = line.split('|')[1]
    if acc not in dic:
        os.system('wget https://www.uniprot.org/uniprot/'+acc+'.fasta')
        for line in open(acc+'.fasta', 'r'):
            if line[0] == '>':
                if 'GN=' in line:
                    gene = line.split('GN=')[1].split(' ')[0]
                    dic[acc] = gene
    l += str(acc) + '\t' + dic[acc] + '\t'
    l += line.replace(' ', '\t')

open('data/GPCRDB_numbering.tsv', 'w').write(l)
