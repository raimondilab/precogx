import os, sys, gzip

ACC2GN = {}
GN2ACC = {}
ACC2ID = {}
for line in gzip.open('data/uniprot_sprot.fasta.gz', 'rt'):
    if line[0] == '>':
        acc = line.split('|')[1]
        id = line.split('|')[2].split(' ')[0]
        ACC2ID[acc] = id
        if 'GN=' in line:
            gene = line.split('GN=')[1].split(' ')[0]
            GN2ACC[gene] = acc
            ACC2GN[acc] = gene

l = ''
for line in open('static/predictor/data_precog/LogRAi_values_final.tsv', 'r'):
    if line[0] == '#':
        header = line.replace('\n', '').split('\t')[1:]
        l += '#GENE\tACC\t' + '\t'.join(header) + '\n'
    else:
        gene = line.split('\t')[0]
        values = line.replace('\n', '').split('\t')[1:]
        #print (values)
        l +=  gene + '\t' + GN2ACC[gene] + '\t' + '\t'.join(values) + '\n'

open('data/shedding.tsv', 'w').write(l)

l = ''
for line in open('static/predictor/data_precog2/dnor_Emax.tsv', 'r'):
    if line[0] == '#':
        header = line.replace('\n', '').split('\t')[5:]
        l += '#GENE\tACC\t' + '\t'.join(header) + '\n'
    else:
        gene = line.split('\t')[2]
        acc = line.split('\t')[1]
        values = line.replace('\n', '').split('\t')[5:]
        #print (values)
        l +=  gene + '\t' + acc + '\t' + '\t'.join(values) + '\n'

open('data/ebbret.tsv', 'w').write(l)

l = ''
for line in open('../precog2/data/IUPHAR_couplings_Marin.tsv', 'r'):
    if line[0] == '#':
        header = line.replace('\n', '').split('\t')[5:]
        l += '#GENE\tACC\tPC\tSC\n'
    else:
        gene = line.split('\t')[1]
        if gene == 'TAAR4P':
            gene = 'TAAR6'
        if gene != '':
            PC = line.split('\t')[4].replace(' family', '').replace('"', '')
            SC = line.split('\t')[5].replace(' family', '').replace('"', '')
            #print (gene, PC, SC)
            l +=  gene + '\t' + GN2ACC[gene] + '\t' + PC + '\t' + SC + '\n'

open('data/iuphar.tsv', 'w').write(l)

l = '#GENE\tACC\tBarr1\tBarr2\n'
for line in open('data/Barr_test_set_STRING_rob_IMEx_HIPPIE.txt', 'r'):
    gene = line.split('\t')[0]
    values = line.replace('\n', '').split('\t')[1:]
    l +=  gene + '\t' + GN2ACC[gene] + '\t' + '\t'.join(values) + '\n'

open('data/string.tsv', 'w').write(l)
