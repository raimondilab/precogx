import os, sys, gzip

ACC2G = {}
for line in open('../../data/GtP_to_UniProt_mapping.tsv', 'r'):
    l = line.replace('"', '')
    if l[0] != '#' and l.split('\t')[0] != 'UniProtKB ID' and l.split('\t')[1] == 'Human':
        gtopdbName = l.split('\t')[2]
        acc = l.split('\t')[0]
        import re
        # as per recommendation from @freylis, compile once only
        CLEANR = re.compile('<.*?>')
        gtopdbName = re.sub(CLEANR, '', gtopdbName)
        ACC2G[acc] = gtopdbName

ACC2GN = {}
GN2ACC = {}
ACC2ID = {}
ID2GN = {}
GN2ID = {}
for line in gzip.open('../../data/uniprot_sprot.fasta.gz', 'rt'):
    if line[0] == '>':
        acc = line.split('|')[1]
        id = line.split('|')[2].split(' ')[0]
        ACC2ID[acc] = id
        if 'GN=' in line:
            gene = line.split('GN=')[1].split(' ')[0]
            GN2ACC[gene] = acc
            ACC2GN[acc] = gene
            ID2GN[id] = gene
            GN2ID[gene] = id

l = ''
for line in open('../predictor/data_precog/LogRAi_values_final.tsv', 'r'):
    if line[0] == '#':
        header = line.replace('\n', '').split('\t')[1:]
        l += '#GENE\tACC\tID\tNAME\t' + '\t'.join(header) + '\n'
    else:
        gene = line.split('\t')[0]
        values = line.replace('\n', '').split('\t')[1:]
        acc = GN2ACC[gene]
        name = ACC2G[acc] if acc in ACC2G else '-'
        #print (values)
        l +=  gene + '\t' + acc + '\t' + GN2ID[gene] + '\t' + name + '\t' + '\t'.join(values) + '\n'

open('../../data/shedding.tsv', 'w').write(l)

l = ''
for line in open('../predictor/data_precog2/dnor_Emax.tsv', 'r'):
    if line[0] == '#':
        header = line.replace('\n', '').split('\t')[5:]
        l += '#GENE\tACC\tID\tNAME\t' + '\t'.join(header) + '\n'
    else:
        gene = line.split('\t')[2]
        values = line.replace('\n', '').split('\t')[5:]
        #print (values)
        if gene == 'CALCR+RAMP3':
            name = ACC2G['P30988-2'] if 'P30988-2' in ACC2G else '-'
            l +=  'CALCR' + '\t' + 'P30988-2' + '\t' + GN2ID['CALCR'] + '\t' + name + '\t' + '\t'.join(values) + '\n'
            name = ACC2G['O60896'] if 'O60896' in ACC2G else '-'
            l +=  'RAMP3' + '\t' + 'O60896' + '\t' + GN2ID['RAMP3'] + '\t' + name + '\t' + '\t'.join(values) + '\n'
        else:
            acc = line.split('\t')[1]
            name = ACC2G[acc] if acc in ACC2G else '-'
            l +=  gene + '\t' + acc + '\t' + GN2ID[gene] + '\t' + name + '\t' + '\t'.join(values) + '\n'

open('../../data/ebbret.tsv', 'w').write(l)

l = ''
for line in open('../../../precog2/data/IUPHAR_couplings_Marin.tsv', 'r'):
    if line[0] == '#':
        header = line.replace('\n', '').split('\t')[5:]
        l += '#GENE\tACC\tID\tNAME\tPC\tSC\n'
    else:
        gene = line.split('\t')[1]
        if gene == 'TAAR4P':
            gene = 'TAAR6'
        if gene != '':
            acc = GN2ACC[gene]
            PC = line.split('\t')[4].replace(' family', '').replace('"', '')
            SC = line.split('\t')[5].replace(' family', '').replace('"', '')
            name = ACC2G[acc] if acc in ACC2G else '-'
            #print (gene, PC, SC)
            l +=  gene + '\t' + acc + '\t' + GN2ID[gene] + '\t' + name + '\t' + PC + '\t' + SC + '\n'

open('../../data/iuphar.tsv', 'w').write(l)

l = '#GENE\tACC\tID\tNAME\tBarr1\tBarr2\n'
for line in open('../../data/Barr_test_set_STRING_rob_IMEx_HIPPIE.txt', 'r'):
    gene = line.split('\t')[0]
    acc = GN2ACC[gene]
    values = line.replace('\n', '').split('\t')[1:]
    name = ACC2G[acc] if acc in ACC2G else '-'
    l +=  gene + '\t' + acc + '\t' + GN2ID[gene] + '\t' + name + '\t '+ '\t'.join(values) + '\n'

open('../../data/string.tsv', 'w').write(l)
