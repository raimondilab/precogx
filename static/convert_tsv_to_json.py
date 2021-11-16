## this is a makeshhift script to convert TSV output of PRECOGx into JSON
import json

data = []
dic = {}
for line in open('OL820//out.tsv', 'r'):
    row = []
    if line[0] != '#':
        row = line.replace('\n', '').split('\t')
        data.append(row)

dic = {'data': data}
with open('OL820/out.json', 'w') as f:
    json.dump(dic, f)
