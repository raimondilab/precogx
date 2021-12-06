import pandas as pd
import numpy, gzip, sys, os
import metrics

pfam_path = 'data/7tm_1.hmm'

'''
##### Create dictionary that contains GPCRs and their STD
def load_aska_std():
	arr = []
	for line in open('/net/home.isilon/ag-russell/bq_gsingh/gpcr/systematic_balanced/LogRAI_v2.txt', 'r'):
		if "SEM" in line.split('\t')[1]:
			for val in line.replace('\n', '').split('\t')[2:]:
				if float(val) != 0.0:
					arr.append(1/float(val))

	mx = max(arr)

	val = numpy.array(val)
	val = val.astype(float, order = 'C')

	dic = {}
	for line in open('/net/home.isilon/ag-russell/bq_gsingh/gpcr/systematic_balanced/LogRAI_v2.txt', 'r'):
		#print line.split('\t')
		if line[0] == '#':
			gprotein_list = line.replace('\n', '').split('\t')[2:]
		elif "Mean" in line.split('\t')[1]:
			gpcr = line.split('\t')[0]
			if '/' in gpcr:
				gpcr = gpcr.split('/')[0]
		elif "SEM" in line.split('\t')[1]:
			dic[gpcr]={}
			for gprotein, val in zip(gprotein_list, line.replace('\n', '').split('\t')[2:]):
				if float(val) != 0.0:
					dic[gpcr][gprotein] = 1/float(val)
				else:
					dic[gpcr][gprotein] = mx

	return dic

##### Create dictionary thats maps Accession to Gene names
def load_uniprot(*control):
	atg={}
	for line in gzip.open('/net/home.isilon/ag-russell/bq_fraimondi/DBs/Uniprot/Uniprot_Genenames.tsv.gz', 'r'):
		if 'HUMAN' in line:
			if atg.has_key(line.split('\t')[0]) == False:
				if len(control) == 0:
					atg[line.split('\t')[0]] = line.split('\t')[2].replace('\n', '').upper()
				elif line.split('\t')[2].replace('\n', '').upper() in control[0]:
					atg[line.split('\t')[0]] = line.split('\t')[2].replace('\n', '').upper()
	print 'Total accessions loaded:', len(atg)
	print atg
	return atg

##### Load Aska's Dataset
def load_aska():
	data = []
	for line in open('/net/home.isilon/ag-russell/bq_fraimondi/PROJECTs/GPCR/GPCR-Coupling/experimental_data/LogRAi_values_final.tsv', 'r'):
		if '#' in line:
			col_name = []
			col_name.append('GPCR')
			col_name += line.replace('\n', '').replace('\r', '').split('\t')[1:]
			col_name = numpy.array(col_name)
			#print col_name
			#sys.exit()
		else:
			row = line.replace('\n', '').replace('\r', '').split()
			data.append(row)
	df = pd.DataFrame(data, columns = col_name)
	return df
'''

def load_bouvier(value):
    data = []
    for line in open('/net/home.isilon/ag-russell/bq_gsingh/gpcr/systematic_balanced/'+value+'.tsv', 'r'):
        if '#' in line:
            col_name = []
            col_name.append('GPCR')
            col_name += line.replace('\n', '').replace('\r', '').replace('\n', '').split('\t')[5:]
            col_name = numpy.array(col_name)
            #print col_name
            #sys.exit()
        else:
            row = []
            row.append(line.split('\t')[2])
            row += line.replace('\n', '').replace('\r', '').replace('\n', '').split('\t')[5:]
            data.append(row)
    df = pd.DataFrame(data, columns = col_name)
    return df


##### Check Aska's Dataset
def check_aska(name, gprotein_list):
	data = []
	for line in open('data/LogRAi_values_final.tsv', 'r'):
		if '#' in line:
			col_name = []
			col_name.append('GPCR')
			col_name += line.replace('\n', '').replace('\r', '').split('\t')[1:]
			col_name = numpy.array(col_name)
		else:
			row = line.replace('\n', '').replace('\r', '').split()
			data.append(row)
	df = pd.DataFrame(data, columns = col_name)
	l = str(name)
	while len(l) <= 20:
		l += ' '
	#print name, df['GPCR'].tolist()
	if name in df['GPCR'].tolist():
		for gprotein in gprotein_list:
			l += '\t' + str(round(float(df[df['GPCR'] == name][gprotein].tolist()[0]), 2))
	else:
		for gprotein in gprotein_list:
			l += '\t-'
	l += '\tASKA_LogRIA\n'

	return l

##### Check Aska's Dataset
def check_bouvier(name, gprotein_list):
	data = []
	for line in open('data/dnor_Emax.tsv', 'r'):
		if '#' in line:
			col_name = []
			col_name.append('GPCR')
			col_name += line.replace('\n', '').replace('\r', '').split('\t')[5:]
			col_name = numpy.array(col_name)
		elif line.split('\t')[3] == 'Class A':
			row = [line.replace('\n', '').replace('\r', '').split('\t')[2]] + line.replace('\n', '').replace('\r', '').split('\t')[5:]
			#print (row)
			data.append(row)
	df = pd.DataFrame(data, columns = col_name)
	l = str(name)
	while len(l) <= 20:
		l += ' '
	#print name, df['GPCR'].tolist()
	if name in df['GPCR'].tolist():
		for gprotein in gprotein_list:
			l += '\t' + str(round(float(df[df['GPCR'] == name][gprotein].tolist()[0]), 2))
	else:
		for gprotein in gprotein_list:
			l += '\t-'
	l += '\tBRET-Assay(dnor_Emax)\n'

	return l

'''
##### Load FASTA
def read_fasta(base, atg):
	obj = {}
	#for line in open('/net/home.isilon/ag-russell/bq_gsingh/gpcr/PF00001_noORs_unipid.fa', 'r'):
	for line in open('/net/home.isilon/ag-russell/bq_gsingh/gpcr/class_A_gpcr.fasta', 'r'):
		if line[0] != '\n' or line[0] != '':
			if line[0] == '>':
				name = line.split('|')[1]
				if atg.has_key(name) == True:
					name = atg[name]

					if obj.has_key(name) == False:
						obj[name] = base(name)
						go = 1
					else:
						go = 0
				else:
					go = 0

			else:
				if go == 1:
					obj[name].add_fasta(line.replace(' ', '').replace('\n', ''))

	print 'Fasta sequences loaded:', len(obj)
	return obj

##### Load FASTA
def load_fasta(atg, *control):
	fasta={}; go = 0

	#for line in open('/net/netfile2/ag-russell/bq_gsingh/uniprot_sprot.fasta', 'r'):
	#for line in open('/net/home.isilon/ag-russell/bq_gsingh/gpcr/PF00001_noORs_unipid.fa', 'r'):
	for line in open('/net/home.isilon/ag-russell/bq_gsingh/gpcr/class_A_gpcr.fasta', 'r'):

		#if line[0] == '>' and 'Homo sapiens' in line:
		if line[0] == '>':
			if atg.has_key(line.split(' ')[0].replace('>', '').split('|')[1].split('|')[0]) == True:
				name = atg[line.split(' ')[0].replace('>', '').split('|')[1].split('|')[0]]
			else:
				continue

			if fasta.has_key(name) == False:
				if len(control) != 0:
					if name in control[0]:
						fasta[name] = ''
						go = 1
				else:
					fasta[name] = ''
					go = 1

			else:
				go = 0

		#elif line[0] == '>' and 'Homo sapiens' not in line:
		#	go = 0

		else:
			if go == 1:
				fasta[name] += line.replace('\n', '')

	print 'Fasta sequences loaded:', len(fasta)
	return fasta


#### Extract other features like TILL, CTLL, KCTL, RCTL, PCTL
def extract_other_features(f, obj, base, atg):
	v = 0; w = 0; proteins = []
	for line in open(f, 'r'):
		#print line
		if line[0:2] == '>>':
			name = line.split('>>')[1].split(' ')[1]
			gene = line.split('>>')[1].replace('\n', '').split('|')[1]
			if atg.has_key(gene) == True:
				gene = atg[gene]
				if obj.has_key(gene) == True:
					if gene not in proteins:
						proteins.append(gene)
						proteins = list(set(proteins))
						v = 1
					else:
						v = 0
				else:
					v = 0
			else:
				v = 0
			continue

		if v == 1:
			if '7tm_1' in line[:25]:
				x = line.split()
			elif name in line[:25]:
				y = line.split()
				obj[gene].ctl = obj[gene].seq[int(y[-1]):].upper()
				obj[gene].ctl_start = int(y[-1]) + 1
				w = 1
				if name == 'GPR132':
					print x
					print y

		if w == 1:
			start = int(x[1])
			count = 0

			for k, (i, j) in enumerate(zip(x[2], y[2])):
				xseq = x[2][:k].replace('-', '').replace('.', '')
				if len(xseq) + start >= 172 and len(xseq) + start <= 204:
					obj[gene].til += j.upper()
				if len(xseq) + start >= 172 and obj[gene].til_start == 0:
					obj[gene].til_start = int(y[1]) + len((y[2])[:k].replace('-', '').replace('.', ''))
				if len(xseq) + start >= 172 and len(xseq) + start <= 204 and obj[gene].til_end >= 0:
					if len((y[2])[:k].replace('-', '').replace('.', '')) > 0:
						obj[gene].til_end = int(y[1]) + len((y[2])[:k+1].replace('-', '').replace('.', '')) - 1

						##if gene == 'GPR132':
						##	print '----'
						##	print y[1], len((y[2])[:k].replace('-', '').replace('.', ''))
						##	print obj[gene].til_start, obj[gene].til_end, k, count+start
						##	print (y[2])[:k+1].replace('-', '').replace('.', '')

				if i!='-' and i!='.':
					count += 1

			w = 0

	print 'Read hmmsearch against 7tm1'
	return obj

#### Extract other features like FILL, SILL
def extract_other_features_extra(file, obj, base, atg):
	v = 0; w = 0; proteins = []
	for line in open(file, 'r'):
		if line[0:2] == '>>':
			name = line.split('>>')[1].split(' ')[1]
			gene = line.split('>>')[1].replace('\n', '').split('|')[1]
			if atg.has_key(gene) == True:
				gene = atg[gene]
				if obj.has_key(gene) == True:
					if gene not in proteins:
						proteins.append(gene)
						proteins = list(set(proteins))
						v = 1
					else:
						v = 0
				else:
					v = 0
			else:
				v = 0
			continue
			continue

		if v == 1:
			if '7tm_1' in line[:25]:
				x = line.split()
			elif name in line[:25]:
				y = line.split()
				if obj[gene].ntl_start == -1:
					obj[gene].ntl = obj[gene].seq[int(y[1]):].upper()
					obj[gene].ntl_start = int(y[1]) - 1
					#print gene, obj[gene].ntl_start
				obj[gene].ctl = obj[gene].seq[int(y[-1]):].upper()
				obj[gene].ctl_start = int(y[-1]) + 1
				w = 1
				if name == 'GPR132':
					print x
					print y

		if w == 1:
			start = int(x[1])
			count = 0
			#start_loop = 120 #48 #94 #12 #172
			#end_loop = 142 #51 #96 #17 #204

			for k, (i, j) in enumerate(zip(x[2], y[2])):
				xseq = x[2][:k].replace('-', '').replace('.', '')

				##TIL
				start_loop = 172
				end_loop = 204
				if len(xseq) + start >= start_loop and len(xseq) + start <= end_loop:
					obj[gene].til += j.upper()
				if len(xseq) + start >= start_loop and obj[gene].til_start == 0:
					obj[gene].til_start = int(y[1]) + len((y[2])[:k].replace('-', '').replace('.', ''))
				if len(xseq) + start >= start_loop and len(xseq) + start <= end_loop and obj[gene].til_end >= 0:
					if len((y[2])[:k].replace('-', '').replace('.', '')) > 0:
						obj[gene].til_end = int(y[1]) + len((y[2])[:k+1].replace('-', '').replace('.', '')) - 1

				##FIL
				start_loop = 12
				end_loop = 17
				if len(xseq) + start >= start_loop and len(xseq) + start <= end_loop:
					obj[gene].fil += j.upper()
				if len(xseq) + start >= start_loop and obj[gene].fil_start == 0:
					obj[gene].fil_start = int(y[1]) + len((y[2])[:k].replace('-', '').replace('.', ''))
				if len(xseq) + start >= start_loop and len(xseq) + start <= end_loop and obj[gene].fil_end >= 0:
					if len((y[2])[:k].replace('-', '').replace('.', '')) > 0:
						obj[gene].fil_end = int(y[1]) + len((y[2])[:k+1].replace('-', '').replace('.', '')) - 1

				##SIL
				start_loop = 94
				end_loop = 96
				if len(xseq) + start >= start_loop and len(xseq) + start <= end_loop:
					obj[gene].sil += j.upper()
				if len(xseq) + start >= start_loop and obj[gene].sil_start == 0:
					obj[gene].sil_start = int(y[1]) + len((y[2])[:k].replace('-', '').replace('.', ''))
				if len(xseq) + start >= start_loop and len(xseq) + start <= end_loop and obj[gene].sil_end >= 0:
					if len((y[2])[:k].replace('-', '').replace('.', '')) > 0:
						obj[gene].sil_end = int(y[1]) + len((y[2])[:k+1].replace('-', '').replace('.', '')) - 1

				##FEL
				start_loop = 48
				end_loop = 51
				if len(xseq) + start >= start_loop and len(xseq) + start <= end_loop:
					obj[gene].fel += j.upper()
				if len(xseq) + start >= start_loop and obj[gene].fel_start == 0:
					obj[gene].fel_start = int(y[1]) + len((y[2])[:k].replace('-', '').replace('.', ''))
				if len(xseq) + start >= start_loop and len(xseq) + start <= end_loop and obj[gene].fel_end >= 0:
					if len((y[2])[:k].replace('-', '').replace('.', '')) > 0:
						obj[gene].fel_end = int(y[1]) + len((y[2])[:k+1].replace('-', '').replace('.', '')) - 1

				##SEL
				start_loop = 120
				end_loop = 142
				if len(xseq) + start >= start_loop and len(xseq) + start <= end_loop:
					obj[gene].sel += j.upper()
				if len(xseq) + start >= start_loop and obj[gene].sel_start == 0:
					obj[gene].sel_start = int(y[1]) + len((y[2])[:k].replace('-', '').replace('.', ''))
				if len(xseq) + start >= start_loop and len(xseq) + start <= end_loop and obj[gene].sel_end >= 0:
					if len((y[2])[:k].replace('-', '').replace('.', '')) > 0:
						obj[gene].sel_end = int(y[1]) + len((y[2])[:k+1].replace('-', '').replace('.', '')) - 1

				if i!='-' and i!='.':
					count += 1

				##TEL
				start_loop = 233
				end_loop = 240
				if len(xseq) + start >= start_loop and len(xseq) + start <= end_loop:
					obj[gene].tel += j.upper()
				if len(xseq) + start >= start_loop and obj[gene].tel_start == 0:
					obj[gene].tel_start = int(y[1]) + len((y[2])[:k].replace('-', '').replace('.', ''))
				if len(xseq) + start >= start_loop and len(xseq) + start <= end_loop and obj[gene].tel_end >= 0:
					if len((y[2])[:k].replace('-', '').replace('.', '')) > 0:
						obj[gene].tel_end = int(y[1]) + len((y[2])[:k+1].replace('-', '').replace('.', '')) - 1

				if i!='-' and i!='.':
					count += 1

			w = 0

	print 'Read hmmsearch against 7tm1'
	return obj
'''

#### Read HMMs
def read_hmm(file):
	x={}; y={}; AA=[]
	for line in open(file, 'r'):
		if 'HMM' in line and 'HMMER' not in line:
			AA = line[13:].split()
		elif '- - -' in line:
			x[line.split()[-5]] = {}
			for i, bit in enumerate(numpy.array(line.split()[1:-5])):
				#x[line.split()[-5]][AA[i]] = 1/float(bit)		## Alignment position with AAs bit scores
				x[line.split()[-5]][AA[i]] = float(bit)		## Alignment position with AAs bit scores
			y[line.split()[0]] = line.split()[-5]	## HMM position to alignment position

	return x, y

'''
#### Read HMMs
def read_hmm_stack(file):
	x={}; y={}; AA=[]
	for line in open(file, 'r'):
		if 'HMM' in line and 'HMMER' not in line:
			AA = line[13:].split()
		elif '- - -' in line:
			x[line.split()[-5]] = {}
			for i, bit in enumerate(numpy.array(line.split()[1:-5])):
				x[line.split()[-5]][AA[i]] = 1/float(bit)		## Alignment position with AAs bit scores
				#x[line.split()[-5]][AA[i]] = float(bit)		## Alignment position with AAs bit scores
			y[line.split()[0]] = line.split()[-5]	## HMM position to alignment position

	return x, y
'''
#### Find max_bits at HMM
def find_max_bits(hmm, f):

	mx = 0
	for aa in hmm[str(f[:-3])]:
		if float(hmm[str(f[:-3])][aa]) > mx:
			mx = float(hmm[str(f[:-3])][aa])

	##mn = 10
	##for aa in hmm[str(f[:-3])]:
	##	if float(hmm[str(f[:-3])][aa]) < mn:
	##		mn = float(hmm[str(f[:-3])][aa])
	##mn = 0

	return mx

def load_iuphar(file):
	d = {}
	for line in open(file, 'r'):
		if line[0] != '#':
			pc = line.split('\t')[-4]
			sc = line.split('\t')[-3]
			gpcr = line.split('\t')[1]

			if d.has_key(gpcr) == False:
				d[gpcr] = {}
				d[gpcr]['pc'] = []
				d[gpcr]['sc'] = []

			pc = pc.replace('family', '').replace(' ', '')
			sc = sc.replace('family', '').replace(' ', '')

			for item in pc.split(','):
				if item != '':
					d[gpcr]['pc'].append(item.split('family')[0].replace(' ', ''))

			for item in sc.split(','):
				if item != '':
					d[gpcr]['sc'].append(item.split('family')[0].replace(' ', ''))

			d[gpcr]['pc'] = list(set(d[gpcr]['pc']))
			d[gpcr]['sc'] = list(set(d[gpcr]['sc']))

	return d


def check_iuphar(name, gprotein_list, d):
	l = str(name)
	while len(l) <= 20:
		l += ' '
	dic = {'Gi1':'Gi/Go', 'Gi2':'Gi/Go', 'GoA':'Gi/Go', 'GoB':'Gi/Go', 'Gz':'Gi/Go', 'Gs':'Gs', 'Gq':'Gq/G11', 'G14':'Gq/G11', 'G15':'Gq/G11', 'G11':'Gq/G11', 'G12':'G12/G13', 'G13':'G12/G13'}
	for gprotein in gprotein_list:
		if d.has_key(name) == True:
			if gprotein in dic:
				if dic[gprotein] in d[name]['pc']:
					l += '\tPC'
				elif dic[gprotein] in d[name]['sc']:
					l += '\tSC'
				else:
					l += '\t-'
			else:
				l += '\t-'
		else:
			l += '\t-'
	l += '\tIUPHAR\n'
	return l

def load_mut_info(file):
	dic = {}
	for line in open(file, 'r'):
		if 'HUMAN' in line.split('\t')[1]:
			if line.split('\t')[2] != '-':
				mut = line.split('\t')[2] + '/' +line.split('\t')[0].split('/')[1]
				if dic.has_key(mut) == False:
					#print mut
					dic[mut] = ' '.join(line.split()[6:])
	return dic

def hmm_search(obj, temp_path):
	l = ''
	for name in obj:
		l += '>' + str(name) + '\n'
		l += str(obj[name].seq) + '\n'
	open(temp_path+'/temp_fasta_file.txt', 'w').write(l)
	os.system('hmmsearch data_precog2/7tm_1.hmm '+temp_path+'/temp_fasta_file.txt > '+temp_path+'/temp_hmm_file.txt')
	return (temp_path+'/temp_hmm_file.txt')

'''
def testing(df, min_max_scaler_all, model):
	col = list(df.columns.values)
	names = df[col[0]].as_matrix()
	X = df[col[1:-1]].as_matrix()
	Y = df[col[-1]].as_matrix()
	if Y.tolist().count(1) != 0 and Y.tolist().count(0) != 0:
		#print Y.tolist().count(1)
		#print Y.tolist().count(0)
		X = min_max_scaler_all.transform(X)
		tp, fp, fn, tn, error = check(Y, model.predict(X))
		mcc, acc, pre, rec, spe, f1m, auc = metrics.main(tp, fp, fn, tn, Y, model.predict(X))
		#print Y
		#print model.predict(X)
		#for name, y, y_prob, y_pred in zip(names, Y, model.predict_proba(X), model.predict(X)):
		#	if name in ['GPR4', 'GPR39', 'S1PR4']:
		#		print name, y, y_prob[1], y_pred
		print str(mcc) +'\t'+ str(acc)+'\t'+ str(pre)+'\t'+ str(rec)+'\t'+ str(spe)+'\t'+ str(f1m)+'\t'+ str(auc) +'\t' + str(Y.tolist().count(1)) + '\t' + str(Y.tolist().count(0))
		#print rec
	else:
		print

def testing_b2(df, min_max_scaler_all, model, weight):
        #print weight
        col = list(df.columns.values)
        X = df[col[1:-1]].as_matrix()
        Y = df[col[-1]].as_matrix()
        if Y.tolist().count(1) != 0 and Y.tolist().count(0) != 0:
                #print Y.tolist().count(1)
                #print Y.tolist().count(0)
                X = min_max_scaler_all.transform(X)
                Y_pred_proba = model.predict_proba(X)[:,1]
                dummy = Y_pred_proba
                dummy[dummy >= weight] = 1
                dummy[dummy < weight] = 0
                Y_pred = dummy
                tp, fp, fn, tn, error = check(Y, Y_pred)
                mcc, acc, pre, rec, spe, f1m, auc = metrics.main(tp, fp, fn, tn, Y, Y_pred)
                #print len(X), len(Y)
                print str(mcc) +'\t'+ str(acc)+'\t'+ str(pre)+'\t'+ str(rec)+'\t'+ str(spe)+'\t'+ str(f1m)+'\t'+ str(auc) +'\t' + str(Y.tolist().count(1)) + '\t' + str(Y.tolist().count(0))
                #print rec
        else:
                print

def check(expected, predicted):
	error = []
	tp=0; tn=0; fp=0; fn=0
	for i in range(0, len(predicted)):
		if expected[i]==predicted[i]:
			error.append(1)
			if expected[i] == 1:
				tp+=1
			else:
				tn+=1
		else:
			error.append(0)
			if expected[i] == 1:
				fn+=1
			else:
				fp+=1

	return tp, fp, fn, tn, error
'''
