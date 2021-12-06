#!/home/marin/.conda/envs/transformer/bin/python
# coding: utf-8

## Libraries to import
import os, sys, time, datetime, argparse, numpy
import pandas as pd
from sklearn.preprocessing import MinMaxScaler
from joblib import dump, load
## User-defined functions to import
sys.path.insert(0, 'data_precog2/')
import load_functions
import analysis
timestamp1 = time.time()

## Set the path to folders
hmm_models_path = 'data_precog2/hmm_models/'
selected_features_path = 'data_precog2/selected_features/'
output_path = 'data_precog2/output/'
feature_files_path = 'data_precog2/feature_files/'
start_ali = 164
end_ali = 735

##########################
## To Parse the input arguments
##########################

path = os.getcwd()
parser = argparse.ArgumentParser(description='PRECOG2: PREdicting COupling probabilities Of G-protein coupled receptors\n(trained on data from Avet et al, 2020)', epilog='Contact: gurdeep.singh[at]bioquant.uni-heidelberg.de')
parser.add_argument('input', help='Path to input file (FASTA format)')
parser.add_argument('uniq_id', help='Unique ID')
parser.add_argument('--hmm', help='Path to input file\'s HMMsearch o/p(against 7tm_1). If absent, this script will generate one for itself using the default settings of HMM')
parser.add_argument('--o', help='Path & name of output file. If absent, the output will be printed on the screen')
parser.add_argument('--hack', help='Path to the DIRECTORY where the hack o/p is to be written. Cannot be empty. Refer to README for more details')
args = parser.parse_args()

fasta_file = args.input
uniq_id = args.uniq_id
if fasta_file[0] != '/':
	fasta_file = path + '/' + args.input
hmm_file = args.hmm
if hmm_file != None:
	if hmm_file[0] != '/':
		hmm_file = path + '/' + hmm_file
out_file = args.o
if out_file != None:
	if out_file[0] != '/':
		out_file = path + '/' + out_file
hack_directory = args.hack
if hack_directory != None:
	if hack_directory[0] != '/':
		hack_directory = path + '/' + hack_directory
##########################

################################
## Base class for all sequences
################################
class base:

	def __init__(self, name, mut):
		self.name = name
		self.seq = ''
		self.til = ''
		self.til_start = 0
		self.til_end = 0
		self.ctl = ''
		self.ctl_start = 0
		self.mutation = {}
		self.mutation['WT'] = {}
		self.mutation['WT']['position'] = {}
		self.add_mut(mut)

	def add_mut(self, mut):
		if mut != '':
			if (mut in self.mutation) == False:
				self.mutation[mut] = {}
				self.mutation[mut]['position'] = {}

	def add_fasta(self, seq):
		self.seq += seq

	def show(self):
		print(self.name)
		print(self.seq)
		print(self.mutation)
		print(self.til)
		print(len(self.til.replace('-', '')))
		print(self.til_start)
		print(self.til_end)
		print(self.ctl)
		print(self.ctl_start)

#################################
## Read the input FASTA sequences
#################################
def read_fasta(fasta_file):
	obj = {}
	for line in open(fasta_file, 'r'):
		if line[0] == '#':
			continue
		if line[0] != '\n' or line[0] != '':
			if line[0] == '>':
				mutation = ''
				name = (line.split('>')[1].replace('\n', '').split())[0]

				if '/' in name:
					mutation = name.split('/')[1]
					name = name.split('/')[0]
					if mutation[-1] == 'a' or mutation[-1] == 'p' or mutation[-1] == 'X' or mutation[0].isalpha() == False or mutation[-1].isalpha() == False or mutation[1:-1].isdigit() == False:
						go = 0
						name = ''
						continue
				else:
					mutation = ''

				if (name in obj) == False:
					obj[name] = base(name, mutation)
					go = 1
				else:
					obj[name].add_mut(mutation)
					go = 0

			else:
				if go == 1:
					obj[name].add_fasta(line.replace(' ', '').replace('\n', ''))
	return obj

obj = read_fasta(fasta_file)
#################################

###############################################
## Read the hmmsearch o/p of input against 7tm1
## or generate one if not provided by the user
###############################################

def extract_other_features(file):
	v = 0; w = 0; proteins = []
	for line in open(file, 'r'):
		if line[0:2] == '>>':
			name = (line.split('>>')[1].replace('\n', '').split())[0]
			if '/' in name:
				name = name.split('/')[0]

			if (name in obj) == True:
				if name not in proteins:
					proteins.append(name)
					proteins = list(set(proteins))
					v = 1
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
				obj[name].ctl = obj[name].seq[int(y[-1]):].upper()
				obj[name].ctl_start = int(y[-1]) + 1
				w = 1

		if w == 1:
			start = int(x[1])
			count = 0

			for k, (i, j) in enumerate(zip(x[2], y[2])):
				xseq = x[2][:k].replace('-', '').replace('.', '')
				if len(xseq) + start >= 172 and len(xseq) + start <= 204:
					obj[name].til += j.upper()
				if len(xseq) + start >= 172 and obj[name].til_start == 0:
					obj[name].til_start = int(y[1]) + len((y[2])[:k].replace('-', '').replace('.', ''))
				if len(xseq) + start >= 172 and len(xseq) + start <= 204 and obj[name].til_end >= 0:
					if len((y[2])[:k].replace('-', '').replace('.', '')) > 0:
						obj[name].til_end = int(y[1]) + len((y[2])[:k+1].replace('-', '').replace('.', '')) - 1
						'''
						if name == 'OR7G2':
							print '----'
							print y[1], len((y[2])[:k].replace('-', '').replace('.', ''))
							print obj[name].til_start, obj[name].til_end, k, count+start
							print (y[2])[:k+1].replace('-', '').replace('.', '')
						'''
				if i!='-' and i!='.':
					count += 1

			w = 0

if os.path.exists(path+'/output/'+uniq_id+'/temp') == False:
    os.system('mkdir '+path+'/output/'+uniq_id+'/temp')

temp_path = path+'/output/'+uniq_id+'/temp/'

## If the file is not provided by user, generate one itself
if hmm_file == None:
	hmm_file = load_functions.hmm_search(obj, temp_path)

extract_other_features(hmm_file)
#print obj['hM3D'].til_start
#print obj['hM3D'].til_end
###########################################################


###################################
## Generate HMM search o/p of input
## against our Gproteins models
###################################

def spinning_cursor():
    while True:
        for cursor in '|/-\\':
            yield cursor

## Create a new fasta file containing only those genes
## that are present in both the fasta and hmmsearch
def new_fasta():
	l = ''
	for name in obj:
		l += '>' + str(name) + '\n' + str(obj[name].seq) + '\n'
	#if os.path.exists('/net/netfile2/ag-russell/bq_gsingh/gpcr/update/aln_psiblast/temp') == False:
	open(temp_path+'/new_fasta_file.txt', 'w').write(l)

new_fasta()
spinner = spinning_cursor()


for files in os.listdir(hmm_models_path):
	if files.endswith('.hmm'):
		#os.system('/net/home.isilon/ag-russell/install/CentOS-7.2.1511-x86_64/bin/hmmsearch '+hmm_models_path+files+' '+'temp/new_fasta_file.txt'+' > temp/'+files.split('.hmm')[0]+'.out')
		os.system('hmmsearch '+hmm_models_path+files+' '+temp_path+'/new_fasta_file.txt'+' > '+temp_path+'/'+files.split('.hmm')[0]+'.out')
		sys.stdout.write(next(spinner))
		sys.stdout.flush()
		time.sleep(0.01)
		sys.stdout.write('\b')

###########################################################

#####################################
## Functions for feature construction
#####################################

def construct_features(name, pos, neg, hmm_pos, hmm_neg, features, row, *mut):
	#df = pd.read_table('/net/netfile2/ag-russell/bq_gsingh/gpcr/update/aln_psiblast/selected_features.txt', sep = '\t', index_col = 0)
	df = pd.read_table(selected_features_path+'/selected_features_all.txt', sep = '\t', index_col = 0)
	map_position = {}
	for x, y in df[['MSA_Pos', 'Domain_Pos']].to_numpy().tolist():
		map_position[str(x)] = str(y.replace('(', '|').replace(')', ''))
	#print map_position

	mutation = ''
	for arg in mut:
		mutation = arg
	#print mutation
	mutation_position = -1
	if mutation != 'WT':
		mutation_position = int(mutation[1:-1])
	#print mutation_position
	#print features
	dom_position = {}

	for f in features:
		remarks = []

		if 'pos' in f:
			if (str(f[:-3]) in pos[name]) == True:
				if pos[name][str(f[:-3])]['position'] == mutation_position:
					AA = mutation[-1]
					dom_position[map_position[str(f[:-3])]] = f
				else:
					AA = pos[name][str(f[:-3])]['aa']

				if AA != '-':
					row.append(1)
				else:
					row.append(0)
			else:
				row.append(0)

			#if f in ['966pos', '967pos'] and mutation in ['T655M', 'WT']:
			#	print f, row[-1], mutation, AA

		elif 'neg' in f:
			if (str(f[:-3]) in neg[name]) == True:
				if neg[name][str(f[:-3])]['position'] == mutation_position:
					AA = mutation[-1]
					dom_position[map_position[str(f[:-3])]] = f
				else:
					AA = neg[name][str(f[:-3])]['aa']

				if AA != '-':
					row.append(1)
				else:
					row.append(0)
			else:
				row.append(0)

		if 'bip' in f:
			if (str(f[:-3]) in pos[name]) == True:
				if pos[name][str(f[:-3])]['position'] == mutation_position:
					AA = mutation[-1]
					dom_position[map_position[str(f[:-3])]] = f
				else:
					AA = pos[name][str(f[:-3])]['aa']

				if AA != '-':
					row.append(float(hmm_pos[str(f[:-3])][AA.upper()]))
				else:
					row.append(float(load_functions.find_max_bits(hmm_pos, f)))
			else:
				row.append(float(load_functions.find_max_bits(hmm_pos, f)))

		elif 'bin' in f:
			if (str(f[:-3]) in neg[name]) == True:
				if neg[name][str(f[:-3])]['position'] == mutation_position:
					AA = mutation[-1]
					dom_position[map_position[str(f[:-3])]] = f
				else:
					AA = neg[name][str(f[:-3])]['aa']

				if AA != '-':
					row.append(float(hmm_neg[str(f[:-3])][AA.upper()]))
				else:
					row.append(float(load_functions.find_max_bits(hmm_neg, f)))
			else:
				row.append(float(load_functions.find_max_bits(hmm_neg, f)))

		elif 'TILL' in f:
			if mutation_position >= obj[name].til_start and mutation_position <= obj[name].til_end:
				sequence = ''
				if obj[name].til != '':
					sequence = list(obj[name].til.replace('-', ''))
					#print sequence, mutation, obj[name].til_start, obj[name].til_end, name, 'til'
					sequence[mutation_position - obj[name].til_start] = mutation[-1]
					sequence = "".join(sequence)
				dom_position['TILL'] = f
			else:
				sequence = obj[name].til
			row.append(len(sequence))

		elif 'CTLL' in f:
			if mutation_position >= obj[name].ctl_start:
				sequence = ''
				if obj[name].ctl != '':
					sequence = list(obj[name].ctl.replace('-', ''))
					#print sequence, mutation, obj[name].ctl_start, name, 'ctl'
					sequence[mutation_position - obj[name].ctl_start] = mutation[-1]
					sequence = "".join(sequence)
				dom_position['CTLL'] = f
			else:
				sequence = obj[name].ctl
			row.append(len(sequence))

		elif '_CTL' in f:
			if mutation_position >= obj[name].ctl_start:
				sequence = ''
				if obj[name].ctl != '':
					sequence = list(obj[name].ctl.replace('-', ''))
					sequence[mutation_position - obj[name].ctl_start] = mutation[-1]
					sequence = "".join(sequence)
					dom_position[f] = f
			else:
				sequence = obj[name].ctl
			row.append(sequence.count(f.split('_')[0]))

		elif '_TIL' in f:
			if mutation_position >= obj[name].til_start and mutation_position <= obj[name].til_end:
				sequence = ''
				if obj[name].til != '':
					sequence = list(obj[name].til.replace('-', ''))
					#print sequence, mutation, obj[name].til_start, obj[name].til_end, name, 'til'
					sequence[mutation_position - obj[name].til_start] = mutation[-1]
					sequence = "".join(sequence)
				dom_position[f] = f
			else:
				sequence = obj[name].til
			row.append(sequence.count(f.split('_')[0]))

		#if mutation != 'WT':
			#print f, remarks
			#sys.exit()
	#print mutation, row
	return row, dom_position

def read_aln(pos, neg, hmm_pos, hmm_neg, l, features, gprotein):
	data = []
	for name in obj:
		if (name in pos) == True and (name in neg) == True:
			for mut in obj[name].mutation:
				if mut == 'WT' or int(mut[1:-1]) <= len(obj[name].seq):
					row = []
					row.append(str(name))
					#row.append(str(mut))
					row, dom_position = construct_features(name, pos, neg, hmm_pos, hmm_neg, features, row, mut)

					for pi in dom_position:
						if (pi in obj[name].mutation[mut]['position']) == False:
							obj[name].mutation[mut]['position'][pi] = {}
							obj[name].mutation[mut]['position'][pi][dom_position[pi]] = []
							obj[name].mutation[mut]['position'][pi][dom_position[pi]].append(str(gprotein))
						else:
							if (str(dom_position[pi]) in obj[name].mutation[mut]['position'][pi]) == False:
								obj[name].mutation[mut]['position'][pi][dom_position[pi]] = []
								obj[name].mutation[mut]['position'][pi][dom_position[pi]].append(str(gprotein))
							else:
								obj[name].mutation[mut]['position'][pi][dom_position[pi]].append(str(gprotein))
								obj[name].mutation[mut]['position'][pi][dom_position[pi]] = list(set(obj[name].mutation[mut]['position'][pi][dom_position[pi]]))
					data.append(row)
	#print data[3]
	#print data[1]
	return numpy.array(data)

def read_gprotein_hmm_out(file, hmm):
	v = 0; w = 0; gpcr = {}
	for line in open(file, 'r'):
		if line[:6] == 'Query:':
                        file_name = line.split()[1]
		if line[0:2] == '>>':
			name = (line.split('>> ')[1].replace('\n', '').split())[0]
			if '/' in name:
				name = name.split('/')[0]
			if (name in obj) == True:
				gpcr[name] = {}
				v = 1
			else:
				name = ''
				v = 0
			continue
		if v == 1:
			#file_name = file.split('.out')[0].split('/')[-1]
			#file_name = file_name.split('_')[0] + '_subali_' + file_name.split('_')[1]
			if file_name in line[:25]:
				x = line.split()
			elif name in line[:25]:
				y = line.split()
				if y[1] != '-' and y[-1] != '-':
					w = 1
		if w == 1:
			start = int(x[1])
			count = 0
			#print x
			#print y
			for k, (i, j) in enumerate(zip(x[2], y[2])):
				if i!='-' and i!='.':
					if (str(count+start) in hmm) == True:
						#if 'GNA12' in file and 'pos' in file:
							#print count+start,
						if int(hmm[str(count+start)]) >= start_ali and int(hmm[str(count+start)]) <= end_ali:
							#if 'GNA12' in file and 'pos' in file:
							#	print hmm[str(count+start)]
							gpcr[name][str(hmm[str(count+start)])] = {}
							gpcr[name][str(hmm[str(count+start)])]['aa'] = j
							if j == '-':
								gpcr[name][str(hmm[str(count+start)])]['position'] = '-'
							else:
								gpcr[name][str(hmm[str(count+start)])]['position'] = int(y[1]) + len((y[2])[:k].replace('-', '').replace('.', ''))
					count += 1
				#print k+int(y[1]),
			#print
			w = 0

	return gpcr

def extract_features(gprotein):
	features = []
	for line in open(feature_files_path+gprotein+'_train.txt', 'r'):
		if line.split('\t')[0] == 'GPCR':
                    features = line.split('\t')[1:-1]
	return features

def extract_model(gprotein):
	for files in os.listdir(output_path):
		if gprotein == files.split('_')[0] and 'model' in files:
			model = load(output_path+files)
			#print files
			break
	return model

def k_fold(file):
	df = pd.read_table(file, lineterminator = '\n', sep = '\t')
	col = list(df.columns.values)
	df[col[1:-1]] = df[col[1:-1]].astype(float)
	min_max_scaler_all = MinMaxScaler()
	min_max_scaler_all.fit_transform(df[col[1:-1]])
	return min_max_scaler_all

#print 'Making predictions for:'
gpcr_list = list(obj.keys())
prediction = {}
gprotein_list = ['GNAI3', 'GNAI1', 'GNAZ', 'GNAO1', 'GNA12', 'GNA13', 'GNAQ', 'GNA14', 'GNA15', 'GNAS', 'GNAL']
gprotein_list = ['Gs','Gi1','Gi2','GoA','GoB','Gz','G12','G13','Gq','G11','G14','G15']
#gprotein_list = ['GNAI1']
dic_ebret = {'G13':'GNA13', 'G12':'GNA12', 'Gs':'GNAS', 'Gq':'GNAQ', 'G11':'GNA11', 'G14':'GNA14', 'G15':'GNA15', 'Gi1':'GNAI1', 'Gi2':'GNAI2', 'Gz':'GNAZ', 'GoA': 'GoA', 'GoB': 'GoB'}
for gprotein in gprotein_list:
    features = extract_features(gprotein)

    hmm_pos, hmm_pos_positions = load_functions.read_hmm(hmm_models_path+gprotein+'_pos.hmm')
    hmm_neg, hmm_neg_positions = load_functions.read_hmm(hmm_models_path+gprotein+'_neg.hmm')
    pos = read_gprotein_hmm_out(temp_path+gprotein+'_pos.out', hmm_pos_positions)
    neg = read_gprotein_hmm_out(temp_path+gprotein+'_neg.out', hmm_neg_positions)
    #print pos['TSHR']['393']
    #sys.exit()
    if hack_directory != None:
        #if gprotein == 'GNA12':
        print('Running HACK option for', gprotein)
        l = analysis.main(pos, neg, hmm_pos, hmm_neg, features, gprotein, list(obj.keys()), obj)
        open(hack_directory+'/'+str(gprotein)+'.txt', 'w').write(l)
        #sys.exit()

    l= 'GPCR'
    for f in features:
        l+='\t' + f
    l+= '\n'
    data = read_aln(pos, neg, hmm_pos, hmm_neg, l, features, gprotein)
    for row in data:
        l += '\t'.join(row) + '\n'

    gprotein = dic_ebret[gprotein]
    open (path+'/output/'+uniq_id+'/ebret/seq_features/'+gprotein+'.txt', 'w').write(l)

sys.exit()
