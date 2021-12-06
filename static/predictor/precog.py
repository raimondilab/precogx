#!/home/marin/.conda/envs/transformer/bin/python

# coding: utf-8
import os, sys, time, datetime, argparse, numpy
path = os.getcwd()
import pandas as pd
from sklearn.preprocessing import MinMaxScaler
from joblib import dump, load

sys.path.insert(0, path+'/data_precog')
import load_functions
import analysis
timestamp1 = time.time()

##########################
## To Parse the arguments
##########################

parser = argparse.ArgumentParser(description='PRECOG (PREdicting COupling probabilities Of G-protein coupled receptors)', epilog='End of help. Contact: gurdeep.singh@bioquant.uni-heidelberg.de')
parser.add_argument('fasta_file', help='path to input file (FASTA formatted); see data/sample.fasta')
parser.add_argument('uniq_id', help='Unique ID')
parser.add_argument('--hmm', help='path to hmmsearch o/p of the input file against 7tm1; if absent, this script will generate one for itself using default settings of HMM')
parser.add_argument('--o', help='path to the output file; if absent, the output will be printed on the screen')
parser.add_argument('--hack', help='path to the directory where the hack output should be stored')
args = parser.parse_args()
fasta_file = args.fasta_file
uniq_id = args.uniq_id
if fasta_file[0] != '/':
	fasta_file = path + '/' + args.fasta_file
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
	hmm_file = load_functions.hmm_search(obj, path, temp_path)

extract_other_features(hmm_file)
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
	open(temp_path+'/new_fasta_file.txt', 'w').write(l)

new_fasta()
spinner = spinning_cursor()

os.chdir(path+'/data_precog/hmm_models/')

for files in os.listdir('.'):
	if files.endswith('.hmm'):
		os.system('hmmsearch '+files+' '+temp_path+'/new_fasta_file.txt'+' >' +temp_path+ '/'+files.split('.hmm')[0]+'.out')
		sys.stdout.write(next(spinner))
		sys.stdout.flush()
		time.sleep(0.01)
		sys.stdout.write('\b')

###########################################################

#####################################
## Functions for feature construction
#####################################

def construct_features(name, pos, neg, hmm_pos, hmm_neg, features, row, *mut):
	df = pd.read_csv(path+'/data_precog/selected_features.txt', sep = '\t', index_col = 0)
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
			file_name = file.split('.out')[0].split('/')[-1]
			file_name = file_name.split('_')[0] + '_subali_' + file_name.split('_')[1]
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
						if int(hmm[str(count+start)]) >= 393 and int(hmm[str(count+start)]) <= 1002:
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
	for line in open(path+'/data_precog/feature_files/'+gprotein+'_train.txt', 'r'):
		if line.split('\t')[0] == 'GPCR':
                    features = line.split('\t')[1:-1]
	return features
#print 'Making predictions for:'
gpcr_list = list(obj.keys())
prediction = {}
gprotein_list = ['GNAI3', 'GNAI1', 'GNAZ', 'GNAO1', 'GNA12', 'GNA13', 'GNAQ', 'GNA14', 'GNA15', 'GNAS', 'GNAL']
for gprotein in gprotein_list:
    features = extract_features(gprotein)

    hmm_pos, hmm_pos_positions = load_functions.read_hmm(path+'/data_precog/hmm_models/'+gprotein+'_pos.hmm')
    hmm_neg, hmm_neg_positions = load_functions.read_hmm(path+'/data_precog/hmm_models/'+gprotein+'_neg.hmm')

    pos = read_gprotein_hmm_out(temp_path+'/'+gprotein+'_pos.out', hmm_pos_positions)
    neg = read_gprotein_hmm_out(temp_path+'/'+gprotein+'_neg.out', hmm_neg_positions)

    if hack_directory != None:
        l = analysis.main(path, pos, neg, hmm_pos, hmm_neg, features, gprotein, list(obj.keys()), obj)
        open(hack_directory+'/'+str(gprotein)+'.txt', 'w').write(l)
            
    l= 'GPCR'
    for f in features:
        l += '\t' + f
            
    l += '\n'
            
    data = read_aln(pos, neg, hmm_pos, hmm_neg, l, features, gprotein)
    for row in data:
        l += '\t'.join(row) + '\n'
    #print (l)
    open (path+'/output/'+uniq_id+'/shed/seq_features/'+gprotein+'.txt', 'w').write(l)

sys.exit()

