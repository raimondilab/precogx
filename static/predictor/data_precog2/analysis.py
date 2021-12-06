import os, sys, time, datetime
import pandas as pd
import numpy
from sklearn.preprocessing import MinMaxScaler
from joblib import dump, load
sys.path.insert(0, 'data/')
import load_functions

def construct_features(name, pos, neg, hmm_pos, hmm_neg, mutation_position, features, row, obj, mut):
	#df = pd.read_table('/net/netfile2/ag-russell/bq_gsingh/gpcr/update/aln_psiblast/selected_features.txt', sep = '\t', index_col = 0)
	#sprint os.getcwd()
	df = pd.read_table('data/selected_features/selected_features_all.txt', sep = '\t', index_col = 0)
	map_position = {}
	for x, y in df[['MSA_Pos', 'Domain_Pos']].as_matrix().tolist():
		map_position[str(x)] = str(y.replace('(', '|').replace(')', ''))
	#print name, mut, mutation_position

	mutation = mut
	dom_position = {}

	for f in features:
		remarks = []

		if 'pos' in f:
			if pos[name].has_key(str(f[:-3])) == True:
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

		elif 'neg' in f:
			if neg[name].has_key(str(f[:-3])) == True:
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
			if pos[name].has_key(str(f[:-3])) == True:
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
			if neg[name].has_key(str(f[:-3])) == True:
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
					#print 'ctl_sequence', sequence, mutation, mutation_position, obj[name].ctl_start, name
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

	#sys.exit()
	return row, dom_position

def read_aln(pos, neg, hmm_pos, hmm_neg, features, gprotein, gpcr, obj):
	data = []
	done_features = []
	AA = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
	if pos.has_key(gpcr) == True and neg.has_key(gpcr) == True:
		for feature in features:
			if 'bip' in feature or 'bin' in feature:
				#print feature
				if feature[:-3] not in done_features:
					if pos[gpcr].has_key(feature[:-3]) == True and neg[gpcr].has_key(feature[:-3]) == True:
						#print feature
						mutation_position = pos[gpcr][feature[:-3]]['position']
						#print gpcr, gprotein, mutation_position, feature, pos[gpcr][feature[:-3]]
						if mutation_position != '-':
							for mut in AA:
								row = []
								row.append(str(gpcr))
								row.append(feature[:-3])
								row.append(str(mutation_position))
								row.append(str(pos[gpcr][feature[:-3]]['aa']))
								row.append(str(mut))
								row, dom_position = construct_features(gpcr, pos, neg, hmm_pos, hmm_neg, mutation_position, features, row, obj, mut)
								'''
								for pi in dom_position:
									if obj[gpcr].mutation[mut]['position'].has_key(pi) == False:
										obj[gpcr].mutation[mut]['position'][pi] = {}
										obj[gpcr].mutation[mut]['position'][pi][dom_position[pi]] = []
										obj[gpcr].mutation[mut]['position'][pi][dom_position[pi]].append(str(gprotein))
									else:
										if obj[gpcr].mutation[mut]['position'][pi].has_key(str(dom_position[pi])) == False:
											obj[gpcr].mutation[mut]['position'][pi][dom_position[pi]] = []
											obj[gpcr].mutation[mut]['position'][pi][dom_position[pi]].append(str(gprotein))
										else:
											obj[gpcr].mutation[mut]['position'][pi][dom_position[pi]].append(str(gprotein))
											obj[gpcr].mutation[mut]['position'][pi][dom_position[pi]] = list(set(obj[gpcr].mutation[mut]['position'][pi][dom_position[pi]]))
								'''
								data.append(row)
					done_features.append(feature[:-3])
	#for row in data:
	#	print row
	#sys.exit()
	#print numpy.array(data[:, 2:])
	return numpy.array(data)

def extract_model(gprotein):
	#for files in os.listdir('/net/netfile2/ag-russell/bq_gsingh/gpcr/update_2/output_VI/'):
	for files in os.listdir('data/output/'):
		if gprotein == files.split('_')[0] and 'model' in files:
			model = load('data/output/'+files)
			break
	return model

def k_fold(file):
	df = pd.read_table(file, lineterminator = '\n', sep = '\t')
	col = list(df.columns.values)
	df[col[1:-1]] = df[col[1:-1]].astype(float)
	min_max_scaler_all = MinMaxScaler()
	min_max_scaler_all.fit_transform(df[col[1:-1]])
	return min_max_scaler_all

def main(pos, neg, hmm_pos, hmm_neg, features, gprotein, gpcr_list, obj):
	l = '#PRECOG2 with Beta-Arrs(v1.0, Jul 2020)\n'
	l += '#'+str(gprotein)+' HACK output\n'
	l += '#gurdeep[dot]singh[at]bioquant[dot]uni-heidelberg[dot]de\n'
	now = datetime.datetime.now()
	l += '#Date: ' + str(now.day) + '-' + str(now.month) + '-' + str(now.year) + '\n'
	l += '#Time: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second) + '\n'
	l += '#Column GPCR: All the input GPCRs\n'
	l += '#Column ALN_POS: Position in MSA used for training the model\n'
	l += '#Column SEQ_POS: Position in the given GPCR\'s sequence\n'
	l += '#Column AA_WT: Amino acid in the Wild Type sequence\n'
	l += '#Column AA_VAR: Amino acid in the Mutated/Variant sequence\n'
	l += '#Column PROB: PRECOG2 probability\n'
	l += '#GPCR\tALN_POS\tSEQ_POS\tAA_WT\tAA_VAR\tPROB\n'
	for gpcr in gpcr_list:
		data = read_aln(pos, neg, hmm_pos, hmm_neg, features, gprotein, gpcr, obj)
		#print data[:, 3:]
		feature_matrix = data[:, 5:]
		model = extract_model(gprotein)
		#min_max = k_fold('/net/netfile2/ag-russell/bq_gsingh/gpcr/update/feature_files/'+str(gprotein)+'_train.txt')
		min_max = k_fold('data/feature_files/'+str(gprotein)+'_train.txt')
		feature_matrix = min_max.transform(numpy.array(feature_matrix))
		Y = model.predict(feature_matrix)
		Y_prob = model.predict_proba(feature_matrix)

		for (name, aln_position, mut_position, aa, mut), y in zip(data[:, :5], Y_prob):
			#obj[name].mutation[mut][gprotein] = round(y[1], 3)
			#print name, aln_position, mut_position, aa, mut, round(y[1], 3)
			l += str(name) + '\t' + str(aln_position) + '\t' + str(mut_position) + '\t' + str(aa) + '\t' + str(mut) + '\t' + str(round(y[1], 3)) + '\n'

	return l
	#print 'Completed.'
