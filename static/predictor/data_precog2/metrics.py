import math
from sklearn.metrics import roc_auc_score as auc_score
from sklearn.metrics import matthews_corrcoef as mcc_score
from sklearn.metrics import accuracy_score as acc_score
from sklearn.metrics import f1_score as f1m_score
from sklearn.metrics import precision_score as pre_score
from sklearn.metrics import recall_score as rec_score
from sklearn.metrics import matthews_corrcoef as mcc_score

def main(tp, fp, fn, tn, exp, pred):
	try:
		mcc = round(mcc_score(exp, pred), 2)
	except ZeroDivisionError:
		mcc = 'NaN'
	try:	
		acc = round(acc_score(exp, pred), 2)
	except ZeroDivisionError:
		acc = 'NaN'
	try:
		pre = round(pre_score(exp, pred), 2)
	except:
		pre = 'NaN'
	try:					
		rec = round(rec_score(exp, pred), 2)
	except ZeroDivisionError:
		rec = 'NaN'
	try:
		spe = round((float(tn)/(tn+fp)), 2)
	except ZeroDivisionError:
		spe = 'NaN'
	try:
		f1m = round(f1m_score(exp, pred), 2)
	except ZeroDivisionError:
		f1m = 'NaN'
	try:
		if len(exp) != list(exp).count(0) and len(exp) != list(exp).count(0):
			auc = round(auc_score(exp, pred), 2)
		else:
			auc = '-'
	except ZeroDivisionError:
		auc = 'NaN'

	return mcc, acc, pre, rec, spe, f1m, auc
