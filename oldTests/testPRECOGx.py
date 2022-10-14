import os, sys
path = os.getcwd()
sys.path.insert(1, path + '/static/predictor/')
import extract
import predict
import precogx

input = 'AVPR2/Y128S\nAVPR2/R137H\nAVPR2/R181C\nAVPR2/R137C'
input = '# G-protein coupled receptor 183\nP32249\n\n# Thyrotropin receptor\nTSHR_HUMAN'
input_file = None
uniq_id, errorCode, flaggedGPCR = precogx.main(15, input, input_file, 'all', os.getcwd())
