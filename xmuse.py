import argparse, sys, os
import muse

description="XMUSE is a pipeline to analyze nucleotide sequence feature using Latent Dirichlet Allocation. You can also call the 'muse' class directly. Please check out https://github.com/xnao25/XMUSE for more details."
parser=argparse.ArgumentParser(description=description)

parser.add_argument("-if","--input_fasta",required=True,help="Path of the input fasta file(s). Use ',' to separate the paths of multiple fasta files.")
parser.add_argument("-it","--input_tag",required=True,help="Tags of the input sequences. Provide the path to a txt file if you have only one fasta input or a ',' separated string for multiple fasta input, where each item corresponds to one fasta.")
parser.add_argument("-m","--mode",default='UWSC',help="Mode to run LDA. Choose from ASC, AWSC, and UWSC. If pre-fitted model is provided, this argument will be ignored.")
parser.add_argument("-k","--kmer",default=4,type=int,help='An integer of k-mer length. If pre-fitted model is provided, this argument will be ignored.')
parser.add_argument("-t","--topic",default=4,type=int,help="An integer of the number of topics. If pre-fitted model is provided, this argument will be ignored.")
parser.add_argument("-sw","--sliding_window",default=0,type=int,help="An interger of the size of sliding window.")
parser.add_argument("-sg","--subgroup",default=1,type=int,help="An integer of subgroup size. Default 1")
parser.add_argument("-s","--startpos",help="An interger of the name/position of the first nucleotide in the sequence.",default=0,type=int)
parser.add_argument("-o","--out_dir",required=True,help="Path to the output directory.")
parser.add_argument("-fm","--fitted_model",default='',help="Path to a pickle file of a muse object. This model will be used to transform the new sequences provided. A pickle file of a muse object will be saved every time after you run this pipeline. Do not use this argument when you want to make a new model.")
#parser.print_help()

def main():
    args=parser.parse_args()
    input_fa=args.input_fasta.split(',') if ',' in args.input_fasta else args.input_fasta
    input_label=args.input_label.split(',') if ',' in args.input_tag else args.input_tag
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)
    muse.run_xmuse(fasta_path=input_fa,label_path=input_label,mode=args.mode,k=args.kmer,topic_num=args.topic,sliding_window=args.sliding_window,subgroup_size=args.subgroup,startpos=args.startpos,modelpath=args.fitted_model,outpath=args.out_dir)

if __name__=='__main__':
    main()