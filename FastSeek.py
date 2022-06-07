import os
import argparse
import re
import shutil

parser = argparse.ArgumentParser(
    description='Seek for target with BLAST.',
    add_help=False,
    usage='\nFastSeek.py --ref [Reference seq] --fna [Input genome] --evalue [Optional, default: 0.001] --ID ['
          'Optional, mark your sequence, default: Target_seq] --protein [Optional, when ref seq was amino acids.] --threads [Optional, threads, default: 6] --silent[Optional, do not print results on screen.]')

parser.add_argument(
    '--ref',
    metavar='[input.gene]',
    required=True,
    type=str)
parser.add_argument(
    '--fna',
    metavar='[input.fna]',
    required=True,
    type=str)
parser.add_argument(
    '--evalue',
    metavar='[Evalue]',
    required=False,
    type=float)
parser.add_argument(
    '--ID',
    metavar='[mark]',
    required=False,
    type=str)
parser.add_argument(
    '--protein',
    action='store_true',
    help='Amino acids',
    required=False)
parser.add_argument(
    '--silent',
    action='store_true',
    help='silent running',
    required=False)
parser.add_argument(
    '--threads',
    metavar='[threads]',
    required=False,
    type=int)

args = parser.parse_args()

if args.evalue:
    evalue = args.evalue
else:
    evalue = 0.001

if args.ID:
    seqID = args.ID
else:
    seqID = 'Target_seq'

if args.threads:
    threads = args.threads
else:
    threads = 6

folder = os.path.exists('temp_blast')
if not folder:
    blastfolder = 'temp_blast'
    os.mkdir(blastfolder)
    print("--Running BLAST--")
else:
    blastfolder = 'temp_blast' + '_1'
    os.mkdir(blastfolder)
    print("Temp folder existed, redirecting to a new folder." + "\n" + "--Running BLAST--")


def blastcmd(db, genome, output, eValue):
    if not args.protein:
        os.system(
            'makeblastdb -in ' + db + ' -dbtype nucl -parse_seqids -out ' + blastfolder + '/' + db +'>/dev/null 2>&1')
        os.system(
            'blastn -query ' + genome + ' -db ' + blastfolder + '/' + db + ' -evalue ' + str(
                eValue) + " -outfmt '6 pident length mismatch gapopen evalue bitscore qstart qend qseq' -num_threads 6 -out " + blastfolder + '/' + output+'>/dev/null 2>&1')
    else:
        os.system(
            'makeblastdb -in ' + db + ' -dbtype prot -parse_seqids -out ' + blastfolder + '/' + db+'>/dev/null 2>&1')
        os.system(
            'tblastn -query ' + genome + ' -db ' + blastfolder + '/' + db + ' -evalue ' + str(
                eValue) + " -outfmt '6 pident length mismatch gapopen evalue bitscore qstart qend qseq' -num_threads 6 -out " + blastfolder + '/' + output+'>/dev/null 2>&1')
    outputfile = str(blastfolder + '/' + output)
    list = []
    matchPattern = re.compile(r'#')
    file = open(outputfile, 'r')
    line = file.readlines()
    for i in line:
        if matchPattern.search(i):
            pass
        else:
            list.append(i)
    file.close()
    os.mkdir(seqID+'Result')
    file = open(seqID+'Result/summary.txt', 'w')
    file.write('Gene\tIdentity\tLength\tMismatch\tGap\tE-value\tScore\tStart\tEnd\n')
    if not args.silent:
        print('Gene\tIdentity\tLength\tMismatch\tGap\tE-value\tScore\tStart\tEnd\n')
    seqorder = 0
    for i in list:
        seqorder = seqorder + 1
        file.write(seqID + '_' + str(seqorder) + '\t' + str(i.split('\t')[0])+'\t'+str(i.split('\t')[1])+'\t'+str(i.split('\t')[2])+'\t'+str(i.split('\t')[3])+'\t'+str(i.split('\t')[4])+'\t'+str(i.split('\t')[5])+'\t'+str(i.split('\t')[6])+'\t'+str(i.split('\t')[7]))
        if not args.silent:
            print(seqID + '_' + str(seqorder) + '\t' + str(i.split('\t')[0])+'\t'+str(i.split('\t')[1])+'\t'+str(i.split('\t')[2])+'\t'+str(i.split('\t')[3])+'\t'+str(i.split('\t')[4])+'\t'+str(i.split('\t')[5])+'\t'+str(i.split('\t')[6])+'\t'+str(i.split('\t')[7]))
    file.close()
    file = open(seqID+'Result/'+seqID+'.fna', 'w')
    fastaorder = 0
    for i in list:
        if int(i.split('\t')[7]) > int(i.split('\t')[6]):
            strand = "forward"
        else:
            strand = 'reverse'
        fastaorder = fastaorder + 1
        file.write('>' + seqID + '_' + str(fastaorder) + ' ' + strand + '\n' + str(i.split('\t')[-1]))
    file.close()


blastcmd(args.ref, args.fna, seqID, evalue)

shutil.rmtree(blastfolder)
