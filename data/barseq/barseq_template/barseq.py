import sys
from Bio import SeqIO
import csv
import copy
import operator


# first argument is sequences
# second is tab-delimited association of barcodes with samples
# Third is tab-delimited association of barcodes with genes
# Fourth is condensed output file
# Fifth is list of identified barcodes for each sample
# Sixth is a suffix for construction of lists of all reads with their UMI and gene for each sample

def hamming2(str1, str2):
	assert len(str1) == len(str2)
	ne = operator.ne
	return sum(map(ne, str1, str2))

samples=[] 
# list of samples in run
samplerbarcodes=[] 
#Rev barcodes corresponding to each sample
genes=[] 
#list of genes
genebarcodes=[]
#list of barcodes corresponding to each gene
badseq=[] 
#variable denoting if we encounter an irregular sequence
samplebarcode1error=[] 
#denotes a barcode not found in samplefbarcodes
flank1error=[] 
#denotes flank sequence which differs from expected before gene barcode
flank2error=[] 
#denotes flank sequence which differs from expected after gene barcode
samplebarcodemismatch=[] 
#denotes if forward and reverse sample barcodes disagree
samplebarcode2error=[] 
#denotes a barcode not found in samplerbarcodes
genebarcodeerror=[] 
#denotes a gene barcode not found in genebarcodes
#initializes arrays

with open(sys.argv[3],'r') as f:
	reader=csv.reader(f,delimiter='\t')
	for gene,genebarcode,in reader:
		genes.append(gene)
		genebarcodes.append(genebarcode)
genedict={key: {} for key in genes}
#makes a dictionary with each key being a different gene.

with open(sys.argv[2],'r') as f:
	reader=csv.reader(f,delimiter='\t')
	for sample,samplerbarcode in reader:
		samples.append(sample)
		samplerbarcodes.append(samplerbarcode)
results={key: copy.deepcopy(genedict) for key in samples}
#makes a dictionary with each key being a sample, with values as the gene dictionary developed above.  Use deepcopy because apparently
#python will keep all sample gene dictionaries the same otherwise.

errorreads=[]
#this will store all the reads that raise an error during processing so you can manually look at them.

results['flank1error']=0
results['flank2error']=0
results['samplebarcode2error']=0
results['genebarcodeerror']=0
errors={}
errors['flank1errors']={}
errors['flank2errors']={}
errors['sample2barcodeerrors']={}
errors['genebarcodeerrors']={}
#initializes result dictionaries
	
fasta_sequences = SeqIO.parse(open(sys.argv[1]),'fastq') #open input fastq sequences 
with open(sys.argv[4],'w') as out_file: #going to put everything in a condensed output file.
	for fasta in fasta_sequences:
		badseq=0
		flank1error=0
		flank2error=0
		samplebarcode2error=0
		genebarcodeerror=0
		#resets error variables
		sequence = str(fasta.seq)
			
		#now we go to the first flank sequence: flank1
		flank1=str(sequence[0:20])
		flank1hamming=hamming2(flank1,'TAGCCATGTCATAGACGTCC')
		# find the hamming distance between the read's flank1 sequence and the sequence it should be
		if flank1hamming > 5:
			flank1error=1
			badseq=1
		# if the hamming distance is geater than 5, then we will consider the read too error-prone to assign reliably.
		
		# now we go to the second flank sequence: flank2
		flank2=str(sequence[27:47])
		flank2hamming=hamming2(flank2,'GTACGTTGATCAAGTCCCGA')
		if flank2hamming > 5:
			flank2error=1
			badseq=1
		# do the same thing as for flank1 above.
		
		#now we go to the second barcode: samplebc2
		samplebc2=str(sequence[47:53])
		samplebc2hamming=[]
		for testbarcode in samplerbarcodes:
			samplebc2hamming.append(hamming2(testbarcode, samplebc2))
		minsamplebc2hamming=min(samplebc2hamming)
		countminsamplebc2hamming=samplebc2hamming.count(minsamplebc2hamming)
		if countminsamplebc2hamming > 1:
			samplebarcode2error=1
			badseq=1
			sampler='unknown'
		elif minsamplebc2hamming > 1:
			samplebarcode2error=1
			badseq=1
			sampler='unknown'
		else:
			sampler=samples[samplebc2hamming.index(minsamplebc2hamming)]
			sample=sampler

		barcode=sequence[53:70] #gets degenerate region
		
		# now we go to the gene barcode
		genebc=str(sequence[20:27])
		genebchamming=[]
		for testbarcode in genebarcodes:
			genebchamming.append(hamming2(testbarcode, genebc))
		mingenebchamming=min(genebchamming)
		countmingenebchamming=genebchamming.count(mingenebchamming)
		if countmingenebchamming > 1:
			genebarcodeerror=1
			badseq=1
		elif mingenebchamming > 1:
			genebarcodeerror=1
			badseq=1
		else:
			gene=genes[genebchamming.index(mingenebchamming)]
		# do the same thing for genebc as for sample1bc and sample2bc above.
		
		errorstring=''.join(map(str,[flank1error, genebarcodeerror, flank2error, samplebarcode2error]))
		#makes a string so one can see what is wrong with a particular read at a glance
		
		if badseq==0 and barcode in results[sample][gene]:
			results[sample][gene][barcode]+=1 #if sequence is legit and degenerate barcode has been seen before, add 1 to running total
		elif badseq==1:
			errorreads.append([sequence, errorstring])
			
			if flank1error==1 and flank1 in errors['flank1errors']:
				results['flank1error']+=1
				errors['flank1errors'][flank1]['count']+=1
			elif flank1error==1:
				results['flank1error']+=1
				errors['flank1errors'][flank1]={}
				errors['flank1errors'][flank1]['count']=1
				errors['flank1errors'][flank1]['hamming']=flank1hamming
			# make running totals of mismatched flank sequences
			
			if flank2error==1 and flank2 in errors['flank2errors']:
				results['flank2error']+=1
				errors['flank2errors'][flank2]['count']+=1
			elif flank2error==1:
				results['flank2error']+=1
				errors['flank2errors'][flank2]={}
				errors['flank2errors'][flank2]['count']=1
				errors['flank2errors'][flank2]['hamming']=flank2hamming
			# make running totals of mismatched flank sequences
			
			if samplebarcode2error==1 and samplebc2 in errors['sample2barcodeerrors']:
				results['samplebarcode2error']+=1
				errors['sample2barcodeerrors'][samplebc2]['count']+=1
			elif samplebarcode2error==1:
				results['samplebarcode2error']+=1
				errors['sample2barcodeerrors'][samplebc2]={}
				errors['sample2barcodeerrors'][samplebc2]['count']=1
				errors['sample2barcodeerrors'][samplebc2]['hamming']=minsamplebc2hamming
			#make running totals of offending sample2 barcodes
			
			if genebarcodeerror==1 and genebc in errors['genebarcodeerrors']:
				results['genebarcodeerror']+=1
				errors['genebarcodeerrors'][genebc]['count']+=1
			elif genebarcodeerror==1:
				results['genebarcodeerror']+=1
				errors['genebarcodeerrors'][genebc]={}
				errors['genebarcodeerrors'][genebc]['count']=1
				errors['genebarcodeerrors'][genebc]['hamming']=mingenebchamming
			#make running totals of offending gene barcodes.
			
		else:
			results[sample][gene][barcode]=1 #otherwise, add a new entry for the new degenerate region you just found.
			
	for data in ['flank1error', 'flank2error', 'samplebarcode2error', 'genebarcodeerror']:
		datavalue=results[data]
		outputstring='%s\t%i\n' % (data, datavalue)
		out_file.write(outputstring)
		del results[data]
	#output total number of errors and delete them from dictionary.
	
	for sample in results:
		for gene in results[sample]:
			totalbcs=len(results[sample][gene])
			outputstring='%s\t%s\t%i\n' % (sample, gene, totalbcs)
			out_file.write(outputstring)
	#output total number of degenerate regions observed for each gene per sample
	
	for error in errors:
		for errortype in errors[error]:
			errorcount=errors[error][errortype]['count']
			errorhamming=errors[error][errortype]['hamming']
			outputstring='%s\t%s\t%s\t%i\t%s\t%i\n' % (error, errortype, 'count', errorcount, 'hamming', errorhamming)
			out_file.write(outputstring)
	#output each offending sequence and the number of times you observed it.
	
with open(sys.argv[5],'w') as out_file: #outputting to larger output file
	for sample in results:
		for gene in results[sample]:			
			for entry in results[sample][gene]:
				entryvalue=results[sample][gene][entry]
				outputstring="%s\t%s\t%s\t%i\n" % (sample, gene, entry, entryvalue)
				out_file.write(outputstring)
	#list each barcode identified and the number of times you identified it.
	for line in errorreads:
		read=line[0]
		errorstring=line[1]
		outputstring="%s\t%s\n" % (read, errorstring)
		out_file.write(outputstring)
	#list each read that raised an error and what was wrong with it.

#outputting UMI file for UMI-tools
for sample in results:
	with open(sample+sys.argv[6],'w') as out_file:
		readnum=1 
		for gene in results[sample]:			
			for UMI in results[sample][gene]:
				UMIcount=results[sample][gene][UMI]
				for i in range(0,UMIcount):
					outputstring="%i_%s\t%s\n" % (readnum, UMI, gene)
					out_file.write(outputstring)
					readnum+=1
