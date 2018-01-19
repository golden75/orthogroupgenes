#!/usr/bin/env python

                ##################################################################################
                ## Neranjan V Perera (2018 January 8)                                           ##
                ## Computational Biology Core                                                   ##
                ## Institute for Systems Genomics                                               ##
                ## University of Connecticut                                                    ##
                ## Copyright                                                                    ##
                ##                                                                              ##
                ## The program will take orthogroup genes and will find the genes of a speies   ##
                ## which belong to a orthogroup                                                 ##
                ##                                                                              ##
		## argument:									##
		## python get_orthogroup_genes.py --orthofile Orthogroups.txt			## 
		## --fasta Funaria.fasta --tag TR -l --outfile sorted_ortholog_genes.fasta	##
                ##################################################################################

import sys,re,argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def info():
        print ("""
                ##################################################################################
                ## Neranjan V Perera (2018 January 8)                                           ##
                ## Computational Biology Core                                                   ##
                ## Institute for Systems Genomics                                               ##
                ## University of Connecticut                                                    ##
                ## Copyright                                                                    ##
                ##                                                                              ##
                ## The program will take orthogroup genes and will find the genes of a speies   ##
                ## which belong to a orthogroup                                                 ##
                ##                                                                              ##
                ##################################################################################
        """)
        return


'''
gets the line containing the gene_Tag and filters out the corrosponding genes 
where only the orthogroups which contains species genes are selected.
[exclusive orthogroup genes for the Tag species will be selected]
input = orthogroup line containing genes
output = gene count of the species 
'''
def orthogroup_genes(ln,sp_tag):
        gene_id_list = ""
        passing_gene_id_list = ""
        tag="^"+sp_tag

	#get the ortholog group line and split the genes into columns
	n_col = len(ln.split())
	word = ln.split()
	
	# if the number of columns are greater than 2 
	# which selects the genes of the TAG species which shears the gene groups with the rest of the species
	# we are removing the orthogroups which have n_col == 1 ; where they are "unassigned genes"
	tag_gene_count = 0
	if(n_col > 2):
		#search for genes which startes with the species-TAG 
		for n in range(1, n_col):
			if(re.search(tag,word[n])):
				tag_gene_count = tag_gene_count + 1
			gene_id_list = gene_id_list + word[n] + "\t"
		
		# Selects the orthogroups with exclusively tag_genes by counting the number of genes in the list is
		# equal to (n_col - 1) ; {n_col-1 is due to removing the orthogroup tag}
		if(tag_gene_count == n_col-1):
			passing_gene_id_list = gene_id_list.strip()
	#return the gene list as a string
	return(passing_gene_id_list)


'''
This will read the gene list and determine the longest and will wirte the FASTA sequence to a file
input = gene list as a string
'''
def select_gene_id(gene_list,f_file,fasta_out,length_flag):
	#need to parse through the fasta file to find the gene in the list
	#then store each gene information in to a TUPLE, with the gene length information
	#then append the TUPLE into a LIST creating a "NESTED-TUPLE"
	#Once the gene information for the ortho-group have been stored:
	#sort the "NESTED-TUPLE" to get the longest gene and write it into a file
	gene_records=[]
	gene_list = gene_list.lstrip(' ')
	gene = gene_list.split("\t")
	for g in gene:
		#reads the fasta file and find the length 
		for seq_record in SeqIO.parse(f_file, "fasta"):
			#match the gene_id to the fasta record;
			#then store it in a tuple
			 if(seq_record.id == g):
				print("match found  gene-id: %s  length: %s" %(seq_record.id,len(seq_record)))
				entry = (len(seq_record),seq_record.seq,seq_record.id,seq_record.description)
				gene_records.append(entry)
	if(gene_records):
		gene_records.sort(reverse=length_flag)
		print("sorted gene: %s  length: %s \n" %(gene_records[0][2],gene_records[0][0]))
		rec = SeqRecord(gene_records[0][1],
			id = gene_records[0][2],
			description = gene_records[0][3])
		SeqIO.write(rec, fasta_out, "fasta")
	return	
	



'''
The main program 
This accepts arguments, where 'orthogroup file' and 'fasta file' is the main required inputs
and as the others are optional, where it will create a out put fasta file with the default name
where it gives the longest genes of a perticular ortho group.

uage: get_orthogroup_genes.py [-h] [--orthofile ORTHO_GROUP_FILE]
                               [--fasta FASTA_FILE] [--tag SPECIES_TAG] [-l]
                               [-s] [--outfile OUTPUT_FILE] [--version]

required arguments:
  --orthofile ORTHO_GROUP_FILE
                        ortho group file name
  --fasta FASTA_FILE    FASTA file name
  --tag SPECIES_TAG     gene identification tag for the species 
			eg: TR is the tag for TR22118|c0_g1_i1_14120|m.21703

optional arguments
  -l                    boolean switch selects the longest gene length. If not
                        used default is set to select the longest
  -s                    boolean switch selects the shortest gene length. If
                        not used default is set to select the longest
  --outfile OUTPUT_FILE
                        output file name of the sorted genes. If not given it
                        will use the default name; sorted_genes.fasta
  --version             show program's version number and exit

output = prints the total gene count
'''

def main():
	#creates a parser object
	parser = argparse.ArgumentParser()
	#create action variables to hold the passed arguments 
	parser.add_argument('--orthofile', action='store', dest='ortho_group_File',
		help='ortho group file name')
	parser.add_argument('--fasta', action='store', dest='fasta_file',
		help='FASTA file name ')
	parser.add_argument('--tag', action='store', dest='species_tag',
		help='gene identification tag for the species eg: TR is the tag for TR22118|c0_g1_i1_14120|m.21703')
	parser.add_argument('-l', action='store_true', default=True,
		dest='gene_length',
		help='boolean switch selects the longest gene length. If not used default is set to select the longest')
	parser.add_argument('-s', action='store_false', default=True,
                dest='gene_length', 
                help='boolean switch selects the shortest gene length. If not used default is set to select the longest')
	parser.add_argument('--outfile', action='store', default='sorted_genes.fasta',
		dest='output_file',
		help='output file name of the sorted genes. If not given it will use the default name; sorted_genes.fasta ') 
	parser.add_argument('--version', action='version', version='%(prog)s 1.0')

	results = parser.parse_args()
	print 'ortho group file name 	=', results.ortho_group_File
	print 'fasta file name		=', results.fasta_file
	print 'species tag		=', results.species_tag
	print 'sorted gene lenght	=', results.gene_length
	print 'output file name 	=', results.output_file
	print '\n\n'

	tag=results.species_tag
	#program checkes whether three main necessary arguments have been recived if not exit
	if(results.ortho_group_File and results.fasta_file and results.species_tag):
		info()
		#open the Orthogroup file to get the genes which belong to each orthogroup
		file = open(results.ortho_group_File , "r")
		##open a file to write the sorted ortholog genes
		out_file = results.output_file
		fasta_out=open(out_file,"w")
		total = 0
		for line in file:
			#
			if(re.search(tag,line)):
				line = line.strip()
				g_list = orthogroup_genes(line,tag)
				if(g_list):
					print("passed gene list :%s" %(g_list))
					#select_gene_id(g_list,sys.argv[2],fasta_out)
					select_gene_id(g_list,results.fasta_file,fasta_out,results.gene_length)
	else:
		print 'Program',sys.argv[0] ,': : error: too few arguments'
		print '''\
                use -h to get the help file
                '''
		




main()

