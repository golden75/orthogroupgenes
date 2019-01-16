**Neranjan Perera**
<h1>get_orthogroupgenes</h1>

The program will selects the grouped ortholog genes for a particular species.		

usage: `get_orthogroup_genes.py [-h] [--orthofile ORTHO_GROUP_FILE] [-fasta FASTA_FILE] [--tag SPECIES_TAG] [-l] [-s]
                                [--outfile OUTPUT_FILE] [--version]`

optional arguments:
 - -h, --help            show this help message and exit
 
optional arguments:
 - -h, --help            show this help message and exit
 
optional arguments:
- -h, --help            show this help message and exit


usage: `get_orthogroup_genes.py [-h] [--orthofile ORTHO_GROUP_FILE]		
                               [--fasta FASTA_FILE] [--tag SPECIES_TAG] [-l]. 
                               [-s] [--outfile OUTPUT_FILE] [--version] `
                               

optional arguments: 
 ... -h, --help            show this help message and exit.
  
  --orthofile ORTHO_GROUP_FILE. 
                        ortho group file name. 
                        
  --fasta FASTA_FILE    FASTA file name. 
  
  --tag SPECIES_TAG     gene identification tag for the species eg: TR is the
                        tag for TR22118|c0_g1_i1_14120|m.21703. 
                        
  -l                    boolean switch selects the longest gene length. If not
                        used default is set to select the longest. 
                        
  -s                    boolean switch selects the shortest gene length. If
                        not used default is set to select the longest. 
                        
  --outfile OUTPUT_FILE
                        output file name of the sorted genes. If not given it
                        will use the default name; sorted_genes.fasta. 
                        
  --version             show program's version number and exit. 

s = "Python syntax highlighting"
print s.

--

```python
import sys,re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import MutableSeq
from array import array



def get_sequence(f_file,gid,out_file):
        f = open(f_file, 'r')
        for seq_record in SeqIO.parse(f,"fasta"):
                if ( gid == seq_record.id):
                        rec = SeqRecord(seq_record.seq,
                                id = seq_record.id,
                                description = "")
                        SeqIO.write(rec, out_file, "fasta")
                        break
        f.close()
        return;


def get_csv_geneID(csv_line):
        ln = csv_line.strip().split(",")
        return ln[1].split("\"")[1].split(":")[1]

def get_seqid(g_id,fasta_out):
        gtf = open(sys.argv[2],"r")
        for gtf_line in gtf:
                gtf_line = gtf_line.split('#', 1)[0]
                if(gtf_line != ""):
                        gln = gtf_line.strip().split("\t")
                        gid = gln[5]
                        if(gid == g_id):
                                print("GeneID:%s\ttableID:\t%s\tprotein-ID:\t%s" %(g_id, gid, gln[7]))
                                get_sequence(sys.argv[3],gln[7],fasta_out)
                                break

```
