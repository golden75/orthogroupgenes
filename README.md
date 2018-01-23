**Neranjan Perera**
<h1>get_orthogroupgenes</h1>

The program will selects the grouped ortholog genes for a particular species.		

usage: `get_orthogroup_genes.py [-h] [--orthofile ORTHO_GROUP_FILE]		
                               [--fasta FASTA_FILE] [--tag SPECIES_TAG] [-l]. 
                               [-s] [--outfile OUTPUT_FILE] [--version] `
                               

optional arguments:.   
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


