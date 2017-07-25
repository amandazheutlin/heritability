import csv, os

# FUNCTIONS
################################################################

# General function to read in a text file, and return a data array (e.g. read trial order)                      
def read_txt(filepath):
    open_file = open(filepath,'r')
    data = open_file.readlines()
    data_out = []
    for line in data:
        line = line.replace('\n','')
        data_out.append(line)
    open_file.close()
    return data_out

# General function to read in a text file, and return a data array (e.g. read trial order)                      
def write_txt(filepath, data):
    print filepath
    open_file = open(filepath,'w')
    for row in data:
        open_file.write(row + '\n')
    open_file.close()
    
def read_csv(filepath):    
    with open(filepath, 'rb') as csvfile:
        data = csv.reader(csvfile, delimiter=',')
        matrix_vals = []
        for row in data:
            matrix_vals.append(row)
    return matrix_vals


# SCRIPT
################################################################

# Define directories/files
out_dir = '/data/swe_gwas/ABZ/immunegenes/comp_h2/GENE_SCRIPTS/'
cov_file = '/data/swe_gwas/ABZ/immunegenes/comp_h2/complement_covar4.csv'
script_types = ['ACE','ADE','AE','CE','DE','E']

# Read in MX script and CSV with covariance values
matrix_vals = read_csv(cov_file)
        
# get the column index for each value
cov_idxs = dict()
for idx,head in enumerate(matrix_vals[0]):
    cov_idxs[head]=idx


# for each model type
for model in script_types:
    mx_script = read_txt(model + 'model.mx')
    # for each gene
    for gene in range(1,len(matrix_vals)):
        gene_name =  matrix_vals[gene][cov_idxs['T1_VAR']].replace('x1','').replace('X1','')
        new_script = []
        script_name = gene_name + '_' + model + 'model_comp4.mx'
		# read the current line of the covariance csv, corresponding to the covariance estimates for a symptom
        cur_cov_values = matrix_vals[gene]
    
    	# read through each line in the mx_script
        for line in mx_script:
            new_line = line
		    
            # for each string to replace
            for cur_key in cov_idxs:
                
                # if the string is in the current line of the script, replace it
 				if cur_key in line: 
					new_line = new_line.replace(cur_key, cur_cov_values[cov_idxs[cur_key]])
            new_script.append(new_line) # append to the new script
        write_txt(out_dir + script_name, new_script)    # write the new script
        mx_cmd = 'mx <' + out_dir + script_name + ' > /data/swe_gwas/ABZ/immunegenes/comp_h2/HERITABILITY_OUT/comp_swedes4/' + gene_name + '_' + model + '_comp4_output.mx'
        os.system(mx_cmd)

    
