import csv, os, operator

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
            matrix_vals.append(row[0].split(','))
    return matrix_vals

# Define directories/files
out_dir = '/data/swe_gwas/ABZ/immunegenes/comp_h2/HERITABILITY_OUT/comp_swedes4/'

# Get Every Output Mx File
output_files = os.listdir(out_dir)

# Because Mx Only Allow Single Characters to Define Arrays, the Output Matrices are Named Funny. Conversion diDict = 
matrix_vals = {'Q':'A_std', 'R':'C_std', 'S':'E_std', 'T':'D_std'}

# Create empty Output dictionaries
ALL_MODELS=dict()
WRITE_ALL_OUT = []
for file in output_files:
    cur_gene = file.split('_')[0]
    # Make Output Field for the Current Gene
    ALL_MODELS[cur_gene] = dict()

# Iterate Through Each Output File
for file in output_files:
    # Get Symptom/Model Info, Read in Data
    print file
    cur_gene = file.split('_')[0]
    cur_model = file.split('_')[1]
    cur_output = read_txt(out_dir + file)
    
    # heritability estimates for this file
    cur_heritability = dict()
    for key in matrix_vals.values():
        cur_heritability[key] = '-'
    
    # Iterate Through Each Line in the Output
    for idx,line in enumerate(cur_output):
        # If the Matrix is Specified in variable 'matrix_vals'
        for mat_letter in matrix_vals.keys():
            if 'MATRIX ' + mat_letter in line:
                if '[=V~*' + matrix_vals[mat_letter].replace('_std','') in cur_output[idx + 2]:
                # if this line corresponds to the output porition of the file
                    cur_est = cur_output[idx+4].split()[1] # the standardized estimate is 4 lines down
                    cur_heritability[matrix_vals[mat_letter]] = cur_est
        
        # Model Fit Estimates
        if 'Chi-squared fit of model >>>>>>>' in line:	
	    chi2 = str(float(line.split('>>>>>>>')[1]))
        if ' Degrees of freedom' in line:
            df = str(float(line.split('>>>>>>>>>>>>>')[1]))
        if 'Probability >>>>>>>>>>>>>>>>>>>>' in line:
            prob = str(float(line.split('>>>>>>>>>>>>>>>>>>>>')[1]))
        if 'Akaike' in line:
            akaike = str(float(line.split('>')[1]))
        if 'RMSEA >>>>>>>>>>>>>>>>>>>>>>>>>>' in line:
            rmsea = str(float(line.split('>>>>>>>>>>>>>>>>>>>>>>>>>>')[1]))
        
        # get unstandardized ACE values
        ACE_letters = {'A':'X','C':'Y','E':'Z'}
        ace_file = cur_gene + '_ACE_comp4_output.mx'
        ace_output = read_txt(out_dir + ace_file)
        # Iterate Through Each Line in the Output
        for idx2,line2 in enumerate(ace_output):
            # If the Matrix is Specified in variable 'matrix_vals'
            for mat_letter in ACE_letters.keys():
                # if this line corresponds to the output portion of the file
                if 'MATRIX ' + mat_letter in line2:
                    if '[=' + ACE_letters[mat_letter] + '*' + ACE_letters[mat_letter] in ace_output[idx2 + 2]:
                        cur_est = ace_output[idx2+4].split()[1] # the standardized estimate is 4 lines down
                        cur_heritability[mat_letter] = cur_est
        
    # put heritability estimates in an ordered string
    herit_out = ''
    for letter in ACE_letters.keys():
	herit_out = herit_out + str(cur_heritability[letter]) + ','
    for letter in ['A_std', 'C_std', 'E_std', 'D_std']:
	herit_out = herit_out + str(cur_heritability[letter]) + ','
    # A C E A_std C_std E_std D_std chi2 dF Prob Akaike RMSEA
    write_out = herit_out + chi2 + ',' + df + ',' + prob + ',' + akaike + ',' + rmsea
    ALL_MODELS[cur_gene][cur_model.replace('.mx','')] = write_out
    WRITE_ALL_OUT.append(cur_gene + ',' + cur_model.replace('.mx','') + ',' + write_out)


# Find the winning model
header_idx = {'A':0,'C':1,'E':2,'A_std':3,'C_std':4,'E_std':5,'D_std':6,'chi2':7,'df':8,'prob':9,'akaike':10,'rmsea':11}
akaike_idx = header_idx['akaike']
all_best_models = ['Gene,BestModel,A,C,E,A_std,C_std,E_std,D_std,chi2,df,prob,akaike,rmsea']
for cur_gene in ALL_MODELS.keys():
    akaike_dict = dict()
    model_names = ALL_MODELS[cur_gene].keys()
    for cur_model in model_names:
        akaike_dict[cur_model] = float(ALL_MODELS[cur_gene][cur_model].split(',')[akaike_idx])
    winning_model = min(akaike_dict.iteritems(), key=operator.itemgetter(1))[0]
    write_out = cur_gene + ',' +  winning_model + ',' + ALL_MODELS[cur_gene][winning_model]
    all_best_models.append(write_out)
write_txt('/data/swe_gwas/ABZ/immunegenes/comp_h2/SUMMARY_OUT/best_fitting_models_comp4.csv', all_best_models)
write_txt('/data/swe_gwas/ABZ/immunegenes/comp_h2/SUMMARY_OUT/all_models_comp4.csv', WRITE_ALL_OUT)

