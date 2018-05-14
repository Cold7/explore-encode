import multiprocessing as mp
import requests, json, collections, time

def download_data(cell):
	n_replicates = 2
	characteristics_dict = {}
	
	exp_open_chrom = ["FAIRE-seq", "ATAC-seq", "DNase-seq", "MNase-seq"]
	exp_meth_dna = ["MRE-seq", "TAB-seq", "WGBS"]
	chipseq = ["ChIP-seq"]
	quant = ["polyA RNA-seq","total RNA-seq"]
	loop = ["ChIA-PET", "Hi-C", "5C"]
	
	########### reading a  list of known tfs in human
	

	
	output = open(cell.replace(" ","_")+".txt","w")
	output.write("computing for "+cell+"\n" )
	characteristics_dict[cell.replace(" ","_")] = {"open_chromatin" : 	{"hg19" : {"G1" : {"treatment" : 0, "genetic_mod" : 0, "none" : 0}, "G2" : {"treatment" : 0, "genetic_mod" : 0, "none" : 0}, "none" : {"treatment" : 0, "genetic_mod" : 0, "none" : 0}}, "GRCh38" : {"G1" : {"treatment" : 0, "genetic_mod" : 0, "none" : 0}, "G2" : {"treatment" : 0, "genetic_mod" : 0, "none" : 0}, "none" : {"treatment" : 0, "genetic_mod" : 0, "none" : 0}}},
													"DNA_meth" : 		{"hg19" : {"G1" : {"treatment" : 0, "genetic_mod" : 0, "none" : 0}, "G2" : {"treatment" : 0, "genetic_mod" : 0, "none" : 0}, "none" : {"treatment" : 0, "genetic_mod" : 0, "none" : 0}}, "GRCh38" : {"G1" : {"treatment" : 0, "genetic_mod" : 0, "none" : 0}, "G2" : {"treatment" : 0, "genetic_mod" : 0, "none" : 0}, "none" : {"treatment" : 0, "genetic_mod" : 0, "none" : 0}}}, 
													"histone_marks" :	{"hg19" : {"G1" : {"treatment" : 0, "genetic_mod" : 0, "none" : 0}, "G2" : {"treatment" : 0, "genetic_mod" : 0, "none" : 0}, "none" : {"treatment" : 0, "genetic_mod" : 0, "none" : 0}}, "GRCh38" : {"G1" : {"treatment" : 0, "genetic_mod" : 0, "none" : 0}, "G2" : {"treatment" : 0, "genetic_mod" : 0, "none" : 0}, "none" : {"treatment" : 0, "genetic_mod" : 0, "none" : 0}}}, 
													"TF_bound" : 		{"hg19" : {"G1" : {"treatment" : 0, "genetic_mod" : 0, "none" : 0}, "G2" : {"treatment" : 0, "genetic_mod" : 0, "none" : 0}, "none" : {"treatment" : 0, "genetic_mod" : 0, "none" : 0}}, "GRCh38" : {"G1" : {"treatment" : 0, "genetic_mod" : 0, "none" : 0}, "G2" : {"treatment" : 0, "genetic_mod" : 0, "none" : 0}, "none" : {"treatment" : 0, "genetic_mod" : 0, "none" : 0}}}, 
													"quantification" : 	{"hg19" : {"G1" : {"treatment" : 0, "genetic_mod" : 0, "none" : 0}, "G2" : {"treatment" : 0, "genetic_mod" : 0, "none" : 0}, "none" : {"treatment" : 0, "genetic_mod" : 0, "none" : 0}}, "GRCh38" : {"G1" : {"treatment" : 0, "genetic_mod" : 0, "none" : 0}, "G2" : {"treatment" : 0, "genetic_mod" : 0, "none" : 0}, "none" : {"treatment" : 0, "genetic_mod" : 0, "none" : 0}}}, 
													"loops" :			{"hg19" : {"G1" : {"treatment" : 0, "genetic_mod" : 0, "none" : 0}, "G2" : {"treatment" : 0, "genetic_mod" : 0, "none" : 0}, "none" : {"treatment" : 0, "genetic_mod" : 0, "none" : 0}}, "GRCh38" : {"G1" : {"treatment" : 0, "genetic_mod" : 0, "none" : 0}, "G2" : {"treatment" : 0, "genetic_mod" : 0, "none" : 0}, "none" : {"treatment" : 0, "genetic_mod" : 0, "none" : 0}}}
													} 	
	
	try:
		URL = "https://www.encodeproject.org/search/?type=Experiment&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&biosample_term_name="+(cell).replace(" ","+")+"&limit=all"
		# GET the object
		response = None
		while response == None:
			try:
				response = requests.get(URL, headers=HEADERS)
				# Extract the JSON response as a python dict
				response_json_dict = response.json()
				
				i=0
				data = collections.defaultdict(dict)
				while True:	
					try:
						if len(response_json_dict["@graph"][i]["replicates"][0]["biological_replicate_number"]) >= n_replicates:
							ID = response_json_dict["@graph"][i]["@id"]
							URL2 = "https://www.encodeproject.org"+ID
							# GET the object
							response2 = None
							while response2 == None:
								try:
									response2 = requests.get(URL2, headers=HEADERS)
									# Extract the JSON response as a python dict
									response_json_dict2 = response2.json()
									assay= "None"
									try:
										assay = response_json_dict2["assay_title"]
									except:
										try:
											assay = response_json_dict2["assay_term_name"]
										except:
											pass
									current_assembly = []
									current_treatment = None
									current_phase = "none"
									current_genetic_modification = None
									status = ""
									#definition for keys to use
									try:
										current_assembly = response_json_dict2["assembly"]
									except:
										pass
									try:
										current_treatment = response_json_dict2["replicates"][0]["library"]["biosample"]["treatments"]
									except:
										pass
									try:
										current_phase = response_json_dict2["replicates"][0]["library"]["biosample"]["phase"]
									except:
										pass
									try:
										current_genetic_modification = response_json_dict2["replicates"][0]["library"]["biosample"]["genetic_modifications"]
									except:
										pass
									try:
										status = response_json_dict2["replicates"][0]["status"]
									except:
										pass

									
									current_exp = False
									if assay in exp_open_chrom:
										current_exp = "open_chromatin"
									elif assay in exp_meth_dna:
										current_exp = "DNA_meth"
									elif assay in chipseq:
										
										#if pipeline is complete and status is released
											#its time to know if it is an TF or a Histone mark
											target = response_json_dict2["target"]["name"].replace("-human","")
											
											current_exp = ""
											if target in  tfs: #if target is a TF
												current_exp = "TF_bound"
											elif target[0] == "H": #maybe it is an histone mark
												current_exp = "histone_marks"
									
									elif assay in quant:
										current_exp = "quantification"
									elif assay in loop:
										current_exp = "loops"
										
									if status == "released" and current_exp!= False and response_json_dict2["assembly"] != [] and (response_json_dict2["internal_status"] == "release ready"): #or response_Sjson_dict2["internal_status"] == "pipeline completed" ):
										for ass in current_assembly:

											current_adding = "none"
											if current_treatment != None:
												current_adding = "treatment"
											elif current_genetic_modification != None:
												current_adding = "genetic_mod"
											#print cell.replace(" ","_"), ass, current_phase, current_adding
											characteristics_dict[cell.replace(" ","_")][current_exp][ass][current_phase][current_adding] += 1
									
								except:
									output.write("error trying to get response from "+URL2+" retrying in 10 seconds")
									time.spleep(10)
					except:
						pass
					i += 1
			except:
				output.write("Some error trying to read a response... waiting 10 second until to retry do get data from "+URL)
				time.spleep(10)
	

	except:
		pass
	output.write("Resume of "+ cell+"\n" )
	output.write(str(json.dumps(characteristics_dict, indent=4, separators=(',', ': '))))
	output.close()
	print "##### end of calculations for "+ cell + "cell line."
	
	
	
global tfs

tfs = []
tf_file = open("TF_list", "r")
for line in tf_file:
	if line[0] != "#" and line != "\n":
		tfs.append(line.split("\t")[0])
			
output = open("review.txt","w")
# Force return from the server in JSON format
HEADERS = {'accept': 'application/json'}

cell_lines = [] #list of final cell that have all of experiments a priori selected


matrix_url = "https://www.encodeproject.org/matrix/?type=Experiment&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&x.limit="
# GET the object
response_matrix = (requests.get(matrix_url, headers=HEADERS)).json()
exp_types = []

#getting exp names
for item in response_matrix["matrix"]["x"]["buckets"]:
	exp_types.append(str(item["key"]))
#for pos in range(len(exp_types)):
#	print pos, exp_types[pos]


#getting experiment quantity
output.write("type\tOpen chromatine\tDNA meth\tchipseq\texpression quantyfication\tchr loops\n")
for bucket in response_matrix["matrix"]["y"]["biosample_type"]["buckets"]:
	data = bucket["biosample_term_name"]["buckets"]
	i = 0
	while True:
		try:
			groups_of_exp = [0,0,0,0,0] # open chromatine, DNA meth, chip-seq, cuantification of gene expression, chromatine loops
			number_of_exp = [0,0,0,0,0]			
			# if there exist at least 1 experiment
			if sum(data[i]["assay_title"]) != 0:
				#chipseq
				if data[i]["assay_title"][0] != 0:
					groups_of_exp[2] = 1
					number_of_exp[2] = data[i]["assay_title"][0]
				
				#dna methylation
				if  data[i]["assay_title"][6] != 0 or data[i]["assay_title"][18] != 0 or data[i]["assay_title"][38] != 0 or data[i]["assay_title"][10] != 0:
					groups_of_exp[1] = 1
					number_of_exp[1] = data[i]["assay_title"][1] + data[i]["assay_title"][6] + data[i]["assay_title"][18] + data[i]["assay_title"][38] + data[i]["assay_title"][10]
				
				#open chromatine
				if data[i]["assay_title"][28] != 0 or data[i]["assay_title"][15] != 0 or data[i]["assay_title"][1] != 0 :
					groups_of_exp[0] = 1
					number_of_exp[0] = data[i]["assay_title"][28] +data[i]["assay_title"][15] + data[i]["assay_title"][1]

				#expression pattern
				if data[i]["assay_title"][3] != 0 or data[i]["assay_title"][7] != 0:
					groups_of_exp[3] = 1
					number_of_exp[3] = data[i]["assay_title"][3] +data[i]["assay_title"][7]

				#chromatine loops
				if data[i]["assay_title"][26] != 0 or data[i]["assay_title"][29] != 0 or data[i]["assay_title"][36] != 0:
					groups_of_exp[4] = 1
					number_of_exp[4] = data[i]["assay_title"][26] +data[i]["assay_title"][29] +  data[i]["assay_title"][36]
				
				
				if sum(groups_of_exp) == len(groups_of_exp):
					cell_lines.append(data[i]["key"])
					output.write(str(data[i]["key"]+"\t"+str(groups_of_exp)+"\t"+str(number_of_exp)+"\n"))
				
		except:
			break
		i += 1
	output.write("########################################\n")

output.write("now compiling some stats about each experiment...\n")
output.close()
##################################################################################################
##
##							Second part of the script.
## Due to the first part we have cell lines with at least one experiment describin one 
## chromatine characteristic, so now is time to describe each line in a dictionary format.
## This format will have this structure
##	{
##		Cell : {
##				characterist 1: {
##								genome_version: {
##												cell_differentiation_state: {
##																			treatment: number, genetic_mod : number, other: number
##	}}}}}}
##												


pool=mp.Pool(processes=int(len(cell_lines))) #for multiprocessing
pool.map(download_data,(cell_lines))	
	
#print "done"
