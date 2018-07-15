import sys
import os
import json
import pprint
from pprint import pprint
from collections import defaultdict		

class ParseMiranda:

	def __init__(self, type, source, organism):
		self.type = type
		self.source = source
		self.organism = organism
		self.my_dir = ""

	def loadJson(self, file):
		'''
			input: file path in current directory
			return: data frame from json
		'''
		df_file = os.path.join(self.my_dir, file)
		with open(df_file) as df_input:    
		    df = json.load(df_input)
		return df

	def loadTsv(self, file, header):
		'''
			input: miranda file in tsv format
			return: data struc containing file input
		'''
		miranda_file = os.path.join(self.my_dir, file)
		miranda_input = open(miranda_file, "r")
		if header == 1:
			next(miranda_input)

		return miranda_input

	def getGeneID(self, transcript_id):
		'''
			input: unique lncRNA transcript ID: LNC_TR_mm10_...
			return: unique lncRNA gene ID: LNC_GE_mm10_... or None
		'''
		for key, value in df_mapping_gts.items():
			if transcript_id in value:
				return key
		return None

	def getMappedToUniqueID(self, lncrna_id):
		'''
			input: noncode lncRNA transcript ID: NONMMU... 
			return: unique lncRNA transcript ID: LNC_TR_mm10_... or None
		'''
		for key, value in df_mapping_ids.items():
			if value["alias"] == lncrna_id:
				return key
		return None

	def getExistingGene(self, gene_id, container):
		'''
			input: unique gene ID: LNC_GE_mm10_...
			return: dictionary entry containing gene ID or None if it does not exist
		'''
		for entry in container:
			if(entry and entry['gene_id'] == gene_id):
				return entry
		return None

	def verifyID(self, lncrna_id):
		'''
			input: gencode gene ID and alias: ENSMUSG00000115958.1|ENSMUST00000193888.5|UTR1
			return: gencode gene ID: ENSMUSG00000115958.1
		'''
		name_array = lncrna_id.split("|") # => ENSMUSG00000003344.14-Upstream-1000
		name_array2 = name_array[0].split("-") # => ENSMUSG00000003344.14
		#print(name_array2[0])
		return name_array2[0]

	def getMirandaInteractions(self, input, output):
		container = []
		for line in input:
			line_array = line.split("\t")
			mirna = line_array[0]
			transcript_id = self.verifyID(line_array[1])
			print(transcript_id)
			container_interactions = {}
			container_interactions['align_score'] = line_array[2]
			container_interactions['energy'] = line_array[3]
			container_interactions['mirna_start'] = line_array[4]
			container_interactions['mirna_end'] = line_array[5]
			container_interactions['lnc_start'] = line_array[6]
			container_interactions['lnc_end'] = line_array[7]
			container_interactions['align_len'] = line_array[8]
			container_interactions['mirna_iden'] = line_array[9]
			container_interactions['lncrna_iden'] = line_array[10]
			container_interactions['mirna_alignment'] = line_array[11]
			container_interactions['alignment'] = line_array[12]
			container_interactions['lncrna_alignment'] = line_array[13].rstrip()
	
			if self.type == 'lncrna':
				transcript_id = self.getMappedToUniqueID(transcript_id)

			if transcript_id:
				container_genes = {}
				container_transcripts = {}
				container_mirnas = {}
				gene_id = None
				gene_id = self.getGeneID(transcript_id)
				if gene_id:
					gene = None
					gene = self.getExistingGene(gene_id, container)
					if gene:	
						if gene['transcript_list']:
							new_transcript = 1
							for transcript in gene['transcript_list']:
								if(transcript and transcript['transcript_id'] == transcript_id):
									new_transcript = 0
									new_micro = 1
									for micro in transcript['interaction_list']:
										if(micro and micro['mirna'] == mirna):
											print("existing gene, existing transcript, existing mirna: " + gene_id + "\t" + transcript_id + "\t" + mirna)
											new_micro = 0
											micro['alignment_list'].append(container_interactions)
									if new_micro == 1:
										print("existing gene, existing transcript, new mirna: " + gene_id + "\t" + transcript_id + "\t" + mirna)
										container_mirnas['mirna'] = mirna
										container_mirnas['alignment_list'] = []
										container_mirnas['alignment_list'].append(container_interactions)
										transcript['interaction_list'].append(container_mirnas)

							if new_transcript == 1:
								print("existing gene, but new transcript, hence new mirna: " + gene_id + "\t" + transcript_id + "\t" + mirna)
								container_mirnas['mirna'] = mirna
								container_mirnas['alignment_list'] = []
								container_mirnas['alignment_list'].append(container_interactions)

								container_transcripts['transcript_id'] = transcript_id
								container_transcripts['interaction_list'] = []
								container_transcripts['interaction_list'].append(container_mirnas)

								gene['transcript_list'].append(container_transcripts)
					else:
						print("new gene, new transcript, new mirna: " + gene_id + "\t" + transcript_id + "\t" + mirna)
						container_genes['gene_id'] = gene_id
						container_genes['gene_type'] = self.type
						container_genes['organism'] = self.organism
						container_genes['source'] = self.source

						container_mirnas['mirna'] = mirna
						container_mirnas['alignment_list'] = []
						container_mirnas['alignment_list'].append(container_interactions)
					
						container_transcripts['transcript_id'] = transcript_id
						container_transcripts['interaction_list'] = []
						container_transcripts['interaction_list'].append(container_mirnas)

						container_genes['transcript_list'] = []
						container_genes['transcript_list'].append(container_transcripts)
					if container_genes:
						container.append(container_genes)
			else:
				print("ERROR: no transcript ID mapping")
		output.write(json.dumps(container, indent=4)) 
	
if __name__ == '__main__':
	input_file   			= sys.argv[1]
	input_source 			= sys.argv[2]
	input_type 				= sys.argv[3]
	input_organism			= sys.argv[4]
	input_genesTranscripts 	= sys.argv[5]
	input_mappingKeys		= sys.argv[6]

	output_file = input_file.replace(".tsv", ".json")
	fout = open(output_file, "w")

	mp = ParseMiranda(input_type, input_source, input_organism)

	df_mapping_ids = mp.loadJson(input_mappingKeys)
	df_mapping_gts = mp.loadJson(input_genesTranscripts)
	fin = mp.loadTsv(input_file, 1)

	mp.getMirandaInteractions(fin, fout)

##MOUSE##
# python3 parseMiranda.py ./mm10/Mirbase_mouse_noncode_filtered.tsv noncode lncrna 2 ./mm10/mm10_genesTranscripts_allDBs_and_pc.json mapping_keys.json
# python3 parseMiranda.py ./mm10/Mirbase_mouse_gencode_lncRNA_filtered_new.tsv gencode lncrna 2 ./mm10/mm10_genesTranscripts_allDBs_and_pc.json mapping_keys.json
# python3 parseMiranda.py ./mm10/Mirbase_mouse_gencode_pc_filtered_new.tsv gencode gene 2 ./mm10/mm10_genesTranscripts_allDBs_and_pc.json mapping_keys.json
# python3 parseMiranda.py ./mm10/Mirbase_mouse_Promotor.tsv all all 2 ./mm10/mm10_genesTranscripts_allDBs_and_pc.json mapping_keys.json

##HUMAN##
# python3 parseMiranda.py ./hg38/Mirbase_human_noncode_filtered.tsv noncode lncrna 1 ./hg38/hg38_genesTranscripts_allDBs_and_pc.json mapping_keys.json
# python3 parseMiranda.py ./hg38/Mirbase_human_lncipedia_filtered.tsv lncipedia lncrna 1 ./hg38/hg38_genesTranscripts_allDBs_and_pc.json mapping_keys.json
# python3 parseMiranda.py ./hg38/Mirbase_human_gencode_lncRNA_filtered_new.tsv gencode lncrna 1 ./hg38/hg38_genesTranscripts_allDBs_and_pc.json mapping_keys.json
# python3 parseMiranda.py ./hg38/Mirbase_human_gencode_pc_filtered_new.tsv gencode gene 1 ./hg38/hg38_genesTranscripts_allDBs_and_pc.json mapping_keys.json
	
