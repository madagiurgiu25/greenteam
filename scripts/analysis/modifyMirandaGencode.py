import sys
import os

def generateModifiedFile(file_in, file_out):
	for line in file_in:
		line_array = line.split("\t")
		mirna = line_array[0]
		align_score = line_array[2]
		energy = line_array[3]
		mirna_start = line_array[4]
		mirna_end = line_array[5]
		lnc_start = line_array[6]
		lnc_end = line_array[7]
		align_len = line_array[8]
		mirna_iden = line_array[9]
		lncrna_iden = line_array[10]
		mirna_alignment = line_array[11]
		alignment = line_array[12]
		lncrna_alignment = line_array[13].rstrip()

		name_array = line_array[1].split("|")
		lncrna = name_array[0]

		file_out.write(mirna + "\t" + lncrna + "\t" + align_score + "\t" + energy + "\t" + mirna_start + "\t" + mirna_end + "\t" + lnc_start + "\t" + lnc_end
			 + "\t" + align_len + "\t" + mirna_iden + "\t" + lncrna_iden + "\t" + mirna_alignment + "\t" + alignment + "\t" + lncrna_alignment + "\n")

if __name__ == '__main__':
	fin = sys.argv[1]
	fout = sys.argv[2]

	f_in = open(fin, "r")
	f_out = open(fout, "w")
	f_out.write(next(f_in))

	generateModifiedFile(f_in, f_out)

# python3 modifyMirandaGencode.py mm10/Mirbase_mouse_gencode_pc_filtered.tsv mm10/Mirbase_mouse_gencode_pc_filtered_new.tsv
# python3 modifyMirandaGencode.py mm10/Mirbase_mouse_Promotor.tsv mm10/Mirbase_mouse_Promotor_new.tsv
# python3 modifyMirandaGencode.py mm10/Mirbase_mouse_Promotor.tsv mm10/Mirbase_mouse_Promotor_new.tsv