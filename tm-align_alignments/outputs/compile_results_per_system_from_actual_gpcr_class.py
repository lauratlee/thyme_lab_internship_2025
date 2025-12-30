#the purpose of this script is to take the key file for gpcrs by class and pull the data for that receptor for its class and make a spreadsheet


#imports
import os,sys

#open a write file to write everything to
write_file = open("pocket_similarity_and_identity_actual_gpcr_class.csv", "w")

#write a header line
write_file.write("human gene,zebrafish gene,actual class,matched,similarity,mismatched,failed alignments,% failed alignments,% identity,% similarity\n")

#read the key file
key_file = open("good_receptor_gpcr_type_key.csv", "r")
for line in key_file.readlines():
	if line.startswith("human gene"):
		continue
	print(line)

	#set defaults for values in case we get nothing
	hgene = line.split(",")[0]
	zgene = line.split(",")[1]
	act_class = line.split(",")[2].split("\n")[0]
	matched = "0"
	similarity = "0"
	mismatched = "0"
	failed_alignments = "0"
	pct_failed_alignments = "0"
	pct_identity = "0"
	pct_similarity = "0"

	#read through the file that corresponds to the class
	for r,d,f in os.walk(os.getcwd()):
		for file in f:
			if file == "[" + act_class + "]human_zebrafish_alignment_summary_2.0_with_similarity.csv":
				#read the file and get the data
				read_file = open(file, "r")
				for line2 in read_file.readlines():
					#try to find a match based on the humand and fis hgenes
					if line2.split(",")[0] == hgene and line2.split(",")[1] == zgene:
						matched = str(line2.split(",")[2])
						similarity = str(line2.split(",")[3])
						mismatched = str(line2.split(",")[4])
						failed_alignments = str(line2.split(",")[5])

						pct_failed_alignments = str(100 * float(failed_alignments) / (float(failed_alignments) + float(mismatched) + float(similarity)))

						pct_identity = str(line2.split(",")[6])
						pct_similarity = str(line2.split(",")[7])

						#write the line
						write_file.write(hgene + "," + zgene + "," + act_class + "," + matched + "," + similarity + "," + mismatched + "," + failed_alignments + "," + pct_failed_alignments + "," + pct_identity + "," + pct_similarity + "\n")