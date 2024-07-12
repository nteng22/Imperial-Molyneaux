#### I ran this locally on my machine using the Terminal. 
### Not necessary to do this on the HPC but it can be done. 

# Essentially need to create a Manifest file with the sample ID, forward and reverse reads.
# Tutorial: https://docs.qiime2.org/2024.2/tutorials/importing/, go to: "“Fastq manifest” formats"

sn -l "path to the fastq files" . # create a symlink to my folder
readlink -f fastq/david/* > files_paths.txt # write out the paths to the fastq files
while IFS= read -r f; do basename "$f" | cut -d'_' -f1 > sampleID.txt; done < files_paths.txt
	# Get the sample ID names of all samples without the other info
uniq sampleID.txt > sampleID_unique.txt # remove the duplicates as the files are R1 and R2. 
wc -l sampleID_unique.txt # see how many sample IDs there are
grep "R1" files_paths.txt > R1_list.txt # pull out absolute path for R1 files
wc -l R1_list.txt # check number of samples
grep "R2" files_paths.txt > R2_list.txt # pull out absolute paths for R2 files
wc -l R2_list.txt # check number of samples
grep -Fvf sampleID_unique.txt R1_list.txt > name_check.txt 
	# Check if sampleID_unique is present in R1 list. Can repeat with R2 but expect they're the same.

paste sampleID_unique.txt R1_list.txt R2_list.txt > Manifest-map.txt # add all columns together
	# first column is sampleID, followed by the path to R1 and path to R2. 
awk -v FS="\t" '{ if (index($2, $1) && index($3, $1)) print $0 }' Manifest-map.txt > newfile
	# Check if each row corresponds to the same sample ID. 
wc -l Manifest-map.txt # need this to be the same, if it isn't suggests one sample went awry. 

awk -vOFS="\t" '$1=$1; BEGIN { str="sample-id forward-absolute-filepath reverse-absolute-filepath"; split(str,arr," "); for(i in arr) printf("%s\t", arr[i]);print}' Manifest-map.txt > NorthUmbria_2405-Map.txt	
	# need to add column heading. 
wc -l NorthUmbria_2405-Map.txt # check +1 number of rows from the previous lists. 