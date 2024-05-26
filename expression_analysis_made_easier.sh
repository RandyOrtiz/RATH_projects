#!/bin/bash

##### Good Merged Assembly to Differential Expression Analysis Results
##### Make sure your folders are organized correctly with raw reads and assembly in each folder, all folders in the same working directory as well
##### Usage: bash expression_analysis_made_easier.sh --go 0006955 --tax 8342
##### Place .sh and .py scripts in working directory


# Paths to software, Change this to reflect your paths
rule_diamond='/home/velox/Documents/Pincho_v01_Master/bin/diamond-linux64/diamond'
rule_dedupe='/home/velox/Documents/Pincho_v01_Master/bin/BBMap_38.86/bbmap/dedupe.sh'
rule_salmon='/home/velox/Documents/Pincho_v01_Master/bin/salmon-1.10.0_linux_x86_64/salmon-latest_linux_x86_64/bin/salmon'

# File names, Change this to reflect the names of your files
rule_assembly='merged_assemblies_busco_mod.fasta'
rule_paired_read_1='R1_001.fastq.gz'
rule_paired_read_2='R2_001.fastq.gz'
rule_unpaired_read='.fastq.gz'
rule_cwd=$(pwd)

### Download Reference Proteome from Uniprot using GO and TAXID
# Function to display usage message
usage() {
    echo "Usage: $0 --go <GO_TERM> --tax <TAX_ID>"
    exit 1
}

# Check if the correct number of arguments is provided
if [ "$#" -ne 4 ]; then
    usage
fi

# Parse command line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --go)
            GO_TERM="$2"
            shift 2
            ;;
        --tax)
            TAX_ID="$2"
            shift 2
            ;;
        *)
            usage
            ;;
    esac
done

# Check if GO_TERM and TAX_ID are set
if [ -z "$GO_TERM" ] || [ -z "$TAX_ID" ]; then
    usage
fi

# Construct the URL
URL1="https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28%28go%3A${GO_TERM}%29+AND+%28taxonomy_id%3A${TAX_ID}%29%29"
URL2="https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession%2Cgo_id&format=tsv&query=%28%28go%3A${GO_TERM}%29+AND+%28taxonomy_id%3A${TAX_ID}%29%29"
URL3="https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession%2Cgo&format=tsv&query=%28%28go%3A${GO_TERM}%29+AND+%28taxonomy_id%3A${TAX_ID}%29%29"
URL4="https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession%2Cgo_p&format=tsv&query=%28%28go%3A${GO_TERM}%29+AND+%28taxonomy_id%3A${TAX_ID}%29%29"
URL5="https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession%2Cgo_f&format=tsv&query=%28%28go%3A${GO_TERM}%29+AND+%28taxonomy_id%3A${TAX_ID}%29%29"

# Define the output file name
rule_reference_proteome="Uniprot_Reference_Proteome_GO${GO_TERM}_TAX${TAX_ID}.fasta.gz"
rule_reference_proteome_GO_table="Uniprot_GO_Table_GO${GO_TERM}_TAX${TAX_ID}.tsv.gz"
rule_reference_proteome_GO_name_table="Uniprot_GO_name_Table_GO${GO_TERM}_TAX${TAX_ID}.tsv.gz"
rule_reference_proteome_GO_bio_table="Uniprot_GO_bio_Table_GO${GO_TERM}_TAX${TAX_ID}.tsv.gz"
rule_reference_proteome_GO_mol_table="Uniprot_GO_mol_Table_GO${GO_TERM}_TAX${TAX_ID}.tsv.gz"












































### Download the uniprot reference proteome
curl -o "$rule_reference_proteome" "$URL1"

# Verify the download was successful
if [ $? -ne 0 ]; then
    echo "Download failed. Please check your internet connection and try again."
    exit 1
fi

echo "Download successful. File saved as $rule_reference_proteome."

# For gzip (.gz) files
gunzip "$rule_reference_proteome"
rule_reference_proteome="Uniprot_Reference_Proteome_GO${GO_TERM}_TAX${TAX_ID}.fasta"


# Replace space with underscores in reference proteome headers
sed -i 's/ /_/g' "$rule_reference_proteome"







### Download the GO table
curl -o "$rule_reference_proteome_GO_table" "$URL2"

# Verify the download was successful
if [ $? -ne 0 ]; then
    echo "Download failed. Please check your internet connection and try again."
    exit 1
fi

echo "Download successful. File saved as $rule_reference_proteome_GO_table."

# For gzip (.gz) files
gunzip "$rule_reference_proteome_GO_table"
rule_reference_proteome_GO_table="Uniprot_GO_Table_GO${GO_TERM}_TAX${TAX_ID}.tsv"


### Modify GO table

temp_file=$(mktemp)

# Create the temporary file and add headers
echo -e "Entry\tGene Ontology IDs" > "$temp_file"

# Read the input file line by line
while IFS=$'\t' read -r entry go_ids; do
    # Skip the header line
    if [ "$entry" == "Entry" ]; then
        continue
    fi
    # Split the GO IDs by semicolon and create new rows
    IFS=';' read -ra go_id_array <<< "$go_ids"
    for go_id in "${go_id_array[@]}"; do
        echo -e "$entry\t$go_id" >> "$temp_file"
    done
done < "$rule_reference_proteome_GO_table"

# Overwrite the input file with the contents of the temporary file
mv "$temp_file" "$rule_reference_proteome_GO_table"
awk -F '\t' 'BEGIN {OFS="\t"} {sub(/^ /, "", $2); print}' $rule_reference_proteome_GO_table > 'temp_file.txt' && mv 'temp_file.txt' $rule_reference_proteome_GO_table





### Download the GO table
curl -o "$rule_reference_proteome_GO_name_table" "$URL3"

# Verify the download was successful
if [ $? -ne 0 ]; then
    echo "Download failed. Please check your internet connection and try again."
    exit 1
fi

echo "Download successful. File saved as $rule_reference_proteome_GO_name_table."

# For gzip (.gz) files
gunzip "$rule_reference_proteome_GO_name_table"
rule_reference_proteome_GO_name_table="Uniprot_GO_name_Table_GO${GO_TERM}_TAX${TAX_ID}.tsv"


### Modify GO table

temp_file=$(mktemp)

# Create the temporary file and add headers
echo -e "Entry\tGene Ontology (GO)" > "$temp_file"

# Read the input file line by line
while IFS=$'\t' read -r entry go_ids; do
    # Skip the header line
    if [ "$entry" == "Entry" ]; then
        continue
    fi
    # Split the GO IDs by semicolon and create new rows
    IFS=';' read -ra go_id_array <<< "$go_ids"
    for go_id in "${go_id_array[@]}"; do
        echo -e "$entry\t$go_id" >> "$temp_file"
    done
done < "$rule_reference_proteome_GO_name_table"

# Overwrite the input file with the contents of the temporary file
mv "$temp_file" "$rule_reference_proteome_GO_name_table"
awk -F '\t' 'BEGIN {OFS="\t"} {sub(/^ /, "", $2); print}' $rule_reference_proteome_GO_name_table > 'temp_file.txt' && mv 'temp_file.txt' $rule_reference_proteome_GO_name_table






### Download the GO table
curl -o "$rule_reference_proteome_GO_bio_table" "$URL4"

# Verify the download was successful
if [ $? -ne 0 ]; then
    echo "Download failed. Please check your internet connection and try again."
    exit 1
fi

echo "Download successful. File saved as $rule_reference_proteome_GO_bio_table."

# For gzip (.gz) files
gunzip "$rule_reference_proteome_GO_bio_table"
rule_reference_proteome_GO_bio_table="Uniprot_GO_bio_Table_GO${GO_TERM}_TAX${TAX_ID}.tsv"


### Modify GO table

temp_file=$(mktemp)

# Create the temporary file and add headers
echo -e "Entry\tGene Ontology (biological process)" > "$temp_file"

# Read the input file line by line
while IFS=$'\t' read -r entry go_ids; do
    # Skip the header line
    if [ "$entry" == "Entry" ]; then
        continue
    fi
    # Split the GO IDs by semicolon and create new rows
    IFS=';' read -ra go_id_array <<< "$go_ids"
    for go_id in "${go_id_array[@]}"; do
        echo -e "$entry\t$go_id" >> "$temp_file"
    done
done < "$rule_reference_proteome_GO_bio_table"

# Overwrite the input file with the contents of the temporary file
mv "$temp_file" "$rule_reference_proteome_GO_bio_table"
awk -F '\t' 'BEGIN {OFS="\t"} {sub(/^ /, "", $2); print}' $rule_reference_proteome_GO_bio_table > 'temp_file.txt' && mv 'temp_file.txt' $rule_reference_proteome_GO_bio_table









### Download the GO table
curl -o "$rule_reference_proteome_GO_mol_table" "$URL5"

# Verify the download was successful
if [ $? -ne 0 ]; then
    echo "Download failed. Please check your internet connection and try again."
    exit 1
fi

echo "Download successful. File saved as $rule_reference_proteome_GO_mol_table."

# For gzip (.gz) files
gunzip "$rule_reference_proteome_GO_mol_table"
rule_reference_proteome_GO_mol_table="Uniprot_GO_mol_Table_GO${GO_TERM}_TAX${TAX_ID}.tsv"


### Modify GO table

temp_file=$(mktemp)

# Create the temporary file and add headers
echo -e "Entry\tGene Ontology (molecular function)" > "$temp_file"

# Read the input file line by line
while IFS=$'\t' read -r entry go_ids; do
    # Skip the header line
    if [ "$entry" == "Entry" ]; then
        continue
    fi
    # Split the GO IDs by semicolon and create new rows
    IFS=';' read -ra go_id_array <<< "$go_ids"
    for go_id in "${go_id_array[@]}"; do
        echo -e "$entry\t$go_id" >> "$temp_file"
    done
done < "$rule_reference_proteome_GO_mol_table"

# Overwrite the input file with the contents of the temporary file
mv "$temp_file" "$rule_reference_proteome_GO_mol_table"
awk -F '\t' 'BEGIN {OFS="\t"} {sub(/^ /, "", $2); print}' $rule_reference_proteome_GO_mol_table > 'temp_file.txt' && mv 'temp_file.txt' $rule_reference_proteome_GO_mol_table




### Setup Query
for dir in */; do
cd $dir

# Annotated file
query=$(echo *$rule_assembly)

# Output FASTA file
sorted_query="sorted_${query}"

# Create a temporary file for sequences with their lengths
temp_file=$(mktemp)

# Read the FASTA file and process sequences
awk '/^>/ {if (seqlen) {print seqlen "\t" seqheader "\t" seq} seqheader=$0; seq=""; seqlen=0; next} {seq=seq""$0; seqlen+=length($0)} END {print seqlen "\t" seqheader "\t" seq}' $query > $temp_file

# Sort by sequence length in descending order and format output
sorted_sequences=$(sort -k1,1nr $temp_file | awk -F'\t' '{print $2"\n"$3}')

# Save the output to a file
echo "$sorted_sequences" > $sorted_query

# Clean up
rm $temp_file

# Print the output file path
echo "Sorted sequences saved to $sorted_query"

### Setup Input
echo $sorted_query
file_count=$(find . -name "*.gz" | wc -l)
if [[ $file_count -eq 2 ]]
then
	read1=$(echo *$rule_paired_read_1)
	read2=$(echo *$rule_paired_read_2)
else
	read1=$(echo *$rule_unpaired_read)
fi
echo $read1
echo $read2

rule_reference_proteome_dir="$rule_cwd/$rule_reference_proteome"

# Diamond Blastx Query against Reference Proteome
$rule_diamond blastx -d $rule_reference_proteome_dir -q $sorted_query -o cc_matches -e 0.0000000001 --outfmt 6 -b5 -c1 -p 24

#diamond will align and list the ranges of the alignments

#Get basenames for file naming
basename_reference=${sorted_query##*/}
basename_reference_noext=${basename_reference%.*}
basename_db=${rule_reference_proteome_dir##*/}
basename_db_noext=${basename_db%.*}

#Use coordinates to trim reconstructions, hence splitting chimeras
#We can use this for atelopus and coqui projects
#Col 1 = SeqID, Col 2 = UniprotID, Col 7 start char, Col 8 end char, Col 3-6,9-12 appended blast data

#Find SeqID (Col 1) and pull ID and sequence from annotation
#create list of OldID
sed -i 's/^/>/g' cc_matches
cat $sorted_query | paste - - > single_line_anno
awk 'NR==FNR{a[$1]=$2;next}{print $0,a[$1]?a[$1]:"NA"}' single_line_anno cc_matches > diamond_fasta 

#remove empty endline first ADD

#trim sequence with seq nums
while read c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 c13; do
echo ${c13} | cut -c${c7}-${c8} 2>/dev/null || echo ${c13} | cut -c${c8}-${c7} 2>/dev/null
done < diamond_fasta > diamond_fasta_col13

#remove old sequences
awk 'NF{NF-=1};1' diamond_fasta > diamond_fasta_col12

#add new sequences
paste diamond_fasta_col12 diamond_fasta_col13 > diamond_fasta_complete_table

#remove first column
cut -d " " -f2- diamond_fasta_complete_table > diamond.fasta

#add >
sed -i 's/^/>/g' diamond.fasta

#convert table to fasta
sed -i -r 's/(.*)\t/\1\n/' diamond.fasta 

#remove dups and linearize
$rule_dedupe in=diamond.fasta out=${basename_db_noext}_${basename_reference_noext}_results.fasta ac=f fastawrap=0

rm cc_matches single_line_anno diamond_fasta diamond_fasta_col13 diamond_fasta_col12 diamond_fasta_complete_table diamond.fasta

sed -i 's/ /__/g' ${basename_db_noext}_${basename_reference_noext}_results.fasta
$rule_salmon index -i transcripts_index -t ${basename_db_noext}_${basename_reference_noext}_results.fasta
if [[ $file_count -eq 2 ]]
then
	$rule_salmon quant -i transcripts_index -l A --validateMappings -o transcripts_quant -p 24 -1 $read1 -2 $read2
else
	$rule_salmon quant -i transcripts_index -l A --validateMappings -o transcripts_quant -p 24 -r $read1
fi

#remove blast details from expression file, quant.sh
awk '{sub(/\__.*$/,"",$1); print $1"\t"$2"\t"$3"\t"$4"\t"$5}' transcripts_quant/quant.sf > ${basename_db_noext}_${basename_reference_noext}_quantmod.sf
rm -r transcripts_index transcripts_quant
cd ..
done




#awk '$5' #remove any row with a 0 for num reads aligned


#before making expression sheet, need db of all ids.
for quant in */*quantmod.sf; do
#prep row titles
cat $quant >> merged_expression_file
done
cut -f1 merged_expression_file > merged_expression_file2
sort -u merged_expression_file2 > merged_expression_file2_uniq
rm merged_expression_file merged_expression_file2





for dir2 in */; do
realdir=$(echo $dir2)
cp merged_expression_file2_uniq ${realdir}merged_expression_file2_uniq
cd $dir2
#prep row data
quant=$(echo *quantmod.sf)
tail -n +2 $quant | sort -k1,1 -k4,4rn | sort -uk1,1 > quant_nodup.sf
cut -f1,4 quant_nodup.sf > quant2.sf
rm quant_nodup.sf
cd ..
done




#use list1 to find matches in list2, append list 2 to list 1
#keep highest tmp value remove dups
for dir3 in */; do
cd $dir3
join -a1 <(sort merged_expression_file2_uniq) <(sort quant2.sf) > tmp.sf
sed -i '/ /! s/$/ 0/' tmp.sf
awk -v OFS="\t" '$1=$1' tmp.sf > table_parts
rm tmp.sf
((counter+=1))
clean_dir3="${dir3%/}"
sed -i "0,/0/s//${clean_dir3}/" table_parts
if [ $counter -ne 1 ]; then
cut -f2 table_parts > tmp2.sf
cp -fr tmp2.sf table_parts
rm tmp2.sf
else
sed -i "s/[|,=:'()#\/+*;\"[\]?<>-]/_/g" table_parts
fi
cd ..
done

find */ -name 'table_parts' -exec paste {} + > DE_table

find . -type f -name 'table_parts' -delete
find . -type f -name 'merged_expression_file2_uniq' -delete
find . -type f -name 'quant2.sf' -delete
find . -type f -name '*quantmod.sf' -delete

#Create metafile


# Define input and output file names
output_file="metafile"

# Create the header for the output file
echo -e "sample\ttreatment\tbatch" > $output_file

# Initialize the batch counter
treatment_type=0
batch_number=0

# Extract the first row of the input file and transpose it to the first column of the output file
awk -F'\t' 'NR==1 { for(i=1; i<=NF; i++) print $i }' 'DE_table' | while read sample; do
  echo -e "$sample\ttreatment_$treatment_type\tbatch_$batch_number" >> $output_file
  ((treatment_type++))
  ((batch_number++))
done

# Remove the 2nd row from the input file and save the result to a temporary file
awk 'NR!=2' $output_file > $temp_file

# Replace the original file with the temporary file
mv $temp_file $output_file

echo "Output file created: $output_file"


### GO Analysis w/ Heatmap
### Create table with expression values, species, and GO terms USING NEW FORMAT GO TABLE(not genes)

rule_reference_proteome="Uniprot_Reference_Proteome_GO${GO_TERM}_TAX${TAX_ID}.fasta"
rule_reference_proteome_GO_table="Uniprot_GO_Table_GO${GO_TERM}_TAX${TAX_ID}.tsv"
rule_reference_proteome_GO_name_table="Uniprot_GO_name_Table_GO${GO_TERM}_TAX${TAX_ID}.tsv"
rule_reference_proteome_GO_bio_table="Uniprot_GO_bio_Table_GO${GO_TERM}_TAX${TAX_ID}.tsv"
rule_reference_proteome_GO_mol_table="Uniprot_GO_mol_Table_GO${GO_TERM}_TAX${TAX_ID}.tsv"

# List of input and output files
rule_go_file_list=($rule_reference_proteome_GO_table $rule_reference_proteome_GO_name_table $rule_reference_proteome_GO_bio_table $rule_reference_proteome_GO_mol_table)
rule_out_file_list=("Exp_Results_Uniprot_GO_Table_GO${GO_TERM}_TAX${TAX_ID}.txt" "Exp_Results_Uniprot_GO_name_Table_GO${GO_TERM}_TAX${TAX_ID}.txt" "Exp_Results_Uniprot_GO_bio_Table_GO${GO_TERM}_TAX${TAX_ID}.txt" "Exp_Results_Uniprot_GO_mol_Table_GO${GO_TERM}_TAX${TAX_ID}.txt")


for i in {0..3}; do
    input_file="${rule_go_file_list[$i]}"
    final_output="${rule_out_file_list[$i]}"
    temp_file=$(mktemp)  # Temporary file to store the modified content

    # Read the input file and create an associative array
    declare -A go_dict
    while IFS=$'\t' read -r entry go_id; do
        if [ "$entry" == "Entry" ]; then
            continue
        fi
        if [ -n "${go_dict[$entry]}" ]; then
            go_dict["$entry"]="${go_dict[$entry]};$go_id"
        else
            go_dict["$entry"]="$go_id"
        fi
    done < "$input_file"

    # Read the file to modify line by line
    {
        read -r header  # Read and keep the header
        echo "$header" > "$temp_file"

        while IFS=$'\t' read -r name value1 value2; do
            modified_line="$name"
            entry=$(echo "$modified_line" | awk -F'|' '{print $2}')
            if [ -n "${go_dict[$entry]}" ]; then
                go_ids="${go_dict[$entry]}"
                IFS=';' read -ra go_id_array <<< "$go_ids"
                for go_id in "${go_id_array[@]}"; do
                    echo -e "$go_id\t$value1\t$value2" >> "$temp_file"
                done
            else
                echo -e "$name\t$value1\t$value2" >> "$temp_file"
            fi
        done
    } < "DE_table"

    # Filter out rows that do not contain any number besides 0 in the second and third columns (excluding the first row and first column)
    awk 'BEGIN{FS=OFS="\t"} NR==1 || $2 ~ /[1-9]/ || $3 ~ /[1-9]/' "$temp_file" > "$temp_file.filtered"

    # Remove rows that do not contain 'GO:' (excluding the header)
    awk 'BEGIN{FS=OFS="\t"} NR==1 || /GO:/' "$temp_file.filtered" > "$temp_file.filtered.go"

    # Sort the file by the first column, excluding the header
    { head -n 1 "$temp_file.filtered.go"; tail -n +2 "$temp_file.filtered.go" | sort -k1,1; } > "$temp_file.sorted"

    # Merge rows with the same header and add their values
    awk 'BEGIN{FS=OFS="\t"} 
    NR==1 {header=$0; next} 
    {
        key=$1;
        for (i=2; i<=NF; i++) {
            sum[key, i] += $i;
        }
        headers=NF;
    }
    END {
        print header;
        for (key in sum) {
            split(key, keys, SUBSEP);
            if (!seen[keys[1]]) {
                seen[keys[1]] = 1;
                row=keys[1];
                for (i=2; i<=headers; i++) {
                    row = row OFS sum[keys[1], i];
                }
                print row;
            }
        }
    }' "$temp_file.sorted" > "$temp_file.merged"

    # Sort the merged output again, excluding the header
    { head -n 1 "$temp_file.merged"; tail -n +2 "$temp_file.merged" | sort -k1,1; } > "$final_output"

    # Clear the associative array for the next iteration
    unset go_dict

    # Remove temporary files
    rm "$temp_file" "$temp_file.filtered" "$temp_file.filtered.go" "$temp_file.sorted" "$temp_file.merged"
done




### String Analysis
### Create Unique list of gene ids to use with string per species























### Generate heatmap

# Ensure all files exist
for file in "${rule_out_file_list[@]}"; do
    if [ ! -f $file ]; then
        echo "$file not found!"
        exit 1
    fi
done

# Run the Python script to generate the heatmaps
python3 create_heatmaps.py "${rule_out_file_list[@]}"

# Check if the heatmaps were created
for file in "${rule_out_file_list[@]}"; do
    image_file="${file/.txt/.png}"
    if [ -f $image_file ]; then
        echo "Heatmap generated successfully: $image_file"
    else
        echo "Failed to generate heatmap for $file."
    fi
done

# Close the terminal
exit









#############################



