#!/bin/bash

input_file="input.csv"  # The file with split Gene Ontology IDs
file_to_modify="file2.csv"  # The file where replacements will be made
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
} < "$file_to_modify"

# Overwrite the original file with the modified content
mv "$temp_file" "$file_to_modify"
