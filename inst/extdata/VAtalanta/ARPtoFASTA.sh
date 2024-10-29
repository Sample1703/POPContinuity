#!/bin/bash

#---------------
#Mathias Currat | 30.07.2024
#script which transform ARP format (diploid SNP) to FASTA format
#---------------

# Define variables
start_pattern="SampleData= {"
end_pattern="}"

# Loop over all .arp files in the current directory
for file in *.arp; do
  # Check if there are any .arp files
  if [ -e "$file" ]; then

    #extract file name to reuse as ouput
    outputfilename=${file/.arp/.fasta}

    #---------------
    #Extract all the data for every samples
    #---------------

    # Perform an action with the file
    echo "Processing file: $file"
    # replace series of space by tabulation
    # sed -r 's/ +/\t/g' "$file" > "tabulated_$file"
    # Use sed to extract text between start_pattern and end_pattern
    sed -n "/$start_pattern/,/$end_pattern/{
      /$start_pattern/b
      /$end_pattern/b
      p
    }" "$file" > "data.1"
    grep -v '^[[:space:]]*$' "data.1" > "data.txt"


    # Initialize line counter
    line_number=1

    # Clear the temporay files if they exist
    > "odd_lines.1"
    > "odd_lines.2"
    > "odd_lines.3"
    > "odd_lines.4"
    > "even_lines.1"
    > "even_lines.2"
    > "even_lines.3"
    > "ID_list.txt"
    > "fastaH1.txt"
    > "fastaH2.txt"
    > $outputfilename

    #---------------
    #Create two files, one with all the odd lines countaining the individual id and one with all even lines without the individual ideas
    #---------------

    # Read the input file line by line
    while IFS= read -r line; do
      # Check if the line number is odd
      if (( line_number % 2 != 0 )); then
        echo "$line" >> "odd_lines.1"
      else
        echo "$line" | sed 's/^\t//' >> "even_lines.1"
      fi
      # Increment the line number
      ((line_number++))
    done < "data.txt"

    # Define the text pattern to add
    text_pattern1=">"

    # Remove the second column (number of individuals) and write to the output file
    awk '{$2=""; sub(/^ /, ""); print}' "odd_lines.1" | sed -r 's/ +/\t/g' > "odd_lines.2"

    # Use sed to add the text pattern before the first column of each line
    sed "s/^/$text_pattern1/" "odd_lines.2" > "odd_lines.3"

    # Define the text pattern to add
     text_pattern2="_H1:"
     text_pattern3="_H2:"

    # Use awk to add the text pattern 2 after the first column of each line in H1
     awk -v pattern="$text_pattern2" '{ $1 = $1pattern } 1' "odd_lines.3" > "odd_lines.4"

    # Read the input file H1 line by line
    while IFS= read -r line; do
    # Read the first column from the line
    ind_ID=$(echo "$line" | awk '{print $1}' | sed 's/_H1/_H2/' )
    echo $ind_ID >> ID_list.1

    #Remove some tabs
    cat "ID_list.1" | sed -r 's/\t//' > "ID_list.txt"
    cat "even_lines.1"  | sed -r 's/^\t//' > "even_lines.2"

    # add individual ID at the beginning of each line in H2
    paste "ID_list.txt" "even_lines.2" > "even_lines.3"

    # Use awk to print each column on a new line
       echo "$line" | awk '{for(i=2;i<=NF;i++) {print $1i-1; if ($i == 0) $i = "a"; else if ($i == 1) $i = "c"; print $i}}'  >> "fastaH1.txt"
    done < "odd_lines.4"

    # Read the input file H2 line by line
    while IFS= read -r line; do
    # Use awk to print each column on a new line
       echo "$line" | awk '{for(i=2;i<=NF;i++) {print $1i-1; if ($i == 0) $i = "a"; else if ($i == 1) $i = "c"; print $i}}'  >> "fastaH2.txt"
    done < "even_lines.3"

    #Merge both fasta files into a single final one, countaining the full dataset.
    cat  "fastaH1.txt" "fastaH2.txt" > $outputfilename

  else
    echo "No .arp files found."
    break
  fi

  # clean unecessary files
  #   rm "tabulated_$file"
     rm data.*
     rm odd_lines.*
     rm even_lines.*
     rm fastaH*.txt
     rm ID_list*
done


