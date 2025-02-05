
import os
import sys
import json
import requests
import argparse
import datetime
import pandas as pd
import session_info

def find_path(path, action, path_type):
    """
    Returns full path of inputs based on current working directory.

        Parameters:
            path (str):   Full or relative path based on current working directory.
                          "./", "../", etc. are acceptable.
            action (str): Either "r" (for inputs) or "w" (for outputs).
            path_type (str): Either "f" (for files) or "d" (for directories).

        Returns:
            abspath (str): Full path of the given input.
    """

    # abspath() removes "/" at the end.
    abspath = os.path.abspath(path).replace("\\", "/")
    if path_type == "d":
        abspath += "/" # Add "/" so that the following will work.

    # In addition to finding the above abspath, "r" and "w" will have other functions.
    # Currently, only "r"+"f", "w"+"d", and "w"+"f" have extra functionality.
    if action == "w":
        # Create parent directories if they don't exist.
        parents = "/".join(abspath.split("/")[:-1])
        if not os.path.isdir(parents):
            print(f"Output directory does not exist. Creating parent(s): {parents}\n", flush=True)
            os.makedirs(parents)
    elif action == "r":
        if path_type == "d":
            # TODO: Add read directory functionality...
            print("Reading directories in does not currently have a functionality", flush=True)
            pass
        elif path_type == "f":
            # Exit if file to be read doesn't exist.
            if not os.path.isfile(abspath):
                print(f"No file found at {abspath}. Exiting.\n", flush=True)
                sys.exit() # check if working            

    return abspath

def fasta_to_df(file):
    """
    Takes in a FASTA-formatted file and returns a pandas DataFrame of the input sequences.

        Parameters:
            file (str): Full path to a FASTA-formatted file

        Returns:
            fasta_df (pandas.DataFrame): pandas DataFrame with short form UniProt ID indexes
                                         columns=[Accession, IDs, Sequence]
    """

    with open(file, "r", encoding="utf-8") as fasta:
        accessions = []
        ids = []
        seqs = []
        for line in fasta:
            if line[0] == ">":
                accessions.append(line.strip())
                # get rid of ">" and remove "|" if they exist (for uniprot entries)
                ids.append(line.strip().split(" ")[0][1:].split("|")[-1])
                seqs.append("")
            else:
                seqs[-1] += line.strip()
    return pd.DataFrame({"Accession": accessions, "IDs": ids, "Sequence": seqs}).set_index("IDs")



def parse_args():
    """
    Takes in command-line arguments and returns an argparse Namespace object.

        Returns:
            arguments (Namespace): Namespace with command-line arguments.
    """

    parser = argparse.ArgumentParser(description="FASTA Splitter")

    parser.add_argument("-n", "--num_files", type=str,help="Number of files after split.")
    parser.add_argument("-f", "--fasta", type=str,
                        help="Full path to sequence file in FASTA format")
    parser.add_argument("-o", "--out_directory", type=str,
                        help="Full path of output directory.")

    return parser.parse_args()

#def write_log(log_path, string):
  
    #with open(log_path, "a", encoding="utf-8") as log:
        #log.write(string)

def main(args):
    """
    Splits a FASTA file into multiple parts.
    """
    
    current_time = str(datetime.datetime.now()).replace("-", "").replace(" ", "_").replace(":", "")
    out_directory = find_path(args.out_directory, action="w", path_type="d")
    script_name = ".".join(__file__.split("/")[-1].split(".")[:-1])
    log_path = find_path(f"{out_directory}/{script_name}.{current_time}.log",
                         action="w", path_type="f")
    # This is a little sketchy, but it redirects all stdout (including session_info) to log_path.
    # There is no other use of log_path here... consider doing something like this in other files...
    # Maybe ask user if stdout should be redirected or not.
    sys.stdout = open(log_path, 'w')

    session_info.show() # Prints to stdout
    print("\n")

    fasta_file = find_path(args.fasta, action="r", path_type="f")
    fasta_filename = fasta_file.split("/")[-1]
    fasta_filename_noext = ".".join(fasta_filename.split(".")[:-1])
    fasta_df = fasta_to_df(fasta_file)

    num_files = int(args.num_files)

    out_directory = find_path(args.out_directory, action="w", path_type="d")
    print("\nInputs:\n", flush=True) # Prints a newline
    print(f"fasta: {fasta_file}\nnum_files:{num_files}", flush=True)
    print(f"out_directory: {out_directory}\n", flush=True)    
    
    
    total_seqs = len(fasta_df)
    base_seqs = int(total_seqs / num_files)
    remainder = total_seqs % num_files
    
    # Initialize out files and save ending indices 
    out_filename = f"{out_directory}/{fasta_filename_noext}_1.txt" # start with 1
    with open(out_filename, "w", encoding="utf-8") as out:
        out.write("")
    files = {1:base_seqs + min(1, remainder)}
    remainder -= min(remainder, 1)
    for file_num in range(2, num_files+1):
        out_filename = f"{out_directory}/{fasta_filename_noext}_{file_num}.txt"
        with open(out_filename, "w", encoding="utf-8") as out:
            out.write("")
        files[file_num] = files[file_num-1] + base_seqs + min(1, remainder)
        remainder -= min(remainder, 1) # If remainder is 0, will subtract 0
    
    file_num = 1
    for i in range(total_seqs):
        if i >= files.get(file_num):
            file_num += 1
        out_filename = f"{out_directory}/{fasta_filename_noext}_{file_num}.txt"
        acc = fasta_df.iloc[i]["Accession"]
        seq = fasta_df.iloc[i]["Sequence"]
        with open(out_filename, "a", encoding="utf-8") as out:
            out.write(f"{acc}\n{seq}\n")
            
if __name__ == "__main__":
    main(parse_args())
