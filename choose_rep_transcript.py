
"""
For tab-delimited tables output by assign_EC_only.py:
Finds coverage for each EC number for each transcript.
Assigns an EC number to each transcript based on EC coverage.
Assigns a transcript to each gene based on EC coverage.
"""

import os
import sys
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

def parse_args():
    """
    Takes in command-line arguments and returns an argparse Namespace object.

        Returns:
            arguments (Namespace): Namespace with command-line arguments.
    """

    parser = argparse.ArgumentParser(description="EC \"Predictor\"")
    # -i and -o must come before tool_name.
    # Ex. -i ./infile.fasta -o ./out_dir annotate

    parser.add_argument("-i", "--infile", type=str,
                        help="Full path to input TSV from assign_EC_only.py")
    parser.add_argument("-o", "--out_directory", type=str,
                        help="Full path of output directory. Must end with \"/\".")

    return parser.parse_args()

#def write_log(log_path, string):
  
    #with open(log_path, "a", encoding="utf-8") as log:
        #log.write(string)

def main(args):
    """
    Chooses a representative transcript and EC number(s) for each gene.
    Multiple EC numbers will be chosen if they appear in BLAST hits the same
    number of times (EC coverage).

        Outputs:
            TSV of genes, transcripts and EC numbers
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
    
    infile_file = find_path(args.infile, action="r", path_type="f")
    infile_filename = infile_file.split("/")[-1]
    infile_df = pd.read_csv(infile_file, sep="\t")
    print(infile_df)
    
    # Order by percent identity in order to break ties
    # between transcripts with the same EC number.
    # Greatest pident value will be the representative transcript.
    # Other ties (unlikely) will be broken by scover, then qcover.
    #infile_df = infile_df.sort_values(by=["gene", "pident", "scover", "qcover"], ascending=False)
    #print(infile_df)
    
    # genes is a dict of dicts of dicts of ints.
    # genes: key=gene, value=transcripts (dict)
    # transcripts: key=transcript, value=ec_counts (dict)
    # ec_counts: key=ec, value=count (int)
    genes = {}
    
    # Listen, every time I want to iterate over a pandas df,
    # I get scolded by strangers on the internet.
    # Future me should read this for the tenth time and decide
    # how much she really cares:
    # https://stackoverflow.com/questions/16476924/how-can-i-iterate-over-rows-in-a-pandas-dataframe
    # https://stackoverflow.com/questions/16476924/how-can-i-iterate-over-rows-in-a-pandas-dataframe/77270285#77270285
    # Here, itertuples seems fine and I only rarely come back to this script.
    # Maybe just don't ever use iterrows.
    for row in infile_df.itertuples(index=False):
        
        gene = row.gene
        transcript = row.transcript
        # Turns string list representation into an actual list.
        ec_nums = row.rec_ec_nums[1:-1].replace("\'","").split(", ")
        pident = row.pident
        qcover = row.qcover
        scover = row.scover

        # Parse EC numbers.
        for ec in ec_nums:
            # Add EC numbers to genes and/or increment their counts.
        
            # I'm doing it this way because if transcripts is saved as
            # Nonetype, then it won't point to a dictionary in genes if
            # I add that dictionary after saving transcripts.
            # In other words, I can only edit transcripts through its parent
            # dictionary (genes) if I initially saved it as a dictionary
            # from genes.
            if not genes.get(gene):
                genes[gene] = {} # Initialize transcripts.
            
            transcripts = genes.get(gene)
        
            if not transcripts.get(transcript):
                transcripts[transcript] = {} # Initialize ec_counts.
        
            ec_counts = transcripts.get(transcript)
            
            if not ec_counts.get(ec):
                ec_counts[ec] = {"count":0,
                                 "sum_pident":0, "avg_pident":0, "max_pident":0,
                                 "sum_qcover":0, "avg_qcover":0, "max_qcover":0,
                                 "sum_scover":0, "avg_scover":0, "max_scover":0} # Initialize ec.
                
            ec_counts[ec]["count"] += 1

            ec_counts[ec]["sum_pident"] += pident
            ec_counts[ec]["avg_pident"] = ec_counts[ec]["sum_pident"] / ec_counts[ec]["count"]
            ec_counts[ec]["max_pident"] = max(ec_counts[ec]["max_pident"], pident)

            ec_counts[ec]["sum_qcover"] += qcover
            ec_counts[ec]["avg_qcover"] = ec_counts[ec]["sum_qcover"] / ec_counts[ec]["count"]
            ec_counts[ec]["max_qcover"] = max(ec_counts[ec]["max_qcover"], qcover)

            ec_counts[ec]["sum_scover"] += scover
            ec_counts[ec]["avg_scover"] = ec_counts[ec]["sum_scover"] / ec_counts[ec]["count"]
            ec_counts[ec]["max_scover"] = max(ec_counts[ec]["max_scover"], scover)            

    out_path = find_path(f"{out_directory}/gene_ECs.{current_time}.txt",
                         action="w", path_type="f")

    with open(out_path, "w", encoding="utf-8") as out:
        out.write("gene\ttranscript\tec_num\tmax_count\tavg_pident\tmax_pident\tavg_qcover\tmax_qcover\tavg_scover\tmax_scover\n")

    # For now, print one line per max transcript/ec combo.
    for gene_name in genes.keys():
        transcripts = genes.get(gene_name)
        
        # dict of (transcript, EC nums) by count
        counts_by_pair = {}
        for transcript_name in transcripts.keys():
            ec_counts = transcripts.get(transcript_name)
            
            for ec_name in ec_counts.keys():
                stats = ec_counts.get(ec_name) # dict
                ec_count = stats["count"]
                avg_pident = stats["avg_pident"]
                
                if not counts_by_pair.get(ec_count):
                    counts_by_pair[ec_count] = []

                counts_by_pair[ec_count].append((transcript_name, ec_name))

        # int of the overall max EC count
        max_count = max(counts_by_pair.keys())

        with open(out_path, "a", encoding="utf-8") as out:
            for pair in counts_by_pair.get(max_count):
                transcript_name = pair[0]
                ec_name = pair[1]
                
                stats = genes.get(gene_name).get(transcript_name).get(ec_name)
                
                avg_pident = stats["avg_pident"]
                max_pident = stats["max_pident"]
                avg_qcover = stats["avg_qcover"]
                max_qcover = stats["max_qcover"]
                avg_scover = stats["avg_scover"]
                max_scover = stats["max_scover"]

                out.write(f"{gene_name}\t{transcript_name}\t{ec_name}\t{max_count}\t{avg_pident}\t{max_pident}\t{avg_qcover}\t{max_qcover}\t{avg_scover}\t{max_scover}\n")

### DECIDE HOW TO CHOOSE TRANSCRIPT IF THERE IS ONLY ONE EC...AVG PIDENT?
### DECIDE HOW TO CHOOSE TRANSCRIPT IF THERE ARE MULTIPLE ECS...BY HAND?

if __name__ == "__main__":
    args = parse_args()
    main(args)