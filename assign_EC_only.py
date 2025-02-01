"""
For BLAST results in which SwissProt is the database, and mRNA transcripts are queries:
Identifies a single top SwissProt hit for each query (if it exists).
Retrieves an EC number for each SwissProt hit.
Outputs a table of transcripts and putative EC numbers.
More processing will be necessary to identify coding sequences in these transcripts.

NOTE: assign_EC.py (also in this folder) was modified a while back to include 
expression information (with diffexp data as input). This version just finds EC 
numbers for transcripts (this makes it easier to do things downstream).
"""

### Start here 2024.09.27 #####
### This is going to be the same as blast_hit_expression, but with an added column for EC number.
### IMPORTANT: THIS IS BACKWARDS FROM blast_hit_expression; transcripts are the query this time.

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

def find_description(json):
    json_description = json.get("proteinDescription")
    rec_name = json_description.get("recommendedName")
    description = rec_name.get("fullName").get("value")
    
    return description

def parse_json(json, ec_source):
    ec_nums = []

    if ec_source == "recommended":
        json_description = json.get("proteinDescription")
        rec_name = json_description.get("recommendedName")
        ec_dicts = rec_name.get("ecNumbers")
        if ec_dicts:
            for ec in ec_dicts:
                ec_nums.append(ec.get("value"))
        print(f"json_description:{json_description}\nrec_name:{rec_name}\nec_nums:{ec_nums}\n\n",
              flush=True)

    elif ec_source == "domain":
        json_description = json.get("proteinDescription")
        includes = json_description.get("includes")
        domains = []
        if includes:
            for domain_dict in includes:
                domains.append(domain_dict.get("recommendedName"))
        for domain in domains:
            ec_dicts = domain.get("ecNumbers")
            if ec_dicts:
                for ec in ec_dicts:
                    ec_nums.append(ec.get("value"))
        print(f"json_description:{json_description}\nincludes:{includes}\ndomains:{domains}\nec_nums:{ec_nums}\n\n",
              flush=True)
    elif ec_source == "rhea":
        comments = json.get("comments")
        if comments:
            #print(comments)
            for comment in comments:
                #print(comment)
                reaction = comment.get("reaction")
                if reaction:
                    ec_nums.append(reaction.get("ecNumber"))
                    #print(ec_nums)
        print(f"comments:{comments}\nec_nums:{ec_nums}\n\n",
              flush=True)
    else:
        print(f"{ec_source} is not a valid ec_source...\n\n",
              flush=True)
    return ec_nums

def get_metadata(mid):
    """
    Takes in a UniProt ID and returns that ID's JSON metadata.

        Parameters:
            mid (str): A UniProt ID.

        Returns:
            response (str): A JSON of the given ID's metadata.
    """

    # TODO: Consider returning only annotations.
    # Returns entire JSON.
    url_metadata = f"https://rest.uniprot.org/uniprotkb/{mid}"
    print(f"Retrieving metadata for {mid}.", flush=True)
    # Set headers to accept JSON.
    headers = {"Accept": "application/json"}

    # For this request specfically, terminate if the first request is not 200.
    # This is to avoid sending requests for non-UniProt accessions.
    current_request = "Metadata retrieval"
    while True:
        try:
            response_metadata = requests.get(url_metadata, headers)
            if response_metadata.status_code != 200:
                print(f"{current_request} status code: " +
                      f"{response_metadata.status_code}.\n" +
                      f"Failed to retrieve metadata for {mid}. Continuing...\n",
                      flush=True)
                return None # Handle Nonetype in script that calls this.
            else:
                break

        except requests.exceptions.ConnectionError as errc:
            print(f"{current_request} caused a connection error: " +
                  f"{errc}. Retrying...\n", flush=True)
        except requests.exceptions.RequestException as err:
            print(f"{current_request} caused exception: {err}. Exiting...\n",
                  flush=True)
            sys.exit()

        # If query fails, try again after 10 seconds.
        time.sleep(10)

    print(f"Metadata retrieved for {mid}.\n", flush=True)

    return response_metadata.json()

def parse_args():
    """
    Takes in command-line arguments and returns an argparse Namespace object.

        Returns:
            arguments (Namespace): Namespace with command-line arguments.
    """

    parser = argparse.ArgumentParser(description="EC \"Predictor\"")
    # -i and -o must come before tool_name.
    # Ex. -i ./infile.fasta -o ./out_dir annotate

    parser.add_argument("-b", "--blast_results", type=str,
                        help="Full path to BLAST results file in outfmt 6")
    #parser.add_argument("-d", "--deseq2_results", type=str,
                        #help="Full path to DESeq2 results file")
    parser.add_argument("-t", "--transcripts", type=str,
                        help="Full path to transcripts file in FASTA format")
    #parser.add_argument("-nfc", "--neg_log2foldchange", type=int,
                        #help="Log2FoldChange threshold for downregulated proteins.",
                        #default=-1)
    #parser.add_argument("-pfc", "--pos_log2foldchange", type=int,
                        #help="Log2FoldChange threshold for upregulated proteins.",
                        #default=1)
    #parser.add_argument("-sit", "--subject_identity_threshold", type=int,
                        #help="Percent cutoff for reporting BLAST hits based on "+
                        #"number of identical residues, divided by full subject length.",
                        #default=50)
    parser.add_argument("-o", "--out_directory", type=str,
                        help="Full path of output directory. Must end with \"/\".")
    
    return parser.parse_args()

#def write_log(log_path, string):
  
    #with open(log_path, "a", encoding="utf-8") as log:
        #log.write(string)

def main(args):
    """
    Identifies transcripts that align well with enzymes in a custom BLAST database.
    Retrieves EVERYTHING. This will take a while. several days of running, probably.

        Outputs:
            TSV of transcript (gene) hits
            FASTAs of transcript hits <- not implemented - too much space used
    """
    
    # This was originally written to filter as well as retrieve EC numbers...don't bother.
    # Retrieve everything and filter later.
    # The below is what I was doing before:
    #For now, uses a 50% identity (based on full subject length) threshold to identify hits.
    #For now, only outputs the top transcript hit for each gene.
    #See how that works.
    
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
    
    # TODO: Make log file include all version information for Python and libraries
    #with open(log_path, "w", encoding="utf-8") as log:
        #log.write(f"{str(session_info.show())}\n")

    blast_results_file = find_path(args.blast_results, action="r", path_type="f")
    blast_results_filename = blast_results_file.split("/")[-1]
    blast_results_df = pd.read_csv(blast_results_file, sep="\t",
                             names=["qseqid", "sseqid", "qlen", "slen", "length", "nident", "evalue"])
    print(blast_results_df)
    
    #deseq2_results_file = find_path(args.deseq2_results, action="r", path_type="f")
    #deseq2_results_filename = deseq2_results_file.split("/")[-1]
    #deseq2_results_df = pd.read_csv(deseq2_results_file, sep="\t",
                                    #header=0,
                                    #names=["gene", "baseMean", "log2FoldChange",
                                           #"lfcSE", "pvalue", "padj"]).set_index("gene")

    #print(deseq2_results_df)
    
    #subjects_file = find_path(args.subjects, action="r", path_type="f")
    #subjects_filename = subjects_file.split("/")[-1]
    #subjects_df = pd.read_csv(subjects_file, sep="\t",
                              #header=0).set_index("BLASTID")

    #print(subjects_df)
    
    transcripts_file = find_path(args.transcripts, action="r", path_type="f")
    transcripts_filename = transcripts_file.split("/")[-1]
    transcripts_df = fasta_to_df(transcripts_file)
    print(transcripts_df)
    
    out_directory = find_path(args.out_directory, action="w", path_type="d")
    print("\nInputs:\n") # Prints a newline
    print(f"blast_results: {blast_results_file}\ntranscripts: {transcripts_file}")
    print(f"out_directory: {out_directory}\n")    
    
    # TODO: parse current_time
    hit_results_file = find_path(f"{out_directory}/EC_results.{current_time}.tsv", action="w", path_type="f")

    #blast_hits_directory = find_path(f"{args.out_directory}/blast_hits/", action="w", path_type="d")

    with open(hit_results_file, "w", encoding="utf-8") as hit_results:
        #hit_results.write("gene\ttranscript\tsubject\tsubject_description\te_value\tperc_sidentity\tbaseMean\tlog2FC\tpadj\n")
        #hit_results.write(f"gene\ttranscript\tsubject\trec_ec_nums\tdom_ec_nums\trhea_ec_nums\te_value\tperc_sidentity\tbaseMean\tlog2FC\tpadj\n")
        #hit_results.write(f"gene\ttranscript\tsubject\tsubject_description\trec_ec_nums\tdom_ec_nums\trhea_ec_nums\te_value\tperc_sidentity\n")
        hit_results.write(f"gene\ttranscript\tsubject\tsubject_description\trec_ec_nums\t" +
                          "dom_ec_nums\trhea_ec_nums\te_value\tperc_sidentity\tperc_qidentity\t" +
                          "qcover\tslen\tqlen_nt\tmax_qlen_aa\talign_len\n")

        # NOTE: I've changed this a lot. There are several iterations of filtering here.
        # NOTE: I have decided I just want this to retrieve everything.
        # NOTE: It will take forever, so I'm going to try to keep a dictionary in memory:
        # NOTE: {subject:{subject_description:sd,rec_ec_nums:ren,dom_ec_nums:den,rhea_ec_nums:rhen}}

        print("Finding EC numbers for all BLAST hits...\n", flush=True)
        print("Be aware that no filtering is being performed here.\n\n",
              flush=True)

        # Dictionary of already retrieved UniProt entries.
        hit_dict = {}
        
        for i in range(len(blast_results_df)):
            subject = blast_results_df.iloc[i]["sseqid"]
            transcript = blast_results_df.iloc[i]["qseqid"]
            gene = ".".join(transcript.split(".")[:-1])
            e_value = blast_results_df.iloc[i]["evalue"]
            nident = int(blast_results_df.iloc[i]["nident"])
            slen = int(blast_results_df.iloc[i]["slen"])
            perc_sidentity = nident / slen * 100
            qlen_nt = int(blast_results_df.iloc[i]["qlen"])
            max_qlen_aa = int(qlen_nt / 3)
            perc_qidentity = nident / max_qlen_aa * 100
            align_len = int(blast_results_df.iloc[i]["length"])
            qcover = align_len / max_qlen_aa * 100

            print("BLAST result:", flush=True)
            print(f"gene: {gene}\ntranscript: {transcript}\nsubject: {subject}\n",
                  flush=True)

            # Save to dictionary or retrieve old entry.
            hit_entry = hit_dict.get(subject)
            if not hit_entry:
                json = get_metadata(subject.split("|")[-1])
        
                if json:
                    subject_description = find_description(json)                
                    rec_ec_nums = parse_json(json, "recommended")
                    dom_ec_nums = parse_json(json, "domain")
                    rhea_ec_nums = parse_json(json, "rhea")

                hit_dict[subject] = {"subject_description":subject_description,
                                     "rec_ec_nums":rec_ec_nums,
                                     "dom_ec_nums":dom_ec_nums,
                                     "rhea_ec_nums":rhea_ec_nums}
            else:
                subject_description = hit_entry.get("subject_description")
                rec_ec_nums = hit_entry.get("rec_ec_nums")
                dom_ec_nums = hit_entry.get("dom_ec_nums")
                rhea_ec_nums = hit_entry.get("rhea_ec_nums")
    
            hit_results.write(f"{gene}\t{transcript}\t{subject}\t{subject_description}\t{rec_ec_nums}\t" +
                              f"{dom_ec_nums}\t{rhea_ec_nums}\t{e_value}\t{perc_sidentity}\t{perc_qidentity}\t" +
                              f"{qcover}\t{slen}\t{qlen_nt}\t{max_qlen_aa}\t{align_len}\n")







        ## {gene1: {features}, gene2: {features}}
        #top_blast_hits = {}

        #print("Finding top BLAST hits...\n\n", flush=True)

        #for i in range(len(blast_results_df)):
            #subject = blast_results_df.iloc[i]["sseqid"]
            #transcript = blast_results_df.iloc[i]["qseqid"]
            #gene = ".".join(transcript.split(".")[:-1])
            #e_value = blast_results_df.iloc[i]["evalue"]
            #perc_sidentity = (int(blast_results_df.iloc[i]["nident"])
                              #/ int(blast_results_df.iloc[i]["slen"])
                              #* 100)

            ##print(f"Query: {gene}\nSubject: {subject}\n\n", flush=True)
            ##try:
                ##baseMean = deseq2_results_df.loc[gene]["baseMean"]
                ##log2FC = deseq2_results_df.loc[gene]["log2FoldChange"]
                ##padj = deseq2_results_df.loc[gene]["padj"]
            ##except KeyError:
                ##baseMean = "-"
                ##log2FC = "-"
                ##padj = "-"
                ##write_log(log_path, f"{gene} not found in DESeq2 results.\n" +
                          ##"Possibly removed due to low genecounts.\n")

            ## TODO: REPLACE THIS WITH A DESCRIPTION FROM JSON
            ##subject_description = subjects_df.loc[subject]["protein_name"]

            ## dict.get(key) hould return None if key does not exist.
            ## There should be nothing else that evaluates to False in values.
            #if top_blast_hits.get(gene) is None:
                ##top_blast_hits[gene] = {"transcript":transcript,
                                        ##"subject":subject,
                                        ###"subject_description":subject_description,
                                        ##"e_value":e_value,
                                        ##"perc_sidentity":perc_sidentity,
                                        ##"baseMean":baseMean,
                                        ##"log2FC":log2FC,
                                        ##"padj":padj}
                #top_blast_hits[gene] = {"transcript":transcript,
                                        #"subject":subject,
                                        ##"subject_description":subject_description,
                                        #"e_value":e_value,
                                        #"perc_sidentity":perc_sidentity}                

            ## Compares top hit with current hit in terms of percent identity, relative to query (gene) length
            #elif top_blast_hits.get(gene).get("perc_sidentity") < perc_sidentity:
                ##top_blast_hits[gene] = {"transcript":transcript,
                                        ##"subject":subject,
                                        ###"subject_description":subject_description,
                                        ##"e_value":e_value,
                                        ##"perc_sidentity":perc_sidentity,
                                        ##"baseMean":baseMean,
                                        ##"log2FC":log2FC,
                                        ##"padj":padj}
                #top_blast_hits[gene] = {"transcript":transcript,
                                        #"subject":subject,
                                        ##"subject_description":subject_description,
                                        #"e_value":e_value,
                                        #"perc_sidentity":perc_sidentity}                

        #print("Finished finding top BLAST hits.\nIdentifying matches and retrieving EC numbers...\n\n",
              #flush=True)
        #for gene in top_blast_hits:
            #print(f"Working on {gene}...\n\n", flush=True)
            #features = top_blast_hits.get(gene)
            #transcript = features["transcript"]
            #subject = features["subject"]
            ##subject_description = features["subject_description"]
            #subject_description = "placeholder"
            #e_value = features["e_value"]
            #perc_sidentity = features["perc_sidentity"]
            ##baseMean = features["baseMean"]
            ##log2FC = features["log2FC"]
            ##padj = features["padj"]

            ## PLACEHOLDER: RETRIEVE EC HERE:
            #json = get_metadata(subject.split("|")[-1])

            #if json:
                #subject_description = find_description(json)                
                #rec_ec_nums = parse_json(json, "recommended")
                #dom_ec_nums = parse_json(json, "domain")
                #rhea_ec_nums = parse_json(json, "rhea")

            #hit_results.write(f"{gene}\t{transcript}\t{subject}\t{subject_description}\t{rec_ec_nums}\t{dom_ec_nums}\t{rhea_ec_nums}\t{e_value}\t{perc_sidentity}\n")
            
            ##if perc_sidentity > args.subject_identity_threshold:
                ##if log2FC != "-":
                    ### Setting a negative pos_log2foldchange and a positive neg_log2foldchange will guarantee all proteins get through.
                    ##if log2FC > args.pos_log2foldchange or log2FC < args.neg_log2foldchange:
                        ###hit_results.write(f"{gene}\t{transcript}\t{subject}\t{subject_description}\t{e_value}\t{perc_sidentity}\t{baseMean}\t{log2FC}\t{padj}\n")
                        ##hit_results.write(f"{gene}\t{transcript}\t{subject}\t{rec_ec_nums}\t{dom_ec_nums}\t{rhea_ec_nums}\t{e_value}\t{perc_sidentity}\t{baseMean}\t{log2FC}\t{padj}\n")

if __name__ == "__main__":
    args = parse_args()
    main(args)
