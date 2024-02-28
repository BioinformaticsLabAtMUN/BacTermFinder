import os
import pandas as pd

def load_fasta_to_df():
    """Load the fasta file into a dictionary"""
    # go through each dir and open files that contain "_100bp.fasta" 
    lisdir = os.listdir(".")
    lisdir = lisdir.sort()
    for dir in os.listdir("."):
        if os.path.isdir(dir):
            for file in os.listdir(dir):
                # file contain 100bp.fasta
                if file.endswith("_100bp_10x_neg.fasta"): # or file.endswith("_100bp.fasta")                           # Change me
                # if file.endswith("_10x_neg.fasta"):                                                               # Change me                            
                    with open(dir + "/" + file, "r") as f:
                        # read the file into a list
                        name = dir + "_" + file.replace(".bed", "").replace(".fasta", "")
                        print(name)
                        lines = f.readlines()
                        # create a dictionary with the key being the header and the value being the sequence
                        fasta_dict = {}
                        for line in lines:
                            if line.startswith(">"):
                                header = line.split(">")[1][:-1]
                                fasta_dict[header] = ""
                            else:
                                fasta_dict[header] += line.strip()
                        # create a dataframe from the dictionary
                        df = pd.DataFrame.from_dict(fasta_dict, orient="index")
                        # rename the columns
                        df.columns = ["sequence"]
                        # reset the index
                        df.reset_index(inplace=True)
                        # rename the columns
                        df.rename(columns={"index": "header"}, inplace=True)
                        # add a column with the directory name
                        df["dir"] = dir
                        name_split = name.split("_")
                        df["specie"] = name_split[2]
                        df["strain"] = name_split[3]
                        df.to_csv(f"all_csv/10x_neg_overlapping_20percent/{name}.csv", columns=None, index=None)                               # Change me
                        # if all neg exists, append to it
                        if os.path.exists(f"all_csv/10x_neg_overlapping_20percent/all_neg.csv"):                                               # Change me
                            df.to_csv(f"all_csv/10x_neg_overlapping_20percent/all_neg.csv", mode="a", header=False, index=False, columns=None) # Change me
                        else:
                            df.to_csv(f"all_csv/10x_neg_overlapping_20percent/all_neg.csv", columns=None, index=None)                          # Change me
                                       

            # with open(dir + "/*_100bp.fasta", "r") as f:
            #     # read the file and split it into a list
            #     fasta = f.read().split(">")
            #     # remove the first element of the list because it's empty
            #     fasta.pop(0)
            #     # create a dictionary with the fasta file
            #     fasta_dict = {}
            #     for i in fasta:
            #         # split each element of the list into a list
            #         i = i.split("\n")
            #         # remove the first element of the list because it's empty
            #         i.pop(0)
            #         # create a key for each element of the list
            #         fasta_dict[i[0]] = i[1]
            #     # create a dataframe from the dictionary
            #     df = pd.DataFrame(fasta_dict.items(), columns=["Name", "Sequence"])
            #     # save the dataframe as a csv file
            #     print(df.head())
            #     # df.to_csv(dir + "/_100bp.csv", index=False)


if __name__ == "__main__":
    load_fasta_to_df()