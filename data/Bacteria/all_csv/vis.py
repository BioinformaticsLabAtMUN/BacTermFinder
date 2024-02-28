# load all of the csv files into a list

import pandas as pd
import os
import matplotlib as plt

def read_all_csv():
    # reads all csv files in the directory and put them in one df
    ls = []
    for dir in os.listdir():
        if dir.endswith(".csv"):
            df_temp = pd.read_csv(dir)
            ls.append(df_temp)
    return ls

def vis_enrichment(seq):
    num_rows = len(seq)
    # create a dictionary with the nucleotides
    nucleotides = {"A": 0, "T": 0, "C": 0, "G": 0}
    # create a dictionary with the nucleotides and their percentages
    nucleotides_perc = {"A": 0, "T": 0, "C": 0, "G": 0}
    for nuc in seq:
        nucleotides[nuc] += 1
    nucleotides_perc["A"] = nucleotides["A"] / num_rows
    nucleotides_perc["T"] = nucleotides["T"] / num_rows
    nucleotides_perc["C"] = nucleotides["C"] / num_rows
    nucleotides_perc["G"] = nucleotides["G"] / num_rows

    return [nucleotides_perc["A"], nucleotides_perc["T"], nucleotides_perc["C"], nucleotides_perc["G"]]

def main():
    ls = read_all_csv()
    # df = pd.concat(ls)
    
    df_A = ls[0]['sequence'].apply(vis_enrichment)
    print(df_A)

    # The df_A is a pd series of lists, i want every every element become the column
    df_A = pd.DataFrame(df_A.tolist(), columns=['A', 'T', 'C', 'G'])
    print(df_A)

    df_A.plot(kind='bar')
    # make the dictionary into a dataframe for every key as a column
    # df_0 = pd.DataFrame.from_dict(df_0, orient='index', columns=['A', 'T', 'C', 'G'] )
    # print(df_0)

if __name__ == "__main__":
    main()