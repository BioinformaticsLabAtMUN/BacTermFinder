# this files reads a genome fasta file
# then create sequences of 100 bb based on sliding window of that got from the sys.argv
# then it generates features for the sequences with FileProcessing module in iLearnPlus
# then loads the trained model and predicts the sequences
# then it writes the results in a file, results have the position of the sequence and the probability and the sequence itself

import os
import pickle
import sys
import time
from multiprocessing import Pool, cpu_count

import numpy as np
import pandas as pd
import tensorflow as tf
from Bio import SeqIO
from Bio.Seq import Seq
from keras.models import load_model
from tqdm import tqdm
import shutil


def extract_sliding_windows(ref_genome_file: str, window_size: int,
                            step_size: int) -> pd.DataFrame:
    """This function extracts sliding windows from a reference genome.

    Args:
        ref_genome_file (str): this is the path to the reference genome file
        window_size (int): this is the size of the sliding window
        step_size (int): this is the step size of the sliding window

    Returns:
        pd.DataFrame: this is a Pandas DataFrame containing the sliding windows
    """
    # Load reference genome into memory
    ref_genome = SeqIO.to_dict(SeqIO.parse(ref_genome_file, "fasta"))

    # Define function to extract windows from a single sequence
    def extract_windows_from_sequence(seq_id):
        seq = ref_genome[seq_id].seq
        seq_len = len(seq)
        windows = []
        u_id = 0
        for i in tqdm(range(0, seq_len - window_size + 1, step_size)):
            u_id += 1
            window_start = i
            window_end = i + window_size
            window_seq = str(seq[window_start:window_end])
            windows.append(
                (u_id, seq_id, window_start, window_end, window_seq, '+'))
            revcomp = str(Seq(window_seq).reverse_complement())
            windows.append(
                (u_id, seq_id, window_start, window_end, revcomp, '-'))
        return windows

    windows = []
    for seq_id in ref_genome.keys():
        seq_windows = extract_windows_from_sequence(seq_id)
        windows.extend(seq_windows)

    # Convert windows to a Pandas DataFrame
    df = pd.DataFrame(
        windows, columns=['u_id', 'seq_id', 'start', 'end', 'seq', 'strand'])

    return df


def df_to_fasta(df: pd.DataFrame, filename: str, train_stat="testing") -> None: #TODO: This could be written to work faster
    """This function writes a Pandas DataFrame to a FASTA file.

    Args:
        df (pd.DataFrame): This is the Pandas DataFrame to be written to a FASTA file.
        filename (str): This is the name of the FASTA file to be written.
        train_stat (str, optional): This is train or test status. It is required by iLearnPlus, but doesn't affect our program functionality. Defaults to "testing".

    Returns:
        None: This function does not return anything.
    """    
    for i in range(len(df)):
        with open(f'{filename}', 'a') as f:
            head = str(df["u_id"][i]) + "_" + str(df["seq_id"][i]) + "_" + str(
                df["start"][i]) + "_" + str(df["end"][i]) + "_" + str(
                    df["strand"][i])
            f.write(f'>{head}|{"-1"}|{train_stat}\n')
            f.write(f'{df["seq"][i]}\n')


def csv_reader_low(path):
    df_test = pd.read_csv(f'{path}', nrows=100)
    float_cols = [c for c in df_test if df_test[c].dtype == "float64"]
    float32_cols = {c: np.float32 for c in float_cols}
    # if bin or PS2 is in the name of the column float32_cols, it is a uint8
    for col in float32_cols.keys():
        if "bin" in col or "PS2" in col:
            float32_cols[col] = np.uint8

    float32_cols['SampleName'] = str
    float32_cols['label'] = bool
    df = pd.read_csv(f'{path}', engine='pyarrow', dtype=float32_cols)
    return df


def join_files(path, i):
    df = csv_reader_low(path + f"/ENAC-{i}-0.csv")
    lens = len(df.columns) - 2
    for file in os.listdir(path):
        if file.endswith(f"{i}.csv") and (not file.endswith(f"ENAC-{i}.csv")):
            # read the file
            df1 = csv_reader_low(path + "/" + file)
            df1 = df1.drop(['label'], axis=1)
            # append the file name to the all column names except the first two columns
            df1.columns = df1.columns[:2].tolist() + [
                file.replace(".csv", '') + "_" + col for col in df1.columns[2:]
            ]

            lens += len(df1.columns) - 2
            # join the files
            col_name = file.replace(".csv", '')
            print("Merging features from bathces", col_name)
            df = pd.merge(
                df,
                df1,
                how='outer',
                on=['SampleName'],
            )
    return df, lens


def col_dropper(path_to_shap, path_to_feature_imp, df):
    shap_importance = pd.read_csv(path_to_shap)
    feature_importances = pd.read_csv(path_to_feature_imp)

    # select 80 percent of feature with quantile to keep
    shap_importance_to_keep = shap_importance[
        shap_importance['feature_importance_vals'] >
        shap_importance['feature_importance_vals'].quantile(0.2)]
    feature_importance_to_keep = feature_importances[
        feature_importances['importance'] >
        feature_importances['importance'].quantile(0.2)]
    # intersection of the two
    features_to_keep = list(
        set(shap_importance_to_keep['col_name']).intersection(
            set(feature_importance_to_keep['feature'])))
    features_to_keep = [
        x for x in features_to_keep if any(c in x for c in [
            'Geary', 'NMBroto', 'PseKNC', 'Z_curve_144bit', 'Z_curve_9bit',
            'ENAC'
        ])
    ]
    # drop the columns that are not in the intersection
    df = df.drop([col for col in df.columns if col not in features_to_keep],
                 axis=1)
    return df


def read_csv_low(file, data_path, input_dim):
    df_test = pd.read_csv(f'{data_path}', nrows=100)
    dtype_cols = [c for c in df_test if df_test[c].dtype == "float64"]

    if file in ['PS2.csv', 'binary.csv']:
        dtype_cols = {c: np.int8 for c in dtype_cols}
    else:
        dtype_cols = {c: np.float32 for c in dtype_cols}

    x = pd.read_csv(
        f'{data_path}',
        dtype=dtype_cols,
        # engine = 'pyarrow',
        #  nrows=80000, #e################################################# comment for production
    )

    sample_names = x['SampleName']
    x.drop(columns=['SampleName', 'label'], inplace=True)

    reshaper_dim = list(input_dim[:])
    reshaper_dim.insert(0, len(x))
    # reshaper_dim.append(1)
    reshaper_dim = tuple(reshaper_dim)
    x = x.values.reshape(reshaper_dim)

    return x, sample_names


if __name__ == '__main__':
    # time it
    start_time = time.time()

    ############################################ Sys inputs ########################################################
    # get the parameters from the command line
    genome_file = sys.argv[1]  # genome file name fasta format
    step_size = int(sys.argv[2])  # step size for sliding window ( aka stride size )
    output_file = sys.argv[3]  # output file name
    batch_size = int(sys.argv[4])  # batch size for iLearnPlus feature generation
    WINDOW_SIZE = 101  # window size for sliding windows, fixed to 101
    ############################################ sliding window generation ##########################################
    # get the sequences
    print("\n")
    print("Sliding windows generation started\n")
    if os.path.exists('Sample.csv'):
        os.remove('Sample.csv')
    df_slide = extract_sliding_windows(genome_file, WINDOW_SIZE, step_size)
    df_slide.to_csv(output_file.split(".")[0] + f'_{genome_file}_sliding_windows.csv',
                    index=False)
    print("\nSliding windows generated")
    ############################################  iLearnPlus ########################################################
    if os.path.exists('df_sample.fasta'):
        os.remove('df_sample.fasta')
    print("\n")
    print("iLearnPlus biological feature generation started")
    # convert the dataframe to fasta format
    df_to_fasta(df_slide.reset_index(drop=True), "df_sample.fasta", "training")

    # remove output folder if it exists
    if os.path.exists('output_sample'):
        shutil.rmtree('output_sample')

    # if os.path.exists('output_sample_merged'):
    #     os.remove('output_sample_merged')

    # generate features
    os.system('python iLearnPlus/util/FileProcessing.py ' +
              'df_sample.fasta' + ' ' + str(batch_size) + ' ' + '16' + ' ' +
              'output_sample')

    files = os.listdir('output_sample')
    
    number_of_batches = len(df_slide) // batch_size
    print("number_of_batches", number_of_batches )
    print("\niLearnPlus biological feature generation finished")
    endtime = time.time()
    print("Feature generation took", endtime - start_time, "s.")

    ############################################ Loading and prediction ########################################################

    # check gpu exists
    if os.system('nvidia-smi') == 0:
        gpu_exists = 1
        print('###### GPU exists, predicting with GPU ########')
    else:
        gpu_exists = 0
        print('###### GPU doesn\'t exists, predicting with CPU  ########')

    input_dim_dict = {
        'ENAC.csv': (97, 4),
        'PS2.csv': (100, 16),
        'NCP.csv': (101, 3),
        'binary.csv': (101, 4),
    }
    

    data_path = 'output_sample/'
    print('\nLoading data')
    # files = os.listdir(data_path)
    for embedding in input_dim_dict.keys():
        # load deep learning model
        print(f"\nLoading Deep Learning Model {embedding} \n")
        model = load_model(f'deep-bactermfinder_3cnn_1d_1cat_reduced_10x_{embedding}_saved_model.h5')

        embedding_wo_csv = embedding.split('.csv')[0]
        out_embed = pd.DataFrame(columns=['SampleName', f'probability_{embedding_wo_csv}'])
        
        for batch in range(number_of_batches + 1):
            dp = f'{data_path}{embedding_wo_csv}-{batch}.csv'
            x, sample_info = read_csv_low(embedding, dp, input_dim_dict[embedding])
            print(f'x shape of {embedding} is: {x.shape}')

            print("Predicting sequences")
            with tf.device('/gpu:0'):
                y_pred = model.predict(x)

            # appending the results
            print("Appending the results")

            sample_info = sample_info.reset_index()
            y_pred = pd.DataFrame(y_pred, )
            y_pred[f'probability_{embedding_wo_csv}'] = y_pred[0]
            
            y_pred = pd.concat([sample_info, y_pred], axis=1)

            out_embed = pd.concat([out_embed, y_pred],ignore_index=True)
            del x
            del y_pred
            del sample_info
            
                 
        # write the results
        out_embed.drop(columns=['index', 0], inplace=True)
        out_embed.to_csv(output_file + '_' + str(embedding_wo_csv) + genome_file + '.csv', index=False)

        del out_embed
        del model

    endtime_pred = time.time()
    print("Predicting sequences took", endtime_pred - endtime, "s.")
    
    # JOIN THE RESULTS, 
    # all of the columns in the rsults are the same exp the probability column
    # so we need to just add that column to the first dataframe
    
    # read the first file
    df = pd.read_csv(output_file + '_' + str(embedding_wo_csv) + genome_file + '.csv')
    # read the rest of the files
    for embedding in input_dim_dict.keys():
        if embedding != 'binary.csv':
            embedding_wo_csv = embedding.split('.csv')[0]
            df1 = pd.read_csv(output_file + '_' + embedding_wo_csv + genome_file +  '.csv')
            df[f'probability_{embedding_wo_csv}'] = df1[f'probability_{embedding_wo_csv}']
            del df1
    # write the results
    df['probability_mean'] = df[ [col for col in df.columns if 'probability' in col] ].mean(axis=1)
    df.to_csv(output_file + genome_file + '_mean.csv', index=False)

    
    ############################################ timing and done ########################################################
    endtime = time.time()
    print("Totally, it took", endtime - start_time,
          f"s to find the terminators of a genome file that has {len(df_slide)} windows")
    print('Done')
