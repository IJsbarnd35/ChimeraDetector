#imports
import pandas as pd
import time
import subprocess
import os


# Call minimap
# Only when used seperately
def create_df(fastq):
    if os.path.exists(fastq):
        #subprocess.call("minimap2 -S b1.fq b1.fq > overlaps.paf")
        #subprocess.call("minimap2 -x ava-ont reads.fa reads.fa > overlaps.paf.")
        df = pd.read_csv(fastq, sep='\t', usecols=[0, 1, 2, 3, 5, 6, 7, 8, 4], header=None)
        df.columns = ['query_id', 'query_length', 'query_start', 'query_end', 'target_id', 'target_length',
                      'target_start', 'target_end', 'strand']
        return df #.to_records(index=False).tolist()
    else:
        print("The file path given does not exist. Please check if you entered the correct path.")

# Read paf file
def find_chimeras(df):
    reads = df["query_name"].unique()
    for read in reads:
        print(read)

# Seperate chimeras from repeats

# find common sequences?


# Execute
if __name__ == '__main__':
    df = create_df("D:/School/stage2/lisan/Graduation/Graduation/overlaps.paf")
    print(df)