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
    
    chimeras = 0
    chimeralist = []
    for read in reads:
        ding = df[(df['query_name'] == read) & (df['target_name'] == read)]
        if len(ding["query_length"].unique()) >= 2:
            # Skip and remove the read if it has multiple sequence lengths, probably a minimap error
            df = df[df.query_name != read]
            df = df[df.target_name != read]
            reads = reads[reads != read]
            continue
        if len(ding) >= 1:
            # Ignore all reads shorter than the threshold
            if ding['query_length'].unique() >= 2000:
                strand = ding['strand'].unique()
                if strand == "-":
                    # check coverage
                    # If 2 different parts, meaning if there is a middle part -> ?
                    for index, match in ding.iterrows():
                        if (match[3] - match[2]) / match[1] >= 0.5:
                            chimeras += 1
                            chimeralist.append(read)
                elif len(strand) == 2:
                    # check if '-' has full coverage
                    for index, match in ding.iterrows():
                        if match[4] == "-" and (match[3] - match[2]) / match[1] >= 0.5:
                            print(read, " is a chimera")
                            chimeras += 1
                            chimeralist.append(read)
    print("chimeras ", chimeras)
    print(chimeralist)

# Seperate chimeras from repeats

# find common sequences?


# Execute
if __name__ == '__main__':
    df = create_df("D:/School/stage2/lisan/Graduation/Graduation/overlaps.paf")
    print(df)