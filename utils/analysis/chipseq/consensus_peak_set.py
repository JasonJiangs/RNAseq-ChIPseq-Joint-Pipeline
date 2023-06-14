import pandas as pd

def read_peak_file(file_path):
    return pd.read_csv(file_path, sep='\t', header=None, names=['chrom', 'start', 'end', 'name', 'score'])

def overlap(a_start, a_end, b_start, b_end):
    return max(0, min(a_end, b_end) - max(a_start, b_start))

def find_consensus_peaks(df, min_overlap=1):
    return df[df['count'] >= min_overlap]

def calculate_percentage_of_consensus_peaks(consensus_peaks, replicate_peaks):
    total_consensus_peaks = len(consensus_peaks)
    overlap_count = 0

    for _, consensus_peak in consensus_peaks.iterrows():
        for _, replicate_peak in replicate_peaks.iterrows():
            if (consensus_peak['chrom'] == replicate_peak['chrom'] and
                    overlap(consensus_peak['start'], consensus_peak['end'],
                            replicate_peak['start'], replicate_peak['end']) > 0):
                overlap_count += 1
                break

    return (overlap_count / total_consensus_peaks) * 100

# Read peak files from multiple replicates
replicate_files = ['../data/chipseq/macs2/LNCaP_DHT_AR_1_PK_summits.bed',
                   '../data/chipseq/macs2/LNCaP_DHT_AR_2_PK_summits.bed']
replicate_data = [read_peak_file(file) for file in replicate_files]

# Merge all peak files into a single DataFrame
merged_peaks = pd.concat(replicate_data).reset_index(drop=True)

# Group by genomic location and count the number of occurrences (i.e., replicates)
merged_peaks['loc'] = merged_peaks['chrom'].astype(str) + ':' + merged_peaks['start'].astype(str) + '-' + merged_peaks['end'].astype(str)
peak_counts = merged_peaks.groupby(['chrom', 'start', 'end', 'loc']).size().reset_index(name='count')

# Define the minimum number of replicates for a peak to be included in the consensus set
min_overlap = len(replicate_data) // 2 + 1

# Create consensus peak set
consensus_peaks = find_consensus_peaks(peak_counts, min_overlap)

# save
consensus_peaks.to_csv('result/consensus_peaks.csv', index=False)

# Calculate the percentage of consensus peaks found in each individual replicate
for i, replicate in enumerate(replicate_data):
    percentage = calculate_percentage_of_consensus_peaks(consensus_peaks, replicate)
    print(f'Replicate {i + 1}: {percentage:.2f}% of consensus peaks')
