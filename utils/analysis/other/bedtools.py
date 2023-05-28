def read_bed_file(filename):
    bed_entries = []
    with open(filename, 'r') as file:
        for line in file:
            split_line = line.strip().split('\t')
            # Extract chromosome, start, and end. Ignore other fields for simplicity
            chromosome, start, end = split_line[0], int(split_line[1]), int(split_line[2])
            bed_entries.append((chromosome, start, end))
    return bed_entries

def write_bed_file(filename, entries):
    with open(filename, 'w') as file:
        for entry in entries:
            file.write('\t'.join(map(str, entry)) + '\n')

def intersect(a_entries, b_entries):
    result = []
    for a_chromosome, a_start, a_end in a_entries:
        for b_chromosome, b_start, b_end in b_entries:
            # If the entries intersect (have the same chromosome and overlapping regions)
            if a_chromosome == b_chromosome and a_start < b_end and b_start < a_end:
                result.append((a_chromosome, max(a_start, b_start), min(a_end, b_end)))
    return result

def subtract(a_entries, b_entries):
    result = []
    for a_chromosome, a_start, a_end in a_entries:
        overlaps = []
        for b_chromosome, b_start, b_end in b_entries:
            if a_chromosome == b_chromosome and a_start < b_end and b_start < a_end:
                overlaps.append((max(a_start, b_start), min(a_end, b_end)))
        if overlaps:
            overlaps.sort()
            prev_end = a_start
            for overlap_start, overlap_end in overlaps:
                if prev_end < overlap_start:
                    result.append((a_chromosome, prev_end, overlap_start))
                prev_end = max(prev_end, overlap_end)
            if prev_end < a_end:
                result.append((a_chromosome, prev_end, a_end))
        else:
            result.append((a_chromosome, a_start, a_end))
    return result

# Load the BED files
summits = read_bed_file('test_data/LNCaP_DHT_AR_1_PK_summits.bed')
promoter_2kb = read_bed_file('test_data/promoters_regions_UD2kb_Mod.bed')
promoter_10kb = read_bed_file('test_data/promoters_regions_UD10kb_Mod.bed')

# Intersect the summits with the promoters
AR_promoter_binding_2kb = intersect(summits, promoter_2kb)
AR_promoter_enhancer_binding_10kb = intersect(summits, promoter_10kb)

# Subtract the 2kb promoter binding from the 10kb promoter binding
AR_enhancer_between_2kb_10kb = subtract(AR_promoter_enhancer_binding_10kb, AR_promoter_binding_2kb)

# Save the results to BED files
write_bed_file('test_result/AR_promoter_binding_UD2kb.bed', AR_promoter_binding_2kb)
write_bed_file('test_result/AR_promoter_enhancer_binding_UD10kb.bed', AR_promoter_enhancer_binding_10kb)
write_bed_file('test_result/AR_enhancer_between_2kb_10kb.bed', AR_enhancer_between_2kb_10kb)
