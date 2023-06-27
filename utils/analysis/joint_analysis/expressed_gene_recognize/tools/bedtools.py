def read_bed_file(filename):
    print('Reading BED file ' + filename)
    bed_entries = []
    with open(filename, 'r') as file:
        for line in file:
            split_line = line.strip().split('\t')
            # Extract chromosome, start, and end. Ignore other fields for simplicity
            chromosome, start, end = split_line[0], int(split_line[1]), int(split_line[2])
            bed_entries.append((chromosome, start, end))
    return bed_entries


def write_bed_file(filename, entries):
    print('Writing BED file ' + filename)
    with open(filename, 'w') as file:
        for entry in entries:
            file.write('\t'.join(map(str, entry)) + '\n')


def intersect(a_entries, b_entries):
    '''
    Find the intersection of two BED files, similar to the bedtools intersect -wa command
    :param a_entries: list of entries in the first BED file, summits
    :param b_entries: list of entries in the second BED file, ranges
    :return: a list of complete ranges from b_entries that has intersection with a_entries
    '''
    print('Intersecting ' + str(len(a_entries)) + ' and ' + str(len(b_entries)) + ' entries')
    # lambda return the larger one of two numbers
    return_larger = lambda x, y: x if x > y else y
    # lambda return the smaller one of two numbers
    return_smaller = lambda x, y: x if x < y else y
    result = []
    for a_chromosome, a_start, a_end in a_entries:
        # print(a_chromosome, a_start, a_end)
        a_larger = return_larger(a_start, a_end)
        a_smaller = return_smaller(a_start, a_end)
        for b_chromosome, b_start, b_end in b_entries:
            b_larger = return_larger(b_start, b_end)
            b_smaller = return_smaller(b_start, b_end)
            # If the entries intersect (have the same chromosome and overlapping regions)
            if a_chromosome == b_chromosome and a_smaller <= b_larger and b_smaller <= a_larger:
                result.append((a_chromosome, a_smaller, a_larger))
    print('Length of intersection is ' + str(len(result)))
    return result


def subtract(a_entries, b_entries):
    '''
    Find the subtraction of two BED files, similar to the bedtools subtract command
    :param a_entries: list of entries in the first BED file
    :param b_entries: list of entries in the second BED file
    :return: list of entries in the first BED file that do not overlap with the second BED file
    '''
    print('Subtracting ' + str(len(a_entries)) + ' and ' + str(len(b_entries)) + ' entries')
    result = []
    for a_chromosome, a_start, a_end in a_entries:
        if_exists = False
        for b_chromosome, b_start, b_end in b_entries:
            # If the entries intersect (have the same chromosome and overlapping regions)
            if a_chromosome == b_chromosome and a_start == b_start and a_end == b_end:
                if_exists = True
                break
        if not if_exists:
            result.append((a_chromosome, a_start, a_end))

    print('Length of subtraction is ' + str(len(result)))
    return result
