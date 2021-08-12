'''
Beagan Nguy
Bi622
Demux assignment

The scripts reads in 4 fastq paired-end files and analyzes them to mapped indexes. 

Usage: python3 twoplex.py 

Optional Parameters:
python3 twoplex.py -i <indexes.txt> -rc <read cut off> -ic <index cut off> -r1 <read1.fastq.gz> -r2 <read2.fastq.gz> -i1 <index1.fastq.gz> -i2 <index2.fastq.gz>

Output:
stats.md - Displays the bins/buckets along with its coressponding percentage compared to the rest of the bins/buckets
52 Fastq Files - Each file is sorted based on it category: Matched(48), Swap(2), Not Good (2)
allfastq.zip - A zip file containing all fastq files

Individual index cutoffs and average read cutoffs are defaulted to 30, which is still good quality with an percent error of 1/1000 = 0.1%.
I did avaerage read cut off because I want this script to do two in one, to both demultiplex and quality control for bad read qualities.

Headers
Headers are appended by the original header, the read end, the first index and second index separated by a '-'. 
Note: The second index is not appended in its reverse complemented form to keep the orginality of the data.
'''
import gzip
from itertools import zip_longest
from zipfile import ZipFile
import sys
import argparse
import matplotlib.pyplot as plt

'''
Arg Parser
'''
def getArgs():
    parser = argparse.ArgumentParser(
        description='Optional arguments: rc, ic, r1, r2, i1, i2'
    )

    parser.add_argument(
        '-i', '-index_lib', help='Input Index Library', type=str, default="/projects/bgmp/shared/2017_sequencing/indexes.txt", required=False
    )
 
    parser.add_argument(
        '-rc', '-read_cut_off', help='Specify read_cutoff', type=int, default=30, required=False
    )

    parser.add_argument(
        '-ic', '-index_cut_off', help='Specify index_cutoff', type=int, default=30, required=False
    )

    parser.add_argument(
        '-r1', '-read_1', help='Specify read 1 file', type=str, default="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz", required=False
    )

    parser.add_argument(
        '-r2', '-read_2', help='Specify read 2 file', type=str, default="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz", required=False
    )

    parser.add_argument(
        '-i1', '-index_1', help='Specify index 1 file', type=str, default="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz", required=False
    )

    parser.add_argument(
        '-i2', '-index_2', help='Specify index 2 file', type=str, default="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz", required=False
    )

    return parser.parse_args()

# Yields each line in the files in parallel
def readFiles(files):
    for lines in zip_longest(*files):
        yield lines

# Returns the reverse complements of a string
def ReverseComplement(seq):
    seq1 = 'ATCGNTAGCNatcgntagcn'
    seq_dict = {seq1[i]:seq1[i+5] for i in range(25) if i < 5 or 10<=i<15}
    rev_comp =  "".join([seq_dict[base] for base in reversed(seq)])
    return rev_comp

# Parse through the index library file and return a dictionary where the key is the index and value is the bin. 
# In addition returns a set of library indexes along with its reverse complemnent sequence.
def readIndexFile(file):
    index_library = set()
    bucket_list = dict()
    with open(file) as fh:
        head = fh.readline()
        for line in fh:
            line = line.rstrip().split('\t')
            bucket_list[line[4]] = [line[3]]
            index_library.add(line[4])
            index_library.add(ReverseComplement(line[4]))
    return bucket_list, index_library

# Convert ascii score to phred score
def convert_phred(letter):
    phred_score = ord(letter) - 33
    return phred_score

# Returns average quality score
def averageQualityScore(q_list):
    sum = 0
    count = 0
    for qscore in q_list:
        sum += convert_phred(qscore)
        count += 1
    average_q1 = sum/count
    return average_q1 

# Write to specified file
# File = The file output address similar to fout
def outputToBucket(head, i1, i2, read, quality, file, read_type):
    index = i1+'-'+i2
    line1 = head + ' : ' +  read_type + ' : ' + index   # Read types: r1, r2, File = File Address
    file.write(line1 + '\n')
    file.write(read + '\n')
    file.write("+" + '\n')
    file.write(quality + '\n')

# Returns the lowest index quality score
def find_min_index(q_list):
    temp_set = set()
    for qscore in q_list:
        temp_set.add(convert_phred(qscore))
    return min(temp_set)

# Peforms demultiplex
class Demultiplexing():
    # Constructor
    def __init__(self, files, bucket_dict, index_libary, read_cut, index_cut):
        self.files = files
        self.bucket_dict = bucket_dict  # Bad Bucket contains low quality reads and indexes
        self.index_library = index_libary
        self.read_cut = read_cut
        self.index_cut = index_cut
        self.organizeIndex()     

    def organizeIndex(self):
        list_of_files = list()      # A list of fastq files names read from the input. Example: read1, read2, index1, index 2
        read_cutoff = self.read_cut
        index_cutoff = self.index_cut
        endcount = 1
        n_set = {'N'}
        head = ''
        i1 = set()
        i2 = set()
        q1 = ''
        q2 = ''
        r1 = ''
        r2 = ''

        # To create a tally of the number of records in each bin/buckets
        bucket_counts = {'swap':0, 'not_good':0}  # Key: bin, Value: [read 1 counts, read 2 counts]
        
        # To be used to tally the total records
        total_counts = 0

        # Open files
        fastq_dict = dict()
        for key, bin in self.bucket_dict.items():       # Key = Index, Value = Bin
            file_address_list = list()
            file_address_list.append(open('{}_r1.fastq'.format(" ".join(bin)),'wt'))
            file_address_list.append(open('{}_r2.fastq'.format(" ".join(bin)),'wt'))
            fastq_dict[key] = file_address_list
            bucket_counts[" ".join(bin)] = 0    # Add more bins to store the counts of reads in each bin

            # Create a list fastq files names to be used downstream for zipping
            list_of_files.append('{}_r1.fastq'.format(" ".join(bin)))
            list_of_files.append('{}_r2.fastq'.format(" ".join(bin)))

        ng_r1 = open('Not_good_r1.fastq', 'w')
        ng_r2 = open('Not_good_r2.fastq', 'w')
        swap_r1 = open('swap_r1.fastq', 'w')
        swap_r2 = open('swap_r2.fastq', 'w')
        
        # Create a list fastq files names to be used downstream for zipping
        list_of_files.append('Not_good_r1.fastq')
        list_of_files.append('Not_good_r2.fastq')
        list_of_files.append('swap_r1.fastq')
        list_of_files.append('swap_r2.fastq')

        for counts, lines in enumerate(readFiles(self.files)): 
            #if endcount == 100:
            #    break
          
           # Assign First header
            if counts == 0:
                head = lines[0].rstrip()
                continue
            
            # Get Quality Scores
            if counts % 4 == 3:
                q1 = lines[0].rstrip()
                q2 = lines[3].rstrip()
                
            
            # Get Reads and Index
            if counts % 4 == 1:
                str_i1 = lines[1].rstrip()
                str_i2 = lines[2].rstrip()
                i1 = {str_i1}
                i2 = {str_i2}
                r1 = lines[0].rstrip()
                r2 = lines[3].rstrip()
                
            # End of Record
            if counts % 4 == 0:
                #print(head, r1, i1, r2, i2, q1, q2)
                qscore_r1 = averageQualityScore(q1)     # Gets the average quality score
                qscore_r2 = averageQualityScore(q2)
                qscore_i1 = find_min_index(str_i1)      # Gets the minimim quality score
                qscore_i2 = find_min_index(str_i2)

                # I did read cut off because I want this script to do two in one, to both demultiplex and quality control.
                if qscore_r1 > read_cutoff and qscore_r2 > read_cutoff and qscore_i1 > index_cutoff and qscore_i2 > index_cutoff:
                    set_i1 = set(i1)
                    set_i2 = set(i2)
                    index_combine = set_i1.union(set_i2)

                    # If there is a "N" in the indexes output to not_good fastq
                    if n_set.issubset(index_combine):
                        # Increment bin count
                        prev_count = bucket_counts['not_good'] + 1
                        bucket_counts['not_good'] = prev_count
                        total_counts += 1

                        # Print to bad bucket, index contains N
                        outputToBucket(head, str_i1, str_i2, r1, q1, ng_r1, "r1")
                        outputToBucket(head, str_i1, str_i2, r2, q2, ng_r2, "r2")

                    else:
                        # Checks to see in index is in the library
                        c1 = i1.intersection(self.index_library)    # Index library includes reverse complements of index
                        c2 = i2.intersection(self.index_library)    # This is why my i2 does not need to be reverse complement in this step

                        # Checks to see if index 1 and index 2 map to a index in the library
                        if c1 and c2:
                            c2 = {ReverseComplement(str_i2)}
                            c3 = c1.intersection(c2)
                            if c3:
                                # Increment bin count <- Keeps a tally of how many times a record is stored into a bin, will be used for stats later
                                bin = " ".join(self.bucket_dict[str_i1])
                                prev_count = bucket_counts[bin] + 1
                                bucket_counts[bin] = prev_count
                                total_counts += 1
                                
                                # Write to good bucket
                                file_r1, file_r2 = [i for i in fastq_dict[str_i1]]      # Get the output bin adress based on index
                                outputToBucket(head, str_i1, str_i2, r1, q1, file_r1, "r1")
                                outputToBucket(head, str_i1, str_i2, r1, q1, file_r2, "r2")
                                #print(c1, c2)
                            else:
                                # Increment bin count
                                prev_count = bucket_counts['swap'] + 1
                                bucket_counts['swap'] = prev_count

                                # Write to swap bucket
                                outputToBucket(head, str_i1, str_i2, r1, q1, swap_r1, "r1")
                                outputToBucket(head, str_i1, str_i2, r2, q2, swap_r2, "r2")
                                #print(c1, c2)
                                pass
                        else:
                            # Increment bin count
                            prev_count = bucket_counts['not_good'] + 1
                            bucket_counts['not_good'] = prev_count
                            total_counts += 1

                            # Write to bad bucket
                            outputToBucket(head, str_i1, str_i2, r1, q1, ng_r1, "r1")
                            outputToBucket(head, str_i1, str_i2, r2, q2, ng_r2, "r2")
                else:
                    # Increment bin count
                    prev_count = bucket_counts['not_good'] + 1
                    bucket_counts['not_good'] = prev_count
                    total_counts += 1

                    #print('To bad bucket')
                    outputToBucket(head, str_i1, str_i2, r1, q1, ng_r1, "r1")
                    outputToBucket(head, str_i1, str_i2, r2, q2, ng_r2, "r2")

                head = lines[0].rstrip()    # The header of the next record
            #endcount +=1

        # Close output files
        for index, file_address_list in fastq_dict.items():
            file_address_list[0].close()
            file_address_list[1].close()
        ng_r1.close()
        ng_r2.close() 
        swap_r1.close() 
        swap_r2.close() 

        # Zip files together
        with ZipFile('all_fastq.gz','w') as zip:
            for filename in list_of_files:
                zip.write(filename)
        
        percentage_bucket_counts = dict()
        # Plot Bar Plot Stats
        for bin, counts in bucket_counts.items():
            percent_counts = (counts/total_counts)*100
            percentage_bucket_counts[bin] = percent_counts
        
        keys = percentage_bucket_counts.keys()
        values = percentage_bucket_counts.values()

        plt.bar(keys, values)   
        plt.title('Percentage of reads from each sample')
        plt.ylabel('Percentage of Reads')
        plt.xlabel('Samples/Bins')
        plt.xticks(rotation = 60)
        plt.tick_params(axis='x', which='major', labelsize=7)
        plt.savefig('Percentage_Reads.png', bbox_inches='tight', dpi=500)
        
        # Print stats to md file
        swap_amount = bucket_counts['swap']
        out_str = 'stats.md'
        fout = open('{}'.format(out_str),'wt')  
        fout.write("Bin\tPercentage\n")      
        for bin, percentage in percentage_bucket_counts.items():
            fout.write("{}\t{:.2f}\n".format(bin, percentage))
        fout.write('Amount of index swapping: {}'.format(swap_amount))
        fout.close()

# Main Function
def main():
    args = getArgs()
    index_filename = args.i
    read_cut = args.rc
    index_cut = args.ic
    read_1 = args.r1
    read_2 = args.r2
    index_1 = args.i1
    index_2 = args.i2
    bucket_dict, index_library = readIndexFile(index_filename)
  
    file_list = [read_1, index_1, index_2, read_2]
    files = [gzip.open('{}'.format(filename), 'rt') for filename in file_list]
    deplex = Demultiplexing(files, bucket_dict, index_library, read_cut, index_cut)
    
if __name__ == "__main__":
	main()
