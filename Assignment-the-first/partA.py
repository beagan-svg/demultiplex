import gzip
from itertools import zip_longest
import sys
import pprint

#r1,i1,i2, r2

# 2,3
def readFiles(files):
    for lines in zip_longest(*files):
        yield lines

# Converts a single character into a phred score
def convert_phred(letter):
    phred_score = ord(letter) - 33
    return (phred_score)

# Some quality scores are empty (q_line)
def updateQuality(q_line, qlist):
    quality_count = 1
    for count, score in enumerate(q_line):
        num = convert_phred(score)
        try:
            qlist[count] += num
        except:
            qlist.append(num)
    return qlist, quality_count

# Average quality list
def averageQualityScore(q_list, num):
    temp_list = list()
    for i in q_list:
        average = i/num
        temp_list.append(average)
    return temp_list

class Demultiplexing():
    def __init__(self, files):
        self.files = files
        self.getIndexAndSeq()

    def getIndexAndSeq(self):
        endcount = 1
        num_qr1 = 0
        num_qr2 = 0
        num_qi1 = 0
        num_qi2 = 0
        qr1 = ''
        qr2 = ''
        qi1 = ''
        qi2 = ''
        qr1_list = list()
        qr2_list = list()
        qi1_list = list()
        qi2_list = list()
	
	#Iterate through the records of all the files
        for counts, lines in enumerate(readFiles(self.files)): 
            # Assign First header

            if endcount == 10000000000000000000000000:
                break
            endcount +=1

            # Get Quality Scores
            if counts % 4 == 3:
                qr1 = lines[0].decode("utf-8").rstrip()
                qi1 = lines[1].decode("utf-8").rstrip()
                qi2 = lines[2].decode("utf-8").rstrip()     
                qr2 = lines[3].decode("utf-8").rstrip()

                #qr1 = lines[0].rstrip()
                #qi1 = lines[1].rstrip()
                #qi2 = lines[2].rstrip()     
                #qr2 = lines[3].rstrip()
               
		# Update the sum list as you go through the records
                qr1_list, q1_count = updateQuality(qr1, qr1_list)
                num_qr1 += q1_count
                qr2_list, q2_count = updateQuality(qr2, qr2_list)
                num_qr2 += q2_count
                qi1_list, q3_count = updateQuality(qi1, qi1_list)
                num_qi1 += q3_count
                qi2_list, q4_count = updateQuality(qi2, qi2_list)
                num_qi2 += q4_count
	
	# Average the Quality Scores
        qr1_list = averageQualityScore(qr1_list, num_qr1)
        qr2_list = averageQualityScore(qr2_list, num_qr2)
        qi1_list = averageQualityScore(qi1_list, num_qi1)
        qi2_list = averageQualityScore(qi2_list, num_qi2)
        print(qr1_list)
        print(qr2_list)
        print(qi1_list)
        print(qi2_list)
            
def main():
    #file_list = ['r1.txt', 'i2.txt', 'i3.txt', 'r4.txt']
    #files = [open('{}'.format(filename)) for filename in file_list]
    file_list = ['1294_S1_L008_R1_001.fastq.gz', '1294_S1_L008_R2_001.fastq.gz', '1294_S1_L008_R3_001.fastq.gz', '1294_S1_L008_R4_001.fastq.gz']
    files = [gzip.open('/projects/bgmp/shared/2017_sequencing/{}'.format(filename)) for filename in file_list]
    deplex = Demultiplexing(files)
    
if __name__ == "__main__":
	main()
