
import os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from matplotlib.ticker import PercentFormatter

from operator import itemgetter


'''
I am creating this file to cehck the MAPQ (mapping quality) in a BAM file. 
This file should produce a distribution of all the MAPQ values. 
'''


file="/home/madhu/datasets/analysis/p3_m6A_RNA_modification_native_RNA_seq/GSM3528749/GSM3528749.bam"
# samtools view "$file" | cut -f5 | sort | uniq -c


stream = os.popen('samtools view {} | cut -f5 | sort | uniq -c'.format(file))
output = stream.read()

# print('The len of output: ', (output[0]))

# print(output[0:4])

# print(output.split('\n')[0].split(' '))
# print(output.split('\n')[1].split(' '))
# print(output.split('\n')[2].split(' '))
# print(output.split('\n')[len(output.split('\n'))-2].split(' '))

## The outputs are: 

'''
(ACONDA) madhu@qgenlabgpu:~$ python /home/madhu/datasets/scripts/chris_scripts/check_MAPQ.py
['', '', '14629', '0']
['', '', '', '', '110', '1']
['', '', '', '', '', '46', '10']
['', '', '', '', '', '32', '9']

'''


# def remove_empty_strings(string):
#     return string != ""
 
# test_list = ["", "GeeksforGeeks", "", "is", "best", ""]
# filtered_list = filter(remove_empty_strings, test_list)
# print(list(filtered_list))  # Output: ['GeeksforGeeks', 'is', 'best']





list_count = []
list_mapQ = []

for a_line in output.split('\n'):
    if a_line == '':
        continue
    s_line = a_line.split(' ')
    while('' in s_line):
        s_line.remove('')
#     print('line: ', s_line)
    count = s_line[0]
    mapQ = s_line[1]
#     print('count: ', count)
#     print('mapQ: ', mapQ)
#     print('\n\n')
    list_count.append(int(count))
    list_mapQ.append(int(mapQ))
    

# print(sorted(tuple((list_count,list_mapQ))))

# x = sorted(tuple((list_count,list_mapQ)))[0]
# y = sorted(tuple((list_count,list_mapQ)))[1]

# print(x)
# print(y)


# print(np.array((list_count,list_mapQ)))

# print(np.array((list_count,list_mapQ)).reshape(len(list_count), 2))

a = np.array(list_mapQ).reshape(len(list_count), 1)
b = np.array(list_count).reshape(len(list_mapQ), 1)

c = np.hstack((a,b))

d = sorted(list(map(list, c)), key = itemgetter(0))

# print(d[0][:])

print([x[0] for x in d])

# print(sorted(list(map(list, c)), key = itemgetter(0)))

# print(sorted(list(c)))

# print(sorted(list(map(list, c)), key=lambda x: x[0]))

# print(np.sort(np.hstack((a,b)), axis=0))

# print(np.hstack((np.array(list_count).reshape(len(list_count), 1), np.array(list_mapQ).reshape(len(list_mapQ), 1))))

# print(np.sort(np.vstack((np.array(list_count).reshape(len(list_count),), np.array(list_mapQ).reshape(len(list_mapQ),)))))

# print(tuple(np.array((list_count,list_mapQ)).reshape(2, len(list_count))))



# # # Creating dataset
# # np.random.seed(23685752)
# # N_points = 10000
# # n_bins = 20
 
# Creating distribution
x = [x[0] for x in d]
y = [x[1] for x in d]

# # # Creating histogram
# # fig, axs = plt.subplots(1, 1,
# #                         figsize =(10, 7),
# #                         tight_layout = True)
 
# # axs.hist(x, bins = 20)


plt.plot(x, y)

# Show plot
plt.show()

name = file.replace('.bam', '.png')

plt.savefig(name)
plt.close()




