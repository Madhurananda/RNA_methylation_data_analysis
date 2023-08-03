
# import pandas as pd
# import openpyxl

# # '''
# # I have tried with the webcrawler, but didn't work. 
# # Now, I have copied it to a Excel file and trying to get the hyperlinks ... 
# # '''



# xls_file = '/home/madhu/datasets/download/p5_decoding_epitranscriptional_native_RNA_seq/original_data/Paper5_datasets.xlsx'

# sh_file = '/home/madhu/datasets/download/p5_decoding_epitranscriptional_native_RNA_seq/original_data/wget_fast5.sh'

# dataframe1 = pd.read_excel(xls_file)
 
# # print(len(dataframe1))




# wb = openpyxl.load_workbook(xls_file)
# ws = wb['Sheet1']

# f = open(sh_file, "w")
# for a_row in range(len(dataframe1)):
#     # This will fail if there is no hyperlink to target
#     try:
#         hLink = ws.cell(row=(a_row+1), column=3).hyperlink.target
#         if 'FAST5' in hLink.upper():
#             ## Add another check depending on the dataset: 
#             if not 'CDNA' in hLink.upper():
#               print('The link is: ', hLink)
#               f.write("wget "+hLink+" &" +"\n" )
#     except:
#         hLink = ws.cell(row=(a_row+1), column=3).value

# f.write("wait\n")
# f.close()


# '''
# This was what I got from Houssem, but it couln't download anything from Javascript websites. 
# '''


import requests
import lxml
from bs4 import BeautifulSoup

'''
I got this script from Houssem. I am using this to download a long list of datasets. 
'''

# url = "https://github.com/nanopore-wgs-consortium/NA12878/blob/master/Genome.md"
url = 'https://www.ebi.ac.uk/ena/browser/view/PRJEB32782'

headers = {
  'User-Agent': 'Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/63.0.3239.132 Safari/537.36 QIHU 360SE'
}

f = requests.get(url, headers = headers)

soup = BeautifulSoup(f.content,'lxml')


reads = soup.find_all("//")

print(reads)

# f = open("wget_fast5.sh", "w")

# for fast5 in reads:
#     f.write("wget "+fast5["href"]+" &" +"\n" )
# f.write("wait")
# f.close()

