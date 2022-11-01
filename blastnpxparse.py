#!/usr/bin/python3

# 2022-10-21 version 8.1: blastnpxparse.py

import sys
import os
import re
import subprocess
import csv
from bs4 import BeautifulSoup
import collections

# command verification 
args = sys.argv

if len(args) !=2:
    sys.exit("\nError\nUsage: python3 [program].py [input file name]\n")

# initial parameters
file_input = args[1]
file_prefix_init = os.path.splitext(os.path.basename(file_input))[0]
file_conversion = '0'
path_input = os.getcwd()
annotation_file = file_prefix_init + '.ann'
ann_file_exist = '0'
BREAD_CRUMB_file = '/home/makawazoe/TXSearchfile/220916_TAX_ID__BREAD_CRUMB_test.tsv'

# file handling
query_taxid =[]
query_organism = []
if os.path.exists(annotation_file) == True:
    ann_file_exist = '1'
    print("\nAnnotation file exists.\n\n\nReading file...\n")
    with open(annotation_file, encoding='utf-8', newline='') as f:
        annotation_tsv = csv.reader(f, delimiter='\t')
        header = next(annotation_tsv)      
        for row in annotation_tsv:
            if 'db_xref' in row:
                p1 = r'taxon\:(.*)'
                taxid_number = re.findall(p1, row[-1])
                query_taxid = query_taxid + taxid_number
            elif 'organism' in row:
                query_organism = query_organism + row[-1:]
            else:
                pass
elif os.path.exists(annotation_file) == False:
    print("\nAnnotation file doesn't exist. Query information will be exracted only from blast result file.\n")
    pass

with open(file_input, "r", encoding="utf-8") as f:
    headline = f.readline()
    if 'xml' in headline:
        print("\nInput file type is XML.\n\n\nReading file ...\n")
        file_prefix = file_prefix_init
        blastxml = f.read()
        if 'blastn' in blastxml:
            mode = 'blastn'
            print("\nblastn output file confirmed.\n")
        elif 'blastp' in blastxml:
            mode = 'blastp'
            print("\nblastp output file confirmed.\n")
        elif 'blastx' in blastxml:
            mode = 'blastx'
            print("\nblastx output file confirmed.\n")
    elif 'Blast4-archive' in headline:
        file_conversion = '1'
        print("\nInput file type is ASN.1.\n\n\nStarting file conversion to XML (Single-file BLAST XML2)\n")
        pass
    else:
        sys.exit("\nError\nInput is neither XML nor ASN.1 file.\nRequirement: blast[n, p, x] output file with run command option '-outfmt 11' (BLAST archive (ASN.1)) or '-outfmt 16' (Single-file BLAST XML2)\n")

if file_conversion == '1':
    outfmt = ['0','16'] 
    for fmt in outfmt:
        file_conv = file_prefix_init + '_' + fmt + os.path.splitext(os.path.basename(file_input))[1]
        file_prefix = file_prefix_init + '_' + fmt
        cmd_list = ['blast_formatter','-archive',file_input,'-out',file_conv,'-outfmt',fmt]
        cmd = map(str, cmd_list)
        blastformatter_run = subprocess.run(cmd)
        print("\nblast_formatter finished.\n\n\nReading file ...\n")
        if fmt == '16':
            with open(file_conv, "r", encoding="utf-8") as f:
                blastxml = f.read()
                if 'blastn' in blastxml:
                    mode = 'blastn'
                    print("\nblastn output file confirmed.\n")
                elif 'blastp' in blastxml:
                    mode = 'blastp'
                    print("\nblastp output file confirmed.\n")
                elif 'blastx' in blastxml:
                    mode = 'blastx'
                    print("\nblastx output file confirmed.\n")
                else:
                    sys.exit("\nError\nInput is not blast[n, p, x] file.\nRequirement: blast[n, p, x] output file with run command option '-outfmt 11' (BLAST archive (ASN.1)) or '-outfmt 16' (Single-file BLAST XML2)\n")
        else:
            pass
elif file_conversion == '0':
    pass

# soup
soup_1 = BeautifulSoup(blastxml, "lxml-xml")
result_queries = soup_1.find_all('BlastOutput2')

# extract and file output
output_file1 = file_prefix + '_' + mode + 'table.tsv'
output_file2 = file_prefix + '_' + mode + 'summary.tsv'
header1 = ["query acc", "query length", "query taxid", "query scientific name", "hit number", "hit acc", "taxid", "scientific name", "hit entry length", "alignment#", "bit-score", "score", "E-value", "identity", "query strand", "hit strand", "query from", "query to", "hit from", "hit to", "alignment length", "gap", "coverage (%)", "identity (%)", "hit title"]
header2 = ["query acc", "query length", "hit number", "hit acc", "hit title", "taxid", "scientific name", "hit entry length", "alignment#", "bit-score", "score", "E-value", "identity", "positive", "query from", "query to", "hit from", "hit to", "alignment length", "gap", "coverage (%)", "identity (%)"]
header3 = ["query acc", "query length", "hit number", "hit acc", "hit title", "taxid", "scientific name", "hit entry length", "alignment#", "bit-score", "score", "E-value", "identity", "positive", "query frame", "query from", "query to", "hit from", "hit to", "alignment length", "gap", "coverage (%)", "identity (%)"]

if mode == 'blastn':
    with open(output_file1, "a", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(header1)
elif mode == 'blastp':
    with open(output_file1, "a", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(header2)
elif mode == 'blastx':
    with open(output_file1, "a", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(header3)

print("\nExtracting ...\n")

for loop_count, result_query in enumerate(result_queries):
    if result_query.find('message') is None:
        pass
    elif result_query.find('message').text == 'No hits found':
        query_acc = result_query.find('query-title').text
        query_summary = str.splitlines(query_acc) + str.splitlines('No hits found')
        with open(output_file2, "a", encoding="utf-8") as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(query_summary)
        continue
    query_acc = result_query.find('query-title').text
    query_len = result_query.find('query-len').text
    if ann_file_exist == '1' and mode == 'blastn':
        query_info = str.splitlines(query_acc) + str.splitlines(query_len) + query_taxid[loop_count:loop_count+1] + query_organism[loop_count:loop_count+1]
    else:
        query_info = str.splitlines(query_acc) + str.splitlines(query_len)
    query_hits = result_query.find_all('Hit')
    if mode == 'blastn':
        hit_sciname_total = [] # blastn specific
        hit_taxid_total = []   # blastn specific
    else:
        hit_title_total = []   # blastp, blastx
    hit_num_total = []
    for query_hit in query_hits:
        hit_num = query_hit.find('num').text
        hit_num_total.append(hit_num)
        if query_hit.find('accession') is None:
            hit_acc = query_hit.find('id').text
        else:
            hit_acc = query_hit.find('accession').text
        if query_hit.find('title') is None:
            hit_title = "null"
            if mode != 'blastn':
                hit_title_total.append(hit_title)
            else:
                pass
        else:
            hit_title_full = query_hit.find('title').text.split(';') # blastp, blastx
            hit_title = hit_title_full[0:1]                          # blastp, blastx
            if mode != 'blastn':
                hit_title_total = hit_title_total + hit_title        # blastp, blastx
            else:
                pass
        if query_hit.find('taxid') is None:
            hit_taxid = "null"
        else:
            hit_taxid = query_hit.find('taxid').text
        if query_hit.find('sciname') is None:
            hit_sciname = "null"
        else:
            hit_sciname = query_hit.find('sciname').text
        if mode == 'blastn':
            hit_sciname_total.append(hit_sciname) # blastn specific
            hit_taxid_total.append(hit_taxid)     # blastn specific
        else:
            pass
        hit_len = query_hit.find('len').text
        if mode == 'blastn':
            hit_info = str.splitlines(hit_num) + str.splitlines(hit_acc) + str.splitlines(hit_taxid) + str.splitlines(hit_sciname) + str.splitlines(hit_len)             # blastn specific
        elif type(hit_title) == list:
            hit_info = str.splitlines(hit_num) + str.splitlines(hit_acc) + hit_title + str.splitlines(hit_taxid) + str.splitlines(hit_sciname) + str.splitlines(hit_len) # blastp, blastx
        elif type(hit_title) == str:
            hit_info = str.splitlines(hit_num) + str.splitlines(hit_acc) + str.splitlines(hit_title) + str.splitlines(hit_taxid) + str.splitlines(hit_sciname) + str.splitlines(hit_len) # blastp, blastn
        hit_hsps = query_hit.find_all('hsps')
        for hit_hsp in hit_hsps:
            alinum = []
            bitscore = []
            score = []
            evalue = []
            ident = []
            query_f = []
            query_t = []
            subj_f = []
            subj_t = []
            alilen = []
            gaps = []
            for hsp_num in hit_hsp.find_all('num'):
                alinum.append(hsp_num.text)
            for hsp_bitscore in hit_hsp.find_all('bit-score'):
                bitscore.append(hsp_bitscore.text)
            for hsp_score in hit_hsp.find_all('score'):
                score.append(hsp_score.text)
            for hsp_eval in hit_hsp.find_all('evalue'):
                evalue.append(hsp_eval.text)
            for hsp_ident in hit_hsp.find_all('identity'):
                ident.append(hsp_ident.text)
            for hsp_query_f in hit_hsp.find_all('query-from'):
                query_f.append(hsp_query_f.text)
            for hsp_query_t in hit_hsp.find_all('query-to'):
                query_t.append(hsp_query_t.text)
            for hsp_subj_f in hit_hsp.find_all('hit-from'):
                subj_f.append(hsp_subj_f.text)
            for hsp_subj_t in hit_hsp.find_all('hit-to'):
                subj_t.append(hsp_subj_t.text)
            for hsp_alilen in hit_hsp.find_all('align-len'):
                alilen.append(hsp_alilen.text)
            for hsp_gaps in hit_hsp.find_all('gaps'):
                gaps.append(hsp_gaps.text)
            if mode == 'blastn':
                query_s = [] # blastn specific
                subj_s = []  # blastn specific
                for hsp_query_s in hit_hsp.find_all('query-strand'): # blastn specific
                    query_s.append(hsp_query_s.text)                 # blastn specific
                for hsp_subj_s in hit_hsp.find_all('hit-strand'):    # blastn specific
                    subj_s.append(hsp_subj_s.text)                   # blastn specific
            elif mode == 'blastp':
                positive = [] # blastp, blastx
                for hsp_positive in hit_hsp.find_all('positive'):    # blastp, blastx
                    positive.append(hsp_positive.text)               # blastp, blastx
            elif mode == 'blastx':
                positive = [] # blastp, blastx
                query_fr = [] # blastx specific
                for hsp_positive in hit_hsp.find_all('positive'):    # blastp, blastx
                    positive.append(hsp_positive.text)               # blastp, blastx
                for hsp_query_fr in hit_hsp.find_all('query-frame'): # blastx specific
                    query_fr.append(hsp_query_fr.text)               # blastx specific
            else:
                pass
        cont_num = len(list(alinum))
        identity = [x / y * 100 for (x, y) in zip(map(int, ident), map(int, alilen))]
        for i in range(cont_num):
            if mode == 'blastx':
                aa_conv_len = int(query_len) // 3                                # blastx specific
                coverage = ((int(alilen[i]) - int(gaps[i])) / aa_conv_len * 100) # blastx specific
            elif mode != 'blastx':
                coverage = ((int(alilen[i]) - int(gaps[i])) / int(query_len) * 100)
            coverage_f = "{:.0f}".format(coverage)
            identity_f = "{:.0f}".format(identity[i])
            if mode == 'blastn' and type(hit_title) == list:
                tsv_out = query_info + hit_info + alinum[i:i+1] + bitscore[i:i+1] + score[i:i+1] + evalue[i:i+1] + ident[i:i+1] + query_s[i:i+1] + subj_s[i:i+1] + query_f[i:i+1] + query_t[i:i+1] + subj_f[i:i+1] + subj_t[i:i+1] + alilen[i:i+1] + gaps[i:i+1] + str.splitlines(coverage_f) + str.splitlines(identity_f) + hit_title
            elif mode == 'blastn' and type(hit_title) == str:
                tsv_out = query_info + hit_info + alinum[i:i+1] + bitscore[i:i+1] + score[i:i+1] + evalue[i:i+1] + ident[i:i+1] + query_s[i:i+1] + subj_s[i:i+1] + query_f[i:i+1] + query_t[i:i+1] + subj_f[i:i+1] + subj_t[i:i+1] + alilen[i:i+1] + gaps[i:i+1] + str.splitlines(coverage_f) + str.splitlines(identity_f) + str.splitlines(hit_title)
            elif mode == 'blastp':
                tsv_out = query_info + hit_info + alinum[i:i+1] + bitscore[i:i+1] + score[i:i+1] + evalue[i:i+1] + ident[i:i+1] + positive[i:i+1] + query_f[i:i+1] + query_t[i:i+1] + subj_f[i:i+1] + subj_t[i:i+1] + alilen[i:i+1] + gaps[i:i+1] + str.splitlines(coverage_f) + str.splitlines(identity_f)
            elif mode == 'blastx':
                tsv_out = query_info + hit_info + alinum[i:i+1] + bitscore[i:i+1] + score[i:i+1] + evalue[i:i+1] + ident[i:i+1] + positive[i:i+1] + query_fr[i:i+1] + query_f[i:i+1] + query_t[i:i+1] + subj_f[i:i+1] + subj_t[i:i+1] + alilen[i:i+1] + gaps[i:i+1] + str.splitlines(coverage_f) + str.splitlines(identity_f)
            with open(output_file1, "a", encoding="utf-8") as f:
                writer = csv.writer(f, delimiter='\t')
                writer.writerow(tsv_out)

    if mode == 'blastn':
        total1 = collections.Counter(hit_taxid_total)   # blastn specific
        if query_taxid[loop_count:loop_count+1] == ["0"]:
            query_summary = str.splitlines(query_acc) + str.splitlines(hit_num_total[-1]) + query_taxid[loop_count:loop_count+1] + query_organism[loop_count:loop_count+1]
            print("\n")
            print("Entry: " + query_acc + "\n")
            print("Number: " + str(loop_count+1) + "\n")
            print("Query organism: " + query_organism[loop_count:loop_count+1] + "\n")
        else:
            query_search_string1 = '"^' + query_taxid[loop_count] + r'\s"'
#            print(query_search_string1)
            grep_cmd_list1 = ["grep"] + str.split(eval(query_search_string1)) + str.split(BREAD_CRUMB_file)
            grep_cmd1 = map(str, grep_cmd_list1)
            grep_run1 = subprocess.run(grep_cmd1, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
#            grep_run1 = subprocess.run(['grep', query_search_string1, BREAD_CRUMB_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            print("\n")
            print("Entry: " + query_acc + "\n")
            print("Number: " + str(loop_count+1) + "\n")
            print(*grep_run1.stdout.split(':::')[1::2])
            query_summary = str.splitlines(query_acc) + str.splitlines(hit_num_total[-1]) + query_taxid[loop_count:loop_count+1] + grep_run1.stdout.split(':::')[1::2]
            print("\n")
            print("Extracting lineage ...\n")
    else:
        query_summary = str.splitlines(query_acc) + str.splitlines(hit_num_total[-1]) + query_taxid[loop_count:loop_count+1]
        total1 = collections.Counter(hit_title_total)   # blastp, blastx

    with open(output_file2, "a", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(query_summary)
        for k, v in total1.items():
            hit_summary = [""]
            if mode != 'blastn':
                hit_summary = hit_summary + str.splitlines(str(v)) + str.splitlines(str(k))
                writer.writerow(hit_summary)
            elif k == 'null':
                hit_summary = hit_summary + str.splitlines(str(v)) + str.splitlines(str(k))
                writer.writerow(hit_summary)
            else:
                hit_search_string1 = '"^' + str(k) + r'\s"'
                grep_cmd_list2 = ["grep"] + str.split(eval(hit_search_string1)) + str.split(BREAD_CRUMB_file)
                grep_cmd2 = map(str, grep_cmd_list2)
                grep_run2 = subprocess.run(grep_cmd2, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                hit_summary = hit_summary + str.splitlines(str(v)) + str.splitlines(str(k)) + grep_run2.stdout.split(':::')[1::2]
                writer.writerow(hit_summary)  

print("\nDone\n")
