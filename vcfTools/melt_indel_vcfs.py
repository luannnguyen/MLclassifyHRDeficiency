#!/usr/bin/env python

debug = 1

#Test donor: DO51959
#Test donor: DO27821

#====== Import packages ======#
## Common
import argparse
import sys
import os
import re
import numpy as np
import pandas as pd
import gzip

## Special
import vcf as pyvcf
from natsort import natsorted,ns
from collections import Counter


#====== Load vcf files (hpc; manifest) ======#
#------ Arguments ------#
parser = argparse.ArgumentParser(description='Merge indel vcfs generated from indelCalling, snowman and svcp')

parser.add_argument('-m', '--manifest_path', help = '*Absolute* path to a text file that lists the input vcf file(s)', required=True)
parser.add_argument('-d', '--donor_tsv_path', help = 'Path to the donor.tsv flat-file database', required=True)
parser.add_argument('-sm1', '--skipmelt1',  action='store_true')
parser.add_argument('-sm2', '--skipmelt2',  action='store_true')

args = parser.parse_args()


#------ Parse manifest ------#
if debug == 1:
    print 'Loading vcfs...'
#manifest_path = '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/merged_vcfs/indels/DO27821/DO27821_indelVcfManifest.tsv'
manifest_path = args.manifest_path
donor_tsv_path = args.donor_tsv_path

vcf_list = []
with open(manifest_path) as manifest:
    vcf_file = manifest.readlines()
    vcf_list.extend(vcf_file)

vcf_files = [vcf.strip('\n') for vcf in vcf_list]
#vcf_files = [os.path.basename(vcf) for vcf in vcf_files]

out_dir = os.path.dirname(manifest_path) # output dir is the dirname of the manifest file

#------ Dictionary of vcfs split by caller ------#
indel_callers = ['indelCalling', 'snowman', 'svcp']
vcfs_by_caller = dict.fromkeys(indel_callers)
for caller in indel_callers:
    indices = [c for c, v in enumerate(vcf_files) if caller in v]
    vcfs_by_caller[caller] = [vcf_files[c] for c in indices]

#------ Get donor id from manifest ------#
DO_re = re.compile(r'DO\d+')
donor_id = DO_re.search(os.path.basename(manifest_path)).group(0)

# #====== Load vcf files (local) ======#
# #------ Arguments ------#
# parser = argparse.ArgumentParser(description='Merge indel vcfs generated from indelCalling, snowman and svcp')
# parser.add_argument('-v', '--vcf_dir', help = 'Directory of the input vcf file(s)', required=True)
# parser.add_argument('-o', '--out_dir', help = 'Output directory', required=True)
# parser.add_argument('-d', '--donor_tsv_path', help = 'Path to the donor.tsv flat-file database', required=True)
#
# args = parser.parse_args()
# vcf_dir = args.vcf_dir
# out_dir = args.out_dir
# donor_tsv_path = args.donor_tsv_path
#
# #------ Paths ------#
# vcf_dir = '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/downloads_subset/indel/head'
# out_dir = '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/downloads_subset/indel/output'
# donor_tsv_path = '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/metadata/donor_pcawg_london.tsv'
#
# vcf_files = [vcf for vcf in os.listdir(vcf_dir) if not vcf.startswith('._') and vcf.endswith('.vcf')]
#
# indel_callers = ['indelCalling', 'snowman', 'svcp']
# vcfs_by_caller = dict.fromkeys(indel_callers)
# for caller in indel_callers:
#     indices = [c for c, v in enumerate(vcf_files) if caller in v]
#     vcfs_by_caller[caller] = [vcf_files[c] for c in indices]

#====== Global functions ======#
if debug == 1:
    print 'Loading global functions...'

#------ VCF reader ------#
def read_vcf(vcf_file):
    return pyvcf.Reader(open(vcf_file, 'r'))
    #return pyvcf.Reader(open(vcf_dir + '/' + vcf_basename, 'r'))


#------ Get vcf file_id and donor_id from donor tsv -----#
df_donor = pd.read_table(donor_tsv_path, sep='\t', usecols=[0,1,2,4,5,6])
def get_vcf_info(vcf, colname):
    # 'File.ID','ICGC.Donor'
    vcf_basename = os.path.basename(vcf)
    file_info = df_donor.loc[df_donor['File.Name'].str.contains(vcf_basename), colname].values[0]
    return file_info

#donor_tsv_path = '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/metadata/donor_pcawg_london.tsv'
#get_vcf_info(vcfs_by_caller['indelCalling'][0], 'ICGC.Donor')

#------ List as strings ------#
def list_as_string(list_object, delim = ','):
    return delim.join(str(element) for element in list_object)

#====== Vcf headers ======#
if debug == 1:
    print 'Loading vcf headers...'

vcf_header_indelCalling='''##fileformat=VCFv4.1
##caller="indelCalling"
##INFO=<ID=FILEID,Number=.,Type=String,Description="File IDs of contributing vcfs">
##FILTER=<ID=GOF,Description="Variant fails goodness-of-fit test.">
##FILTER=<ID=badReads,Description="Variant supported only by reads with low quality bases close to variant position, and not present on both strands.">
##FILTER=<ID=alleleBias,Description="Variant frequency is lower than expected for het">
##FILTER=<ID=Q20,Description="Variant quality is below 20.">
##FILTER=<ID=HapScore,Description="Too many haplotypes are supported by the data in this region.">
##FILTER=<ID=MQ,Description="Root-mean-square mapping quality across calling region is low.">
##FILTER=<ID=strandBias,Description="Variant fails strand-bias filter">
##FILTER=<ID=SC,Description="Variants fail sequence-context filter. Surrounding sequence is low-complexity">
##FILTER=<ID=QD,Description="Variants fail quality/depth filter.">
##FILTER=<ID=ALTC,Description="Alternative reads in control and other filter not PASS">
##FILTER=<ID=VAF,Description="Variant allele frequency in tumor < 10% and other filter not PASS">
##FILTER=<ID=VAFC,Description="Variant allele frequency in tumor < 5% or variant allele frequency in control > 5%">
##FILTER=<ID=QUAL,Description="Quality of entry too low and/or low coverage in region">
##FILTER=<ID=ALTT,Description="Less than three variant reads in tumor">
##FILTER=<ID=GTQ,Description="Quality for genotypes below thresholds">
##FILTER=<ID=GTQFRT,Description="Quality for genotypes below thresholds and variant allele frequency in tumor < 10%">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AF,Number=1,Type=Float,Description="Mean tumor ALT allele frequency across replicate vcfs">
##SAMPLE=<ID=TUMOR,SampleName=tumor_NA,Individual=NA,Description="Tumor">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TUMOR
'''

vcf_header_snowman='''##fileformat=VCFv4.1
##caller="snowman"
##INFO=<ID=FILEID,Number=.,Type=String,Description="File IDs of contributing vcfs">
##FILTER=<ID=PASS,Description="3+ split reads, 60 contig MAPQ">
##FILTER=<ID=REPEAT,Description="3+ split reads, 60 contig MAPQ">
##FILTER=<ID=WEAKCIGARMATCH,Description="For indels <= 5 bp, require 8+ split reads or 4+ and 3+ cigar matches">
##FILTER=<ID=WEAKASSEMBLY,Description="4 or fewer split reads">
##FILTER=<ID=LOWMAPQ,Description="Assembly contig has less than MAPQ 60">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AF,Number=1,Type=Float,Description="Mean tumor ALT allele frequency across replicate vcfs">
##SAMPLE=<ID=TUMOR,SampleName=tumor_NA,Individual=NA,Description="Tumor">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TUMOR
'''

vcf_header_svcp='''##fileformat=VCFv4.1
##caller="indelCalling"
##INFO=<ID=FILEID,Number=.,Type=String,Description="File IDs of contributing vcfs">
##FILTER=<ID=F004,Description="Tum medium read depth strand bias check: Calls In 8% Reads Bt Depth 10 And 200 (inclusive)">
##FILTER=<ID=F005,Description="Tum high read depth strand bias check: Calls In 4% Reads > Depth 200">
##FILTER=<ID=F006,Description="Small call excessive repeat check: Fail if Length <= 4 and Repeats > 9">
##FILTER=<ID=F010,Description="Variant must not exist within the Unmatched Normal Panel">
##FILTER=<ID=F012,Description="Germline: When length < 11 and depth > 9, fail if the variant is seen in both 20% of normal reads AND 20% of tumour reads in either pindel or bwa">
##FILTER=<ID=F018,Description="Sufficient Depth: Pass if depth > 10">
##FILTER=<ID=F015,Description="No normal calls">
##FILTER=<ID=F016,Description="Verify indel condition">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AF,Number=1,Type=Float,Description="Mean tumor ALT allele frequency across replicate vcfs">
##SAMPLE=<ID=TUMOR,SampleName=tumor_NA,Individual=NA,Description="Tumor">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TUMOR
'''

vcf_header_final = '''##fileformat=VCFv4.1
##INFO=<ID=NCALLERS,Number=1,Type=Integer,Description="Number of callers reporting variant">
##INFO=<ID=FILEID_indelCalling,Number=.,Type=String,Description="File IDs of contributing vcfs produced by indelCalling">
##INFO=<ID=FILEID_snowman,Number=.,Type=String,Description="File IDs of contributing vcfs produced by snowman">
##INFO=<ID=FILEID_svcp,Number=.,Type=String,Description="File IDs of contributing vcfs produced by svcp">
##FILTER=<ID=PASS,Description="PASS present once in FILTER">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AF,Number=1,Type=Float,Description="Tumor ALT allele frequency">
##SAMPLE=<ID=TUMOR,SampleName=tumor_NA,Individual=NA,Description="Tumor">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TUMOR
'''

#%%%%%%%%%%%%%%%%%%%%%%%% Melt #1 %%%%%%%%%%%%%%%%%%%%%%%%#
if not args.skipmelt1:

    #====== Melt #1 functions ======#
    #------ Get relevant vcf fields by caller as dictionary ------#
    def merge_vcfs_by_caller(dict_vcfs_by_caller, caller):
        # Dictionary structure:
        # coord_key: CHROM, POS, REF, ALT
        # FILTER: [+], INFO: [FileID], TUMOR:[GT,AF]
        d_merged_records = {}

        #dict_vcfs_by_caller = vcfs_by_caller
        #caller = 'svcp'
        for vcf in dict_vcfs_by_caller[caller]:
            #vcf_basename = vcfs_by_caller[caller][0]

            vcf_reader = read_vcf(vcf)

            file_id = get_vcf_info(vcf, 'File.ID')

            for record in vcf_reader:
                #...... Get coordinates ......#
                coordinates = [record.CHROM, record.POS, record.REF, str(record.ALT[0])]
                coord_key = list_as_string(coordinates)

                if debug == 1:
                    print 'Populating dictionary of variants (melt #1) | File ID: %s | Position: %s' %(file_id, coord_key)

                #...... Common metadata ......#
                if coord_key not in d_merged_records:
                    # FILTER[+], INFO[FileID], TUMOR[GT,AF]
                    d_merged_records[coord_key] = [[], [], []]

                ## 0. FILTER: [+]; PASS is not reported by pyvcf

                if len(record.FILTER) == 0:
                    d_merged_records[coord_key][0].append('PASS')
                else:
                    d_merged_records[coord_key][0].append(list_as_string(record.FILTER))

                ## 1. INFO: [FileID]
                d_merged_records[coord_key][1].append(file_id)

                #...... Caller specific meta data ......#
                if caller == 'indelCalling':

                    ## 2. TUMOR: [GT, AF]
                    genotype = record.samples[1]['GT']
                    AF = round( float(record.samples[1]['NV']) / float(record.samples[1]['NR']) ,4)

                    d_merged_records[coord_key][2].append( [genotype, AF] )

                elif caller == 'snowman':
                    ## 2. TUMOR: [GT, AF]
                    genotype = './.' # Snowman has no tumor genotype data; this part of script can be optimized
                    AF = round( float(record.INFO['TFRAC'] ) ,4)

                    d_merged_records[coord_key][2].append( [genotype, AF] )

                elif caller == 'svcp':
                    ## 2. TUMOR: [GT, AF]
                    genotype = './.' # svcp reports ./. for all genotypes; this part of script can be optimized

                    alt_reads = float(record.samples[1]['PP']) + float(record.samples[1]['NP'])
                    total_reads = float(record.samples[1]['PR']) + float(record.samples[1]['NR'])

                    #AF = round(alt_reads / total_reads, 4)
                    #AF = float(0)
                    if total_reads == 0:
                        AF = float(0) # prevent divide by zero
                    else:
                        AF = round(alt_reads/total_reads, 4)

                    d_merged_records[coord_key][2].append( [genotype, AF] )

        ## Output preliminary dictionary
        return d_merged_records

    # merged_indelCalling_vcfs = merge_vcfs_by_caller(vcfs_by_caller, 'indelCalling')
    # merged_indelCalling_vcfs

    #------ Write vcf lines; merge replicate vcfs by caller ------#
    def dict_record_as_line(k, v, caller):
        ## coordinates
        l_coords = k.split(',')
        CHROM = l_coords[0]
        POS = l_coords[1]
        REF = l_coords[2]
        ALT = l_coords[3]

        ## metadata
        ID = '.'
        QUAL = '.'
        FILTER = list_as_string(v[0], delim=';') #FIlTER
        INFO = 'FILEID=' + list_as_string(v[1]) # Convert FileID lists to strings
        FORMAT = 'GT:AF'

        if caller == 'indelCalling':
            genotype_most_freq = Counter([l[0] for l in v[2]]).most_common(1)[0][0]
        else:
            genotype_most_freq = './.'

        AF_mean = round(np.mean([l[1] for l in v[2]]), 4)

        TUMOR = genotype_most_freq + ':' + str(AF_mean)

        ## print lines
        return '\t'.join([CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,TUMOR])

    #====== Execute merge and write #1 ======#
    #indel_callers = ['indelCalling'] #debugging
    for caller in indel_callers:

        # DO_re = re.compile(r'DO\d+')
        # donor_id = DO_re.search(os.path.basename(manifest_path)).group(0)

        melt1_vcf_path = '%s/%s%s%s%s' %(out_dir, donor_id, '_', caller, '_merged.vcf.gz')
        #melt1_vcf_path = '%s/%s/%s%s%s%s' %(out_dir, donor_id, donor_id, '_', caller, '_merged.vcf')

        with gzip.open(melt1_vcf_path, 'wb', 9) as melt1_vcf:
        #with open(melt1_vcf_path, 'w') as melt1_vcf:

            ## Write header
            vcf_header = eval('%s%s' %('vcf_header_', caller))
            melt1_vcf.write(vcf_header)

            ## Write body
            dict_merged_vcfs_by_caller = merge_vcfs_by_caller(vcfs_by_caller, caller)

            ## Progress counter
            total_variants = len(dict_merged_vcfs_by_caller)
            progress_counter = 0

            if debug == 1:
                print 'Sorting dictionary (melt #1)...'

            for k, v in natsorted(dict_merged_vcfs_by_caller.iteritems()):
                progress_counter = progress_counter + 1

                if debug == 1:
                    print 'Writing lines (melt #1) | Progress: %s/%s | Position: %s' %(progress_counter, total_variants, k)

                line = dict_record_as_line(k, v, caller) + '\n'
                #print line
                melt1_vcf.write(line)

    del dict_merged_vcfs_by_caller  # remove large dict object from memory

    #====== Write 'done' file if script completes ======#
    if debug == 1:
        print 'Writing done_melt1.txt'
    done_melt1_txt = open('%s/%s' %(out_dir, 'done_melt1.txt'), 'w')
    done_melt1_txt.close()



    # #====== Test read vcf ======#
    # test_vcf = '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/downloads_subset/indel/output/DO51959/DO51959_indelCalling_merged.vcf.gz'
    #
    # test_vcf_reader = pyvcf.Reader(open(test_vcf, 'r'))
    # #test_record = next(test_vcf_reader)
    #
    # for i in test_vcf_reader:
    #     print i


#%%%%%%%%%%%%%%%%%%%%%%%% Melt #2 %%%%%%%%%%%%%%%%%%%%%%%%#
if not args.skipmelt2:

    #====== Load merged vcfs (by caller) ======#
    #manifest_path = '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/merged_vcfs/indels/DO51959/DO51959_indelVcfManifest.tsv'
    #out_dir = os.path.dirname(manifest_path)

    merged_vcf_files = [ vcf for vcf in os.listdir(out_dir) if vcf.endswith('_merged.vcf.gz') ] # not vcf.startswith('._') and
    merged_vcf_files = [ out_dir+'/'+vcf for vcf in merged_vcf_files ]

    #test_vcf_reader = read_vcf(vcf_files_merged_by_caller[0])


    #====== Melt #2 functions ======#
    #------ Read vcfs merged by caller into dictionary ------#
    #indel_callers = ['indelCalling', 'snowman', 'svcp']
    def merge_vcfs_across_callers(vcf_files_merged_by_caller):
        d_merged_records = {}

        for vcf in vcf_files_merged_by_caller:
            caller = [caller for caller in indel_callers if caller in vcf][0]

            vcf_reader = read_vcf(vcf)

            for record in vcf_reader:
                ## Get coordinates
                coordinates = [record.CHROM, record.POS, record.REF, str(record.ALT[0])]
                coord_key = list_as_string(coordinates)

                if debug == 1:
                    print 'Populating dictionary of variants (melt #2) | Position: %s' %coord_key

                ## Add coordinates to dictionary if not exist
                if coord_key not in d_merged_records:
                    # FILTER[+], INFO[FileID], TUMOR[GT,AF]
                    d_merged_records[coord_key] = ['', [], []]

                ## 0. FILTER: [+]; PASS is not reported by pyvcf
                if len(record.FILTER) == 0:
                    d_merged_records[coord_key][0] += 'PASS'
                else:
                    d_merged_records[coord_key][0] += list_as_string(record.FILTER, ';')

                ## 1. INFO: [FileID]
                INFO = 'FILEID_' + caller + '=' + list_as_string(record.INFO['FILEID'])
                d_merged_records[coord_key][1].append(INFO)
                #d_merged_records[coord_key][1].append(record.INFO['FILEID'])

                ## 2. TUMOR: [GT, AF]
                genotype = record.samples[0]['GT']
                AF = record.samples[0]['AF']
                d_merged_records[coord_key][2].append([genotype, AF])

        return d_merged_records


    #------ Write vcf lines; merge replicate vcfs by caller ------#
    def dict_record_as_line2(k, v):
        l_coords = k.split(',')
        CHROM = l_coords[0]
        POS = l_coords[1]
        REF = l_coords[2]
        ALT = l_coords[3]

        ## FILTER: PASS if 'PASS' in FILTER for at least 1 vcf
        if 'PASS' in v[0]:
            FILTER = 'PASS'
        else:
            FILTER = 'FAIL'

        ## INFO:
        FILEID_string = list_as_string(v[1],';')
        NCALLERS = sum([ caller in FILEID_string for caller in indel_callers ])
        INFO = 'NCALLERS=' + str(NCALLERS) + ';' + FILEID_string

        ## TUMOR:
        GT_not_null = []
        AF_l = []
        regex = re.compile(r'\d/\d')
        for l in v[2]:
            # GT
            GT_match = regex.search(l[0])
            if GT_match:
                GT_not_null.append(GT_match.group(0))

            # AF
            AF_l.append(l[1])

        if len(GT_not_null) != 0:
            GT = Counter(GT_not_null).most_common()[0][0]
        else:
            GT = './.'

        AF = round( np.mean(AF_l), 4)

        TUMOR = GT + ':' + str(AF)

        ID = '.'
        QUAL = '.'
        FORMAT = 'GT:AF'

        return '\t'.join([CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, TUMOR])

    #====== Execute merge and write #2 ======#
    dict_merged_vcfs_across_callers = merge_vcfs_across_callers(merged_vcf_files)

    # DO_re = re.compile(r'DO\d+')
    # donor_id = DO_re.search(os.path.basename(manifest_path)).group(0)

    melt2_vcf_path = '%s/%s%s' %(out_dir, donor_id, '_indel.vcf.gz')

    with gzip.open(melt2_vcf_path, 'wb', 9) as melt2_vcf:
        ## Write header
        melt2_vcf.write(vcf_header_final)

        ## Write body
        # Progress counter
        total_variants = len(dict_merged_vcfs_across_callers)
        progress_counter = 0

        if debug == 1:
            print 'Sorting dictionary (melt #2)...'

        for k,v in natsorted(dict_merged_vcfs_across_callers.iteritems()):
            progress_counter = progress_counter + 1
            if debug == 1:
                print 'Writing lines (melt #2) | Progress: %s/%s | Position: %s' % (progress_counter, total_variants, k)

            line = dict_record_as_line2(k, v) + '\n'
            melt2_vcf.write(line)

    # ====== Write 'done' file if script completes ======#
    if debug == 1:
        print 'Writing done_melt2.txt'
    done_melt2_txt = open('%s/%s' % (out_dir, 'done_melt2.txt'), 'w')
    done_melt2_txt.close()

    # #====== Test read vcf ======#
    # test_vcf = '/Users/lnguyen//hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/merged_vcfs/indels/DO51959/DO51959_indel.vcf.gz'
    #
    # test_vcf_reader = pyvcf.Reader(open(test_vcf, 'r'))
    # #test_record = next(test_vcf_reader)
    #
    # for i in test_vcf_reader:
    #     print i

