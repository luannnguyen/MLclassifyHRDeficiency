#!/usr/bin/env python

debug = 1

#Test donor: DO51959

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
manifest_path = args.manifest_path
donor_tsv_path = args.donor_tsv_path

#manifest_path = '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/downloads_subset/snv2/DO51959_snvVcfManifest.tsv'
#donor_tsv_path= '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/metadata/donor_pcawg_london.tsv'

vcf_list = []
with open(manifest_path) as manifest:
    vcf_file = manifest.readlines()
    vcf_list.extend(vcf_file)

vcf_files = [vcf.strip('\n') for vcf in vcf_list]
#vcf_files = [os.path.basename(vcf) for vcf in vcf_files]

out_dir = os.path.dirname(manifest_path) # output dir is the dirname of the manifest file

#------ Dictionary of vcfs split by caller ------#
snv_callers = ['mutect','snvCalling','svcp']
vcfs_by_caller = dict.fromkeys(snv_callers)
for caller in snv_callers:
    indices = [c for c, v in enumerate(vcf_files) if caller in v]
    vcfs_by_caller[caller] = [vcf_files[c] for c in indices]

#------ Get donor id from manifest ------#
DO_re = re.compile(r'DO\d+')
donor_id = DO_re.search(os.path.basename(manifest_path)).group(0)


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

vcf_header_mutect ='''##fileformat=VCFv4.1,
##caller=mutect,
##INFO=<ID=FILEID,Number=.,Type=String,Description="File IDs of contributing vcfs">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AF,Number=1,Type=Float,Description="Mean tumor ALT allele frequency across replicate vcfs">
##SAMPLE=<ID=TUMOR,SampleName=tumor_NA,Individual=NA,Description="Tumor">"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TUMOR
'''

vcf_header_snvCalling ='''##fileformat=VCFv4.1,
##caller=snvCalling,
##FILTER=<ID=RE,Description="variant in UCSC_27Sept2013_RepeatMasker.bed.gz region and/or SimpleTandemRepeats_chr.bed.gz region, downloaded from UCSC genome browser and/or variant in segmental duplication region, annotated by annovar">
##FILTER=<ID=BL,Description="variant in DAC-Blacklist from ENCODE or in DUKE_EXCLUDED list, both downloaded from UCSC genome browser">
##FILTER=<ID=DP,Description="<= 5 reads total at position in tumor">
##FILTER=<ID=SB,Description="Strand bias of reads with mutant allele = zero reads on one strand">
##FILTER=<ID=TAC,Description="less than 6 reads in Tumor at position">
##FILTER=<ID=dbSNP,Description="variant in dbSNP135">
##FILTER=<ID=DB,Description="variant in 1000Genomes, ALL.wgs.phase1_integrated_calls.20101123.snps_chr.vcf.gz or dbSNP">
##FILTER=<ID=HSDEPTH,Description="variant in HiSeqDepthTop10Pct_chr.bed.gz region, downloaded from UCSC genome browser">
##FILTER=<ID=MAP,Description="variant overlaps a region from wgEncodeCrgMapabilityAlign100mer.bedGraph.gz:::--breakPointMode --aEndOffset=1 with a value below 0.5, punishment increases with a decreasing mapability">
##FILTER=<ID=SBAF,Description="Strand bias of reads with mutant allele = zero reads on one strand and variant allele frequency below 0.1">
##FILTER=<ID=FRQ,Description="variant allele frequency below 0.05">
##FILTER=<ID=TAR,Description="Only one alternative read in Tumor at position">
##FILTER=<ID=UNCLEAR,Description="Classification is unclear">
##FILTER=<ID=DPHIGH,Description="Too many reads mapped in control at this region">
##FILTER=<ID=DPLOWC,Description="Only 5 or less reads in control">
##FILTER=<ID=1PS,Description="Only two alternative reads, one on each strand">
##FILTER=<ID=ALTC,Description="Alternative reads in control">
##FILTER=<ID=ALTCFR,Description="Alternative reads in control and tumor allele frequency below 0.3">
##FILTER=<ID=FRC,Description="Variant allele frequency below 0.3 in germline call">
##FILTER=<ID=YALT,Description="Variant on Y chromosome with low allele frequency">
##FILTER=<ID=VAF,Description="Variant allele frequency in tumor < 5 times allele frequency in control">
##FILTER=<ID=BI,Description="Bias towards a PCR strand or sequencing strand">
##INFO=<ID=FILEID,Number=.,Type=String,Description="File IDs of contributing vcfs">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AF,Number=1,Type=Float,Description="Mean tumor ALT allele frequency across replicate vcfs">
##SAMPLE=<ID=TUMOR,SampleName=tumor_NA,Individual=NA,Description="Tumor">"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TUMOR
'''

vcf_header_svcp ='''##fileformat=VCFv4.1,
##caller=svcp,
##FILTER=<ID=DTH,Description="Less than 1/3 mutant alleles were >= 25 base quality">
##FILTER=<ID=RP,Description="Coverage was less than 8 and no mutant alleles were found in the first 2/3 of a read (shifted 0.08 from the start and extended 0.08 more than 2/3 of the read length)">
##FILTER=<ID=MN,Description="More than 0.05 of mutant alleles that were >= 15 base quality found in the matched normal">
##FILTER=<ID=PT,Description="Mutant alleles all on one direction of read (1rd allowed on opposite strand) and in second half of the read. Second half of read contains the motif GGC[AT]G in sequenced orientation and the mean base quality of all bases after the motif was less than 20">
##FILTER=<ID=MQ,Description="Mean mapping quality of the mutant allele reads was < 21">
##FILTER=<ID=SR,Description="Position falls within a simple repeat using the supplied bed file">
##FILTER=<ID=CR,Description="Position falls within a centromeric repeat using the supplied bed file">
##FILTER=<ID=PH,Description="Mutant reads were on one strand (permitted proportion on other strand: 0.04), and mean mutant base quality was less than 21">
##FILTER=<ID=HSD,Description="Position falls within a high sequencing depth region using the supplied bed file">
##FILTER=<ID=GI,Description="Position falls within a germline indel using the supplied bed file">
##FILTER=<ID=VUM,Description="Position has >= 3 mutant allele present in at least 1 percent unmatched normal samples in the unmatched VCF.">
##FILTER=<ID=SE,Description="Coverage is >= 10 on each strand but mutant allele is only present on one strand">
##FILTER=<ID=MNP,Description="Tumour sample mutant allele proportion - normal sample mutant allele proportion < 0.2">
##INFO=<ID=FILEID,Number=.,Type=String,Description="File IDs of contributing vcfs">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AF,Number=1,Type=Float,Description="Mean tumor ALT allele frequency across replicate vcfs">
##SAMPLE=<ID=TUMOR,SampleName=tumor_NA,Individual=NA,Description="Tumor">"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TUMOR
'''

vcf_header_final = '''##fileformat=VCFv4.1
##INFO=<ID=NCALLERS,Number=1,Type=Integer,Description="Number of callers reporting variant">
##INFO=<ID=FILEID_mutect,Number=.,Type=String,Description="File IDs of contributing vcfs produced by mutect">
##INFO=<ID=FILEID_snvCalling,Number=.,Type=String,Description="File IDs of contributing vcfs produced by snvCalling">
##INFO=<ID=FILEID_svcp,Number=.,Type=String,Description="File IDs of contributing vcfs produced by svcp">
##FILTER=<ID=FAIL,Description="PASS not present in FILTER of any vcfs">
##FILTER=<ID=PASS1,Description="PASS present in FILTER of 1 vcf">
##FILTER=<ID=PASS2,Description="PASS present in FILTER of 2 vcfs">
##FILTER=<ID=PASS3,Description="PASS present in FILTER of 3 vcfs">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AF,Number=1,Type=Float,Description="Tumor ALT allele frequency">
##SAMPLE=<ID=TUMOR,SampleName=tumor_NA,Individual=NA,Description="Tumor">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TUMOR
'''

#%%%%%%%%%%%%%%%%%%%%%%%% Melt #1 %%%%%%%%%%%%%%%%%%%%%%%%#
if not args.skipmelt1:

    # ====== Melt #1 functions ======#
    #------ Get relevant vcf fields by caller as dictionary ------#
    def merge_vcfs_by_caller(dict_vcfs_by_caller, caller):
        # Dictionary structure:
        # coord_key: CHROM, POS, REF, ALT
        # FILTER: [+], INFO: [FileID], TUMOR:[GT,AF]
        d_merged_records = {}

        # dict_vcfs_by_caller = vcfs_by_caller
        # caller = 'mutect'
        for vcf in dict_vcfs_by_caller[caller]:
            # vcf = vcfs_by_caller[caller][0]

            vcf_reader = read_vcf(vcf)
            #print record.INFO

            file_id = get_vcf_info(vcf, 'File.ID')

            for record in vcf_reader:
                # ...... Get coordinates ......#
                coordinates = [record.CHROM, record.POS, record.REF, str(record.ALT[0])]
                coord_key = list_as_string(coordinates)

                if debug == 1:
                    print 'Populating dictionary of variants (melt #1) | File ID: %s | Position: %s' % (file_id, coord_key)

                # ...... Common metadata ......#
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
                if caller == 'mutect':

                    ## 2. TUMOR: [GT, AF]
                    genotype = './.'
                    AF = record.INFO['tumor_f'][0]

                    d_merged_records[coord_key][2].append([genotype, AF])

                elif caller == 'snvCalling':
                    ## 2. TUMOR: [GT, AF]
                    genotype = record.samples[1]['GT']
                    AF = float(record.INFO['AF'][1])

                    d_merged_records[coord_key][2].append([genotype, AF])

                elif caller == 'svcp':
                    ## 2. TUMOR: [GT, AF]
                    genotype = record.samples[1]['GT']
                    genotype = genotype.replace('|','/')

                    AF = float(record.samples[1]['PM'])

                    d_merged_records[coord_key][2].append([genotype, AF])

        ## Output preliminary dictionary
        return d_merged_records

    # dict_merged_records = merge_vcfs_by_caller(vcfs_by_caller,'mutect')
    #
    # for k,v in dict_merged_records.iteritems():
    #     # if len(v[1])>1:
    #     #     print c,v
    #
    #     l_coords = k.split(',')
    #     CHROM = l_coords[0]
    #     POS = l_coords[1]
    #     REF = l_coords[2]
    #     ALT = l_coords[3]
    #
    #     ## metadata
    #     ID = '.'
    #     QUAL = '.'
    #     FILTER = list_as_string(v[0], delim=';') #FIlTER
    #     INFO = 'FILEID=' + list_as_string(v[1]) # Convert FileID lists to strings
    #     FORMAT = 'GT:AF'
    #
    #     if caller == 'mutect':
    #         genotype_most_freq = './.'
    #     else:
    #         genotype_most_freq = Counter([l[0] for l in v[2]]).most_common(1)[0][0]
    #
    #     #print [float(l[1]) for l in v[2]]
    #
    #     AF_mean = round(np.mean([float(l[1]) for l in v[2]]),4)
    #
    #     TUMOR = genotype_most_freq + ':' + str(AF_mean)
    #
    #     ## print lines
    #     print '\t'.join([CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,TUMOR])



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

        if caller == 'mutect':
            genotype_most_freq = './.'
        else:
            genotype_most_freq = Counter([l[0] for l in v[2]]).most_common(1)[0][0]

        AF_mean = round(np.mean([float(l[1]) for l in v[2]]),4)

        TUMOR = genotype_most_freq + ':' + str(AF_mean)

        ## print lines
        return '\t'.join([CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,TUMOR])

    # for k,v in d_merged_records.iteritems():
    #     print dict_record_as_line(k,v,'svcp')

    #====== Execute merge and write #1 ======#
    #snv_callers = ['svcp'] #debugging
    for caller in snv_callers:

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


#%%%%%%%%%%%%%%%%%%%%%%%% Melt #2 %%%%%%%%%%%%%%%%%%%%%%%%#
if not args.skipmelt2:

    #====== Load merged vcfs (by caller) ======#
    #manifest_path = '/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/merged_vcfs/snv/DO1000/DO1000_indelVcfManifest.tsv'
    #out_dir = os.path.dirname(manifest_path)

    merged_vcf_files = [ vcf for vcf in os.listdir(out_dir) if vcf.endswith('_merged.vcf.gz') ] # not vcf.startswith('._') and
    #merged_vcf_files = [ vcf for vcf in os.listdir(out_dir) if vcf.endswith('_merged.vcf') ] # not vcf.startswith('._') and
    merged_vcf_files = [ out_dir+'/'+vcf for vcf in merged_vcf_files ]

    #test_vcf_reader = read_vcf(vcf_files_merged_by_caller[0])


    #====== Melt #2 functions ======#
    #------ Read vcfs merged by caller into dictionary ------#
    #snv_callers = ['mutect','snvCalling','svcp']
    def merge_vcfs_across_callers(vcf_files_merged_by_caller):
        d_merged_records = {}

        for vcf in vcf_files_merged_by_caller:
            caller = [caller for caller in snv_callers if caller in vcf][0]

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
                    d_merged_records[coord_key] = [[], [], []]

                ## 0. FILTER: [+]; PASS is not reported by pyvcf
                # if len(record.FILTER) == 0:
                #     d_merged_records[coord_key][0] += 'PASS'
                # else:
                #     d_merged_records[coord_key][0] += list_as_string(record.FILTER, ';')

                if len(record.FILTER) == 0:
                    d_merged_records[coord_key][0].append('PASS')
                else:
                    d_merged_records[coord_key][0].append(list_as_string(record.FILTER))

                ## 1. INFO: [FileID]
                INFO = 'FILEID_' + caller + '=' + list_as_string(record.INFO['FILEID'])
                d_merged_records[coord_key][1].append(INFO)
                #d_merged_records[coord_key][1].append(record.INFO['FILEID'])

                ## 2. TUMOR: [GT, AF]
                genotype = record.samples[0]['GT']
                AF = record.samples[0]['AF']
                d_merged_records[coord_key][2].append([genotype, AF])

                #print  record.FILTER

        return d_merged_records

    # merge_vcfs_across_callers(merged_vcf_files)
    # test_d_merged_records = merge_vcfs_across_callers(merged_vcf_files)
    #


    #------ Write vcf lines; merge replicate vcfs by caller ------#
    def dict_record_as_line2(k, v):
        l_coords = k.split(',')
        CHROM = l_coords[0]
        POS = l_coords[1]
        REF = l_coords[2]
        ALT = l_coords[3]

        ## FILTER: PASS if 'PASS' in FILTER for at least 1 vcf

        # if 'PASS' in v[0]:
        #     FILTER = 'PASS'
        # else:
        #     FILTER = 'FAIL'

        pass_count = v[0].count('PASS')
        if pass_count > 0:
            FILTER = 'PASS' + str(pass_count)
        else:
            FILTER = 'FAIL'

        ## INFO:
        FILEID_string = list_as_string(v[1],';')
        NCALLERS = sum([ caller in FILEID_string for caller in snv_callers ])
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

    # for k,v in test_d_merged_records.iteritems():
    #     line = dict_record_as_line2(k, v) + '\n'
    #     print line
    #     pass_count = v[0].count('PASS')
    #     if pass_count > 0:
    #         FILTER = 'PASS' + pass_count
    #     else:
    #         FILTER = 'FAIL'

    #====== Execute merge and write #2 ======#
    dict_merged_vcfs_across_callers = merge_vcfs_across_callers(merged_vcf_files)

    # DO_re = re.compile(r'DO\d+')
    # donor_id = DO_re.search(os.path.basename(manifest_path)).group(0)

    melt2_vcf_path = '%s/%s%s' %(out_dir, donor_id, '_snv.vcf.gz')

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
            #print line
            melt2_vcf.write(line)

    # ====== Write 'done' file if script completes ======#
    if debug == 1:
        print 'Writing done_melt2.txt'
    done_melt2_txt = open('%s/%s' % (out_dir, 'done_melt2.txt'), 'w')
    done_melt2_txt.close()

