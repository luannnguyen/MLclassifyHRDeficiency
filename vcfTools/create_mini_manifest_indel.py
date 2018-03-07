
import sys
import os
import re
import numpy as np
import pandas as pd

donor_tsv_path = '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/metadata/donor_pcawg_london.tsv'

df_donor = pd.read_table(donor_tsv_path, sep='\t', usecols=['Analysis.ID','File.Name','Variant.Type','Caller','ICGC.Donor'])
df_donor_indel = df_donor.loc[df_donor['Variant.Type'] == 'INDEL', ]

donor_ids = list(set(list(df_donor['ICGC.Donor'])))


merged_vcf_path_local = '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/merged_vcfs/indels'
downloads_vcf_path = '/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/downloads'
for donor in donor_ids:
    out_path = merged_vcf_path_local + '/' + donor

    # if not os.path.exists(out_path):
    #     os.makedirs(out_path)

    df = df_donor_indel.loc[ df_donor_indel['ICGC.Donor'] == donor ]
    analysis_ids = list(df['Analysis.ID'].values)
    file_names = list(df['File.Name'].values)

    vcf_paths = [downloads_vcf_path + '/' + i + '/' + j for i,j in zip(analysis_ids, file_names)]

    print 'Creating manifest for %s' %donor
    with open(out_path + '/' + donor + '_indelVcfManifest.tsv', 'w') as manifest_tsv:
        for i in vcf_paths:
            manifest_tsv.write(i + '\n')
