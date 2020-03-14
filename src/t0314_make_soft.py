# -*- coding: utf-8 -*-
from path import Path
import pymisca.ext as pyext
def WORKDIR():return Path('./rawfile.workdir').makedirs_p()
validate_fastq_py = ('/home/feng/envs/0726-polyq/src/validate_fastq.py')
def _readData(fn,**kw):
    return pyext.readData(Path(__file__).dirname()/fn,encoding='utf8',**kw)

mcurr = _readData(
        '/home/feng/envs/upGeo/results/0424-database-raw/mcurr.csv',
        guess_index=0)
#m#curr['
mcurr['FULL_PATH'] = mcurr['FULL_PATH'].astype(str)
mcurr = mcurr.loc[~mcurr['FULL_PATH'].str.contains("Raw_data/184R_Q_reseq181_combined")]
_rawMeta = mcurr
def rawMeta(): return _rawMeta
def DEBUG():return _DEBUG ##runtime

if 1:

    def df_mappedData_full():
        df = _readData('/home/feng/envs/upGeo/results/0407-database-mapped/mcurr.csv')
        dfc = df
        #     dfc = df.loc[df['DATAACC'].isin(_accs)]
        dfc = dfc.loc[dfc['SIZE']>2.]

        blacklist = '''
        ### 184RS19 mistaken
        /home/feng/writable/teamkj/__backup/mapped-data/RNA-seq/Mapped_data/184R_Q/S19/star_out/181R-Q-851-elf3-ZT8-22C_S3.fastq.gz/
        ### 184RS22 mistaken
        /home/feng/writable/teamkj/__backup/mapped-data/RNA-seq/Mapped_data/184R_Q/S22/star_out/181R-Q-854-4195-ZT8-27C_S6.fastq.gz/
        ### 184RS27 mistaken 
        /home/feng/writable/teamkj/__backup/mapped-data/RNA-seq/Mapped_data/184R_Q/S27/star_out/181R-Q-859-4195-ZT12-22C_S11.fastq.gz/
        ### 176C duplication
        /home/feng/writable/teamkj/__backup/mapped-data/ChIP-seq/Mapped_data/176C/176C/
        #### name_sorted 
        /home/feng/writable/teamkj/__backup/mapped-data/RNA-seq/Mapped_data/139R_Noemie/176C_BURLAT_220316_rep1_N-40660957/tophat_results/
        /home/feng/writable/teamkj/__backup/mapped-data/RNA-seq/Mapped_data/139R_Noemie/176C_BURLAT_220316_rep1_N-40660957/tophat_results/176C_BURLAT_220316_rep1_S11_trimmo_ensembl_nomixed_unstranded/176C_BURLAT_220316_rep1_S11_trimmo_paired_2_10_5_1_tophat_ensembl_TAIR10_nomixed_unstranded_sorted_rmdup_picard_name_sorted.bam
        /home/feng/writable/teamkj/__backup/mapped-data/ChIP-seq/Mapped_data/182C/S24/35SELF3-27C-HA_S24_peaks.narrowPeak
        /home/feng/writable/teamkj/__backup/mapped-data/ChIP-seq/Mapped_data/182C/S24/35SELF3-27C-HA_S24_RPKM.bw
        Mapped_data/pipeline_rnaseq/199R
    '''.replace(' ','').strip().splitlines()
        
        for k in blacklist:
            dfc = dfc.loc[~dfc["FULL_PATH"].str.contains(k)]
        dfc = dfc.loc[~dfc['FILEACC'].str.endswith("orig.bam")]
        dfc = dfc.loc[~dfc['FILEACC'].str.endswith("unmapped.bam")]
        dfc = dfc.loc[~dfc['FILEACC'].str.endswith("name_sorted.bam")] #### name_sorted.bam is for htseq-count
        dfc = dfc.loc[~(~dfc['FILEACC'].str.endswith("RPKM.bw") & dfc['FILEACC'].str.endswith(".bw"))]
        # dfc = dfc.loc[~(~dfc['FILEACC'].str.endswith("RPKM.bw") & dfc['FILEACC'].str.endswith(".bw"))]
        # dfc.to_csv('df_mappedData_full.csv')
        return dfc

    import sys
    def df_mappedData_rnaseq():
        _accs = DATA_ACC_RNASEQ()
        dfc = df_mappedData_full()
        dfc = dfc.loc[dfc['DATAACC'].isin(_accs)]
        #     assert (dfc.groupby("EXT").apply(len)==len(set(_accs))).all()
        #dfc = dfc.reset_index(drop=0).pivot_table(index='DATAACC',columns='EXT',values='FULL_PATH',aggfunc=lambda x:x)    
        dfc = dfc.reset_index(drop=0).pivot_table(index='DATAACC',columns='EXT',values='FULL_PATH',
    aggfunc=lambda x:[sys.stdout.write(x.to_csv()) if DEBUG() else [],x][1])    
        return dfc

    def DATA_ACC_RNASEQ():
        if not DEBUG():
            s ='143RS1  143RS2  143RS3  143RS4  143RS5  143RS6  144RS10 144RS11 144RS12 144RS13 144RS14 144RS15 144RS16 144RS17 144RS18 144RS19 144RS20 144RS21 148RS1  148RS2  148RS3  148RS4  148RS5  148RS6  148RS7  148RS8  149RS1  149RS10 149RS11 149RS12 149RS13 149RS16 149RS17 149RS18 149RS19 149RS2  149RS20 149RS21 149RS22 149RS23 149RS24 149RS25 149RS26 149RS27 149RS28 149RS29 149RS3  149RS30 149RS4  149RS5  149RS6  149RS7  149RS8  149RS9  150RS1  150RS2  150RS3  150RS4  150RS5  150RS6  150RS7  150RS8  169RS17 169RS18 169RS19 169RS20 169RS21 169RS22 169RS23 169RS24 169RS25 169RS26 169RS27 169RS28 169RS29 169RS30 169RS31 169RS32 169RS33 169RS34 193RS1  193RS10 193RS11 193RS12 193RS13 193RS14 193RS2  193RS3  193RS4  193RS5  193RS6  193RS7  193RS8  193RS9  196RS1  196RS10 196RS11 196RS12 196RS13 196RS14 196RS15 196RS16 196RS2  196RS3  196RS4  196RS5  196RS6  196RS7  196RS8  196RS9  199RS1  199RS10 199RS11 199RS12 199RS13 199RS14 199RS15 199RS16 199RS2  199RS3  199RS4  199RS5  199RS6  199RS7  199RS8  199RS9  201RS1  201RS10 201RS11 201RS12 201RS13 201RS14 201RS15 201RS16 201RS17 201RS18 201RS19 201RS2  201RS20 201RS21 201RS22 201RS23 201RS3  201RS4  201RS5  201RS6  201RS7  201RS8  201RS9'
            s = s.split()
        else:
            s = ['143RS1','201RS9']
        return s


# df= (df_mappedData_rnaseq())
# import pdb;pdb.set_trace();

def template_common():
    template = u'''
#### filepath: {{sample.fname}}

^SAMPLE =  {{sample.data_acc}}

!Sample_molecule = {{sample.library_molecule}}
!Sample_library_strategy = {{sample.library_strategy}}
!Sample_instrument_model = NextSeq 500

!Sample_title = {{ sample.data_acc}}
!Sample_type = SRA
!Sample_source_name = {{sample.tissue}}
!Sample_organism = {{sample.species}}

###[CHECK-GENOTYPE-and-ANITBODY!]
!Sample_characteristics = genotype: {{sample.genotype}}
!Sample_characteristics = cultivar: {{sample.cultivar}}
!Sample_characteristics = tissue: {{sample.tissue}}
!Sample_characteristics = photoperiod: {{sample.photoperiod}}
!Sample_characteristics = temperature: {{sample.temperature}}
!Sample_characteristics = replicate: {{sample.replicate}}
!Sample_genome_build = {{sample.genome}}


### growth protocol
{{sample.growth_protocol}}
###[CHECK ANTIBODY!!!]
{{sample.extraction_protocol}}

###[KEEP-THIS-SECTION!]
{{sample.processing_protocol}}
{{sample.rawfile_detail}}
'''.strip()
    # template_commone = '123'
    # template = '\n'.join([x.strip() for x in template.splitlines()])
    return template

from numpy import nan

def processing_rnaseq_brachy():
    s = '''
!Sample_data_processing = Adapters were trimmed off from raw reads with Trimmomatic with argument "ILLUMINACLIP:$FA_ADAPTER:6:30:10 LEADING:3 TRAILING:3 MINLEN:36 SLIDINGWINDOW:4:15". 
!Sample_data_processing = Raw reads were aligned with Hisat2 with arguments "--no-mixed --rna-strandness RF --dta --fr" to produce a SAM file.
!Sample_data_processing = Duplicate reads were removed with Picard using default
setting
!Sample_data_processing = Alignments in SAM file were assembled into transcripts abundances with stringtie with argument "--rf".


!Sample_data_processing = Supplementary_files_format_and_content: *.bam: HISAT2-aligned and picard deduplicated genomic alignment .
!Sample_data_processing = Supplementary_files_format_and_content: .stringtie.count: TSV table containing abundance of transcripts with bed-formatted coordinates. 
!Sample_data_processing = Supplementary_files_format_and_content: _htseq_count.ct: TSV table containing htseq-counted transcripts.
!Sample_supplementary_file_1 = {{sample.file_bam}}
!Sample_supplementary_file_2 = {{sample.file_count}}
!Sample_supplementary_file_3 = {{sample.file_ct}}

'''
    return s

def meta_df_rnaseq():
    # mcurr0 = _readData('/home/feng/static/results/0318-makeRNA-polyQ/mcurr0.csv',guess_index=0)
    mcurr0 = _readData('0630-submit-brachy-rnaseq.csv',guess_index=0)
    mcurr0.columns = mcurr0.columns.str.lower()
    mcurr0 = mcurr0.rename(columns={
        'dataacc':'data_acc',
        'gtype':'genotype',
        'temp':'temperature',

    })
    const_dict = {
        "tissue":"fully expanded leaves",
        'library_strategy':'RNA-Seq',
        'library_molecule':'total RNA',
        'species':"txid15368", ### arabidopsis
        "cultivar":"Bd21-3",
        "genome":'Bd21-v3.1',
        "replicate":"single",
        "fname":nan,
        # "light":nan,
        # "fragmentation_method":nan,
        # "chip_antibody":nan,
        'growth_protocol':
        u'!Sample_growth_protocol = Seeds were imbibed in distilled water at 4 ºC '
        u'for two days before sowing. Plants were grown in 5 parts John Innes #2,'
        u' 3 parts peat, 1 parts silver sand, 3 parts course vermiculite,'
        u' Osmocote 2.7 g/L. All plants were grown in growth cabinets in '
        u'{{sample.photoperiod}} with constant temperature {{sample.temperature}},'
        u' 65 % humidity '
        u'and 350 μmol m-2 s-1 PPFD (Photosynthetic Photon Flux Density). Samples '
        u'were collected after {{sample.age}} at {{sample.ztime}} by snap-freezing in liquid nitrogen.',

        "extraction_protocol":
        u"!Sample_extract_protocol = "
        u"Qiagen RNeasy Mini Kit (74104) was used to extract RNA. "
        u"RNA quality and integrity were assessed on the Agilent 2200 TapeStation system. "
        u"\n"
        u'!Sample_library_construction_protocol = ' 
        u"Library preparation was performed with 1µg total RNA using the NEBNext® Ultra™ "
        u"Directional RNA Library Prep Kit for Illumina® (E7420L). The libraries were "
        u"sequenced on a NextSeq500 (Illumina) running a final pooled library. Each pool "
        u"contained 24 to 30 samples and was sequenced using NextSeq® 500/550 High Output"
        u" Kit v2 (150 cycles) TG-160-2002 on a NextSeq500 (Illumina). "
        ,

        # "processing_protocol": ',
        "processing_protocol":"{{sample.processing_protocol}}",
        "rawfile_detail":"{{sample.rawfile_detail}}"
    }

    for key,v in const_dict.items():
        mcurr0[key]= v

        
    # print mcurr0.iloc[:1].to_csv(sep='\t',index=0)
    # if DEBUG():
    #     dfc = mcurr0.iloc[:1]
    # else:
    dfc = mcurr0.loc[mcurr0["data_acc"].isin(DATA_ACC_RNASEQ())]
    assert dfc.data_acc.is_unique
    return dfc


def sample_rnaseq_processing_protocol(sample):
    from pymisca.events import CopyEvent,LinkEvent
    from pymisca.ext import f as _f
    import os
    
    OUTDIR = WORKDIR() / sample['data_acc']/ 'supp'
#     sample.data_acc_control = _get_data_acc_control(sample)
    rec = df_mappedData_rnaseq().loc[sample['data_acc']]

    for attrName,key in [
        ('file_bam','bam'),
        ('file_bw','bw'),
        ('file_count','count'),
        ('file_ct','ct'),
        # ('file_count','txt'),
#         ('file_npk','narrowPeak')
    ]:
        
        fullname= rec[key]
        sample[attrName+'_orig'] = fullname
        if pyext.pd.isnull(fullname) or (fullname ==[]):
            sample[attrName] = 'NA'
        else:
            basename = os.path.basename(fullname)
            sample[attrName] = LinkEvent(
                    CopyEvent(fullname,OUTDIR / basename).dest,
                    WORKDIR()/_f("ftp/{sample.data_acc}.supp.{basename}"),1
                ).dest.relpath(WORKDIR() / 'ftp')    
    if sample['file_ct'] == 'NA':
        assert sample['file_count'] != 'NA'
        template = u'''
    !Sample_data_processing = Adapters were trimmed off from raw reads with Trimmomatic with argument "ILLUMINACLIP:$FA_ADAPTER:6:30:10 LEADING:3 TRAILING:3 MINLEN:36 SLIDINGWINDOW:4:15". 
    !Sample_data_processing = Raw reads were aligned with Hisat2 with arguments "--no-mixed --rna-strandness RF --dta --fr" to produce a SAM file.
    !Sample_data_processing = Duplicate reads were removed with Picard using default
    setting
    !Sample_data_processing = Alignments in SAM file were assembled into transcripts abundances with stringtie with argument "--rf".


    !Sample_data_processing = Supplementary_files_format_and_content: .bam: HISAT2-aligned and picard deduplicated genomic alignment .
    !Sample_data_processing = Supplementary_files_format_and_content: .stringtie.count: TSV table containing abundance of transcripts with bed-formatted coordinates. 
    !Sample_data_processing = Supplementary_files_format_and_content: RPKM.bw:  RPKM normalised bigwig files
    !Sample_supplementary_file_1 = {{sample.file_bam}}
    !Sample_supplementary_file_2 = {{sample.file_count}}
    !Sample_supplementary_file_3 = {{sample.file_bw}}
    '''
    else:
        assert sample['file_ct'] != 'NA'
        template = u'''
    !Sample_data_processing = Adapters were trimmed off using Trimmomatic with "ILLUMINACLIP:$FA_ADAPTER:2:10:5:1"
    !Sample_data_processing = The trimmed reads were aligned using Tophat with "--max-multihits --library-type fr-firststrand --no-mixed"
    !Sample_data_processing = Duplicate reads were removed with Picard using default setting
    !Sample_data_processing = Alignments in SAM file were assembled into transcripts abundances using htseq-count with "-r name -s no -f bam -t exon -i gene_id" against the GTF annotation"

    !Sample_data_processing = Supplementary_files_format_and_content: *.bam: Tophat-aligned and picard deduplicated genomic alignment .
    !Sample_data_processing = Supplementary_files_format_and_content: _htseq_count.ct: TSV table containing htseq-counted transcript abundances.
    !Sample_supplementary_file_1 = {{sample.file_bam}}
    !Sample_supplementary_file_2 = {{sample.file_ct}}
    '''

    res = pyext.jf2(template,)
    return res


def clean(s):
    out = []
    for x in s.splitlines():
        x = x.strip()
        if not x or x.startswith('#'):
            continue
        out.append(x)
    return '\n'.join(out)

def sample_get_rawfile_detail(sample,WORKDIR=WORKDIR, validate_fastq_py=validate_fastq_py):

    # return 
    from pymisca.events import LinkEvent
    from pymisca.ext import f as _f
    node = pyext.file__asModule(validate_fastq_py)
    node.rawMeta = rawMeta()
    node.DATA_ACC = sample['data_acc']
    node.WORKDIR = WORKDIR()
#     pyext.path.Path('/home/feng/envs/0726-polyq/WORKDIR.submit/').realpath()
#         node.WORKDIR = pyext.path.Path('/home/feng/envs/0830-polyq/WORKDIR/').realpath()
    node.valid_fastq()
    sample.rawfile_nodes = nodes = node.combined_valid_fastq()['OUTPUT_NODES']
    sample.rawfile_files_orig = [x['OUTPUT_FILE'] for x in nodes]
#     sample.rawfile_files = [x['OUTPUT_FILE'].relpath(WORKDIR()) for x in nodes]

    #### Relinking because GEO needs a flat directory tree
    sample.rawfile_files = [
        LinkEvent(
            x['OUTPUT_FILE'],
            WORKDIR()/"ftp"/_f('{sample.data_acc}.{x["OUTPUT_FILE"].basename()}'),
            1,).dest.relpath(WORKDIR()/"ftp") for x in nodes]
    
#     print nodes[0]._data.keys()
    sample.rawfile_checksums = [x['FILE_MD5']['MD5_HEX'] for x in nodes]
    sample.rawfile_readlengths = [ '75' for x in nodes]
    sample.rawfile_is_paired = 'paired-end' if len(nodes) > 1 else 'single'
    template = u'''
!Sample_raw_file_name = {{','.join(sample.rawfile_files)}}
!Sample_raw_file_type = fastq
!Sample_raw_file_checksum = {{','.join(sample.rawfile_checksums)}}
!Sample_raw_file_read_length = {{','.join(sample.rawfile_readlengths)}}
!Sample_raw_file_single_or_paired-end = {{sample.rawfile_is_paired}}
!Sample_raw_file_instrument_model = NextSeq 500
'''
    
    return pyext.jf2(template)



import io

from attrdict import AttrDict
from jinja2 import Template,StrictUndefined
import pymisca.ext as pyext
def main():
    temp = template_common()
    samples = meta_df_rnaseq()
    df = samples.fillna('NA')
    samples = pyext.df__iterdict( df)
    samples = list(samples)
    samples = [AttrDict(x) for x in samples]
    # Samples
    from pprint import pprint
#    pprint(samples)

    # print(samples[0])
    for sample in samples:
        sample['rawfile_detail'] = sample_get_rawfile_detail(sample)    
        sample['processing_protocol'] = sample_rnaseq_processing_protocol(sample)
        # processing_rnaseq_brachy(),
        # for sample in samples:
        pprint(sample)
        # pprint(samples)
        res = Template(temp, undefined=StrictUndefined).render(**locals())
        res = Template(res, undefined=StrictUndefined).render(**locals())
        # res = Template(res, undefined=StrictUndefined).render(**locals())
        # res =   pyext.jf2(temp)
        # break
        # res = pyext.jf2(res)
        res = clean(res)
        fn = (WORKDIR()+'.soft_files').makedirs_p()/sample['data_acc']+'.soft.txt'
        with io.open( fn,'w',encoding='utf8') as f:
            f.write(res)
            # print(res)
        # break
    # print(template_common())
    pass



if __name__ == '__main__':
    import sys
    if '--real' in sys.argv:
        _DEBUG = 0
    else:
        _DEBUG = 1
    main()
