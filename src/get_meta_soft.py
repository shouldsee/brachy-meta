# -*- coding: utf-8 -*-
'''
https://github.com/ENCODE-DCC/ucscGb/blob/master/ucscGb/externalData/geo/soft.py
'''
import pymisca.ext as pyext
import attrdict

from pymisca import shell

from delay_and_cache import delay_and_cache as _dac
from delay_and_cache import StaticClass,_ctf,_fkey

def DEBUG(): return 0
def WORKDIR():return pyext.path.Path('$HOME/envs/0830-polyq/WORKDIR').expand().realpath()    

def WORKDIR():
    return pyext.path.Path('$PWD/get_meta_soft.work').expand().realpath()
    # HOME/envs/0726-polyq/WORKDIR.submit/').expand().realpath()
    return pyext.path.Path('$HOME/envs/0726-polyq/WORKDIR.submit/').expand().realpath()

#### Overwrite file-reading stream
from util import _get_file
_readData = lambda x,**kw: pyext.readData(_get_file(x),**kw)
# def DATA_ACC_RNASEQ():
#     lst =[
#         _readData("/home/feng/envs/0726-polyq/src/BOXPLOT_1_OX-lines.tsv").dropna().index.tolist(),
#         _readData("/home/feng/envs/0726-polyq/src/BOXPLOT2_mutants.tsv").dropna().index.tolist(),
#         _readData("/home/feng/envs/0726-polyq/src/BOXPLOT3_Qlines.tsv").dropna().index.tolist(),
#     ]
#     lst = pyext.stringList__flatten(lst)
#     lst = list(set(lst))
#     return lst

# def ACC2BNAME():
#     return


MCURR0 =  _readData('upGeo/0312-meta-copy/meta_rna.tsv',guess_index=0)
def meta_all():
    pd = pyext.pd
#     buf = pd.DataFrame([])
    buf = []
    mcurr0 = _readData('/home/feng/static/results/0318-makeRNA-polyQ/mcurr0.csv',guess_index=0)
    mcurr0 = mcurr0.drop(columns='RunID')
#     print(mcurr0['RunID']==mcurr0['runID']).mean()
    mcurr0.columns = mcurr0.columns.str.lower()
    mcurr0['fname'] = mcurr0['fname_']
    buf.append(mcurr0)
    
    mcurr0 = MCURR0   
    # mcurr0 = _readData('/home/feng/meta/meta_chip.tsv',guess_index=0)
    mcurr0.columns = mcurr0.columns.str.lower()
    buf.append(mcurr0)
    buf = pd.concat(buf,axis=0)
    return buf

def df_figureRegistry():
    pd = pyext.pd
    dfc = pyext.pd.DataFrame( BUFFER_CHIPSEQ() + BUFFER_RNASEQ(),
                             columns=['DATA_ACC','FIGURE_ACC','DATE_ADDED','POSITION_IN_GRAPH']
                            )
    
    ser = dfc.groupby('DATA_ACC').apply(lambda _df:':'.join(["[%s]"%x for x in _df['FIGURE_ACC']])).to_frame('ALL_FIGURE_ACC')
    # dfc = pyext.pd.concat([dfc,ser],axis=1)
    dfc = pd.merge(dfc,ser,left_on='DATA_ACC',right_index=True)    
    
    dfc['BASENAME'] = meta_all().set_index('dataacc')['bname'].reindex(dfc['DATA_ACC']).tolist()
    dfc = dfc.sort_values(['FIGURE_ACC','POSITION_IN_GRAPH'])
    return dfc

import itertools
def BUFFER_RNASEQ():
    lst =[
        ("/home/feng/envs/0726-polyq/src/BOXPLOT_1_OX-lines.tsv",'20190822'),
#         .dropna().index.tolist(),
        ("/home/feng/envs/0726-polyq/src/BOXPLOT2_mutants.tsv",'20190822'),
#         .dropna().index.tolist(),
        ("/home/feng/envs/0726-polyq/src/BOXPLOT3_Qlines.tsv",'20190822'),
        
        ("/home/feng/envs/0726-polyq/src/FIG-S5_HEATMAP2_genotypes.tsv", "20190906"),
        ("/home/feng/envs/0726-polyq/src/FIG-S6_HEATMAP1_Q-lines.tsv","20190906" ),
        ("/home/feng/envs/0726-polyq/src/FIG-S7_HEATMAP3.tsv", "20190906"),#         .dropna().index.tolist(),
    ]
    buf = []
    for FNAME,DATE in lst:
        figureName = pyext.os.path.basename(FNAME)
        df = pyext.readData(FNAME).dropna()
#         df['POSITION_IN_GRAPH'] = -1
        buf += list(zip(
                 df.index,
                 itertools.repeat(figureName),
                 itertools.repeat(DATE),
                 df.reset_index().index + 1,
#                 dfPOSITION_IN_GRAPH＇］,
        ))

#     lst = pyext.stringList__flatten(lst)
#     lst = list(set(lst))
    return buf

def DATA_ACC_RNASEQ():
    s ='143RS1  143RS2  143RS3  143RS4  143RS5  143RS6  144RS10 144RS11 144RS12 144RS13 144RS14 144RS15 144RS16 144RS17 144RS18 144RS19 144RS20 144RS21 148RS1  148RS2  148RS3  148RS4  148RS5  148RS6  148RS7  148RS8  149RS1  149RS10 149RS11 149RS12 149RS13 149RS16 149RS17 149RS18 149RS19 149RS2  149RS20 149RS21 149RS22 149RS23 149RS24 149RS25 149RS26 149RS27 149RS28 149RS29 149RS3  149RS30 149RS4  149RS5  149RS6  149RS7  149RS8  149RS9  150RS1  150RS2  150RS3  150RS4  150RS5  150RS6  150RS7  150RS8  169RS17 169RS18 169RS19 169RS20 169RS21 169RS22 169RS23 169RS24 169RS25 169RS26 169RS27 169RS28 169RS29 169RS30 169RS31 169RS32 169RS33 169RS34 193RS1  193RS10 193RS11 193RS12 193RS13 193RS14 193RS2  193RS3  193RS4  193RS5  193RS6  193RS7  193RS8  193RS9  196RS1  196RS10 196RS11 196RS12 196RS13 196RS14 196RS15 196RS16 196RS2  196RS3  196RS4  196RS5  196RS6  196RS7  196RS8  196RS9  199RS1  199RS10 199RS11 199RS12 199RS13 199RS14 199RS15 199RS16 199RS2  199RS3  199RS4  199RS5  199RS6  199RS7  199RS8  199RS9  201RS1  201RS10 201RS11 201RS12 201RS13 201RS14 201RS15 201RS16 201RS17 201RS18 201RS19 201RS2  201RS20 201RS21 201RS22 201RS23 201RS3  201RS4  201RS5  201RS6  201RS7  201RS8  201RS9'
    s = s.split()
    return s
    # return list(set(zip(*BUFFER_RNASEQ())[0]))


def DATA_ACC_CHIPSEQ():
    return list(set(zip(*BUFFER_CHIPSEQ())[0]))

def BUFFER_CHIPSEQ():
    '''
    See: 0726-polyq/src/make_chipseq_pileups.py
    '''
    buf = [
        
        ("192CS17","fig-2c","20190822"),
        ("192CS18","fig-2c","20190822"),
        
        ("192CS1","fig-2d","20190822"),
        ("192CS2","fig-2d","20190822"),
        ("192CS3","fig-2d","20190822"),
        ("192CS4","fig-2d","20190822"),
        
        ("189CS10","fig-2e","20190905"),
        ("189CS11","fig-2e","20190905"),
        ("189CS16","fig-2e","20190905"),
        ("189CS17","fig-2e","20190905"),
        
#         "192CS10",
#         "192CS11",
#         "192CS16",
#         "192CS17",
    ]
    
#     for FNAME in lst:
    # if 1:
    #     FNAME=  INPUTDIR() /'src/0726-figure-meta.tsv'
    #     figureName = pyext.os.path.basename(FNAME)
    #     figureName = 'figS4E_0905'
    #     df = pyext.readData(FNAME)[['figS4E_0905']].dropna()
    #     df.columns = ['POSITION_IN_GRAPH']
    #     buf += list(zip(
    #              df.index,
    #              itertools.repeat(figureName),
    #              itertools.repeat('20190905'),
    #              df['POSITION_IN_GRAPH']
    #     ))

    return buf

@_dac
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
    !Sample_characteristics = photoperiod: {{sample.light}}
    !Sample_characteristics = temperature: {{sample.temperature}}
    !Sample_characteristics = replicate: {{sample.replicate}}
    !Sample_characteristics = fragmentation method: {{sample.fragmentation_method}}
    !Sample_characteristics = chip antibody: {{sample.chip_antibody}}
    !Sample_genome_build = {{sample.genome}}
    
    !Sample_growth_protocol = {{sample.growth_protocol}}
    ###[CHECK ANTIBODY!!!]
    {{sample.extraction_protocol}}
    
    ###[KEEP-THIS-SECTION!]
    {{sample.processing_protocol}}

    {{sample.rawfile_detail}}
    '''
    template = '123'
    # template = '\n'.join([x.strip() for x in template.splitlines()])
    return template


from numpy import nan
def meta_df_rnaseq():
    # mcurr0 = _readData('/home/feng/static/results/0318-makeRNA-polyQ/mcurr0.csv',guess_index=0)
    mcurr0 = _readData('./0630-submit-brachy-rnaseq.csv',guess_index=0)
    mcurr0.columns = mcurr0.columns.str.lower()
    mcurr0 = mcurr0.rename(columns={
        'dataacc':'data_acc',
        'gtype':'genotype',
        'temp':'temperature',

    })
    const_dict = {
        'tissue':'seedling',
        'library_strategy':'RNA-Seq',
        'library_molecule':'total RNA',
        'species':"txid3702", ### arabidopsis
        "cultivar":"Col-0",
        "genome":"TAIR10",
        "replicate":"single",
        "fname":nan,
        "light":nan,
        "fragmentation_method":nan,
        "chip_antibody":nan,
        "growth_protocol":u"Seeds were sterilised and sown on ½ X Murashige and Skoog-agar (MS-agar) plates at pH 5.7. Seedlings were grown at 22˚C and 170 µmol/m2/s light with a 16-hour long-day photoperiod for 8-10 days and then treated at ZT1 or ZT2.5 with 30 µM cycloheximide in deionised water at 22˚C or mock controls (0.1% v/v DMSO) at 22˚C or 4˚C.  ",
        "extraction_protocol": u'''
    !Sample_extract_protocol = Samples were harvested by flash-freezing in liquid nitrogen. Total RNA was extracted from 30-40 mg of ground seedling tissue using the MagMAX™-96 Total RNA Isolation Kit (Thermo Fisher Scientific)  according to the manufacturer's instructions. RNA quality and integrity was assessed on the Agilent 2200 TapeStation.
    !Sample_library_construction_protocol = Library preparation was performed with 1 ug of high integrity total RNA using the Next Ultra Directional RNA library preparation kit (NEB) according to the manufacturer's instructions. The libraries were sequenced on a NextSeq 500 (Illumina) using paired-end sequencing with a NextSeq 500/550 High Output v2 kit (150 cycles).
    ''',
        "processing_protocol":"{{sample.processing_protocol}}",
        "rawfile_detail":"{{sample.rawfile_detail}}"
    }
    for key,v in const_dict.items():
        mcurr0[key]= v
        
    # print mcurr0.iloc[:1].to_csv(sep='\t',index=0)
    if DEBUG():
        dfc = mcurr0.iloc[:1]
    else:
        dfc = mcurr0.loc[mcurr0["data_acc"].isin(DATA_ACC_RNASEQ())]
    assert dfc.data_acc.is_unique
    return dfc

def meta_df_chipseq():
    mcurr0 = MCURR0
    # mcurr0 = _readD   ata('/home/feng/meta/meta_chip.tsv',guess_index=0)
#     mcurr0 = _readData('/home/feng/static/results/0318-makeRNA-polyQ/mcurr0.csv',guess_index=0)
    mcurr0.columns = mcurr0.columns.str.lower()
    mcurr0 = mcurr0.rename(columns={
        'dataacc':'data_acc',
        'gtype':'genotype',
        'temp':'temperature',
    })
    const_dict = {
        'tissue':'seedling',
        'library_strategy':'ChIP-Seq',
        'library_molecule':'genomic DNA',
        'species':"txid3702", ### arabidopsis
        "cultivar":"Col-0",
        "genome":"TAIR10",
        "replicate":"single",
        
        "genotype":nan,
        "age":nan,
        "ztime":nan,
        "light":nan,
        "temperature":nan,
        "fragmentation_method":"sonication",
        "chip_antibody":"{{chip_antibody}}",
        "growth_protocol":u"Seeds were sterilised and sown on ½ X Murashige and Skoog-agar (MS-agar) plates at pH 5.7. Seedlings were grown at 22˚C and 170 µmol/m2/s light with a 16-hour long-day photoperiod for 8-10 days and then treated at ZT1 or ZT2.5 with 30 µM cycloheximide in deionised water at 22˚C or mock controls (0.1% v/v DMSO) at 22˚C or 4˚C.  ",
        "extraction_protocol": u'''
    !Sample_extract_protocol = 3g seedlings for each set were fixed under vacuum for 20 min in 1xPBS  (10 mM PO43−, 137 mM NaCl, and 2.7 mM KCl) containing 1% Formaldehyde (F8775 SIGMA). The reaction was quenched by adding glycine to a final concentration of 62 mM. Chromatin immunoprecipitation (ChIP) was performed as described (12), with the exception that 100 µl FLAG M2 agarose affinity gel antibody was used (A2220 SIGMA-Aldrich) per sample. 
    !Sample_library_construction_protocol = Sequencing libraries were prepared using TruSeq ChIP Sample Preparation Kit (Illumina  IP-202-1024) and samples sequenced on NextSeq500 machine from Illumina using NextSeq® 500/550 High Output Kit v2 (75 cycles) TG-160-2005
             '''
        ,
        "processing_protocol":"{{sample.processing_protocol}}",
        "rawfile_detail":"{{sample.rawfile_detail}}"
    }
    for key,v in const_dict.items():
        mcurr0[key]= v
    # print mcurr0.iloc[:1].to_csv(sep='\t',index=0)
    if DEBUG():
        dfc = mcurr0.iloc[:1]
    else:
        dfc = mcurr0.loc[mcurr0["data_acc"].isin(DATA_ACC_CHIPSEQ())]
    assert dfc.data_acc.is_unique
    return dfc

import pymisca.ext as pyext
from pymisca import shell
import path



INPUTDIR  = lambda res = pyext.path.Path('.').realpath():res
# def sample_get_title(sample):
def sample_init_full():
    its = [
        pyext.df__iterdict( meta_df_chipseq().fillna("NA")),
        pyext.df__iterdict( meta_df_rnaseq().fillna("NA")),
                           ]
    it = sum(map(list,its),[])
                               
                           
# #     it = pyext.df__iterdict( meta_df_rnaseq().fillna("NA"))
#     it =( list(pyext.df__iterdict( meta_df_chipseq().fillna("NA")))
#     +  list(pyext.df__iterdict( meta_df_rnaseq().fillna("NA")) ))
    def _worker(sample):
        sample = attrdict.AttrDict(sample)
        sample.title = "_".join([sample[k] for k in "data_acc,age,tissue,genotype,ztime,temperature".split(",")])
        return sample
    
    return map(_worker, it)
def get_soft_text():
#     OUTPUT_FILE = WORKDIR() / _fkey() ("get_soft_text.tar.gz")
#     dfc = meta_df_rnaseq()
    # BUFFER_ALL = BUFFER_RNASEQ() + BUFFER_CHIPSEQ()

    FORCE = 1

    OUTDIR = WORKDIR()/_fkey()
    shell.real__dir(dirname= OUTDIR)
    
    if not FORCE and shell.file__notEmpty(OUTDIR + '.tar.gz'):
        pass
    else:
        # df_figureRegistry().to_csv(OUTDIR / 'figureRegistry.csv',index=0)
        def _worker(sample):
            sample = attrdict.AttrDict(sample)
            sample.title = "_".join([sample[k] for k in "data_acc,age,tissue,genotype,ztime,temperature".split(",")])
            # print(template_common())
            res  = res= pyext.jf2( template_common())
            res = '\n'.join([x.strip() for x in res.splitlines()])
            sample.soft_text = res        
            pyext.printlines([sample.soft_text], OUTDIR / pyext.f("{sample.data_acc}.soft.txt"))
#             %sample.data_acc)

            
        [_worker(sample) for sample in pyext.df__iterdict( meta_df_rnaseq().fillna("NA"))]
        [_worker(sample) for sample in pyext.df__iterdict( meta_df_chipseq().fillna("NA"))]
        with pyext.TempDirScope(getTempDirName=lambda: OUTDIR/'_temp',
                                force=1) as tdir:
            tdir = tdir.d
#         tdir = OUTDIR / '_temp'
            shell.real__dir(dirname= tdir)
#             shell.shellexec(' '.join(["tar" ,"-C",tdir,"-xvzf", INPUTDIR()/"./src/polyq-0830.get_soft_text_jaehoon.tar.gz",]))
#             shell.shellexec(' '.join(["tar" ,"-C",tdir,"-xvzf", INPUTDIR()/"./src/polyq-0905.get_soft_text_Jaehoon0906.tar.gz",]))
    
#             shell.shellexec(' '.join(["tar" ,"-C",tdir,"-xvzf", INPUTDIR()/"./src/polyq-0907.katja.get_soft_text.tar.gz",]))
    
            shell.shellexec(' '.join(["unzip","-d",tdir, "-o",INPUTDIR()/"./src/polyq-0907-jaehoon.get_soft_text.zip", ]))
            for fname in tdir.glob('*.soft.txt'):
                DATA_ACC = fname.basename().replace('.soft.txt','',1)
    #             in df_figureRegistry()['DATA_ACC'].tolist()
                assert DATA_ACC in df_figureRegistry()['DATA_ACC'].tolist(),(DATA_ACC,  "JH0905")
                print('[JH0906]',fname.basename(),)
                pyext.file__link(fname, OUTDIR/fname.basename(),force=1)

            
        CMD = [
            "cd",OUTDIR,
            "&&tar","-cvzf",OUTDIR.realpath()+'.tar.gz',"*",
        ]
        CMD = ' '.join(pyext.stringList__flatten(CMD))
        res = shell.shellexec(CMD)
#         pd.DataFrame()
        
@_dac
def rawMeta():
    mcurr = _readData(
        '/home/feng/envs/upGeo/results/0424-database-raw/mcurr.csv',
        guess_index=0)
    mcurr = mcurr.loc[~mcurr['FULL_PATH'].str.contains("Raw_data/184R_Q_reseq181_combined")]
    return mcurr
#     mcurr = _readData(
#         '/home/feng/envs/upGeo/results/0424-database-raw/mcurr.csv',
#         guess_index=0)
#     mcurr = mcurr.loc[~mcurr['FULL_PATH'].str.contains("Raw_data/184R_Q_reseq181_combined")]
#     mcurr = mcurr.loc[~mcurr["FULL_PATH"].isin(['/home/feng/writable/teampw/__backup/syno3/raw-data/PW_HiSeq_data/RNA-seq/Raw_data/184R_Q_reseq181_combined/181RESEQ_Q_870/181R-Q-870-4196-ZT12-27C-adult_S22_L005_R1_001.fastq.gz'],
#        ['/home/feng/writable/teampw/__backup/syno3/raw-data/PW_HiSeq_data/RNA-seq/Raw_data/184R_Q_reseq181_combined/181RESEQ_Q_870/181R-Q-870-4196-ZT12-27C-adult_S22_L006_R1_001.fastq.gz'],
#        ['/home/feng/writable/teampw/__backup/syno3/raw-data/PW_HiSeq_data/RNA-seq/Raw_data/184R_Q_reseq181_combined/181RESEQ_Q_870/181R-Q-870-4196-ZT12-27C-adult_S22_L007_R1_001.fastq.gz'],
#        ['/home/feng/writable/teampw/__backup/syno3/raw-data/PW_HiSeq_data/RNA-seq/Raw_data/184R_Q_reseq181_combined/181RESEQ_Q_870/181R-Q-870-4196-ZT12-27C-adult_S22_L008_R1_001.fastq.gz']
#                                              )]
#     rawMeta = mcurr    
# get_soft_text(meta_df_rnaseq())
# get_soft_text(meta_df_chipseq())

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
'''.replace(' ','').strip().splitlines()
    
    for k in blacklist:
        dfc = dfc.loc[~dfc.index.str.contains(k)]
    dfc = dfc.loc[~dfc['FILEACC'].str.endswith("orig.bam")]
    dfc = dfc.loc[~dfc['FILEACC'].str.endswith("unmapped.bam")]
    dfc = dfc.loc[~(~dfc['FILEACC'].str.endswith("RPKM.bw") & dfc['FILEACC'].str.endswith(".bw"))]
    
    return dfc

@_dac
def df_mappedData_chipseq():
    _accs = DATA_ACC_CHIPSEQ()
    dfc = df_mappedData_full()
    dfc = dfc.loc[dfc['DATAACC'].isin(_accs)]
    #     assert (dfc.groupby("EXT").apply(len)==len(set(_accs))).all()
    dfc = dfc.reset_index(drop=0).pivot_table(index='DATAACC',columns='EXT',values='FULL_PATH',aggfunc=lambda x:x)    
    return dfc

@_dac
def df_mappedData_rnaseq():
    _accs = DATA_ACC_RNASEQ()
    dfc = df_mappedData_full()
    dfc = dfc.loc[dfc['DATAACC'].isin(_accs)]
    #     assert (dfc.groupby("EXT").apply(len)==len(set(_accs))).all()
    dfc = dfc.reset_index(drop=0).pivot_table(index='DATAACC',columns='EXT',values='FULL_PATH',aggfunc=lambda x:x)    
    return dfc


def sample_rnaseq_processing_protocol(sample):
    from pymisca.events import CopyEvent,LinkEvent
    from pymisca.ext import f as _f
    import os
    
    OUTDIR = WORKDIR() / sample['data_acc']/ 'supp'
#     sample.data_acc_control = _get_data_acc_control(sample)
    rec = df_mappedData_rnaseq().loc[sample['data_acc']]

    for attrName,key in [
        ('file_count','txt'),
        ('file_bam','bam'),
#         ('file_npk','narrowPeak')
    ]:
        
        fullname= rec[key]
        sample[attrName+'_orig'] = fullname
        if pyext.pd.isnull(fullname):
            sample[attrName] = 'NA'
        else:
            basename = os.path.basename(fullname)
            sample[attrName] = LinkEvent(
                CopyEvent(fullname,
                OUTDIR / basename
                ).dest,
               WORKDIR()/_f("ftp/{sample.data_acc}.supp.{basename}"),
               1
            ).dest.relpath(WORKDIR() / 'ftp')    
    
        
#     sample.file_count = CopyEvent( rec['txt'], 
#                                OUTDIR / os.path.basename(rec['txt'])
#                               ).dest.relpath(WORKDIR())
#     sample.file_bam = CopyEvent(rec['bam'],
#                                 OUTDIR / os.path.basename(rec['bam'])
#                                ).dest.relpath(WORKDIR())
    template = u'''
!Sample_data_processing = Raw fastq were uploaded to Bluebee Genomics Platform and analysed with \
Quantseq FWD analysis pipeline (https://www.lexogen.com/quantseq-data-analysis/). Briefly, the raw reads \
were trimmed with BBDuk and aligned to GTF-annotated genome with STAR. 

!Sample_data_processing = Supplementary_files_format_and_content: *.txt: TSV table containing abundance of transcripts by STAR. 
!Sample_data_processing = Supplementary_files_format_and_content: *.bam: Genomic alignment by STAR.

!Sample_supplementary_file_1 = {{sample.file_bam}}
!Sample_supplementary_file_2 = {{sample.file_count}}
    '''
    res = pyext.jf2(template,)
    return res

def sample_get_rawfile_detail(sample):
    from pymisca.events import LinkEvent
    from pymisca.ext import f as _f
    node = pyext.file__asModule('/home/feng/envs/0726-polyq/src/validate_fastq.py')
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

def sample_template_find_curated(sample):
    sample['template_curated']  = res = u''.join( pyext.readData( WORKDIR() /'get_soft_text'/ sample['data_acc'] + '.soft.txt','it'))
    return res
def sample_get_expt_type(sample):
    c = "C" in sample['data_acc'].split('S')[0]
    r = "R" in sample['data_acc'].split('S')[0]
    if c and not r:
        return "CHIPSEQ"
    elif not c and r:
        return "RNASEQ"
    elif c and r:
        assert 0,(sample['data_acc'],)
    else:
        assert 0,(sample['data_acc'],)
#         return
        
def sample_template_finalise(sample):
    sample_template_find_curated(sample)
    sample.expt_type = sample_get_expt_type(sample)
    if sample['expt_type'] == 'CHIPSEQ':
        sample['processing_protocol']= sample_chipseq_processing_protocol(sample)
    elif sample['expt_type'] == 'RNASEQ':
        sample['processing_protocol']= sample_rnaseq_processing_protocol(sample)
        
    sample['rawfile_detail'] = sample_get_rawfile_detail(sample)    
    sample['template_final'] = res = pyext.jf2(sample['template_curated'])
    return res


# ns.sample_get_rawfile_detail = sample_get_rawfile_detail

def sample_chipseq_processing_protocol(sample):
    from pymisca.events import CopyEvent,LinkEvent
    from pymisca.ext import f as _f
    import os
    def _get_data_acc_control(sample):
        buf = '''
        198C,195CS13
        176C,176CS21
        182C,176CS21
        189C,176CS21
        192C,192CS19
        '''.replace(' ','')
        mapper= dict([x.strip().split(',') for x in buf.splitlines() if x.strip()])
        res = mapper.get(sample.data_acc.split('S')[0])
        return res
    
    OUTDIR = WORKDIR() / sample['data_acc']/ 'supp'
    sample.data_acc_control = _get_data_acc_control(sample)
    rec = df_mappedData_chipseq().loc[sample['data_acc']]

#     fullname= rec['bw']
#     basename = os.path.basename(fullname)
#     sample.file_bw = LinkEvent( CopyEvent( fullname,
#                                           OUTDIR / basename,
#                                          ).dest,
#                                WORKDIR()/_f("ftp/{sample.data_acc}.supp.{basename}"),
#                               1,).dest.relpath(WORKDIR()/ 'ftp')
    for attrName,key in [
        ('file_bw','bw'),
        ('file_bam','bam'),
        ('file_npk','narrowPeak')]:
        
        fullname= rec[key]
        sample[attrName+'_orig'] = fullname
        if pyext.pd.isnull(fullname):
            sample[attrName] = 'NA'
        else:
            basename = os.path.basename(fullname)
            sample[attrName] = LinkEvent(
                CopyEvent(fullname,
                OUTDIR / basename
                ).dest,
               WORKDIR()/_f("ftp/{sample.data_acc}.supp.{basename}"),
               1
            ).dest.relpath(WORKDIR() / 'ftp')
    #     if not len(rec['narrowPeak']):
#     if pyext.pd.isnull(rec['narrowPeak']):
#         sample.file_npk = 'NA'
#     else:
#         fullname= rec['narrowPeak']
#         basename = os.path.basename(fullname)
#         sample.file_npk = LinkEvent(
#             CopyEvent(fullname,
#                       OUTDIR /  basename,
#                       ).dest,
#            WORKDIR()/_f("ftp/{sample.data_acc}.supp.{basename}"),
#            1 
#         ).dest.relpath(WORKDIR()/'ftp')
        
    template = u'''
!Sample_data_processing = Adapters were trimmed off from raw reads with Trimmomatic with argument "ILLUMINACLIP:$FA_ADAPTER:6:30:10 LEADING:3 TRAILING:3 MINLEN:36 SLIDINGWINDOW:4:15". \
Raw reads were mapped to the genome "TAIR10" with Bowtie2 under argument:"--no-mixed --no-discordant --no-unal -k2". Any read that mapped to more than one genomic location was discarded. \
PCR duplicate reads were removed with Picard using default setting.
!Sample_data_processing = Genomic binding profile was quantified in RPKM (Reads Per Kilobase per Million mapped reads) using a bin-size of 10bp. "deeptools.bamCoverage" is used.
!Sample_data_processing = For each treated ChIP-Seq library, peaks were called against a control {{sample.data_acc_control}} using MACS2 with argument "--keep-dup 1 -p 0.1".
!Sample_data_processing = Supplementary_files_format_and_content: *.bam: Genomic alignements that were sorted, deduplicated and filtered for uniq-mapped reads.

!Sample_data_processing = Supplementary_files_format_and_content: *_RPKM.bw: RPKM-normalised bigwig track at 10bp resolution
!Sample_data_processing = Supplementary_files_format_and_content: *.narrowPeak: containing MACS2-called peaks. as described
!Sample_supplementary_file_1 = {{sample.file_bam}}
!Sample_supplementary_file_2 = {{sample.file_bw}}
!Sample_supplementary_file_3 = {{sample.file_npk}}
    '''
    res = pyext.jf2(template,)
    return res

def sample_get_depth(sample,
                     OUTPUT_FILE,
                     FORCE=0):
    if not FORCE and pyext.file__notEmpty(OUTPUT_FILE):
        pass
    else:
        d = pyext._DICT_CLASS()
        d['DATA_ACC']=sample['data_acc']
        CMD = ["gzip","-d<",sample['rawfile_files_orig'][0],'|wc','-l']
        res = pyext.shellexec(' '.join(CMD))
        d['READ_COUNT_RAW'] = int(res.strip())//4
        CMD = ["cat",sample["file_bam_orig"],
               "|samtools","view","-F0x4","-F0x100","-c",
              ]
        res = pyext.shellexec(' '.join(CMD))        
        d['READ_COUNT_UNIQ_MAPPED'] = int(res.strip())
        with open(OUTPUT_FILE,'w') as f:
            pyext.json.dump(d, f,indent=4)
    return pyext.readData(OUTPUT_FILE)

def dump_sequencing_info(sample):    
    template = '''
sample_accession={sample.data_acc}
number_of_reads={pyext.size__humanReadable(read_count_raw,suffix='')}
number_of_uniq_mapped_reads={pyext.size__humanReadable(read_count_uniq_mapped,suffix='')}
rawfile_readlengths={sample.rawfile_readlengths[0]}
rawfile_is_paired={sample.rawfile_is_paired}

    '''.strip()
    return pyext.template__format(template,sample)
def dict_sanitise(sample):
    for k in sample.keys():
        if k.isupper():
            sample[str(k.lower())] = sample.pop(k)
    return sample

def dump_peak_quality(sample):
    if not pyext.file__notEmpty(sample['npkfile']):
        print (sample['data_acc'],sample['npkfile'])
        template= '''
        sample_accession={sample.data_acc}
        number_of_peaks_below_5%FDR=NA
        number_of_peaks_above_5fold_enrichment=NA
        '''.strip().replace('\t','').replace(' ','')
    else:
        df = _readData(sample['npkfile'],columns=pyext.columns.bed,header=None)

        template= '''
        sample_accession={sample.data_acc}
        number_of_peaks_below_5%FDR={((df['neglogQval'] > -pyext.np.log10(0.05) )).sum()}
        number_of_peaks_above_5fold_enrichment={(df['FC'] > 5 ).sum()}
        '''.strip().replace('\t','').replace(' ','')
    
    return pyext.template__format(template,sample)

if __name__ == '__main__':
    get_soft_text()

#     mcurr = _readData(
#         '/home/feng/envs/upGeo/results/0424-database-raw/mcurr.csv',
#         guess_index=0)
#     mcurr = mcurr.loc[~mcurr['FULL_PATH'].str.contains("Raw_data/184R_Q_reseq181_combined")]

#     mcurr = mcurr.loc[~mcurr["FULL_PATH"].isin(['/home/feng/writable/teampw/__backup/syno3/raw-data/PW_HiSeq_data/RNA-seq/Raw_data/184R_Q_reseq181_combined/181RESEQ_Q_870/181R-Q-870-4196-ZT12-27C-adult_S22_L005_R1_001.fastq.gz'],
#        ['/home/feng/writable/teampw/__backup/syno3/raw-data/PW_HiSeq_data/RNA-seq/Raw_data/184R_Q_reseq181_combined/181RESEQ_Q_870/181R-Q-870-4196-ZT12-27C-adult_S22_L006_R1_001.fastq.gz'],
#        ['/home/feng/writable/teampw/__backup/syno3/raw-data/PW_HiSeq_data/RNA-seq/Raw_data/184R_Q_reseq181_combined/181RESEQ_Q_870/181R-Q-870-4196-ZT12-27C-adult_S22_L007_R1_001.fastq.gz'],
#        ['/home/feng/writable/teampw/__backup/syno3/raw-data/PW_HiSeq_data/RNA-seq/Raw_data/184R_Q_reseq181_combined/181RESEQ_Q_870/181R-Q-870-4196-ZT12-27C-adult_S22_L008_R1_001.fastq.gz']
#                                              )]
#     rawMeta = mcurr

    mappedMeta = _readData(
        '/home/feng/work/results/0407-database-mapped/mcurr.csv',
                                    guess_index=0)
    for DATA_ACC in DATA_ACC_RNASEQ() + DATA_ACC_CHIPSEQ():
        print DATA_ACC
        node = pyext.file__asModule('/home/feng/envs/0726-polyq/src/validate_fastq.py')
        node.rawMeta = rawMeta()
        node.DATA_ACC = DATA_ACC
    #     '150RS1'
        node.WORKDIR = pyext.path.Path('/home/feng/envs/0726-polyq/WORKDIR.submit/').realpath()
#         node.WORKDIR = pyext.path.Path('/home/feng/envs/0830-polyq/WORKDIR/').realpath()
        node.valid_fastq()
        node.combined_valid_fastq()
#         break
        
        
    print("[DONE]")