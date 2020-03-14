import path

import pymisca.ext as pyext
from pymisca.ext import df__format,df__paste0,shellexec
import pandas as pd

from delay_and_cache import delay_and_cache as _dac
from delay_and_cache import StaticClass,_ctf,_fkey

import synotil.ptn
# import 

from pymisca.atto_jobs_list.file_md5 import file_md5
from pymisca.atto_jobs_list.fastqFiles_combine import fastqFiles_combine


def rawMeta(): raise ValueError("Must be pd.DataFrame()")
def WORKDIR(): raise ValueError("Must be pd.DataFrame()")
    
def valid_fastq(
  DATA_ACC,
  rawMeta,
               ):
    assert DATA_ACC is not None
    def rawFile__validateChunk(dfc,idKeys = None):
        '''Validate the df_raw
    '''
        if idKeys is None:
            idKeys = ['runID','sampleID']
        idKeys = list(idKeys)
        dfc = dfc.sort_values(idKeys + ['read', 'chunk'])
        gp  = dfc.groupby(idKeys + ['read'])
        for (key,df) in gp:
            assert len(df) == 4,key
        if not dfc['ext'].iloc[0].startswith('.'):
            dfc['ext'] = dfc['ext'].map(lambda x:'.%s'%x)
        dfc['fnameCombined'] = pyext.df__paste0(dfc,
                                                  idKeys + [
                                                      'read',
                                                      'ext'
                                                  ],
                                                  sep='_',
                                                  )
        dfc['fnameCombinedSize'] = 0
        return dfc

    rawCurr = rawMeta().loc[rawMeta()["DATAACC"]==DATA_ACC]
    # rawCurr = rawMeta.query(query)
    rawCurr = pd.concat([rawCurr,rawCurr['BASENAME'].str.extract(synotil.ptn.baseSpace)],axis=1)
    rawCurr['fname'] = rawCurr['FULL_PATH']

    rawCurr = rawCurr.sort_values(['DATAACC','BASENAME','SIZE'],ascending=False
                                 ).groupby(['BASENAME','DATAACC']).first()
    rawCurr = rawCurr.sort_values(['BASENAME'])

    DF_LIST = list(x for x in rawCurr.groupby(['RUN_ID','SAMPLE_ID','read','ext']))
    for (_,_,read,ext),df in DF_LIST:
        assert int(read) in [1,2],(DATA_ACC,read,)
        assert ext in ["fastq","fastq.gz"],(DATA_ACC,ext,)
        assert len(df) in [1,4],pyext.ppJson((
            DATA_ACC(),len(df),
            df[['FULL_PATH',]].values))
#         assert df.read.un
    return dict(
        DATA_ACC=DATA_ACC,
        DF_LIST = DF_LIST,
             )
    
def combined_valid_fastq(DATA_ACC,valid_fastq, WORKDIR):
    it = valid_fastq()
    res = []
    for k,df in valid_fastq()['DF_LIST']:
        res += [fastqFiles_combine({"INPUT_FILE_LIST":df.fname,
                            'OUTDIR':WORKDIR()/DATA_ACC()/_fkey(),
                                   })]
                
#                                    DATA_ACC()/})]
    assert len(res) in [1,2], (DATA_ACC(), res)
    return dict(
        DATA_ACC=DATA_ACC(),
        OUTPUT_NODES= res,
             )    
#     return res