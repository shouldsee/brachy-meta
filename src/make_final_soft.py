# -*- coding: utf-8 -*-
import pymisca.ext as pyext
ns = pyext.file__asModule('./get_meta_soft.py')

# dfc = node.df_figureRegistry()

# import pymisca.ext as pyext
# node = pyext.file__asModule('src/get_meta_soft.py')
# node.WORKDIR = lambda: pyext.path.Path('./WORKDIR.submit/').realpath()
ns.get_soft_text()
pyext.MDFile('./WORKDIR.submit/get_soft_text.tar.gz')



WORKDIR = ns.WORKDIR
# OUTDIR = WORKDIR() / "final_soft"
OUTDIR = WORKDIR() / "get_soft_text"
(WORKDIR()/'ftp').makedirs_p()
pyext.real__dir(dirname = OUTDIR)
_samples = ns.sample_init_full()


for sample in _samples:
    print '[template_finalise]',sample['data_acc'],'...'
    # try:
    if 1:
        ns.sample_template_find_curated(sample)
        #     continue
        ns.sample_template_finalise(sample)
        res = sample['template_final']
        res = '\n'.join([x.strip() for x in res.splitlines()])
        sample.soft_text = res        
        pyext.printlines([sample.soft_text], OUTDIR / pyext.f("{sample.data_acc}.autofilled.soft.txt"))
    try:
    	pass
    except Exception as e:
        print('FAILED')
        print(str(e))
        

template=  u'''
^SERIES = 0829-polyq
!Series_title = RNA-Seq and ChIP-Seq profiling of ELF3, an prion-like domain-containig in ELF3 that functions as a
thermosensor in Arabidopsis.
!Series_summary = Temperature is a major environmental variable governing plant growth and
development. ELF3 contains a polyglutamine (polyQ) repeat 8–10, embedded within a predicted prion domain (PrD). We find the length of the polyQ repeat correlates with thermal responsiveness. Plants from hotter climates appear to have lost the PrD domain, and these versions of ELF3 are stable at high temperature and lack thermal responsiveness. ELF3 temperature sensitivity is also modulated by the levels of ELF4, indicating that ELF4 can stabilise ELF3 function. This RNA-Seq dataset provides evidence for the hypothetical ELF3 function of temperature sensing .
!Series_overall design = Single samples were taken at each time point. RNA-Seqs and ChIP-Seqs were performed for different genotypes at different temperature and objective time. 
!Series_contributor = Jaehoon Jung
!Series_contributor = Katja, Jaeger
!Series_contributor = Feng, Geng
{% for sample in _samples %}
!Series_sample_id = {{sample.data_acc}}{% endfor %}
'''
res = pyext.jf2(template)
pyext.printlines([res], OUTDIR / pyext.f("SERIES.soft.txt"))

from pymisca import shell

OFNAME = "0830-polyq-submit.soft"
shell.shellexec(" ".join([
    "cd",OUTDIR,
    "&&cat","SERIES.soft.txt","*autofill*",
                          "|grep","-v","^#",
                          "|grep","-v","^$",
                          "|sed",u"'s/^\xEF\xBB\xBF//g'",
                          ">" + OUTDIR / OFNAME,
#                           "|tee","0830-polyq-submit.soft",
]),silent=1)
pyext.file__link(OUTDIR/OFNAME, WORKDIR()/"ftp"/OFNAME, force=1)
                
CMD = [
    "cd",OUTDIR,
    "&&tar","-cvzf",OUTDIR.realpath()+'.tar.gz',"*",
]
CMD = ' '.join(pyext.stringList__flatten(CMD))
res = shell.shellexec(CMD)
# pyext.printlines([pyext.jf2(template)], WORKDIR()/ 'get_soft_text' / '0829-polyq.SERIES.soft.txt')
#     break