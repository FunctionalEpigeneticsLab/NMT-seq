#!/usr/bin/env python
import os
import re
import gzip
import argparse
import subprocess

DNAlib = "NMT_default"
chroms = ['chr{}'.format(x) for x in list(range(1,23))+['X', 'Y', 'M']]

def sagg(a, b):
    return(a+b)

def map_sagg(atup, ltup):
    ltup = tuple(map(int, ltup))
    stup = tuple(map(sagg, atup, ltup))
    return stup

def run_stats(workdir, outdir, qcdir, mapdir, plotMbias, Rscript):
    workdir = os.path.abspath(workdir)
    print('Working dir: '+workdir)
    SampleIDs_list = os.listdir(os.path.join(workdir, mapdir))
    SampleIDs = [s for s in SampleIDs_list if 'GC' in s and 'downsampled' not in s]
    outdir = os.path.join(workdir, outdir)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    fastqc_dir = os.path.join(workdir, qcdir)
    DNAalign_workdir = os.path.join(workdir, mapdir)
    outfh = os.path.join(outdir, 'DNAalign.batch.stats.tsv')

    chroms = ['chr{}'.format(x) for x in list(range(1,23))+['X', 'Y', 'M']]
    
    if plotMbias!='n':
        if not os.path.isfile(plotMbias):
            print('Warning: plotMbias script path does not exist!')

    with open(outfh, 'w') as fo:
        #fo.write('Sample\tRaw_readpairs\tProcessed_seq\tUniq_map\tUniq_map_rate\tNo_map_rate\tMulti_map_rate\tDiscarded_seq\tTopS_map\tComp_topS_map\tComp_botS_map\tBotS_map\tMethylated_CpG_rate\tMethylated_CHG_rate\tMethylated_CHH_rate\tMethylated_Cunknown_rate\tDuplication_rate\tFiltered_fragments\tFiltered_Methylated_CpG_rate\tFiltered_Methylated_CHG_rate\tFiltered_Methylated_CHH_rate\tFiltered_CpG_positions\tFiltered_GpC_positions\tGCG_mrate\tGCHH_mrate\tGCHG_mrate\tCCG_mrate\tCCHH_mrate\tCCHG_mrate\tHCG_mrate\tHCHH_mrate\tHCHG_mrate\tHCH_mrate\n')
        fo.write('Sample\tRaw_readpairs\tProcessed_seq\tUniq_map\tUniq_map_rate\tNo_map_rate\tMulti_map_rate\tDiscarded_seq\tTopS_map\tComp_topS_map\tComp_botS_map\tBotS_map\tMethylated_CpG_rate\tMethylated_CHG_rate\tMethylated_CHH_rate\tMethylated_Cunknown_rate\tDuplication_rate\tFiltered_fragments\tFiltered_Methylated_CpG_rate\tFiltered_Methylated_CHG_rate\tFiltered_Methylated_CHH_rate\n')
        for sampleid in SampleIDs:
            sample_sub_workdir = os.path.join(DNAalign_workdir, sampleid)
            (anaseq, uniqmap, nomap, multimap, discardseq, chrMcount, topmap, comptopmap, compbotmap, botmap, mCpG, mCHG, mCHH, mCunknown, umCpG, umCHG, umCHH, umCunknown) = tuple([0]*18)
            (duppairs, filterfrag, dedupmCpG, dedupmCHG, dedupmCHH, dedupumCpG, dedupumCHG, dedupumCGG) = tuple([0]*8)
            #(nomeCpGnum, nomeGpCnum) = (0, 0)
            (uniqmaprate, nomaprate, multimaprate, discardrate, mCpGrate, mCHGrate, mCHHrate, mCunknownrate, duprate, dedupmCpGrate, dedupmCHGrate, dedupmCHHrate) = tuple([-1]*12)
            (GmCG, GumCG, GmCHH, GumCHH, GmCHG, GumCHG, CmCG, CumCG, CmCHH, CumCHH, CmCHG, CumCHG, HmCG, HumCG, HmCHH, HumCHH, HmCHG, HumCHG, GCG_mrate, GCHH_mrate, GCHG_mrate, CCG_mrate, CCHH_mrate, CCHG_mrate, HCG_mrate, HCHH_mrate, HCHG_mrate, HCH_mrate) = tuple([0]*28)
            procline0 = '******\nProcessing sample: %s\n******' % (sampleid)
            print(procline0)

            zip_file = os.path.join(fastqc_dir, sampleid+'_R1_001_fastqc.zip')
            procline00 = 'Processing fastqc output: %s' % (zip_file)
            print(procline00)
            command = f"unzip -p {zip_file} */fastqc_data.txt | grep 'Total Sequences' | cut -f 2"
            process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()
            totalreadpairs = -1
            if process.returncode == 0:
                rawreadpairs = int(stdout.decode().strip())
            else:
                print(f"Error processing {zip_file}: {stderr.decode()}")

            for samfile in os.listdir(sample_sub_workdir):
                samfile = '%s/%s' % (sample_sub_workdir, samfile)
                #if samfile.endswith("merge.deduplicated.filtered.sorted.bam"):
                    #procline1 = 'Count reads mapped to chrM from bam: %s' % (samfile)
                    #print(procline1)
                    #chrMcountcmd = '%s view -f 3 %s chrM | wc -l' % (samtools, samfile)
                    #subprocess.call(chrMcountcmd, shell=True)
                    #chrMcount = subprocess.run(chrMcountcmd, shell=True,stdout=subprocess.PIPE)
                    #chrMcount = chrMcount.stdout.decode()
                    #chrMcount = int(chrMcount)

                if samfile.endswith("PE_report.txt") or samfile.endswith("SE_report.txt"):
                    procline2 = 'Processing Bismark lane/read-end report file: %s' % (samfile)
                    print(procline2)
                    (lanaseq, luniqmap, lnomap, lmultimap, ldiscardseq, ltopmap, lcomptopmap, lcompbotmap, lbotmap, lmCpG, lmCHG, lmCHH, lmCunknown, lumCpG, lumCHG, lumCHH, lumCunknown) = tuple([0]*17)

                    with open(samfile, 'r') as fh0:
                        for l0 in fh0:
                            if 'analysed in total' in l0:
                                lanaseq = l0.strip().split('\t')[1]
                            elif 'unique best hit' in l0:
                                luniqmap = l0.strip().split('\t')[1]
                            elif 'no alignments under any condition' in l0:
                                lnomap = l0.strip().split('\t')[1]
                            elif 'did not map uniquely' in l0:
                                lmultimap = l0.strip().split('\t')[1]
                            elif 'discarded because genomic sequence could not be extracted' in l0:
                                ldiscardseq = l0.strip().split('\t')[1]
                            elif l0.startswith('CT/GA/CT:') or l0.startswith('CT/CT:'):
                                ltopmap = l0.strip().split('\t')[1]
                            elif l0.startswith('GA/CT/CT:') or l0.startswith('GA/CT:'):
                                lcomptopmap = l0.strip().split('\t')[1]
                            elif l0.startswith('GA/CT/GA:') or l0.startswith('GA/GA:'):
                                lcompbotmap = l0.strip().split('\t')[1]
                            elif l0.startswith('CT/GA/GA:') or l0.startswith('CT/GA:'):
                                lbotmap = l0.strip().split('\t')[1]
                            elif 'Total methylated C\'s in CpG context' in l0:
                                lmCpG = l0.strip().split('\t')[1]
                            elif 'Total methylated C\'s in CHG context' in l0:
                                lmCHG = l0.strip().split('\t')[1]
                            elif 'Total methylated C\'s in CHH context' in l0:
                                lmCHH = l0.strip().split('\t')[1]
                            elif 'Total methylated C\'s in Unknown context' in l0:
                                lmCunknown = l0.strip().split('\t')[1]
                            elif 'Total unmethylated C\'s in CpG context' in l0:
                                lumCpG = l0.strip().split('\t')[1]
                            elif 'Total unmethylated C\'s in CHG context' in l0:
                                lumCHG = l0.strip().split('\t')[1]
                            elif 'Total unmethylated C\'s in CHH context' in l0:
                                lumCHH = l0.strip().split('\t')[1]
                            elif 'Total unmethylated C\'s in Unknown context' in l0:
                                lumCunknown = l0.strip().split('\t')[1]


                    linput = (lanaseq, luniqmap, lnomap, lmultimap, ldiscardseq, ltopmap, lcomptopmap, lcompbotmap, lbotmap, lmCpG, lmCHG, lmCHH, lmCunknown, lumCpG, lumCHG, lumCHH, lumCunknown)
                    tinput = (anaseq, uniqmap, nomap, multimap, discardseq, topmap, comptopmap, compbotmap, botmap, mCpG, mCHG, mCHH, mCunknown, umCpG, umCHG, umCHH, umCunknown)
                    (anaseq, uniqmap, nomap, multimap, discardseq, topmap, comptopmap, compbotmap, botmap, mCpG, mCHG, mCHH, mCunknown, umCpG, umCHG, umCHH, umCunknown) = map_sagg(tinput, linput)

                elif samfile.endswith("deduplication_report.txt"):
                    procline3 = 'Processing Bismark duplication summary file: %s' % (samfile)
                    print(procline3)
                    # in case of processing pairend in single-end form, aggregate read data
                    rduppairs = 0

                    with open(samfile, 'r') as fh1:
                        for l1 in fh1:
                            if 'Total number duplicated alignments removed' in l1:
                                rduppairs = l1.strip().rsplit(' ',1)[0]
                                rduppairs = rduppairs.split('\t')[1]
                                rduppairs = int(rduppairs)

                    duppairs += rduppairs

                elif samfile.endswith("_splitting_report.txt"):
                    procline4 = 'Processing Bismark C count report file: %s' % (samfile)
                    print(procline4)

                    with open(samfile, 'r') as fh2:
                        for l2 in fh2:
                            if 'lines in total' in l2:
                                filterfrag = l2.strip().split(' ')[1]
                            elif 'Total methylated C\'s in CpG context' in l2:
                                dedupmCpG = l2.strip().split('\t')[1]
                            elif 'Total methylated C\'s in CHG context' in l2:
                                dedupmCHG = l2.strip().split('\t')[1]
                            elif 'Total methylated C\'s in CHH context' in l2:
                                dedupmCHH = l2.strip().split('\t')[1]
                            elif 'Total C to T conversions in CpG context' in l2:
                                dedupumCpG = l2.strip().split('\t')[1]
                            elif 'Total C to T conversions in CHG context' in l2:
                                dedupumCHG = l2.strip().split('\t')[1]
                            elif 'Total C to T conversions in CHH context' in l2:
                                dedupumCHH = l2.strip().split('\t')[1]
                    (filterfrag, dedupmCpG, dedupmCHG, dedupmCHH, dedupumCpG, dedupumCHG, dedupumCHH) = tuple(map(int, (filterfrag, dedupmCpG, dedupmCHG, dedupmCHH, dedupumCpG, dedupumCHG, dedupumCHH)))

                #For NOMe-seq output, only count number of sites on autosome, sex and mitochondria chromosomes
                elif samfile.endswith(".NOMe.CpG.cov.filtered.gz"):
                    procline5 = 'Processing Bismark NOMe CpG coverage filtered file: %s' % (samfile)
                    print(procline5)

                    with gzip.open(samfile, 'r') as fh3:
                        for l3 in fh3:
                            l3 = l3.decode('utf-8')
                            curchrinfo = l3.strip().split('\t')[0]
                            if curchrinfo in chroms:
                                nomeCpGnum += 1
                            else:
                                continue

                elif samfile.endswith(".NOMe.GpC.cov.filtered.gz"):
                    procline6 = 'Processing Bismark NOMe GpC coverage filtered file: %s' % (samfile)
                    print(procline6)

                    with gzip.open(samfile, 'r') as fh4:
                        for l4 in fh4:
                            l4 = l4.decode('utf-8')
                            curchrinfo = l4.strip().split('\t')[0]
                            if curchrinfo in chroms:
                                nomeGpCnum += 1
                            else:
                                continue

                elif samfile.endswith(".M-bias.txt") and os.path.isfile(plotMbias):
                    if plotMbias!='n' and os.path.isfile(plotMbias):
                        procline7 = 'Plotting M-bias: %s' % (samfile)
                        print(procline7)
                        if os.path.isfile(plotMbias):
                            rscriptcmd = '%s %s %s %s %s' % (Rscript, plotMbias, samfile, sampleid, DNAlib)
                            subprocess.call(rscriptcmd, shell=True)
                        else:
                           print('Proceed without plotting M-bias...')

                elif samfile.endswith(".cytosine_context_summary.txt"):
                    procline8 = 'Processing Bismark cytosine context summary: %s' % (samfile)
                    print(procline8)
                    (GmCG, GumCG, GmCHH, GumCHH, GmCHG, GumCHG, CmCG, CumCG, CmCHH, CumCHH, CmCHG, CumCHG, HmCG, HumCG, HmCHH, HumCHH, HmCHG, HumCHG, GCG_mrate, GCHH_mrate, GCHG_mrate, CCG_mrate, CCHH_mrate, CCHG_mrate, HCG_mrate, HCHH_mrate, HCHG_mrate, HCH_mrate) = tuple([0]*28)

                    with open(samfile, 'r') as fh5:
                        for l5 in fh5:
                            if not l5.startswith('upstream'):
                                (upbase, trinuc, fournuc, cntmc, cntumc, permc) = l5.strip().split('\t')
                                (cntmc, cntumc) = tuple(map(int, (cntmc, cntumc)))

                                if (upbase == "G"):
                                    if (trinuc == "CGA" or trinuc == "CGT" or trinuc == "CGG" or trinuc == "CGC"):
                                        GmCG += cntmc
                                        GumCG += cntumc
                                    elif (trinuc == "CAG" or trinuc == "CTG" or trinuc == "CCG"):
                                        GmCHG += cntmc
                                        GumCHG += cntumc
                                    else:
                                        GmCHH += cntmc
                                        GumCHH += cntumc
                                elif (upbase == "C"):
                                    if (trinuc == "CGA" or trinuc == "CGT" or trinuc == "CGG" or trinuc == "CGC"):
                                        CmCG += cntmc
                                        CumCG += cntumc
                                    elif (trinuc == "CAG" or trinuc == "CTG" or trinuc == "CCG"):
                                        CmCHG += cntmc
                                        CumCHG += cntumc
                                    else:
                                        CmCHH += cntmc
                                        CumCHH += cntumc
                                elif (upbase == "A" or upbase == "T"):
                                    if (trinuc[0:2] == "CG"):
                                        HmCG += cntmc
                                        HumCG += cntumc
                                    elif (trinuc == "CAG" or trinuc == "CTG" or trinuc == "CCG"):
                                        HmCHG += cntmc
                                        HumCHG += cntumc
                                    else:
                                        HmCHH += cntmc
                                        HumCHH += cntumc
                    GCG_mrate = '{0:.2%}'.format(GmCG/(GmCG+GumCG))
                    GCHH_mrate = '{0:.2%}'.format(GmCHH/(GmCHH+GumCHH))
                    GCHG_mrate = '{0:.2%}'.format(GmCHG/(GmCHG+GumCHG))
                    CCG_mrate = '{0:.2%}'.format(CmCG/(CmCG+CumCG))
                    CCHH_mrate = '{0:.2%}'.format(CmCHH/(CmCHH+CumCHH))
                    CCHG_mrate = '{0:.2%}'.format(CmCHG/(CmCHG+CumCHG))
                    HCG_mrate = '{0:.2%}'.format(HmCG/(HmCG+HumCG))
                    HCHH_mrate = '{0:.2%}'.format(HmCHH/(HmCHH+HumCHH))
                    HCHG_mrate = '{0:.2%}'.format(HmCHG/(HmCHG+HumCHG))
                    HCH_mrate = '{0:.2%}'.format((HmCHH+HmCHG)/(HmCHH+HumCHH+HmCHG+HumCHG))
                else:
                    continue

            if (anaseq != 0):
                uniqmaprate = '{0:.2%}'.format(uniqmap/anaseq)
                nomaprate = '{0:.2%}'.format(nomap/anaseq)
                multimaprate = '{0:.2%}'.format(multimap/anaseq)
                discardrate = '{0:.2%}'.format(discardseq/anaseq)
                duprate = '{0:.2%}'.format(duppairs/(anaseq-nomap))

            if (mCpG != 0 and umCpG != 0):
                mCpGrate = '{0:.2%}'.format(mCpG/(mCpG+umCpG))

            if (mCHG != 0 and umCHG != 0):
                mCHGrate = '{0:.2%}'.format(mCHG/(mCHG+umCHG))

            if (mCHH != 0 and umCHH != 0):
                mCHHrate = '{0:.2%}'.format(mCHH/(mCHH+umCHH))

            if (mCunknown != 0 and umCunknown != 0):
                mCunknownrate = '{0:.2%}'.format(mCunknown/(mCunknown+umCunknown))

            if (dedupmCpG != 0 and dedupumCpG != 0):
                dedupmCpGrate = '{0:.2%}'.format(dedupmCpG/(dedupmCpG+dedupumCpG))

            if (dedupmCHG != 0 and dedupumCHG != 0):
                dedupmCHGrate = '{0:.2%}'.format(dedupmCHG/(dedupmCHG+dedupumCHG))

            if (dedupmCHH != 0 and dedupumCHH != 0):
                dedupmCHHrate = '{0:.2%}'.format(dedupmCHH/(dedupmCHH+dedupumCHH))

            anaseq = '{:,}'.format(anaseq)
            uniqmap = '{:,}'.format(uniqmap)
            topmap = '{:,}'.format(topmap)
            comptopmap = '{:,}'.format(comptopmap)
            compbotmap = '{:,}'.format(compbotmap)
            botmap = '{:,}'.format(botmap)
            nomeCpGnum = '{:,}'.format(nomeCpGnum)
            nomeGpCnum = '{:,}'.format(nomeGpCnum)
            filterfrag = '{:,}'.format(filterfrag)

            #(rawreadpairs, anaseq, uniqmap, discardseq, topmap, comptopmap, compbotmap, botmap, filterfrag, nomeCpGnum, nomeGpCnum) = tuple(map(lambda x: int(x.replace(',', '')) if isinstance(x, str) else x, (rawreadpairs, anaseq, uniqmap, discardseq, topmap, comptopmap, compbotmap, botmap, filterfrag, nomeCpGnum, nomeGpCnum)))

            metrics_ori = (sampleid, rawreadpairs, anaseq, uniqmap, uniqmaprate, nomaprate, multimaprate, discardseq, topmap, comptopmap, compbotmap, botmap, mCpGrate, mCHGrate, mCHHrate, mCunknownrate, duprate, filterfrag, dedupmCpGrate, dedupmCHGrate, dedupmCHHrate, nomeCpGnum, nomeGpCnum, GCG_mrate, GCHH_mrate, GCHG_mrate, CCG_mrate, CCHH_mrate, CCHG_mrate, HCG_mrate, HCHH_mrate, HCHG_mrate, HCH_mrate)

            metrics = tuple(
                    element.replace(',', '').replace('%', '') if isinstance(element, str) else element
                    for element in metrics_ori
                    )
            pfline = '\t'.join(['%s'] * len(metrics)) + '\n' % tuple(metrics)
            fo.write(pfline)
    
    print("Summary file generated!")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--workdir', type=str, default='.', help='Path to the root working directory')
    parser.add_argument('--qcdir', type=str, default='fastqc/fastqc_raw', help='Name of the fastqc directory')
    parser.add_argument('--mapdir', type=str, default='trim', help='Name of the mapping directory')
    parser.add_argument('--outdir', type=str, default='summary', help='Name of the output directory')
    parser.add_argument('--plotMbias', type=str, default='n', help='Plot M-bias (provide path to the script) or not (n)')
    parser.add_argument('--Rscript', type=str, default='Rscript', help='Path to Rscript program')

    args = parser.parse_args()
    workdir = args.workdir
    qcdir = args.qcdir
    mapdir = args.mapdir
    outdir = args.outdir
    plotMbias = args.plotMbias
    Rscript = args.Rscript
    run_stats(workdir, outdir, qcdir, mapdir, plotMbias, Rscript)
