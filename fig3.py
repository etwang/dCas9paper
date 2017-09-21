import os, sys, subprocess, operator
import matplotlib
matplotlib.use('Agg')
from pylab import *
import miso_utils
import seaborn as sns 
from stats import pearsonr
from adjustText import adjust_text
from numpy import random


# Compare sigmoid fits to DM05 vs DM04.
def compareTAtoXia(tibsigmoid_f, mono_f, out_f, min1=1, min2=1, txt_f='none', filterlist='none'):
    if filterlist is not 'none':
        filtered = []   
        for line in open(filterlist):
            event, gene = line.strip().split()
            filtered.append(event)
    
    eventToGene = {}
    eventToDelta1 = readTibSigmoid(tibsigmoid_f, min1)
    eventToDelta2, eventToGene = readMonotonic(mono_f, min2)

    data = []
    genes = []
    events = []
    for event in eventToDelta1:
        if event in eventToDelta2:
            if filterlist == 'none' or event in filtered:
                psi11, psi12, delta1 = eventToDelta1[event]     #dmpsi, ctpsi, dm-ct
                psi21, psi22, delta2 = eventToDelta2[event]     # 
                data.append([psi11, psi12, psi21, psi22, delta1, delta2])
                events.append(event)
                if event in eventToGene:
                    genes.append(eventToGene[event])
                else:
                    genes.append('n/a')
    data = array(data)
    print data.shape

    figure(figsize=(2.2, 2.1))
    sns.set_style('white')

    scatter(data[:, 4], data[:, 5], s=4, color='r', alpha=.8)
    cc = corrcoef(data)
    print cc 
    if txt_f is not False:
        txt = open(txt_f, 'w')
    texts = []
    for x in range(data.shape[0]):
        if txt_f is not False:
            if data[x, 4] * data[x, 5] > 0:
                txt.write(events[x] + "\t" + genes[x] + "\n")
        if abs(data[x, 4]) > .25 and abs(data[x, 5]) > .25 and data[x, 4] * data[x, 5] > 0:
            texts.append(text(data[x, 4], data[x, 5], genes[x], fontsize=6, color='k',\
                zorder=1))
    if txt_f is not False:
        txt.close()
    cc, p = pearsonr(data[:, 4].tolist(), data[:, 5].tolist())
    text(1, -.9, '$Pearson R$: %.2f'%(cc), fontsize=6, ha='right')

    xlabel('$\Delta\Psi$, Tibialis Anterior\n(DM1 - Unaffected)', fontsize=8, ma='center')
    ylabel('$\Delta\Psi$, Primary myoblasts\n(DM1 - Unaffected)', fontsize=8, ma='center')
    xlim(-1, 1)
    ylim(-1, 1)
    axhline(y=0, linestyle='--', color='#CCCCCC', dashes=(4, 2), lw=.5)
    axvline(x=0, linestyle='--', color='#CCCCCC', dashes=(4, 2), lw=.5)
    xticks(linspace(-1, 1, 5), fontsize=8)
    yticks(linspace(-1, 1, 5), fontsize=8)

    sns.despine()
    adjust_text(texts, data[:, 4], data[:, 5], arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    subplots_adjust(left=.3, bottom=.2)
    savefig(out_f)


# Compare sigmoid fits of TA to rescue.
def compareTAtoBF(tibsigmoid_f, bf_f, out_f, min1=1, min2=1, txt_f='none', filterlist='none'):

    if filterlist is not 'none':
        filtered = []   
        for line in open(filterlist):
            event, gene = line.strip().split()
            filtered.append(event)
 
    eventToGene = {}
    eventToDelta1 = readTibSigmoid(tibsigmoid_f, min1)
    eventToDelta2, eventToGene = readBF(bf_f, min2)

    data = []
    genes = []
    events = []
    for event in eventToDelta1:
        if event in eventToDelta2:
            if filterlist == 'none' or event in filtered:
                psi11, psi12, delta1 = eventToDelta1[event]     #dmpsi, ctpsi, dm-ct
                psi21, psi22, delta2, bf = eventToDelta2[event] #CAG, empty, psi1-psi2 
                if abs(psi11 - psi22) < .33:
                    data.append([psi11, psi12, psi21, psi22, delta1, delta2, bf])
                    events.append(event)
                    if event in eventToGene:
                        genes.append(eventToGene[event])
                    else:
                        genes.append('n/a')
    data = array(data)
    print data.shape

    figure(figsize=(2.2, 2.1))
    sns.set_style('white')
    scatter(data[:, 4], data[:, 5], s=4, color='g', alpha=.8)
    cc, p = pearsonr(data[:, 4].tolist(), data[:, 5].tolist())
    print cc, p
    if txt_f is not 'none':
        txt = open(txt_f, 'w')
        txt.write("#Event\tGene\tTAdpsi\tdpsi\tBF\n")
    texts = []
    for x in range(data.shape[0]):
        print genes[x], data[x, :]
        if abs(data[x, 4]) >= .25 and abs(data[x, 5]) > .25 and data[x, 4] * data[x, 5] > 0 and \
            data[x, 6] > 5:
            texts.append(text(data[x, 4], data[x, 5], genes[x], fontsize=6, color='k',\
                zorder=1))
        if txt_f is not 'none':
            #if data[x, 4] * data[x, 5] > 0:
            txt.write("\t".join(map(str, [events[x], genes[x], data[x, 4], data[x, 5], data[x, 6]])) + "\n")
    if txt_f is not 'none':
        txt.close()

    #adjust_text(texts, force_text=0.05)
    xlabel('$\Delta\Psi$, Tibialis Anterior\n(DM1 - Unaffected)', fontsize=8, ma='center')
    ylabel('$\Delta\Psi$, DM1 myoblasts\n(dSaCas9$_{empty}$ - dSaCas9$_{CAG}$)',\
        fontsize=8, ma='center')
    xlim(-1, 1)
    ylim(-1, 1)
    xticks(linspace(-1, 1, 5), fontsize=8)
    yticks(linspace(-1, 1, 5), fontsize=8)
    cc, p = pearsonr(data[:, 4].tolist(), data[:, 5].tolist())
    text(1, -.9, '$Pearson R$: %.2f'%(cc), fontsize=6, ha='right')

    axhline(y=0, linestyle='--', color='#CCCCCC', dashes=(4, 2), lw=.5)
    axvline(x=0, linestyle='--', color='#CCCCCC', dashes=(4, 2), lw=.5)

    sns.despine()
    adjust_text(texts, data[:, 4], data[:, 5], arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    subplots_adjust(left=.3, bottom=.2)
    savefig(out_f)


# Compare sigmoid fits of TA to rescue - uses new data from 73017.
def compareTAtoMonotonic_73017(tibsigmoid_f, mono_f, out_f, min1=1, min2=1, txt_f='none', filterlist='none'):

    if filterlist is not 'none':
        filtered = []   
        for line in open(filterlist):
            event, gene = line.strip().split()
            filtered.append(event)
 
    eventToGene = {}
    eventToDelta1 = readTibSigmoid(tibsigmoid_f, min1)
    eventToMono, eventToGene = readMonotonic(mono_f, min2)

    data = []
    genes = []
    events = []
    for event in eventToDelta1:
        if event in eventToMono:
            if filterlist == 'none' or event in filtered:
                psi11, psi12, delta1 = eventToDelta1[event]     #dmpsi, ctpsi, dm-ct
                psi21, psi22, delta2 = eventToMono[event]       #empty, CAG, psi1-psi2  
                #if abs(psi11 - psi21) < .33:        # only compare those events with similar baseline in DM
                data.append([psi11, psi12, psi22, psi21, delta1, -delta2])
                events.append(event)
                if event in eventToGene:
                    genes.append(eventToGene[event])
                else:
                    genes.append('n/a')
    data = array(data)
    print data.shape

    figure(figsize=(2.2, 2.1))
    sns.set_style('white')
    scatter(data[:, 4], data[:, 5], s=4, color='g', alpha=.8)
    cc, p = pearsonr(data[:, 4].tolist(), data[:, 5].tolist())
    print cc, p
    if txt_f is not 'none':
        txt = open(txt_f, 'w')
        txt.write("#Event\tGene\tTAdpsi\tdpsi\tBF\n")
    texts = []
    for x in range(data.shape[0]):
        print genes[x], data[x, :]
        if abs(data[x, 4]) >= .25 and abs(data[x, 5]) > .25 and data[x, 4] * data[x, 5] > 0:
            texts.append(text(data[x, 4], data[x, 5], genes[x], fontsize=6, color='k',\
                zorder=1))
        if txt_f is not 'none':
            #if data[x, 4] * data[x, 5] > 0:
            txt.write("\t".join(map(str, [events[x], genes[x], data[x, 4], data[x, 5]])) + "\n")
    if txt_f is not 'none':
        txt.close()

    #adjust_text(texts, force_text=0.05)
    xlabel('$\Delta\Psi$, Tibialis Anterior\n(DM1 - Unaffected)', fontsize=8, ma='center')
    ylabel('$\Delta\Psi$, DM1 myoblasts\n(dSaCas9$_{empty}$ - dSaCas9$_{CAG}$)',\
        fontsize=8, ma='center')
    xlim(-.8, .8)
    ylim(-.8, .8)
    xticks(linspace(-.8, .8, 5), fontsize=8)
    yticks(linspace(-.8, .8, 5), fontsize=8)
    cc, p = pearsonr(data[:, 4].tolist(), data[:, 5].tolist())
    text(1, -.9, '$Pearson R$: %.2f'%(cc), fontsize=6, ha='right')

    axhline(y=0, linestyle='--', color='#CCCCCC', dashes=(4, 2), lw=.5)
    axvline(x=0, linestyle='--', color='#CCCCCC', dashes=(4, 2), lw=.5)

    sns.despine()
    adjust_text(texts, data[:, 4], data[:, 5], arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    subplots_adjust(left=.3, bottom=.2)
    savefig(out_f)


# Compare sigmoid fits to DM05 vs DM04.
def TAheatmap(tibpsi_f, tibgroups_f, sigmoid_f, bf_f, out_f, filterlist):

    filtered = []   
    for line in open(filterlist):
        event, gene = line.strip().split()
        filtered.append(event)
   
    eventToEC50 = {} 
    for line in open(sigmoid_f):
        if not line.startswith("#"):
            vals = line.strip().split("\t")
            event, gene, ec50, slope, ec50err, slopeerr, fit, ctpsi, dmpsi =\
                vals
            dmpsi = float(dmpsi)
            ctpsi = float(ctpsi)
            eventToEC50[event] = float(ec50)
 
    eventToGene = {}
    eventToPsi, mbnls = readTibPsi(tibpsi_f, tibgroups_f)
    eventToDelta2, eventToGene = readBF(bf_f, '0,0')

    data = []
    genes = []
    events = []
    for event in filtered:
        if event in eventToDelta2 and event in eventToPsi:
            tapsi, ctlpsi = eventToPsi[event]
            ordered = zip(tapsi, mbnls)
            ordered.sort(key=operator.itemgetter(1))
            tapsi, orderedmbnls = zip(*ordered)
            tapsi = array(tapsi)
            psi2, psi1, delta = eventToDelta2[event]
            if abs(mean(ctlpsi) - psi1) < .3:
                if mean(ctlpsi) - mean(tapsi[:10]) < 0:
                    tapsi = 1 - tapsi
                    psi1 = 1 - psi1
                    psi2 = 1 - psi2

                normed = (tapsi - min(tapsi)) / (max(tapsi) - min(tapsi))
                npsi1 = (psi1 - min(tapsi)) / (max(tapsi) - min(tapsi))
                npsi2 = (psi2 - min(tapsi)) / (max(tapsi) - min(tapsi))

                #print eventToGene[event], mean(ctlpsi), psi2, tapsi, mbnls
                data.append([normed, eventToEC50[event], npsi1, npsi2]) 
                gene = eventToGene[event]
                if gene.endswith(","):
                    gene = gene[:-1]
                genes.append(gene)

    data.sort(key=operator.itemgetter(1))
    mat = [x[0] for x in data] 
    mat = array(mat)
    start = [x[2] * mat.shape[1] for x in data]
    end = [x[3] * mat.shape[1] for x in data]
    
    print mat.shape
    pcolor(mat, cmap=get_cmap('YlOrRd'))

    for i in range(len(start)):
        arrow(start[i], i + .5, end[i] - start[i], 0, fc="k", ec="k",\
            head_width=0.25, head_length=0.2)
        print start[i], end[i]

    yticks(arange(mat.shape[0]) + .5, genes, fontsize=8)
    xlim(0, mat.shape[1])
    ylim(0, mat.shape[0])
    savefig(out_f) 


def plotSigmoid(tibpsi_f, groups_f, tibsigmoid_f, consolidated_f, event, out_f):
  
    cmd = 'head -1 ' + consolidated_f 
    pid = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    header = pid.stdout.readline().split("\t")[1:]
    idx = [x for x in range(len(header)) if header[x].endswith("_mean")]
    samples = [header[x][:-5] for x in idx]
            
    cmd = 'grep ' + event + ' ' + consolidated_f 
    pid = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    psi = map(float, pid.stdout.readline().split("\t")[1:-3])
    psi = [psi[x] for x in idx]
    
    cmd = 'grep ' + event + ' ' + tibsigmoid_f
    pid = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    siginfo = pid.stdout.readline().split("\t")
    ec50, slope, psimin, psimax, ec50err, slopeerr, psiminerr, psimaxerr  = map(float, siginfo[2:10])
    
    ctsamples = [] 

    for line in open(groups_f):
        sample, status, mbnl = line.strip().split()
        if status != 'DM1_Tibialis':
            ctsamples.append(sample)
    
    cmd = 'head -1 ' + tibpsi_f 
    pid = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    tibheader = pid.stdout.readline().split("\t")[1:]
    ctidx = array([x for x in range(len(tibheader)) if tibheader[x] in ctsamples])
    dmidx = array([x for x in range(len(tibheader)) if tibheader[x] not in ctsamples])

    eventToPsi, mbnls = readTibPsi(tibpsi_f, groups_f)
    tapsi, ctlpsi = eventToPsi[event.strip()]

    # Plot sigmoid fit
    x = linspace(0, 1, 20)
    y = psimin + (psimax - psimin) / (1 + exp(-slope * (x - ec50))) 

    # Plot 95% conf interval
    psiminvec = random.normal(psimin, psiminerr, 100)
    psimaxvec = random.normal(psimax, psimaxerr, 100)
    ec50vec = random.normal(ec50, ec50err, 100)
    slopevec = random.normal(slope, slopeerr, 100)
 
    print len(mbnls), len(tapsi) 
    sns.set_style("white")
    figure(figsize=(3, 2))

    scatter(mbnls[ctidx], tapsi[ctidx], s=10, color='g', edgecolor='k', zorder=3)
    scatter(mbnls[dmidx], tapsi[dmidx], s=10, color='k', zorder=3)
    #scatter(mbnls[ctidx], tapsi[ctidx], s=7, color='#5F7EC6', zorder=3)
    #scatter(mbnls[dmidx], tapsi[dmidx], s=7, color='#E97E25', zorder=3)

    plot(x, y, 'k--', dashes=(4, 2), lw=1)
    yvec = []
    for i in range(100):
        yvec.append(psiminvec[i] + (psimaxvec[i] - psiminvec[i]) / (1 + exp(-slopevec[i] * (x - ec50vec[i]))))

    yvec = array(yvec)
    yerrlo = []
    yerrhi = []
    for i in range(yvec.shape[1]):
        yerrlo.append(percentile(yvec[:, i], 5))
        yerrhi.append(percentile(yvec[:, i], 95))
    fill_between(x, yerrlo, yerrhi, color='#EEEEEE')

    xlim(0, 1)
    ylim(0, 1)
    xticks(linspace(0, 1, 5), fontsize=8)
    yticks(linspace(0, 1, 5), fontsize=8)
    ylabel('$\Psi$, MBNL1 exon 5', fontsize=8)
    xlabel('[MBNL]$_{inferred}$', fontsize=8)

    ctlpsi = []
    dmpsi = []
    for i in range(len(samples)):
        if samples[i].startswith("DM04"):
            print psi[i]
            ctlpsi.append(psi[i])
            axhline(psi[i], color='#264A9C', lw=.5, zorder=8)
        else:
            print psi[i]
            dmpsi.append(psi[i])
            axhline(psi[i], color='#FFAE6B', lw=.5, zorder=8)
    fill_between([0, 1], min(ctlpsi), max(ctlpsi), color='#264A9C', alpha=.2)
    fill_between([0, 1], min(dmpsi), max(dmpsi), color='#FFAE6B', alpha=.2)
    

    sns.despine()
    subplots_adjust(left=.2, bottom=.2)
    savefig(out_f)


# Read a monotonic file. Return dictionary of event->delta psi for
# events where maxz >= minz.
def readMonotonic(monotonic_f, minvals):

    minz, mind = map(float, minvals.split(","))
    eventToMono = {}
    eventToGene = {}
    for line in open(monotonic_f):
        if not line.startswith("#"):
            event, gene, trueval, m, s, psi1, psi2, delta, z, maxz = \
                line.strip().split("\t")
            if "," not in psi1:
                if float(maxz) >= minz and abs(float(delta)) >= mind:
                    eventToMono[event] = [float(psi2), float(psi1), float(delta)]
                    if gene.endswith(","):
                        gene = gene[:-1]
                    eventToGene[event] = gene
    return eventToMono, eventToGene
   
def readTibSigmoid(tibsigmoid_f, maxfit):
  
    maxfit = float(maxfit) 
    eventToDelta = {}
    for line in open(tibsigmoid_f):
        if not line.startswith("#"):
            vals = line.strip().split("\t")
            event, gene, ec50, slope, psimin, psimax,\
                ec50err, slopeerr, psiminerr, psimaxerr, fit, ctpsi, dmpsi =\
                vals
            dmpsi = float(dmpsi)
            ctpsi = float(ctpsi)
            if float(fit) <= maxfit: 
                eventToDelta[event] = [dmpsi, ctpsi, dmpsi - ctpsi]
    return eventToDelta

def readBF(bf_f, minvals):
  
    minbf, mind = map(float, minvals.split(","))
    eventToSummary, header = miso_utils.getSummary(bf_f)
    eventToDelta = {}
    eventToGene = {}
    for event in eventToSummary:
        bf = float(eventToSummary[event]['bayes factor'][0])
        psi1 = float(eventToSummary[event]['sample 1 mean'][0])
        psi2 = float(eventToSummary[event]['sample 2 mean'][0])
        delta = float(eventToSummary[event]['delta psi'][0])
        if bf >= minbf and abs(float(delta)) >= mind:
            eventToDelta[event] = [psi1, psi2, -delta, bf]
            gene = eventToSummary[event]['symb']
            if gene.endswith(","):
                gene = gene[:-1]
            eventToGene[event] = gene
    return eventToDelta, eventToGene


def readTibPsi(tibpsi_f, groups_f):
     
    dmsamples = [] 
    ctsamples = [] 
    sampleToMBNL = {}
    for line in open(groups_f):
        sample, status, mbnl = line.strip().split()
        mbnl = float(mbnl)
        sampleToMBNL[sample] = mbnl
        if status == 'DM1_Tibialis':
            dmsamples.append(sample)
        else:
            ctsamples.append(sample)
   
    eventToPsi = {}      
    for line in open(tibpsi_f):
        if line.startswith("#"):
            header = line.strip().split("\t")[1:]
            mbnls = array([sampleToMBNL[x] for x in header])
        else:
            vals = line.strip().split("\t")
            event = vals[0]
            try:
                psi = array(map(float, vals[1:]))
                ctlpsi = [psi[x] for x in ctidx]
                eventToPsi[event] = [psi, ctlpsi]
            except:
                pass
    return eventToPsi, mbnls 

  
# Fit and plot sigmoid curves.
def plotSigmoidsFull(consolidated_f, order_f, xiamonotonic_f, fits_f, out_f):

    sns.set_style('white')
    sampleToMBNL = {}
    sampleToDisease = {}
    for line in open(order_f):
        sample, group, mbnl = line.strip().split()
        sampleToMBNL[sample] = float(mbnl)
        sampleToDisease[sample] = group

    eventToMono, eventToGene = readMonotonic(xiamonotonic_f, '1,.1')

    eventToPsi = {}
    eventToGene = {}
    mbnls = []
    statuses = []
    for line in open(consolidated_f):
        if line.startswith("#"):
            header = line.strip().split("\t")

            # Get the mean psi value indexes
            meanidx = [x for x in range(len(header)) if header[x].endswith("_mean") \
                and header[x][:-5] in sampleToMBNL]

            # Match up the MBNL values
            for x in meanidx:
                mbnls.append(sampleToMBNL[header[x][:-5]])
                statuses.append(sampleToDisease[header[x][:-5]])
            ctidx = array([x for x in range(len(meanidx)) if statuses[x] == 'CT_Tibialis'])
            dmidx = array([x for x in range(len(meanidx)) if statuses[x] != 'CT_Tibialis'])
        else:
            vals = line.strip().split("\t")
            event = vals[0]
            psi = [vals[x] for x in meanidx]
            gene = vals[-2]
            if ',' in gene:
                gene = gene.split(",")[0]

            # two isoform cases
            if ',' not in "\t".join(psi):
                nonNA = [x for x in range(len(psi)) if psi[x] != 'n/a']
                goodpsi = [float(psi[x]) for x in nonNA]
                goodmbnl = [mbnls[x] for x in nonNA]
                controlpsi = [float(psi[x]) for x in nonNA \
                    if statuses[x] == 'CT_Tibialis']
                eventToPsi[event] = [goodpsi, goodmbnl, controlpsi]
                eventToGene[event] = gene

            # multi isoform cases
            else: 
                niso = len(psi[0].split(","))
                multipsi = [x.split(",") for x in psi]
                for i in range(niso):
                    nonNA = [x for x in range(len(multipsi)) if \
                        multipsi[x][0] != 'n/a']
                    goodpsi = [float(multipsi[x][i]) for x in nonNA]
                    goodmbnl = [mbnls[x] for x in nonNA]
                    controlpsi = [float(multipsi[x][i]) for x in nonNA \
                        if statuses[x] == 'CT_Tibialis']
                    eventToPsi[event + "." + str(i)] = [goodpsi, goodmbnl, controlpsi]
                    eventToGene[event + "." + str(i)] = gene
    print len(eventToPsi), 'events'
  
    events = [[e, eventToGene[e]] for e in eventToPsi]
    events.sort(key=operator.itemgetter(1))
    events = [e[0] for e in events]

    from matplotlib.backends.backend_pdf import PdfPages

    fig = figure(figsize=(8.5, 11))
    n = 0 
    pdf = PdfPages(out_f)
    nrows = 7
    ncols = 6
    for line in open(fits_f):
        if not line.startswith("#"):
            vals = line.strip().split("\t")
            event, gene = vals[:2]
            ec50, slope, psimin, psimax, ec50err, slopeerr, psiminerr, psimaxerr,\
                fit, ctpsi, dmpsi = map(float, vals[2:])
            if fit < 1.3:
                psi, mbnls, controlpsi = eventToPsi[event] 
                psi = array(psi)
                mbnls = array(mbnls)

                if n == nrows * ncols:
                    sns.despine()
                    suptitle('Figure S4', fontsize=10, ha='left')
                    subplots_adjust(wspace=.5, hspace=1.6)
                    pdf.savefig(fig)
                    close()
                    fig = figure(figsize=(8.5, 11))
                    n = 0
                
                ax = subplot(nrows, ncols, n + 1)

                # Plot points
                scatter(mbnls, psi, s=1, color='k')

                # Plot fit
                x = linspace(0, 1, 20)
                plot(x, sigmoidfull(x, ec50, slope, psimin, psimax),\
                    lw=1, color='#CCCCCC', zorder=-1)

                if event in eventToMono:
                    psi1, psi2, dpsi = eventToMono[event]
                    if (mean(psi) - mean(controlpsi)) * dpsi > 0:
                        axhline(eventToMono[event][0], color='#FF0000', alpha=.8, lw=.5)
                        axhline(eventToMono[event][1], color='#264A9C', alpha=.8, lw=.5)

                # Plot error
                #psiminvec = random.normal(psimin, psiminerr, 100)
                #psimaxvec = random.normal(psimax, psimaxerr, 100)
                #ec50vec = random.normal(ec50, ec50err, 100)
                #slopevec = random.normal(slope, slopeerr, 100)
                #yvec = []
                #for j in range(100):
                #    yvec.append(psiminvec[j] + (psimaxvec[j] - psiminvec[j]) / (1 + exp(-slopevec[j] * (x - ec50vec[j]))))
                #yvec = array(yvec)
                #yerrlo = []
                #yerrhi = []
                #for j in range(yvec.shape[1]):
                #    yerrlo.append(percentile(yvec[:, j], 5))
                #    yerrhi.append(percentile(yvec[:, j], 95))
                #fill_between(x, yerrlo, yerrhi, color='#EEEEEE')

                xlim(-.1, 1.1)
                ylim(-.1, 1.1)
                if n >= ncols * (nrows - 1):
                    xlabel('[MBNL]$_{inferred}$', fontsize=6)
                if n%ncols == 0:
                    ylabel('$\Psi$', fontsize=6)
                xticks(linspace(0, 1, 3), fontsize=6)
                yticks(linspace(0, 1, 5), fontsize=6)
                title(gene, fontsize=6)

                ttl = ax.title
                ttl.set_position([.5, 1.8])

                exons = event.split("@")
                text(.5, 2.1, "\n".join(exons), va='top', fontsize=4, ha='center')
                text(.5, 1.5, 'EC50: %.1f, Slope: %.1f\n$\Psi_{min}$: %.2f, $\Psi_{max}$: %.2f, Fit: %.3f'%\
                    (ec50, slope, psimin, psimax, fit),\
                    va='top', ha='center', fontsize=4)
                print n
                n += 1

    sns.despine()
    subplots_adjust(wspace=.5, hspace=1.6)
    suptitle('Data S1', fontsize=10, ha='left')
    pdf.savefig(fig)
    close()
    pdf.close()


# Sigmoid helper function
def sigmoidfull(x, x0, k, psimin, psimax):
     y = psimin + (psimax - psimin) / (1 + np.exp(-k*(x-x0)))
     return y

## Fig 3C
#plotSigmoid('/home/eric.t.wang/scratch/DMseq/miso/tib_psi.txt',\
#    '/home/eric.t.wang/scratch/DMseq/miso/groups_DM1',\
#    '/home/eric.t.wang/scratch/DMseq/miso/sigmoidsfull.txt',\
#    '/home/eric.t.wang/scratch/dCas9/miso/nonUTRevents_consolidated_group1.txt',
#    'chr3:152163071:152163328:+@chr3:152164493:152164546:+@chr3:152165409:152165562:+',\
#    '/home/eric.t.wang/scratch/dCas9/miso/mbnl1_sigmoid.pdf')
## Fig 3D
#compareTAtoXia('/home/eric.t.wang/scratch/DMseq/miso/sigmoidsfull.txt',\
#    '/home/eric.t.wang/scratch/dCas9/miso/nonUTRevents_monotonic_group1.txt',\
#    '/home/eric.t.wang/scratch/dCas9/miso/TAsigmoidfull.vs.group1.pdf',\
#    1.3, '1,.1',\
#    '/home/eric.t.wang/scratch/dCas9/miso/TAsigmoidfull.vs.group1.overlap_1.3_1')
#
## Fig 3F
#compareTAtoBF('/home/eric.t.wang/scratch/DMseq/miso/sigmoidsfull.txt',\
#    '/home/eric.t.wang/scratch/dCas9/miso/nonUTRevents_summaries/DM05_diff_CAG_vs_DM05_diff_empty.miso_bf',\
#    '/home/eric.t.wang/scratch/dCas9/miso/original.pdf',\
#    1.3, '0,0',\
#    '/home/eric.t.wang/scratch/dCas9/miso/TAsigmoidfull.vs.CAGvsEmpty.events.txt',\
#    '/home/eric.t.wang/scratch/dCas9/miso/TAsigmoidfull.vs.group1.overlap_1.3_1')
#
# Data S1 
#plotSigmoidsFull('/home/eric.t.wang//scratch/DMseq/miso/tibialis_consolidated_nonUTR_summaries.txt',\
#    '/home/eric.t.wang/scratch/DMseq/miso/groups_DM1',\
#    '/home/eric.t.wang/scratch/dCas9/miso/nonUTRevents_monotonic_group1.txt',\
#    '/home/eric.t.wang/scratch/DMseq/miso/sigmoidsfull.txt',\
#    '/home/eric.t.wang/scratch/dCas9/miso/DataS1.pdf')



