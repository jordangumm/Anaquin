import re
import os
import sys
import math
import copy
import json
import datetime
import tempfile
import functools
import subprocess
import collections
from operator import add

D_AXIS_SIZE  = 8
D_TITLE_SIZE = 8

def run(cmd):
    out = subprocess.check_output(cmd, shell=True)
    return out.decode("utf-8") 

def runR(file, out, width='NA', height='NA', o = {}, panelS=6):
    if not os.path.exists(file):
        return None
    name = tempfile.NamedTemporaryFile().name
    with open(name, 'w') as tmp:
        envs  = "Sys.setenv(family='Helvetica')\n"
        envs += "Sys.setenv(report='Anaquin')\n"
        for i in o:
            envs = envs + str("\nSys.setenv({0}={1})".format(i, o[i]))
        if out is not None:
            envs = envs + '\nSys.setenv(panelF="' + os.path.abspath(out) + '")'        
        if panelS is not None:
            envs = envs + '\nSys.setenv(panelS="' + str(panelS) + '")'
            cmd = str("{0}\npdf(NULL)\nsource('{1}')\n".format(envs, os.path.abspath(file)))
        else:
            cmd = str("{0}\npdf(NULL)\nsource('{1}')\npdf(NULL)\nggsave('{2}', width={3}, height={4})\n".format(envs, os.path.abspath(file), os.path.abspath(out), width, height))
        tmp.write(cmd)

    try:
        r = run("Rscript --vanilla " + name)
    except subprocess.CalledProcessError:
        return None
    return r

def isFloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False

#
# Convert a TSV file into HTML
#

def tsv2HTML(file, keys, names=None, prec=None, showHead=True):
    x = " <thead><tr>\n"
    for i in range(len(keys)):
        key = keys[i] if names is None else names[i]
        if showHead:
            x = x + "<th class='tHead'>" + key + "</th>\n"
    x += "</tr></thead>\n"
    tsv = readTSV(file)
    for i in range(len(tsv[keys[0]])):
        x += "<tbody><tr>\n"
        for j in range(len(keys)):
            key = keys[j]
            v = tsv[key][i]
            if isFloat(v):
                if prec is None:
                    x += "<td class='tCell'>" + toFloat(v, 2) + "</td>\n"
                else:
                    x += "<td class='tCell'>" + toFloat(v, prec[j]) + "</td>\n"
            else:
                x += "<td class='tCell'>" + str(v) + "</td>\n"            
        x += "</tr></tbody>\n"
    return x

def RLadTable(file, cols):
    s = tempfile.NamedTemporaryFile(delete=False).name
    o = tempfile.NamedTemporaryFile(delete=False).name        
    with open(s, "w") as w:
        w.write("x <- read.csv('" + file + "', sep='\\t')\n")
        w.write("x <- x[grepl('_LD_', x$Name),]; x$Copy <- NaN\n") # Only synthetic ladders
        w.write("x[grepl('_A', x$Name),]$Copy <- 1; x[grepl('_B', x$Name),]$Copy <- 2; x[grepl('_C', x$Name),]$Copy <- 4; x[grepl('_D', x$Name),]$Copy <- 8\n")
        w.write("stopifnot(sum(is.nan(x$Copy)) == 0);x <- x[,c('Copy','Med')]\n")
        w.write("x <- x[x$Med != '-',]; x$Med <- as.numeric(as.character(x$Med))\n")
        w.write("m <- aggregate(.~Copy, x, mean); s <- aggregate(.~Copy, x, sd); c <- 100*(s$Med/m$Med)\n")
        w.write("q0 <- aggregate(.~Copy, x, min); q25 <- aggregate(.~Copy, x, function(x) { quantile(x, 0.25) });q50 <- aggregate(.~Copy, x, function(x) { quantile(x, 0.50) });q75 <- aggregate(.~Copy, x, function(x) { quantile(x, 0.75) });q100 <- aggregate(.~Copy, x, max)\n")
        w.write("x <- data.frame(Copy=m$Copy, Mean=m$Med, SD=s$Med, CV=c, Q0=q0$Med, Q25=q25$Med, Q50=q50$Med, Q75=q75$Med, Q100=q100$Med, Ratio=NaN);\n")
        w.write("f = function(i,j,x) { if (i %in% x$Copy && j %in% x$Copy) { x[x$Copy==j,]$Ratio <- x[x$Copy==j,]$Q50 / x[x$Copy==i,]$Q50; x } }\n")
        w.write("x <- f(4,8,x);x <- f(2,4,x);x <- f(1,2,x)\n")
        w.write("x$Copy <- paste(x$Copy, 'cn', sep='')\n")
        w.write("x <- x[, c(" + str(cols).replace("[", "").replace("]", "") + ")]\n")
        w.write("write.table(x, file='" + o + "', sep='\\t', quote=F, row.names=F);print(mean(x$Ratio, na.rm=TRUE))")
    return (s, o)

#
# Sort a data file in R
#

def RSort(file, sort, cols):
    s = tempfile.NamedTemporaryFile(delete=False).name
    o = tempfile.NamedTemporaryFile(delete=False).name        
    with open(s, "w") as w:
        w.write("x <- read.table('" + file  + "', sep='\\t', header=T)\n")
        w.write("x <- x[,colnames(x) %in% " + cols + "]\n")
        w.write("x <- x[with(x, order(" + sort + ")),]\n")
        w.write("write.table(x,'" + o + "',row.names=F,quote=F,sep='\\t')")
    return (s, o)

#
# Filter for certain columns
#

def RFilter(file, col):
    s = tempfile.NamedTemporaryFile(delete=False).name
    o = tempfile.NamedTemporaryFile(delete=False).name        
    with open(s, "w") as w:
        w.write("x <- read.table('" + file  + "', sep='\\t', header=T)\n")
        w.write("x <- x[x$" + col + "!='-',]\n")
        w.write("write.table(x,'" + o + "',row.names=F,quote=F,sep='\\t')")
    return (s, o)

#
# Aggregation with mean and CV in R
#

def RMeanCV(file, cols):
    s = tempfile.NamedTemporaryFile(delete=False).name
    o = tempfile.NamedTemporaryFile(delete=False).name  
    with open(s, "w") as w:
        w.write("x <- read.table('" + file  + "', sep='\\t', header=T, stringsAsFactors=F)\n")
        w.write("x <- x[,colnames(x) %in% " + cols + "]\n")
        w.write("x[x=='-'] <- NA\n")
        w.write("x <- x[complete.cases(x),]\n")
        w.write("x[] <- lapply(x, as.numeric)\n")
        w.write("m <- aggregate(.~Mix, x, mean)\n") # Mean for all sequin by mixture level
        w.write("s <- aggregate(.~Mix, x, sd)\n")   # SD for all sequin by mixture level
        w.write("stopifnot(all(m[,1] == s[,1]))\n")
        w.write("x <- data.frame(Name=m[,1], Mean=m[,2], CV=s[,2]/m[,2])\n")
        w.write("x <- x[with(x, order(Name)),]\n")
        w.write("write.table(x,'" + o + "',row.names=F,quote=F,sep='\\t')")
    return (s, o)

def isMiss(x):
    return x == '-' or x is None or x == "NA"

def median(x):
    x = [i for i in x if not isMiss(i)]    
    sortedLst = sorted(x)
    lstLen = len(x)
    index = (lstLen - 1) // 2
    if (lstLen % 2):
        return sortedLst[index]
    else:
        return (sortedLst[index] + sortedLst[index + 1])/2.0

def myPercent(data, percent):
    def f(N, percent, key=lambda x:x):
        if not N:
            return None
        k = (len(N)-1) * percent
        f = math.floor(k)
        c = math.ceil(k)
        if f == c:
            return key(N[int(k)])
        d0 = key(N[int(f)]) * (c-k)
        d1 = key(N[int(c)]) * (k-f)        
        return d0+d1
    return functools.partial(f, percent=percent/100.0)(sorted(data))

def myMean(x):
    x = [i for i in x if not isMiss(i)]
    n = len(x)
    if n < 1:
        raise ValueError('mean requires at least one data point')
    return sum(x)/n # in Python 2 use sum(data)/float(n)

def mySD(data, ddof=0):
    def _ss(data):
        c = myMean(data)
        ss = sum((x-c)**2 for x in data)
        return ss

    n = len(data)
    if n < 2:
        raise ValueError('variance requires at least two data points')
    ss = _ss(data)
    pvar = ss/(n-ddof)
    return pvar**0.5

def knownKeys():
    return [ 'Total_Total', 'Mutation_SNP', 'Mutation_Indel', "Genotype_Homozygous", "Genotype_Heterozygous", 'Genotype_Somatic' ]

def attr(x):
    return x.split('_')[0]

def grp(x):
    return x.split('_')[1]

def formatAttr(x):
    return grp(x) + ' (' + attr(x) + ')'

# Columns to ignore while dumping
def dumpIgnores():
    return ['GCcontent', 'GeneContext', 'MobileElement', 'SimpleRepeat', 'QSI', 'QSS']

def sortKeys(x):
    r1 = [ "Total", "Mutation", "Genotype", "GCcontent", "GeneContext", "SimpleRepeat", "MobileElement" ]
    r2 = [ "SNP", "Indel",                                        \
           "ATrich", "GCrich", "CodingRegion", "NoncodingRegion", \
           "DNA", "LINE", "LTR",         \
           "Mono", "Tri", "Di", "Quad",  \
           "Heterozygous", "Homozygous" ]
    
    def cmp(x, y):
        x1 = x.split("(")[1].replace(")", "")
        x2 = x.split(" ")[0]
        y1 = y.split("(")[1].replace(")", "")
        y2 = y.split(" ")[0]
        x1 = r1.index(x1) if x1 in r1 else 100
        y1 = r1.index(y1) if y1 in r1 else 100
        x2 = r2.index(x2) if x2 in r2 else 100
        y2 = r2.index(y2) if y2 in r2 else 100

        if x1 > y1:
            return 1
        elif x1 < y1:
            return -1
        elif x2 > y2:
            return 1
        elif x2 < y2:
            return -1
        
        a1 = attr(x)
        a2 = attr(y)
        if a1 != a2:
            return 1 if a1 > a2 else -1
        else:
            g1 = grp(x)
            g2 = grp(y)            
            return 1 if g1 > g2 else -1        
    x.sort(cmp=cmp)
    return x

def isAttrStr(x):
    return '_' in x and len(x.split('_')) == 2

# Convert a string to a fixed-preicison floating string representation
def toFloat(x, r=2):
    return '-' if isMiss(x) else str(('{0:.' + str(r) + 'f}').format(float(x)))

# Convert a string to float and accordingly adjust for missing values
def str2Float(x, na20=True):
    return 0 if isMiss(x) and na20 else float(x)

def tryFloat(x):
    try:
        return float(int(x))
    except Exception:
        pass    
    try:
        return float(x)
    except Exception:
        return str(x)
    
def add(x):
    return sum([float(i) for i in filter(lambda x: not isMiss(x), x)])

def med(x, r=0):
    x = [i for i in x if not isMiss(i)]

    if len(x) == 0:
        return str('-')
    else:
        x = [float(i) for i in x]        
        if r == 0:
            med = str(int(round(myPercent(x, 50), r)))
        else:
            f = '%.' + str(r) + 'f'
            med = str(f % myPercent(x, 50))

        if len(x) > 1:
            if r == 0:
                return med + ' +- ' + str(int(round(mySD(x, ddof=1), r)))                
            else:
                return med + ' +- ' + str(round(mySD(x, ddof=1), r))                
        else:
            return  med + ' +- 0'            

def parseFBED(file):
    x = {}
    with open(file) as f:
        for line in f:
            toks = line.strip().split('\t')
            if len(toks) == 0:
                continue            
            assert(len(toks) >= 4)
            
            tmp = toks[3].split('_')
            key = tmp[0] # Eg: GCContent
            grp = tmp[1] # Eg: GCrich
            key = key + '_' + grp
            
            if not key in x:
                x[key] = []
            x[key].append({ "name":toks[0], 'start':int(toks[1]), 'end':int(toks[2]), 'val':toks[3] })
    return x

def parseTSV(file, only = None, skip = None, keys = None):
    assert only is None or type(only) is dict

    # The data that we read from a TSV file
    def data():
        return { 'ql':[], 'dp':[], 'll':[], 'ef':[], 'of':[], 'nr':[], 'nv':[], 'tr':[], 'tv':[], 'l':[], 'n':[], "md":[] }

    c = {}
    x = { 'Total_Total':data(), 'Mutation_SNP':data(), 'Mutation_Indel':data(), "Genotype_Homozygous":data(), "Genotype_Heterozygous":data(), 'Genotype_Somatic':data() }

    with open(file) as f:
        for line in f:
            toks = line.strip().split('\t')
            if toks[0] == "Name":
                for i in range(0, len(toks)):
                    c[toks[i]] = int(i)
                continue

            # Stay for conditions
            if only is not None:
                if 'Label' in only and not toks[c['Label']] in only['Label']:
                    continue
                elif 'ExpFreq' in only and toks[c['ExpFreq']] != '-' and not float(toks[c['ExpFreq']]) in only['ExpFreq']:
                    continue
                if 'Genotype' in only and not toks[c['Genotype']] in only['Genotype']:
                    continue
            
            # Skip for conditions
            if skip is not None:
                b = False
                for i in skip:
                    if i in c and toks[c[i]] == skip[i]:
                        b = True
                        break
                if b:
                    continue

            if "nan" in toks:
                print("Warning: NAN detected")
                continue
     
            #
            # List of metrics, not all available
            #

            l = toks[c['Chrom']] + ':' + toks[c['Position']]

            # Median 
            md = toFloat(toks[c['Med']]) if 'Med' in c else None

            def depth():
                if "Depth" in c:
                    return toFloat(toks[c["Depth"]])
                elif "Ref.Depth" in c and "Var.Depth" in c:
                    return str2Float(toks[c["Ref.Depth"]]) + str2Float(toks[c["Var.Depth"]])
                else:
                    return None

            dp = depth()
            
            # Expected allele frequency
            ef = float(toks[c['ExpFreq']]) if 'ExpFreq' in c and toks[c['ExpFreq']] != '-' else None
            
            # Length of the sequin
            ll = toFloat(toks[c['Size']]) if 'Size' in c else None

            assert(ll != '0')
            
            def oAF():
                if 'Obs.Freq (Tumor)' in c:
                    k = 'Obs.Freq (Tumor)'
                elif 'ObsFreq' in c:
                    k = 'ObsFreq'
                else:
                    k = None
                return toFloat(toks[c[k]]) if k in c else None
            
            # Observed allele frequency
            of = oAF()

            # Quality score
            ql = toFloat(toks[c['Qual']]) if 'Qual' in c else None

            if (ql is None or ql == '-') and 'QSI' in c and toks[c['QSI']] != '-':
                ql = toFloat(toks[c['QSI']])
            elif (ql is None or ql == '-') and 'QSS' in c and toks[c['QSS']] != '-':
                ql = toFloat(toks[c['QSS']])

            nr = toFloat(toks[c['Ref.Depth (Normal)']]) if 'Ref.Depth (Normal)' in c else None
            nv = toFloat(toks[c['Var.Depth (Normal)']]) if 'Var.Depth (Normal)' in c else None
            tr = toFloat(toks[c['Ref.Depth (Tumor)']])  if 'Ref.Depth (Tumor)'  in c else None
            tv = toFloat(toks[c['Var.Depth (Tumor)']])  if 'Var.Depth (Tumor)'  in c else None

            if nr is None and 'RefCount' in c:
                nr = toFloat(toks[c['RefCount']]) # Reuse "normal" for single sample
            if nv is None and 'VarCount' in c:
                nv = toFloat(toks[c['VarCount']]) # Reuse "normal" for single sample

            if of != '-' and tr is not None and tv is not None:
                of_ = float(tv) / (float(tr) + float(tv))
                assert(abs(float(of) - float(of)) <= 0.1)

            def add(key):
                if key not in x:
                    x[key] = data()
                x[key]['n'].append(toks[c["Name"]])
                x[key]['l'].append(l)
                x[key]['dp'].append(dp)
                x[key]['ll'].append(ll)
                x[key]['ef'].append(ef)
                x[key]['of'].append(of)
                x[key]['nr'].append(nr)
                x[key]['nv'].append(nv)
                x[key]['tr'].append(tr)
                x[key]['tv'].append(tv)
                x[key]['ql'].append(ql)
                x[key]['md'].append(md)

            add('Total_Total')

            if toks[c['Type']] == 'SNP':
                add('Mutation_SNP')
            else:
                add('Mutation_Indel')

            if 'Genotype' in c and toks[c['Genotype']] == 'Homozygous':
                add("Genotype_Homozygous")
            elif 'Genotype' in c and toks[c['Genotype']] == 'Heterozygous':
                add("Genotype_Heterozygous")
            else:
                add('Genotype_Somatic')

            for key in keys:
                if toks[c[attr(key)]].split('_')[0] == grp(key):
                    add(key) # Eg: GCcontent_GCrich

    return x

def parseKM(file, only = None, skip = None, keys = None):
    assert only is None or type(only) is dict

    def data():
        return { 'cn':[], 'seq':[], 'll':[] }

    c = {}
    x = { 'Total_Total':data() }
    
    with open(file) as f:
        for line in f:
            toks = line.strip().split('\t')            
            if toks[0] == 'Sequence':
                for i in range(0, len(toks)):
                    c[toks[i]] = int(i)
                continue
            if only is not None:
                if 'Label' in only and not toks[c['Label']] in only['Label']:
                    continue
            if skip is not None:
                b = False
                for i in skip:
                    if i in c and toks[c[i]] in skip[i]:
                        b = True
                        break
                if b:
                    continue
            
            def add(key):
                if key not in x:
                    x[key] = data()
                x[key]['cn'].append(toks[c['Count']])
                x[key]['seq'].append(toks[c['Sequence']])
                x[key]['ll'].append(toks[c['Chrom']] + ':' + toks[c['Start']] + '-' + toks[c['End']])

            add('Total_Total')

            for key in keys:
                if toks[c[key]] != '-':
                    add(toks[c[key]])

    return x

def parseKTSV(data, src):
    name = data["name"]
    base = data["base"] + os.sep + name   
    repo = data["base"] + os.sep + "report_files/" + name 
    bed  = parseFBED(base + "_features.bed")
    ger = parseTSV(src, {'Label':'Germline'}, None, list(bed.keys()))
    som = parseTSV(src, {'Label':'Somatic'},  None, list(bed.keys()))

    def parseGR():
        g = {}
        keys = list(ger.keys())
        
        for i in range(0, len(keys)):
            key = keys[i]
            g[key] = { 'kn':'0', 'n':'-', 'tp':'-', 'fn':'-', 'sn':'-', 'rc':'-', 'vc':'-' }

            if key in ger:
                grp = ger[key]        
                obs = grp['of']        
                ntp = len([i for i in obs if i is not None and i != '-'])
                nfn = len(obs) - ntp

                g[key]['n']  = len(grp['of'])
                g[key]['tp'] = ntp
                g[key]['fn'] = nfn
                g[key]['md'] = med(grp['md'])
                g[key]['nr'] = med(grp['nr'])
                g[key]['nv'] = med(grp['nv'])
                
                if (ntp + nfn) != 0:
                    g[key]['sn'] = toFloat(float(ntp) / (ntp + nfn))
        return g
        
    def parseSO():    
        x = {}
        x['ntp'] = len([i for i in som['Total_Total']['of'] if i is not None and i != '-'])        
        x['nfn'] = len(som['Total_Total']['of']) - x['ntp']
        x['af']  = sorted(list(set([float(i) for i in som['Total_Total']['ef'] if i != '-' and i != None])), reverse=True)
    
        for af in x['af']:
            af  = str(af)
            grp = parseTSV(src, {'Genotype':'Somatic', 'ExpFreq':[float(af)]}, {'Label':'All'}, list(bed.keys()))
            x[af] = calcAttrPerf(grp, bed, ["Genotype_Homozygous", "Genotype_Heterozygous"])            
            x[af]['tp'] = len([i for i in grp['Total_Total']['of'] if i is not None and i != '-'])
            x[af]['fn'] = len(grp['Total_Total']['of']) - x[af]['tp']
            x[af]['tp'] = str(x[af]['tp'])
            x[af]['fn'] = str(x[af]['fn'])
            x[af]['nr'] = med(grp['Total_Total']['nr'])
            x[af]['nv'] = med(grp['Total_Total']['nv'])
            
        return x

        #
        # Each copy ladder has multiple sequin. Compute the median of those numbers.
        #
        #  TODO: The numbers should have been among all k-mers for all sequin in a copy.
        #
        
        cn = {}
        cp = []
        for i in x:
            min = toFloat(median([j["q0"]  for j in x[i]]))
            q25 = toFloat(median([j["q25"] for j in x[i]]))
            q50 = toFloat(median([j["q50"] for j in x[i]]))
            q75 = toFloat(median([j["q75"] for j in x[i]]))
            max = toFloat(median([j["q100"] for j in x[i]]))
            mu = toFloat(median([j["mu"] for j in x[i]]))
            sd = toFloat(median([j["sd"] for j in x[i]]))
            cv = toFloat(median([j["cv"] for j in x[i]]))
            cn[i] = float(q50)
            cp.append({ "cn":str(i) + "n", "cv":cv, "sd":sd, "mu":mu, "q0":min, "q25":q25, "q50":q50, "q75":q75, "q100":max })

        for i in cp:
            thisCN = int(i["cn"].replace("n", ""))
            lastCN = int(i["cn"].replace("n", "")) - 1            
            i["rt"] = "-" if not lastCN in cn else str(cn[thisCN] / cn[lastCN])

        # Mean ratio
        mr = toFloat(myMean([float(i["rt"]) for i in cp if i["rt"] != "-"]))
        
        return (cp, mr)

    return { "ger":parseGR(), "som":parseSO() }

def sStats(data):
    name = data["name"]
    base = data["base"] + os.sep + name

    stats = {}
    with open(base + '_summary.stats', 'r') as r:
        for line in r:
            line = line.strip().replace('%', '').replace('(', '').replace(')', '')
            toks = [i for i in line.split(' ') if len(i) > 0]
            
            if len(toks) > 1:
                last = toks[len(toks) - 1]            
                if 'Sequin regions' in line:
                    stats['regs'] = last
                elif 'Sequin size' in line:
                    stats['size'] = last

    return stats

def calcAttrPerf(tsv, bed, ignores=None):
    x = { 'keys':knownKeys() + list(bed.keys()) }    
    if ignores is not None:
        for i in ignores:
            x['keys'].remove(i)

    for key in x['keys']:
        def count(x):
            # Use size for counting in the len() function
            return len(x['ll'])

        null = key not in tsv or len(tsv[key]['n']) == 0
                    
        x[key] = {}
        x[key]['n']  = count(tsv[key]) if not null else '-'
        x[key]['of'] = med(tsv[key]['of'], 4) if not null else '-'

        nr = [float(i) for i in tsv[key]['nr'] if i != "-"] if not null else None
        nv = [float(i) for i in tsv[key]['nv'] if i != "-"] if not null else None

        x[key]['nd'] = med([sum(i) for i in zip(nr, nv)] if not null else '-')
        x[key]['tr'] = med(tsv[key]['tr']) if not null else '-'
        x[key]['tv'] = med(tsv[key]['tv']) if not null else '-'
        x[key]['ql'] = med(tsv[key]['ql']) if not null else '-'
   
    return x

def parseSTSV(data):
    name = data["name"]
    base = data["base"] + os.sep + name
    bed  = parseFBED(data["base"] + os.sep + 'genome_files/' + name + '_features.bed')
    tsv  = parseTSV(base + "_sequin.tsv", {'Label':['TP','FP','FN']}, None, list(bed.keys()))
    tsv.pop("Genotype_Somatic", None)

    efs = set([float(i) for i in tsv['Total_Total']['ef'] if i != '-' and i != None])    
    assert(len(efs) > 0)

    tsv['ef'] = ''
    for af in sorted(set(efs), reverse=True):
        all = parseTSV(base + "_sequin.tsv", {'Label':['TP','FN'], 'ExpFreq':[float(af)]}, None, list(bed.keys()))
        tp  = parseTSV(base + "_sequin.tsv", {'Label':['TP'], 'ExpFreq':[float(af)]}, None, list(bed.keys()))
        fn  = parseTSV(base + "_sequin.tsv", {'Label':['FN'], 'ExpFreq':[float(af)]}, None, list(bed.keys()))

        # Make sure we have something
        assert(len(all['Total_Total']['of']) > 0)

        # Only detected measurements
        obs = [i for i in all['Total_Total']['of'] if i != '-']

        af = str(af)
        tsv['ef'] += (',' if len(tsv['ef']) > 0 else '') + af

        tsv[af] = {}
        tsv[af]['tp'] = len(tp['Total_Total']['l'])
        tsv[af]['fn'] = len(fn['Total_Total']['l'])

        if len(obs) == 0:
            tsv[af]['ob'] = '-'
            tsv[af]['dn']  = '-'
            tsv[af]['dr']  = '-'
            tsv[af]['dv']  = '-'
            tsv[af]['ql']  = '-'
        else:
            # That might include non-detected sequin
            obs = all['Total_Total']['of']

            tsv[af]['ob'] = med(obs, 4)
            tsv[af]['ql'] = med(tp['Total_Total']['ql'])
            
            nr = [float(i) for i in tp['Total_Total']['nr']]
            nv = [float(i) for i in tp['Total_Total']['nv']]

            tsv[af]['dn'] = med([sum(x) for x in zip(nr, nv)])
            tsv[af]['dr'] = med([float(i) for i in tp['Total_Total']['tr']])
            tsv[af]['dv'] = med([float(i) for i in tp['Total_Total']['tv']])

    tsv['fp'] = calcAttrPerf(parseTSV(base + "_sequin.tsv", {'Label':['FP']}, \
                             None, list(bed.keys())), bed, ["Genotype_Homozygous", "Genotype_Heterozygous", "Genotype_Somatic"])
    tsv['sv'] = calcAttrPerf(parseTSV(base + "_sequin.tsv", {'Label':['SV']}, \
                             None, list(bed.keys())), bed, ["Genotype_Homozygous", "Genotype_Heterozygous", "Genotype_Somatic"])
    return tsv

#
# Parse summary statistics and copy the key-value pairs
#

def parseStats(data):
    name = data["name"]
    base = data["base"] + os.sep + name
    stats = {}
    with open(base + '_summary.stats', 'r') as r:
        for line in r:
            toks = [i.strip() for i in line.strip().split(':') if len(i) > 0]
            if len(toks) > 1:
                last = toks[len(toks) - 1]
                stats[' '.join(toks[0:len(toks)-1]).replace(':', '')] = last
            elif len(toks) == 1:
                stats[toks[0]] = ""
    return stats

def parseGTSV(data):
    name = data["name"]
    base = data["base"] + os.sep + name   
    repo = data["base"] + os.sep + 'report_files/' + name 
    geno = data["base"] + os.sep + 'genome_files/' + name 

    bed = parseFBED(geno + '_features.bed')
    tsv = parseTSV(base + "_sequin.tsv", {'Label':['TP', 'FP', 'FN']}, None, list(bed.keys()))
    gtp = parseTSV(base + "_sequin.tsv", {'Label':['TP']}, None, list(bed.keys()))
    gfp = parseTSV(base + "_sequin.tsv", {'Label':['FP']}, None, list(bed.keys()))
    gfn = parseTSV(base + "_sequin.tsv", {'Label':['FN']}, None, list(bed.keys()))
    gsm = parseTSV(base + "_sequin.tsv", {'Label':['SV']}, None, list(bed.keys()))

    def count(x):
        # Use 'll' for simply counting in the len() function
        return len(x['ll'])

    x = { 'strs':knownKeys() + list(bed.keys()) }
    x['strs'].remove('Genotype_Somatic')
    
    for i in range(0, len(x['strs'])):
        key = x['strs'][i]
        
        ll = int(add(tsv[key]['ll'])) if key in tsv else '-' # Default implementation for sequin size
        
        tp = count(gtp[key]) if key in gtp else 0
        fp = count(gfp[key]) if key in gfp else 0
        fn = count(gfn[key]) if key in gfn else 0
        sn = float(tp) / (tp + fn) if tp + fn > 0 else '-'
        pc = float(tp) / (tp + fp) if tp + fp > 0 else '-'        
        f1 = toFloat(2.0 * ((pc * sn) / (pc + sn))) if pc != '-' and sn != '-' and (sn + pc) != 0 else '-'

        x[key] = {}
        x[key]['size'] = str(ll)
        x[key]['N']   = count(tsv[key]) if key in tsv else '-'
        x[key]['DP']  = med(tsv[key]['dp']) if key in tsv else '-'
        x[key]['TP']  = tp
        x[key]['FP']  = fp
        x[key]['FN']  = fn
        x[key]['SN']  = toFloat(sn)
        x[key]['PC']  = toFloat(pc)
        x[key]['FDR'] = toFloat(fp / (float(ll) / 1000)) if fp != 0 and ll > 0 else '-'
        x[key]['TQ']  = med(gtp[key]['ql']) if key in gtp else '-'
        x[key]['S']   = count(gsm[key]) if key in gsm else '-'
        x[key]['SQ']  = med(gsm[key]['ql'] if key in gsm else '-')

    return x

def replace(data, t, x, y = None):
    if t not in data:
        raise Exception(t + ' not found')
    if isinstance(y, float):
        y = round(y, 2)
    if x is None:
        return data.replace(t, '-')
    elif y is None:
        return data.replace(t, str(x))
    return data.replace(t, str(x[y]))

def loadHTML(file, css):
    with open(css, "r") as c:    
        with open(file, 'r') as f:
            return f.read().replace('__Date__', datetime.datetime.now().strftime("%B %d, %Y %I:%M%p")).replace("__CSS__", "<style>\n" + c.read() + "\n</style>")

def writeHTML(dst, name, data, css = None):
    x = ""
    if css is not None:
        for i in css:
            x += ("." + i + " { " + css[i] + ";}\n")
    data = data.replace("%CSS%", x)
    
    file = dst + os.sep + name + '_report.html'
    with open(file, 'w') as f:
        f.write(data)
    return data

def dumpTSV(file, ignores, custom=None):
    c  = {}
    c2 = {} # indices to headers
    head = ''
    lines = []
    with open(file) as f:
        for line in f:
            toks = line.strip().split('\t')
            if "Name" in toks[0]:
                head += '<tr>'
                for i in range(0, len(toks)):
                    c[toks[i]] = int(i)
                    c2[int(i)] = toks[i]
                    if toks[i] not in ignores:
                        head += ('<th class="ghead">' + toks[i] + '</th>')
                head += '</tr>\n'
                continue
            elif toks[c['Label']] != 'TP' and toks[c['Label']] != 'FN':
                continue
            for i in ignores:
                if i in c:
                    toks[c[i]] = None
            lines.append({ 'chr':c['Chrom'], 'pos':c['Position'], 'toks':toks })

    lines = sorted(lines, key = lambda x: (x['chr'], x['pos']))

    html = ''
    for line in lines:
        html += '<tr>'
        for i in range(0,len(line['toks'])):
            if line['toks'][i] != None:
                v = line['toks'][i]                
                if custom is not None and c2[i] in custom:
                    v = toFloat(v, custom[c2[i]])
                html += ('<td>' + v + '</td>')
        html += '</tr>\n'

    return (head, html)

def dumpStats(s, html):
    for key in s:
        if ('__' + key + '__') in html:
            html = replace(html, '__' + key + '__', s[key])
    return html

def germline(data, out, html, css):
    name = data["name"]
    base = data["base"] + os.sep + name   
    repo = data["base"] + os.sep + 'report_files/' + name 
    geno = data["base"] + os.sep + 'genome_files/' + name 
    gsum = parseStats(data)
    html = dumpStats(gsum, loadHTML(html, css))

    gt   = parseGTSV(data)
    gtab = ''
    tmp  = "<tr><td class='gfirst'>__T__</td>\n\
                <td class='gcell'>__N__</td>\n\
                <td class='gcell blue'>__TP__</td>\n\
                <td class='gcell'>__FP__</td>\n\
                <td class='gcell'>__FN__</td>\n\
                <td class='gcell'>__SN__</td>\n\
                <td class='gcell'>__PC__</td>\n\
                <td class='gcell'>__FDR__</td>\n\
                <td class='gcell'>__DP__</td>\n\
                <td class='gcell green'>__S__</td></tr>"
    
    hasSamp = gsum["Input file (second)"] != '-'
    
    x = ''
    m = {}
    for i in list(gt.keys()):
        if i != "strs":
            m[formatAttr(i)] = i

    for key in sortKeys(m.keys()):
        key = m[key]
        tab = copy.copy(tmp)

        if not isAttrStr(key):
            continue
        
        tab = replace(tab, '__T__',  formatAttr(key))
        tab = replace(tab, '__N__',   gt[key]['N'])
        tab = replace(tab, '__DP__',  gt[key]['DP'])
        tab = replace(tab, '__TP__',  gt[key]['TP'])
        tab = replace(tab, '__FP__',  gt[key]['FP'])
        tab = replace(tab, '__FN__',  gt[key]['FN'])
        tab = replace(tab, '__SN__',  gt[key]['SN'])
        tab = replace(tab, '__PC__',  gt[key]['PC'])
        tab = replace(tab, '__FDR__', gt[key]['FDR'])        
        if "S" in gt[key]:
            tab = replace(tab, "__S__", gt[key]["S"])
        else:
            tab = replace(tab, "__S__", "")
        gtab += tab

    html = replace(html, '__T1__', gtab)

    (head, tab) = dumpTSV(base + "_sequin.tsv", ['Size'] + dumpIgnores(), { 'Ref.Depth':0, 'ExpFreq':1, 'Var.Depth':0, 'Qual':0, 'Obs.Freq':4 })
    html = replace(html, '__HEAD__',  head)
    html = replace(html, '__T2__', tab)
    
    writeHTML(out, name, html, { "notCombined": "display:none" } if ("Reference regions (decoy)" in gsum) else { "combined": "display:none" })

def somatic(data, out, html, css):
    name = data["name"]
    base = data["base"] + os.sep + name   
    repo = data["base"] + os.sep + 'report_files/' + name 
    geno = data["base"] + os.sep + 'genome_files/' + name 
    tsv  = parseSTSV(data)
    ssum = parseStats(data)
    html = dumpStats(ssum, loadHTML(html, css))

    def parseAUC(x):
        p = [i for i in x.split('\n') if ' | ' in i and 'AUC' not in i]
        AUC = {}
        for i in p:
            toks = i.split(' | ')
            AUC[float(toks[0].replace('|', ''))] = str(round(float(toks[1].replace('|', '')), 2))
        return AUC

    runR(geno + '_ladder.R', repo + '_F1.png', 5, 5, { 'axis.size':8, 'title.size':11, 'legend.direction':'"vertical"', 'legend.position':'"none"' })
    html = replace(html, '__S_F1__', 'report_files/' + name + '_F1.png')
    
    runR(geno + '_qualFilter.R', repo + '_F2.png', 5, 5, { 'axis.size':8, 'title.size':11 })            
    html = replace(html, '__S_F2__', 'report_files/' + name + '_F2.png')

    AUC = parseAUC(runR(geno + '_ROC.R', repo + '_F3.png', 5, 5, \
                    { 'axis.size':8, 'title.size':11, 'legend.title.size':8, 'legend.text.size':8 }))
    html = replace(html, '__S_F3__', 'report_files/' + name + '_F3.png')

    STable = ''
    tmp = '<tr><td>__AF__</td>\n\
               <td>__OBS__</td>\n\
               <td>__TP__</td>\n\
               <td>__FN__</td>\n\
               <td>__DR__</td>\n\
               <td>__DV__</td>\n\
               <td>__DN__</td>\n\
               <td>__QL__</td>\n\
               <td>__AUC__</td></tr>\n'

    for af in tsv['ef'].split(','):
        tab = tmp
        tab = tab.replace('__AF__',  str('{0:.3f}'.format(float(af))))
        tab = tab.replace('__OBS__', tsv[af]['ob'])
        tab = tab.replace('__TP__',  str(tsv[af]['tp']))
        tab = tab.replace('__FN__',  str(tsv[af]['fn']))
        tab = tab.replace('__DR__',  str(tsv[af]['dr']))
        tab = tab.replace('__DV__',  str(tsv[af]['dv']))
        tab = tab.replace('__DN__',  str(tsv[af]['dn']))
        tab = tab.replace('__QL__',  str(tsv[af]['ql']))
        if af in AUC:
            tab = tab.replace('__AUC__', str(AUC[af]))
        elif round(float(af), 6) in AUC:
            tab = tab.replace('__AUC__', str(AUC[round(float(af), 6)]))            
        else:
            tab = tab.replace('__AUC__', '-')
        STable += tab

    html = replace(html, '__T1__', STable)
    
    (head, tab) = dumpTSV(base + "_sequin.tsv", ['Size'] + dumpIgnores())
    html = replace(html, '__HEAD__', head)
    html = replace(html, '__SEQS__', tab)
    
    tmp = '<tr><td>__T__</td>\n\
               <td>__N__</td>\n\
               <td>__OF__</td>\n\
               <td>__TR__</td>\n\
               <td>__TV__</td>\n\
               <td>__ND__</td>\n\
               <td>__QL__</td></tr\n'

    def fill(att, hack=False):
        x = ''
        m = {}
        for i in att["keys"]:
            m[formatAttr(i)] = i
            
        for key in sortKeys(m.keys()):
            # TODO: Hack for chrQ straification
            if hack and (not "Total" in key and not "Mutation" in key):
                continue
            
            key = m[key]
            tab = copy.copy(tmp)
            tab = replace(tab, '__T__',  formatAttr(key))
            tab = replace(tab, '__N__',  att[key]['n'])
            tab = replace(tab, '__OF__', att[key]['of'])
            tab = replace(tab, '__TR__', att[key]['tr'])
            tab = replace(tab, '__TV__', att[key]['tv'])
            tab = replace(tab, '__ND__', att[key]['nd'])
            tab = replace(tab, '__QL__', att[key]['ql'])
            x += tab
        return x
    
    html = replace(html, "__T2__", fill(tsv['fp']))
    html = replace(html, "__T3__", fill(tsv['sv'], hack=True))    

    writeHTML(out, name, html, { "notCombined": "display:none" } if ("Reference regions (decoy)" in ssum) else { "combined": "display:none" })

#
# Generate a table for synthetic ladder
#

def splitSNTable(sn):
    x = ""
    tmp  = "<tr><td class='gfirst'>__CN__</td>\n\
                <td class='gcell'>__MU__</td>\n\
                <td class='gcell'>__CV__</td>\n\
                <td class='gcell'>__SD__</td>\n\
                <td class='gcell'>__RT__</td>\n\
                <td class='gcell'>__00__</td>\n\
                <td class='gcell'>__25__</td>\n\
                <td class='gcell'>__50__</td>\n\
                <td class='gcell'>__75__</td>\n\
                <td class='gcell'>__100__</td></tr>"

    for i in sn:
        tab = copy.copy(tmp)
        tab = replace(tab, "__CN__",  i["cn"])	
        tab = replace(tab, "__CV__",  i["cv"])	
        tab = replace(tab, "__SD__",  i["sd"])
        tab = replace(tab, "__MU__",  i["mu"])
        tab = replace(tab, "__RT__",  i["rt"])
        tab = replace(tab, "__00__",  i["q0"])
        tab = replace(tab, "__25__",  i["q25"])	
        tab = replace(tab, "__50__",  i["q50"])	
        tab = replace(tab, "__75__",  i["q75"])	
        tab = replace(tab, "__100__", i["q100"])
        x += tab

    return x

def readTSV(file):
    assert(os.path.exists(file))
    c = {} # Column names
    x = {} # Data
    with open(file) as f:
        for line in f:
            toks = line.strip().split('\t')
            if len(c) == 0:
                for i in range(0, len(toks)):
                    c[i] = toks[i]
                    x[toks[i]] = []
                continue
            for i in range(0, len(toks)):
                x[c[i]].append(tryFloat(toks[i]))
    return x

def parseROut(x):
    x = x.replace("[1]", "").split("\n")
    return [str(i.strip()) for i in x if i.strip() != ""]

#
# Fix and validate for all "split" reports
#

def validReport(html):
    # Second input file not provided?
    if "__User sequence file (second)__" in html:
        html = html.replace("__User sequence file (second)__", "-")
    return html

def meta(data, out, html, css):
    name = data["name"]
    base = data["base"] + os.sep + name
    stat = parseStats(data) 
    html = dumpStats(stat, loadHTML(html, css))

    try:
        # Do we have calibrated sequin TSV?
        tsv = (base + "_sequin_calibrated.tsv") if stat["Sequin Calibration"] != "NA" else (base + "_sequin.tsv")

        if not os.path.exists(tsv):
            raise Exception(tsv + " not found")
        
        (s1, r1) = RFilter(tsv, "Mix")
        (s2, r2) = RMeanCV(r1, "c('Mix', 'Count')")

        run("Rscript " + s1)
        run("Rscript " + s2)
        html = replace(html, "__T2__", tsv2HTML(r2, ["Name", "Mean", "CV"], ["Expected Abundance (Normalised)", "Measured Abundance (Mean read Count)", "CV"]))
    except:
        html = replace(html, "__T2__", "")
    
    # Mean Ratio for synthetic ladder
    try:
        # Do we have calibrated ladder TSV?
        tsv = (base + "_ladder_calibrated.tsv") if stat["Ladder Calibration"] != "NA" else (base + "_ladder.tsv")

        if not os.path.exists(tsv):
            raise Exception(tsv + " not found")

        cols = ["Copy", "Mean", "Ratio", "SD", "CV", "Q0", "Q25", "Q50", "Q75", "Q100"]
        (s3, r3) = RLadTable(tsv, cols)

        html = html.replace("__MR__", toFloat(run("Rscript " + s3).replace("[1] ", "").strip()))
        html = replace(html, "__T1__", tsv2HTML(r3, cols, cols))    
    except:
        html = html.replace("__MR__", "-")
        html = replace(html, "__T1__", "Synthetic table not shown. Not enough synthetic sequin detected.")    

    try:
        x = parseROut(runR(base + "_abundance.R", out + "/F4.png", 5, 5, { 'axis.title.y.r':0, 'axis.text':5, 'axis.size':D_AXIS_SIZE, 'title.size':D_TITLE_SIZE, 'legend.direction':'"vertical"', 'legend.position':'"none"' }))
        assert(len(x) == 2) # Slope and then R2       
        html = replace(html, "__Slope_2__", toFloat(x[0], 2)) # Calibrated slope for abundance ladder
        html = replace(html, "__R2_2__", toFloat(x[1], 2))    # Calibrated R2 for abundance ladder
        html = replace(html, "__F4__", 'F4.png')
    except:
        html = replace(html, "__F4__", "")
        html = replace(html, "__R2_2__", "-")
        html = replace(html, "__Slope_2__", "-")
    
    try:
        x = parseROut(runR(base + "_ladderCopy.R", out + "/F1.png", 5, 5, { 'annotate.text':2.5, 'axis.title.y.r':0, 'axis.text':5, 'axis.size':D_AXIS_SIZE, 'title.size':D_TITLE_SIZE, 'legend.direction':'"vertical"', 'legend.position':'"none"' }))
        assert(len(x) == 2) # Slope and then R2       
        html = replace(html, "__Slope_1__", toFloat(x[0], 2)) # Calibrated slope for synthetic ladder
        html = replace(html, "__R2_1__", toFloat(x[1], 2))    # Calibrated R2 for synthetic ladder
        html = replace(html, "__F1__", 'F1.png')
    except:
        html = replace(html, "__F1__", "")
        html = replace(html, "__R2_1__", "-")
        html = replace(html, "__Slope_1__", "-")

    if runR(base + "_ladderDensity.R", out + "/F2.png", 5, 5, { 'axis.title.y.r':0, 'axis.text':5, 'axis.size':D_AXIS_SIZE, 'title.size':D_TITLE_SIZE, 'legend.direction':'"vertical"', 'legend.position':'"none"' }) is None:
        html = replace(html, "__F2__", "")
    else:
        html = replace(html, "__F2__", 'F2.png')

    writeHTML(out, name, validReport(html))

def rna(data, out, html, css):
    name = data["name"]
    base = data["base"] + os.sep + name
    stat = parseStats(data) 
    html = dumpStats(stat, loadHTML(html, css))

    def RGAgg(src):
        R = tempfile.NamedTemporaryFile(delete=False).name
        r = tempfile.NamedTemporaryFile(delete=False).name
        with open(R, "w") as w:
            w.write("library(Anaquin)\n")
            w.write("x <- read.table('" + src + "', sep='\\t', header=T)\n")
            w.write("x[,3:ncol(x)] <- sapply(x[,3:ncol(x)], as.numeric.factor); x[is.na(x)] <- 0\n")
            w.write("x <- aggregate(.~Gene, x[,c(-1)], sum)\n")
            w.write("x$Mix <- round(x$Mix, 4)\n") # Rounding off floating errors
            w.write("write.table(x,'" + r + "',row.names=F,quote=F,sep='\\t')\n")
        return (R, r)

    # Do we have calibrated sequin TSV?
    tsv = (base + "_sequin_calibrated.tsv") if stat["Sequin Calibration"] != "NA" else (base + "_sequin.tsv")

    (s1, r1) = RGAgg(tsv)
    (s3, r3) = RMeanCV(r1, "c('Mix', 'Med')")
    (s4, r4) = RSort(tsv, "Mix", "c('Mix', 'Med')")

    run("Rscript " + s1) # Gene-level aggregation to r1
    run("Rscript " + s3) # MeanCV at the gene level
    run("Rscript " + s4) # Sort isoform results

    html = replace(html, "__T1__", tsv2HTML(r3, ["Name", "Mean", "CV"], ["Expected Sequin Gene Expression (attomol/ul)", "Observed Expression (median TPM)", "Coefficient Variation"]))
    
    try:
        x = parseROut(runR(base + "_isoform.R", out + "/F1.png", 5, 5, { 'annotate.text':2.5, 'axis.title.y.r':0, 'axis.text':5, 'axis.size':D_AXIS_SIZE, 'title.size':D_TITLE_SIZE, 'legend.direction':'"vertical"', 'legend.position':'"none"' }))
        assert(len(x) == 2) # Slope and then R2       
        html = replace(html, "__Slope_1__", toFloat(x[0], 2))
        html = replace(html, "__R2_1__", toFloat(x[1], 2)) 
        html = replace(html, "__F1__", 'F1.png')
    except:
        html = replace(html, "__F1__", "")
        html = replace(html, "__R2_1__", "-")
        html = replace(html, "__Slope_1__", "-")
    
    try:
        x = parseROut(runR(base + "_gene.R", out + "/F2.png", 5, 5, { 'annotate.text':2.5, 'axis.title.y.r':0, 'axis.text':5, 'axis.size':D_AXIS_SIZE, 'title.size':D_TITLE_SIZE, 'legend.direction':'"vertical"', 'legend.position':'"none"' }))
        assert(len(x) == 2) # Slope and then R2       
        html = replace(html, "__Slope_2__", toFloat(x[0], 2))
        html = replace(html, "__R2_2__", toFloat(x[1], 2)) 
        html = replace(html, "__F2__", 'F2.png')
    except:
        html = replace(html, "__F2__", "")
        html = replace(html, "__R2_2__", "-")
        html = replace(html, "__Slope_2__", "-")
    
    writeHTML(out, name, validReport(html))

def genome(data, out, html, css, src=None):
    name = data["name"]
    base = data["base"] + os.sep + name   
    repo = data["base"] + os.sep + 'report_files/' + name 
    geno = data["base"] + os.sep + 'genome_files/' + name 
    html = dumpStats(parseStats(data), loadHTML(html, css))
    
    if src is None:
        src = base + "_sequin.tsv"

    tsv = parseKTSV(data, src)

    try:
        cols = ["Copy", "Mean", "Ratio", "SD", "CV", "Q0", "Q25", "Q50", "Q75", "Q100"]
        (s1, r1) = RLadTable(src, cols)

        # Mean Ratio for synthetic ladder
        html = html.replace("__MR__", toFloat(run("Rscript " + s1).replace("[1] ", "").strip()))
        html = replace(html, "__T1__", tsv2HTML(r1, cols, cols, showHead=False))
    except:
        html = html.replace("__MR__", "")
        html = html.replace("__T1__", "")

    stab = ""
    tmp  = "<tr><td>__EXP__</td>\n\
                <td>__OBS__</td>\n\
                <td>__TP__</td>\n\
                <td>__FN__</td>\n\
                <td>__NR__</td>\n\
                <td>__NV__</td></tr>\n"
    s = tsv['som']

    for key in tsv['som']["af"]:
        key = str(key)
        tab = copy.copy(tmp)
        tab = tmp.replace('__EXP__', key)
        tab = tab.replace('__OBS__', s[key]["Total_Total"]['of'])
        tab = tab.replace('__TP__',  s[key]['tp'])
        tab = tab.replace('__FN__',  s[key]['fn'])
        tab = tab.replace('__NR__',  s[key]['nr'])
        tab = tab.replace('__NV__',  s[key]['nv'])
        stab += tab
    html = replace(html, '__T2__', stab)

    try:
        x = parseROut(runR(base + "_ladderCopy.R", out + "/F1.png", 5, 5, { 'annotate.text':2.5, 'axis.title.y.r':0, 'axis.text':5, 'axis.size':D_AXIS_SIZE, 'title.size':D_TITLE_SIZE, 'legend.direction':'"vertical"', 'legend.position':'"none"' }))
        assert(len(x) == 2) # Slope and then R2       
        html = replace(html, "__Slope_1__", toFloat(x[0], 2)) # Calibrated slope for synthetic ladder
        html = replace(html, "__R2_1__", toFloat(x[1], 2))    # Calibrated R2 for synthetic ladder
        html = replace(html, "__F1__", 'F1.png')
    except:
        html = replace(html, "__F1__", "")
        html = replace(html, "__R2_1__", "-")
        html = replace(html, "__Slope_1__", "-")

    if runR(base + '_ladderDensity.R', out + '/F2.png', 5, 5, { 'axis.title.y.r':0, 'axis.text':5, 'axis.size':6, 'title.size':8, 'legend.direction':'"vertical"', 'legend.position':'"none"' }) is None:
        html = replace(html, "__F2__", '')
    else:
        html = replace(html, "__F2__", 'F2.png')

    if runR(base + '_ladderVariation.R', out + '/F3.png', 5, 5, { 'axis.title.y.r':0, 'axis.text':5, 'axis.size':6, 'title.size':8, 'legend.direction':'"vertical"', 'legend.position':'"none"' }) is None:
        html = replace(html, "__F3__", '')
    else:
        html = replace(html, "__F3__", 'F3.png')

    try:
        x = parseROut(runR(base + '_somatic.R', out + '/F4.png', 5, 5, { 'axis.title.y.r':0, 'axis.text':5, 'axis.size':6, 'title.size':8, 'legend.direction':'"vertical"', 'legend.position':'"none"' }))
        assert(len(x) == 2) # Slope and then R2
        html = replace(html, "__F4__", 'F4.png')        
        html = replace(html, "__R2_4__", toFloat(x[1], 2)) 
        html = replace(html, "__Slope_4__", toFloat(x[0], 2))
    except:
        html = replace(html, "__F4__", "")
        html = replace(html, "__R2_4__", "-")
        html = replace(html, "__Slope_4__", "-")

    return writeHTML(out, name, validReport(html))

def calibrate(data, out, html, css):
    html = genome(data, out, html, css, data["base"] + os.sep + "calibrate_calibrated_sequin.tsv")
    html = html.replace("PARTITION SUMMARY", "PARTITION SUMMARY (pre-calibration)")
    html = html.replace("SOMATIC VARIANTS", "SOMATIC VARIANTS (post-calibration)")
    writeHTML(out, data["name"], validReport(html)) # Replace the HTML written by genome()
    return html

#
# Eg: python scripts/report.py genome output output split scripts/template/gSplit.html scripts/template/style.css
#     python scripts/report.py rna output output rna scripts/template/rSplit.html scripts/template/style.css
#     python scripts/report.py meta output output meta scripts/template/mSplit.html scripts/template/style.css
#     python scripts/report.py calibrate output output calibrate scripts/template/gSplit.html scripts/template/style.css
#     python scripts/report.py germline output output germline scripts/template/germline.html scripts/template/style.css
#     python scripts/report.py somatic output output somatic scripts/template/somatic.html scripts/template/style.css
#

if __name__ == '__main__':    
    mode = sys.argv[1]
    out  = sys.argv[2]

    if not os.path.exists(out):
        raise Exception(out + " not existed")
    if mode == "germline":
        germline({ "base": os.path.abspath(sys.argv[3]), "name": sys.argv[4] }, out, sys.argv[5], sys.argv[6])
    elif mode == "somatic":
        somatic({ "base": os.path.abspath(sys.argv[3]), "name": sys.argv[4] }, out, sys.argv[5], sys.argv[6])
    elif mode == "calibrate":
        calibrate({ "base": os.path.abspath(sys.argv[3]), "name": sys.argv[4] }, out, sys.argv[5], sys.argv[6])
    elif mode == "genome":
        genome({ "base": os.path.abspath(sys.argv[3]), "name": sys.argv[4] }, out, sys.argv[5], sys.argv[6])
    elif mode == "rna":
        rna({ "base": os.path.abspath(sys.argv[3]), "name": sys.argv[4] }, out, sys.argv[5], sys.argv[6])
    elif mode == "meta":
        meta({ "base": os.path.abspath(sys.argv[3]), "name": sys.argv[4] }, out, sys.argv[5], sys.argv[6])
    else:
        raise Exception('Unknown ' + str(mode))
#<<@@@@>>
