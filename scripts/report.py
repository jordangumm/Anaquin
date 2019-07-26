import os
import sys
import datetime
import tempfile
import functools
import subprocess

D_AXIS_SIZE  = 8
D_TITLE_SIZE = 8

def run(cmd):
    x = os.popen(cmd).read()
    return "" if x is None else x

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
        return run("Rscript --vanilla " + name)
    except subprocess.CalledProcessError:
        return None

def parseROut(x):
    x = x.replace("[1]", "").split("\n")
    return [str(i.strip()) for i in x if i.strip() != ""]

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

def germline(base, out, html, css):
    html = loadHTML(html, css)
    html = copyStats(base + "/germline_summary.txt", html)
    html = replace(html, "__T1__", copyTable(base + "/germline_variant_table.tsv"))
    html = replace(html, "__T2__", copyTable(base + "/germline_sequin.tsv", ["NAME", "CHROM", "POSITION", "LABEL", "GENOTYPE", "TYPE", "EXP_FREQ", "REF_DEPTH", "VAR_DEPTH", "OBS_FREQ", "QUAL"]))
    writeHTML(out, "germline", validReport(html))

def somatic(base, out, html, css):
    html = loadHTML(html, css)
    html = copyStats(base + "/somatic_summary.txt", html)
    html = replace(html, "__T1__", copyTable(base + "/somatic_variant_table.tsv"))
    html = replace(html, "__T2__", copyTable(base + "/somatic_sequin.tsv"))

    runR(base + '/report_files/somatic_ladder.R', out + "/F1.png", 5, 5, { 'axis.size':8, 'title.size':11, 'legend.direction':'"vertical"', 'legend.position':'"none"' })
    html = replace(html, '__F1__', "F1.png")
    
    runR(base + '/report_files/somatic_qualFilter.R', out + "/F2.png", 5, 5, { 'axis.size':8, 'title.size':11 })            
    html = replace(html, '__F2__', "F2.png")

    runR(base + '/report_files/somatic_ROC.R', out + "/F3.png", 5, 5, { 'axis.size':8, 'title.size':11 })            
    html = replace(html, '__F3__', "F3.png")

    writeHTML(out, "somatic", validReport(html))

def validReport(html):
    # Second input file not provided?
    if "__User sequence file (second)__" in html:
        html = html.replace("__User sequence file (second)__", "-")
    return html

def copyTable(file, cols = None, tdStyle="width:100px"):
    isHead = True
    html = "<table style='margin-top:20px'>"
    cs = []
    with open(file) as f:
        for line in f:
            toks = line.strip().split('\t')
            if isHead:
                isHead = False
                html += "<thead><tr>"
                for i in range(0, len(toks)):
                    if cols is None or toks[i] in cols:
                        cs.append(i)
                        html += ('<th class="ghead">' + toks[i] + "</th>")
                html += "</tr></thead><tbody>\n"
                continue
            html += "<tr>"
            for i in range(len(toks)):
                if i in cs:
                    html += ("<td style='" + tdStyle + "'>" + toks[i] + "</td>")
            html += "</tr>\n"
    return html + "</tbody></table>"
    
def copyStats(file, html):
    block = ""
    with open(file, "r") as r:
        for line in r:
            line = line.strip()
            if len(line) == 0:
                continue
            elif not ":" in line:
                block = line
                continue
            else:
                toks = [i.strip() for i in line.strip().split(':') if len(i) > 0]
                assert(len(toks) >= 1)

                if len(toks) > 1:
                    key = ' '.join(toks[0:len(toks)-1]).replace(':', '')
                    val = toks[len(toks) - 1]
                elif len(toks) == 1:
                    key = toks[0]
                    val = ""
                if ("__" + key + "__") in html:
                    html = replace(html, "__" + key + "__", val)
                elif block != "" and ("__" + block + "-" + key + "__") in html:
                    html = replace(html, "__" + block + "-" + key + "__", val)      
    return html

def meta(base, out, html, css):
    html = loadHTML(html, css)
    html = copyStats(base + "/meta_summary.txt", html)
    html = replace(html, "__T1__", copyTable(base + "/meta_ladder_table.tsv"))
    html = replace(html, "__T2__", copyTable(base + "/meta_sequin_table.tsv"))

    if runR(base + "/report_files/meta_ladderCopy.R", out + "/F1.png", 5, 5, { 'annotate.text':2.5, 'axis.title.y.r':0, 'axis.text':5, 'axis.size':D_AXIS_SIZE, 'title.size':D_TITLE_SIZE, 'legend.direction':'"vertical"', 'legend.position':'"none"' }) is not None:
        html = replace(html, "__F1__", "F1.png")
    else:
        html = replace(html, "__F1__", "")

    if runR(base + "/report_files/meta_ladderDensity.R", out + "/F2.png", 5, 5, { 'axis.title.y.r':0, 'axis.text':5, 'axis.size':D_AXIS_SIZE, 'title.size':D_TITLE_SIZE, 'legend.direction':'"vertical"', 'legend.position':'"none"' }) is not None:
        html = replace(html, "__F2__", "F2.png")
    else:
        html = replace(html, "__F2__", "")

    if runR(base + "/report_files/meta_abundance.R", out + "/F3.png", 5, 5, { 'axis.title.y.r':0, 'axis.text':5, 'axis.size':D_AXIS_SIZE, 'title.size':D_TITLE_SIZE, 'legend.direction':'"vertical"', 'legend.position':'"none"' }) is not None:
        html = replace(html, "__F3__", "F3.png")
    else:
        html = replace(html, "__F3__", "")

    writeHTML(out, "meta", validReport(html))

def rna(base, out, html, css):
    html = loadHTML(html, css)
    html = copyStats(base + "/rna_summary.txt", html)
    html = replace(html, "__T1__", copyTable(base + "/rna_sequin_gene_table.tsv"))

    if runR(base + "/report_files/rna_isoform.R", out + "/F1.png", 5, 5, { 'annotate.text':2.5, 'axis.title.y.r':0, 'axis.text':5, 'axis.size':D_AXIS_SIZE, 'title.size':D_TITLE_SIZE, 'legend.direction':'"vertical"', 'legend.position':'"none"' }):
        html = replace(html, "__F1__", "F1.png")
    else:
        html = replace(html, "__F1__", "")

    if runR(base + "/report_files/rna_gene.R", out + "/F2.png", 5, 5, { 'annotate.text':2.5, 'axis.title.y.r':0, 'axis.text':5, 'axis.size':D_AXIS_SIZE, 'title.size':D_TITLE_SIZE, 'legend.direction':'"vertical"', 'legend.position':'"none"' }):
        html = replace(html, "__F2__", "F2.png")
    else:
        html = replace(html, "__F2__", "")
    
    writeHTML(out, "rna", validReport(html))    

def genome(base, out, html, css, src=None):
    html = loadHTML(html, css)
    html = copyStats(base + "/split_summary.txt", html)
    html = replace(html, "__T1__", copyTable(base + "/split_ladder_table.tsv"))
    html = replace(html, "__T2__", copyTable(base + "/split_sequin_table.tsv"))    
    
    if runR(base + "/report_files/split_ladderCopy.R", out + "/F1.png", 5, 5, { 'annotate.text':2.5, 'axis.title.y.r':0, 'axis.text':5, 'axis.size':D_AXIS_SIZE, 'title.size':D_TITLE_SIZE, 'legend.direction':'"vertical"', 'legend.position':'"none"' }) is not None:
        html = replace(html, "__F1__", "F1.png")
    else:
        html = replace(html, "__F1__", "")

    if runR(base + "/report_files/split_ladderDensity.R", out + "/F2.png", 5, 5, { 'axis.title.y.r':0, 'axis.text':5, 'axis.size':D_AXIS_SIZE, 'title.size':D_TITLE_SIZE, 'legend.direction':'"vertical"', 'legend.position':'"none"' }) is not None:
        html = replace(html, "__F2__", "F2.png")
    else:
        html = replace(html, "__F2__", "")

    if runR(base + "/report_files/split_somatic.R", out + "/F4.png", 5, 5, { 'axis.title.y.r':0, 'axis.text':5, 'axis.size':D_AXIS_SIZE, 'title.size':D_TITLE_SIZE, 'legend.direction':'"vertical"', 'legend.position':'"none"' }) is not None:
        html = replace(html, "__F4__", "F4.png")
    else:
        html = replace(html, "__F4__", "")
    
    return writeHTML(out, "split", validReport(html))

def calibrate(base, out, html, css):
    html = loadHTML(html, css)
    html = copyStats(base + "/calibrate_summary.txt", html)
    html = replace(html, "__T2__", copyTable(base + "/calibrate_sequin_table.tsv"))

    if runR(base + "/report_files/calibrate_somatic.R", out + "/F4.png", 5, 5, { 'axis.title.y.r':0, 'axis.text':5, 'axis.size':D_AXIS_SIZE, 'title.size':D_TITLE_SIZE, 'legend.direction':'"vertical"', 'legend.position':'"none"' }) is not None:
        html = replace(html, "__F4__", "F4.png")
    else:
        html = replace(html, "__F4__", "")

    if runR(base + "/report_files/calibrate_insert.R", out + "/F5.png", 5, 5, { 'axis.title.y.r':0, 'axis.text':5, 'axis.size':D_AXIS_SIZE, 'title.size':D_TITLE_SIZE, 'legend.direction':'"vertical"', 'legend.position':'"none"' }) is not None:
        html = replace(html, "__F5__", "F5.png")
    else:
        html = replace(html, "__F5__", "")

    writeHTML(out, "calibrate", validReport(html))

#
# Eg: python scripts/report.py genome output output/report_files scripts/template/gSplit.html scripts/template/style.css
#     python scripts/report.py rna output output/report_files scripts/template/rSplit.html scripts/template/style.css
#     python scripts/report.py meta output output/report_files scripts/template/mSplit.html scripts/template/style.css
#     python scripts/report.py calibrate output/report_files output scripts/template/calibrate.html scripts/template/style.css
#     python scripts/report.py germline output/report_files output scripts/template/germline.html scripts/template/style.css
#     python scripts/report.py somatic output/report_files output scripts/template/somatic.html scripts/template/style.css
#

if __name__ == '__main__':    
    mode = sys.argv[1]
    base = os.path.abspath(sys.argv[2])
    out  = os.path.abspath(sys.argv[3])

    assert(os.path.exists(base))
    if not os.path.exists(out):
        os.makedirs(out)

    if mode == "germline":
        germline(base, out, sys.argv[4], sys.argv[5])
    elif mode == "somatic":
        somatic(base, out, sys.argv[4], sys.argv[5])
    elif mode == "calibrate":
        calibrate(base, out, sys.argv[4], sys.argv[5])
    elif mode == "genome":
        genome(base, out, sys.argv[4], sys.argv[5])
    elif mode == "rna":
        rna(base, out, sys.argv[4], sys.argv[5])
    elif mode == "meta":
        meta(base, out, sys.argv[4], sys.argv[5])
    else:
        raise Exception('Unknown ' + str(mode))
#<<@@@@>>