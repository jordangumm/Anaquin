#include <fstream>
#include "kseq.h"
#include "ProcessReads.h"
#include "tools/samtools.hpp"
#include "parsers/parser_bam.hpp"

int __Combined__;
std::vector<const Anaquin::ReadName *> __HackBAMDecoy2__;

// Combined reads
Anaquin::AlignedReads __CombinedBAM1__;

// Combined BAM header
bam_hdr_t *__CombinedBAMHeader__;

int ProcessReads(KmerIndex& index, const ProgramOptions& opt, MinCollector& tc)
{
    int limit = 1048576;
    std::vector<std::pair<const char*, int>> seqs;
    seqs.reserve(limit/50);

    //SequenceReader SR(opt);

    // need to receive an index map
    std::ios_base::sync_with_stdio(false);

    //int tlencount = (opt.fld == 0.0) ? 10000 : 0;
    size_t numreads = 0;
    size_t nummapped = 0;

    if (opt.pseudobam) {
        index.writePseudoBamHeader(std::cout);
    }

    MasterProcessor MP(index, opt, tc);
    MP.processReads();
    numreads = MP.numreads;
    nummapped = MP.nummapped;

    // write output to outdir
    if (opt.write_index) {
        std::string outfile = opt.output + "/counts.txt";
        std::ofstream of;
        of.open(outfile.c_str(), std::ios::out);
        tc.write(of);
        of.close();
    }

    return numreads;
}

void MasterProcessor::processReads() {
    // start worker threads
    if (!opt.batch_mode) {
        std::vector<std::thread> workers;
        for (int i = 0; i < opt.threads; i++) {
            workers.emplace_back(std::thread(ReadProcessor(index,opt,tc,*this,i)));
        }

        // let the workers do their thing
        for (int i = 0; i < opt.threads; i++) {
            workers[i].join(); //wait for them to finish
        }

        // now handle the modification of the mincollector
        for (auto &t : newECcount) {
            if (t.second <= 0) {
                continue;
            }
            int ec = tc.increaseCount(t.first); // modifies the ecmap

            if (ec != -1 && t.second > 1) {
                tc.counts[ec] += (t.second-1);
            }
        }
    } else {
        std::vector<std::thread> workers;
        int num_ids = opt.batch_ids.size();
        int id =0;
        while (id < num_ids) {
            // TODO: put in thread pool
            workers.clear();
            int nt = std::min(opt.threads, (num_ids - id));
            for (int i = 0; i < nt; i++,id++) {
                workers.emplace_back(std::thread(ReadProcessor(index, opt, tc, *this, id)));
            }

            for (int i = 0; i < nt; i++) {
                workers[i].join();
            }

            if (opt.umi) {
                // process the regular EC umi now
                for (int i = 0; i < nt; i++) {
                    int l_id = id - nt + i;
                    auto &umis = batchUmis[l_id];
                    std::sort(umis.begin(), umis.end());
                    size_t sz = umis.size();
                    nummapped += sz;
                    if (sz > 0) {
                        ++batchCounts[l_id][umis[0].first];
                    }
                    for (auto j = 1u; j < sz; j++) {
                        if (umis[j-1] != umis[j]) {
                            ++batchCounts[l_id][umis[j].first];
                        }
                    }
                    umis.clear();
                }
            }
        }

        int num_newEcs = 0;
        if (!opt.umi) {
            // for each cell
            for (int id = 0; id < num_ids; id++) {
                // for each new ec
                for (auto &t : newBatchECcount[id]) {
                    // add the ec and count number of new ecs
                    if (t.second <= 0) {
                        continue;
                    }
                    int ecsize = index.ecmap.size();
                    int ec = tc.increaseCount(t.first);
                    if (ec != -1 && ecsize < index.ecmap.size()) {
                        num_newEcs++;
                    }
                }
            }
            // for each cell
            for (int id = 0; id < num_ids; id++) {
                auto& c = batchCounts[id];
                c.resize(c.size() + num_newEcs,0);
                // for each new ec
                for (auto &t : newBatchECcount[id]) {
                    // count the ec
                    if (t.second <= 0) {
                        continue;
                    }
                    int ec = tc.findEC(t.first);
                    assert(ec != -1);
                    ++c[ec];
                }
            }
        } else {
            // UMI case
            // for each cell
            for (int id = 0; id < num_ids; id++) {
                // for each new ec
                for (auto &t : newBatchECumis[id]) {
                    // add the new ec
                    int ecsize = index.ecmap.size();
                    int ec = tc.increaseCount(t.first);
                    if (ec != -1 && ecsize < index.ecmap.size()) {
                        num_newEcs++;
                    }
                }
            }
            // for each cell
            for (int id = 0; id < num_ids; id++) {
                auto& c = batchCounts[id];
                c.resize(c.size() + num_newEcs,0);
                std::vector<std::pair<int, std::string>> umis;
                umis.reserve(newBatchECumis[id].size());
                // for each new ec
                for (auto &t : newBatchECumis[id]) {
                    // record the ec,umi
                    int ec = tc.findEC(t.first);
                    umis.push_back({ec, std::move(t.second)});
                }
                // find unique umis per ec
                std::sort(umis.begin(), umis.end());
                size_t sz = umis.size();
                if (sz > 0) {
                    ++batchCounts[id][umis[0].first];
                }
                for (auto j = 1u; j < sz; j++) {
                    if (umis[j-1] != umis[j]) {
                        ++batchCounts[id][umis[j].first];
                    }
                }
                for (auto x : c) {
                    num_umi += x;
                }
            }
        }
    }
}

void MasterProcessor::update(const std::vector<int>& c, const std::vector<std::vector<int> > &newEcs,
                             std::vector<std::pair<int, std::string>>& ec_umi, std::vector<std::pair<std::vector<int>, std::string>> &new_ec_umi,
                             int n, std::vector<int>& flens, std::vector<int> &bias, int id) {
    // acquire the writer lock
    std::lock_guard<std::mutex> lock(this->writer_lock);

    if (!opt.batch_mode) {
        for (auto i = 0u; i < c.size(); i++) {
            tc.counts[i] += c[i]; // add up ec counts
            nummapped += c[i];
        }
    } else {
        if (!opt.umi) {
            for (auto i = 0u; i < c.size(); i++) {
                batchCounts[id][i] += c[i];
                nummapped += c[i];
            }
        } else {
            for (auto &t : ec_umi) {
                batchUmis[id].push_back(std::move(t));
            }
        }
    }

    if (!opt.batch_mode) {
        for(auto &u : newEcs) {
            ++newECcount[u];
        }
    } else {
        if (!opt.umi) {
            for(auto &u : newEcs) {
                ++newBatchECcount[id][u];
            }
        } else {
            for (auto &u : new_ec_umi) {
                newBatchECumis[id].push_back(std::move(u));
            }
        }
    }
    if (!opt.umi) {
        nummapped += newEcs.size();
    } else {
        nummapped += new_ec_umi.size();
    }

    if (!flens.empty()) {
        int local_tlencount = 0;
        for (int i = 0; i < flens.size(); i++) {
            tc.flens[i] += flens[i];
            local_tlencount += flens[i];
        }
        tlencount += local_tlencount;
    }

    if (!bias.empty()) {
        int local_biasCount = 0;
        for (int i = 0; i < bias.size(); i++) {
            tc.bias5[i] += bias[i];
            local_biasCount += bias[i];
        }
        biasCount += local_biasCount;
    }

    numreads += n;
    // releases the lock
}

void MasterProcessor::outputFusion(const std::stringstream &o) {
    std::string os = o.str();
    if (!os.empty()) {
        std::lock_guard<std::mutex> lock(this->writer_lock);
        ofusion << os << "\n";
    }
}

ReadProcessor::ReadProcessor(const KmerIndex& index, const ProgramOptions& opt, const MinCollector& tc, MasterProcessor& mp, int _id) :
    paired(!opt.single_end), tc(tc), index(index), mp(mp), id(_id) {
    // initialize buffer
    bufsize = 1ULL<<23;
    buffer = new char[bufsize];

    if (opt.batch_mode) {
        assert(id != -1);
        batchSR.files = opt.batch_files[id];
        if (opt.umi) {
            batchSR.umi_files = {opt.umi_files[id]};
        }
        batchSR.paired = !opt.single_end;
    }

    seqs.reserve(bufsize/50);
    if (opt.umi) {
        umis.reserve(bufsize/50);
    }
    newEcs.reserve(1000);
    counts.reserve((int) (tc.counts.size() * 1.25));
    clear();
}

ReadProcessor::ReadProcessor(ReadProcessor && o) :
    paired(o.paired),
    tc(o.tc),
    index(o.index),
    mp(o.mp),
    id(o.id),
    bufsize(o.bufsize),
    numreads(o.numreads),
    seqs(std::move(o.seqs)),
    names(std::move(o.names)),
    quals(std::move(o.quals)),
    umis(std::move(o.umis)),
    newEcs(std::move(o.newEcs)),
    flens(std::move(o.flens)),
    bias5(std::move(o.bias5)),
    batchSR(std::move(o.batchSR)),
    counts(std::move(o.counts)) {
    buffer = o.buffer;
    o.buffer = nullptr;
    o.bufsize = 0;
}

ReadProcessor::~ReadProcessor() {
    if (buffer != nullptr) {
        delete[] buffer;
        buffer = nullptr;
    }
}

void ReadProcessor::operator()() {
    while (true) {
        // grab the reader lock
        if (mp.opt.batch_mode) {
            if (batchSR.empty()) {
                return;
            } else {
                batchSR.fetchSequences(buffer, bufsize, seqs, names, quals, bams, rcs, umis, false);
            }
        } else {
            std::lock_guard<std::mutex> lock(mp.reader_lock);
            if (mp.SR.empty()) {
                // nothing to do
                return;
            } else {
                // get new sequences
                mp.SR.fetchSequences(buffer, bufsize, seqs, names, quals, bams, rcs, umis, mp.opt.pseudobam || mp.opt.fusion);
            }
            // release the reader lock
        }

        // process our sequences
        processBuffer();

        // update the results, MP acquires the lock
        mp.update(counts, newEcs, ec_umi, new_ec_umi, paired ? seqs.size()/2 : seqs.size(), flens, bias5, id);
        clear();
    }
}

static unsigned __i__ = 0;

void ReadProcessor::processBuffer()
{
    const char* s1 = nullptr;
    const char* s2 = nullptr;
    const char* r1 = nullptr;
    const char* r2 = nullptr;
    const char* q1 = nullptr;
    const char* q2 = nullptr;
    void * b1 = nullptr;
    void * b2 = nullptr;
    bool rc1, rc2;

    int l1,l2;

    for (auto i = 0u; i < seqs.size(); i++)
    {
        s1 = seqs[i].first;
        l1 = seqs[i].second;
        r1 = names[i].first;
        q1 = quals[i].first;
        b1 = bams.empty() ? nullptr : bams[i];
        rc1 = rcs.empty() ? false : rcs[i];

        i++;
        s2 = seqs[i].first;
        l2 = seqs[i].second;
        r2 = names[i].first;
        q2 = quals[i].first;
        b2 = bams.empty() ? nullptr : bams[i];
        rc2 = rcs.empty() ? false : rcs[i];

        numreads++;
        __i__++;
        
        if ((__i__ % 10000) == 0)
        {
            extern void wait(const std::string &);
            wait(std::to_string(__i__));
        }

        extern void KPartition(unsigned, const char *, const char *, const char *, void *, bool, const char *, const char *, const char *, void *, bool);
        KPartition(id, r1, s1, q1, b1, rc1, r2, s2, q2, b2, rc2);
    }
}

void ReadProcessor::clear() {
    numreads=0;
    memset(buffer,0,bufsize);
    newEcs.clear();
    counts.clear();
    counts.resize(tc.counts.size(),0);
    ec_umi.clear();
    new_ec_umi.clear();
}

/** -- sequence reader -- **/
SequenceReader::~SequenceReader() {
    if (fp1) {
        gzclose(fp1);
    }
    if (paired && fp2) {
        gzclose(fp2);
    }

    kseq_destroy(seq1);
    if (paired) {
        kseq_destroy(seq2);
    }

    // check if umi stream is open, then close
}

// Returns true if there is more left to read from the files
bool SequenceReader::fetchSequences(char *buf,
                                    const int limit,
                                    std::vector<std::pair<const char *, int> > &seqs,
                                    std::vector<std::pair<const char *, int> > &names,
                                    std::vector<std::pair<const char *, int> > &quals,
                                    std::vector<void *> &bams,
                                    std::vector<bool>   &rcs,
                                    std::vector<std::string> &umis,
                                    bool full)
{
    std::string line;
    std::string umi;

    seqs.clear();
    umis.clear();
    bams.clear();
    rcs.clear();

    if (full) {
        names.clear();
        quals.clear();
    }

    bool usingUMIfiles = !umi_files.empty();
    //int umis_read = 0;

    int bufpos = 0;
    int pad = (paired) ? 2 : 1;

    while (true)
    {
        if (!state) { // should we open a file
            if (current_file >= files.size()) {
                // nothing left
                return false;
            } else {
                // close the current file
                if(fp1) {
                    gzclose(fp1);
                }
                if (paired && fp2) {
                    gzclose(fp2);
                }
                // close current umi file
                if (usingUMIfiles) {
                    // read up the rest of the files
                    f_umi->close();
                }

                switch (__Combined__)
                {
                    case 0:
                    {
                        // open the next one
                        fp1 = gzopen(files[current_file].c_str(),"r");
                        seq1 = kseq_init(fp1);
                        l1 = kseq_read(seq1);
                        if (paired) {
                            current_file++;
                            fp2 = gzopen(files[current_file].c_str(),"r");
                            seq2 = kseq_init(fp2);
                            l2 = kseq_read(seq2);
                        }
                        if (usingUMIfiles) {
                            // open new umi file
                            f_umi->open(umi_files[current_file]);
                        }
                        break;
                    }
                        
                    case 1: { break; }
                }
                
                state = true;
            }
        }

        switch (__Combined__)
        {
            case 1:
            {
                for (auto i = 0; i < limit; i++)
                {
                    if (_bamI >= __HackBAMDecoy2__.size())
                    {
                        current_file++; // move to next file
                        state = false;  // haven't opened file yet
                        break;
                    }
                    else
                    {
                        const auto &n = *(__HackBAMDecoy2__[_bamI]);
                        const auto &x = __CombinedBAM1__[n];
                        
                        if (x.first && x.second)
                        {
                            rcs.push_back(x.first->rc);
                            rcs.push_back(x.second->rc);
                            bams.emplace_back(x.first->b);
                            bams.emplace_back(x.second->b);
                            seqs.emplace_back (nullptr, 0);
                            seqs.emplace_back (nullptr, 0);
                            quals.emplace_back(nullptr, 0);
                            quals.emplace_back(nullptr, 0);
                            names.emplace_back(n.c_str(), n.size());
                            names.emplace_back(n.c_str(), n.size());
                        }
                    }
                    
                    // Read next record
                    _bamI++;
                }
                
                if (_bamI < __HackBAMDecoy2__.size())
                {
                    return true; // Read it next time
                }

                break;
            }
                
            case 0:
            {
                // the file is open and we have read into seq1 and seq2
                if (l1 > 0 && (!paired || l2 > 0)) {
                    int bufadd = l1 + l2 + pad;
                    // fits into the buffer
                    if (full) {
                        nl1 = seq1->name.l;
                        if (paired) {
                            nl2 = seq2->name.l;
                        }
                        bufadd += (l1+l2) + pad + (nl1+nl2)+ pad;
                    }
                    if (bufpos+bufadd< limit) {
                        char *p1 = buf+bufpos;
                        memcpy(p1, seq1->seq.s, l1+1);
                        bufpos += l1+1;
                        seqs.emplace_back(p1,l1);
                        
                        if (usingUMIfiles) {
                            std::stringstream ss;
                            std::getline(*f_umi, line);
                            ss.str(line);
                            ss >> umi;
                            umis.emplace_back(std::move(umi));
                        }
                        if (full) {
                            p1 = buf+bufpos;
                            memcpy(p1, seq1->qual.s,l1+1);
                            bufpos += l1+1;
                            quals.emplace_back(p1,l1);
                            p1 = buf+bufpos;
                            memcpy(p1, seq1->name.s,nl1+1);
                            bufpos += nl1+1;
                            names.emplace_back(p1,nl1);
                        }
                        
                        if (paired) {
                            char *p2 = buf+bufpos;
                            memcpy(p2, seq2->seq.s,l2+1);
                            bufpos += l2+1;
                            seqs.emplace_back(p2,l2);
                            if (full) {
                                p2 = buf+bufpos;
                                memcpy(p2,seq2->qual.s,l2+1);
                                bufpos += l2+1;
                                quals.emplace_back(p2,l2);
                                p2 = buf + bufpos;
                                memcpy(p2,seq2->name.s,nl2+1);
                                bufpos += nl2+1;
                                names.emplace_back(p2,nl2);
                            }
                        }
                    } else {
                        return true; // read it next time
                    }
                    
                    // read for the next one
                    l1 = kseq_read(seq1);
                    if (paired) {
                        l2 = kseq_read(seq2);
                    }
                } else {
                    current_file++; // move to next file
                    state = false; // haven't opened file yet
                }
                
                break;
            }
        }
    }
}

bool SequenceReader::empty() {
    return (!state && current_file >= files.size());
}

SequenceReader::SequenceReader(SequenceReader&& o) :
    fp1(o.fp1),
    fp2(o.fp2),
    seq1(o.seq1),
    seq2(o.seq2),
    l1(o.l1),
    l2(o.l2),
    nl1(o.nl1),
    nl2(o.nl2),
    paired(o.paired),
    files(std::move(o.files)),
    umi_files(std::move(o.umi_files)),
    f_umi(std::move(o.f_umi)),
    current_file(o.current_file),
    state(o.state) {
    o.fp1 = nullptr;
    o.fp2 = nullptr;
    o.seq1 = nullptr;
    o.seq2 = nullptr;
    o.state = false;
}
