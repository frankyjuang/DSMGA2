/***************************************************************************
 *   Copyright (C) 2015 Tian-Li Yu and Shih-Huan Hsu                       *
 *   tianliyu@ntu.edu.tw                                                   *
 ***************************************************************************/

#include <list>
#include <vector>
#include <algorithm>
#include <iterator>

#include <iostream>
#include <cstdio>
#include "chromosome.h"
#include "dsmga2.h"
#include "fastcounting.h"
#include "statistics.h"


using namespace std;


DSMGA2::DSMGA2 (int n_ell, int n_nInitial, int n_maxGen, int n_maxFe, int fffff) {


    previousFitnessMean = -INF;
    ell = n_ell;                    // length of chromosome
    nCurrent = (n_nInitial/2)*2;    // population size          // has to be even

    Chromosome::function = (Chromosome::Function)fffff;
    Chromosome::nfe = 0;            // number of function evaluation
    Chromosome::lsnfe = 0;
    Chromosome::hitnfe = 0;
    Chromosome::hit = false;

    selectionPressure = 2;
    maxGen = n_maxGen;
    maxFe = n_maxFe;

    graph.init(ell);

    bestIndex = -1;
    masks = new list<int>[ell];
    selectionIndex = new int[nCurrent];
    orderN = new int[nCurrent];
    orderELL = new int[ell];
    population = new Chromosome[nCurrent];
    fastCounting = new FastCounting[ell];

    for (int i = 0; i < ell; i++)
        fastCounting[i].init(nCurrent);


    pHash.clear();
    for (int i=0; i<nCurrent; ++i) {
        population[i].initR(ell);
        double f = population[i].getFitness();
        pHash[population[i].getKey()] = f;
    }

    if (GHC) {
        for (int i=0; i < nCurrent; i++)
            population[i].GHC();
    }
}


DSMGA2::~DSMGA2 () {
    delete []masks;
    delete []orderN;
    delete []orderELL;
    delete []selectionIndex;
    delete []population;
    delete []fastCounting;
}



bool DSMGA2::isSteadyState () {

    if (stFitness.getNumber () <= 0)
        return false;

    if (previousFitnessMean < stFitness.getMean ()) {
        previousFitnessMean = stFitness.getMean () + EPSILON;
        return false;
    }

    return true;
}



int DSMGA2::doIt (bool output) {
    generation = 0;

    while (!shouldTerminate ()) {
        oneRun (output);
    }
    return generation;
}


void DSMGA2::oneRun (bool output) {

    if (CACHE)
        Chromosome::cache.clear();

    mixing();


    double max = -INF;
    stFitness.reset ();

    for (int i = 0; i < nCurrent; ++i) {
        double fitness = population[i].getFitness();
        if (fitness > max) {
            max = fitness;
            bestIndex = i;
        }
        stFitness.record (fitness);

    }

    if (output)
        showStatistics ();

    ++generation;
}


bool DSMGA2::shouldTerminate () {
    bool
    termination = false;

    if (maxFe != -1) {
        if (Chromosome::nfe > maxFe)
            termination = true;
    }

    if (maxGen != -1) {
        if (generation > maxGen)
            termination = true;
    }


    if (population[0].getMaxFitness() <= stFitness.getMax() )
        termination = true;


    if (stFitness.getMax() - EPSILON <= stFitness.getMean() )
        termination = true;

    return termination;

}


bool DSMGA2::foundOptima () {
    return (stFitness.getMax() > population[0].getMaxFitness());
}


void DSMGA2::showStatistics () {

    printf ("Gen:%d  Fitness:(Max/Mean/Min):%f/%f/%f \n",
            generation, stFitness.getMax (), stFitness.getMean (),
            stFitness.getMin ());
    fflush(NULL);
}



void DSMGA2::buildFastCounting() {

    if (SELECTION) {
        for (int i = 0; i < nCurrent; i++)
            for (int j = 0; j < ell; j++) {
                fastCounting[j].setVal(i, population[selectionIndex[i]].getVal(j));
            }

    } else {
        for (int i = 0; i < nCurrent; i++)
            for (int j = 0; j < ell; j++) {
                fastCounting[j].setVal(i, population[i].getVal(j));
            }
    }

}

int DSMGA2::countOne(int x) const {

    int n = 0;

    for (int i=0; i<fastCounting[0].lengthLong; ++i) {
        unsigned long val = 0;

        val = fastCounting[x].gene[i];

        n += myBD.countOne(val);
    }

    return n;
}


int DSMGA2::countXOR(int x, int y) const {

    int n = 0;

    for (int i=0; i<fastCounting[0].lengthLong; ++i) {
        unsigned long val = 0;


        val = fastCounting[x].gene[i];

        val ^= fastCounting[y].gene[i];

        n += myBD.countOne(val);
    }

    return n;
}

int DSMGA2::countAND3(int x, int y, int z) const {

    int n = 0;

    for (int i = 0; i < fastCounting[0].lengthLong; ++i) {
        unsigned long val = 0;

        val = fastCounting[x].gene[i];
        cout << val << endl;
        val &= fastCounting[y].gene[i];
        cout << val << endl;
        val &= fastCounting[z].gene[i];
        cout << val << endl;

        n += myBD.countOne(val);
    }

    return n;
}


void DSMGA2::restrictedMixing(Chromosome& ch) {

    int r = myRand.uniformInt(0, ell-1);

    list<int> mask = masks[r];

    size_t size = findSize(ch, mask);       // return longest available mask length
    if (size > (size_t)ell/2)
        size = ell/2;

    // prune mask to exactly size
    while (mask.size() > size)
        mask.pop_back();


    bool taken = restrictedMixing(ch, mask);

    EQ = true;
    if (taken) {

        genOrderN();

        for (int i=0; i<nCurrent; ++i) {

            if (EQ)
                backMixingE(ch, mask, population[orderN[i]]);
            else
                backMixing(ch, mask, population[orderN[i]]);
        }
    }

}

void DSMGA2::backMixing(Chromosome& source, list<int>& mask, Chromosome& des) {

    Chromosome trial(ell);
    trial = des;
    for (list<int>::iterator it = mask.begin(); it != mask.end(); ++it)
        trial.setVal(*it, source.getVal(*it));

    if (trial.getFitness() > des.getFitness()) {
        pHash.erase(des.getKey());
        pHash[trial.getKey()] = trial.getFitness();
        des = trial;
        return;
    }

}

void DSMGA2::backMixingE(Chromosome& source, list<int>& mask, Chromosome& des) {

    Chromosome trial(ell);
    trial = des;
    for (list<int>::iterator it = mask.begin(); it != mask.end(); ++it)
        trial.setVal(*it, source.getVal(*it));

    if (trial.getFitness() > des.getFitness()) {
        pHash.erase(des.getKey());
        pHash[trial.getKey()] = trial.getFitness();

        EQ = false;
        des = trial;
        return;
    }

    if (trial.getFitness() >= des.getFitness()) {
        pHash.erase(des.getKey());
        pHash[trial.getKey()] = trial.getFitness();

        des = trial;
        return;
    }

}

bool DSMGA2::restrictedMixing(Chromosome& ch, list<int>& mask) {

    bool taken = false;
    size_t lastUB = 0;

    for (size_t ub = 1; ub <= mask.size(); ++ub) {

        size_t size = 1;
        Chromosome trial(ell);
        trial = ch;

        for (list<int>::iterator it = mask.begin(); it != mask.end(); ++it) {

            trial.flip(*it);

            ++size;
            if (size > ub) break;
        }

        if (isInP(trial)) continue;     // is in population (pHash)


        if (trial.getFitness() >= ch.getFitness()) {
            pHash.erase(ch.getKey());
            pHash[trial.getKey()] = trial.getFitness();

            taken = true;
            ch = trial;
        }

        if (taken) {
            lastUB = ub;
            break;
        }
    }

    if (lastUB != 0) {
        while (mask.size() > lastUB)
            mask.pop_back();
    }

    return taken;

}

size_t DSMGA2::findSize(Chromosome& ch, list<int>& mask) const {

    DLLA candidate(nCurrent);
    for (int i=0; i<nCurrent; ++i)
        candidate.insert(i);

    size_t size = 0;
    for (list<int>::iterator it = mask.begin(); it != mask.end(); ++it) {

        int allele = ch.getVal(*it);

        for (DLLA::iterator it2 = candidate.begin(); it2 != candidate.end(); ++it2) {
            if (population[*it2].getVal(*it) == allele)
                candidate.erase(*it2);

            if (candidate.isEmpty())
                break;
        }

        if (candidate.isEmpty())
            break;

        ++size;
    }

    return size;


}

size_t DSMGA2::findSize(Chromosome& ch, list<int>& mask, Chromosome& ch2) const {

    size_t size = 0;
    for (list<int>::iterator it = mask.begin(); it != mask.end(); ++it) {
        if (ch.getVal(*it) == ch2.getVal(*it)) break;
        ++size;
    }
    return size;
}

void DSMGA2::mixing() {

    if (SELECTION)
        selection();

    //* really learn model
    buildFastCounting();
    buildGraph();
    // ================================================
    for (int i=0; i<ell; ++i) {
        findClique(i, masks[i]);        // generate connection priority for each bit
    }

    int repeat = (ell>50)? ell/50: 1;

    for (int k=0; k<repeat; ++k) {

        genOrderN();
        for (int i=0; i<nCurrent; ++i) {
            restrictedMixing(population[orderN[i]]);
            if (Chromosome::hit) break;
        }
        if (Chromosome::hit) break;
    }


}

inline bool DSMGA2::isInP(const Chromosome& ch) const {

    unordered_map<unsigned long, double>::const_iterator it = pHash.find(ch.getKey());
    return (it != pHash.end());
}

inline void DSMGA2::genOrderN() {
    myRand.uniformArray(orderN, nCurrent, 0, nCurrent-1);
}

inline void DSMGA2::genOrderELL() {
    myRand.uniformArray(orderELL, ell, 0, ell-1);
}

void DSMGA2::buildGraph() {

    int *one = new int [ell];
    char* key = new char[10];

    for (int i = 0; i < ell; ++i) {
        one[i] = countOne(i);
    }

    for (int i=0; i<ell; ++i) {

        for (int j=i+1; j<ell; ++j) {

            int n00x, n01x, n10x, n11x;
            int nijX = countXOR(i, j);

            n11x = (one[i] + one[j] - nijX) / 2;
            n10x = one[i] - n11x;
            n01x = one[j] - n11x;
            n00x = nCurrent - n01x - n10x - n11x;

            double p00x = (double)n00x/(double)nCurrent;
            double p01x = (double)n01x/(double)nCurrent;
            double p10x = (double)n10x/(double)nCurrent;
            double p11x = (double)n11x/(double)nCurrent;

            double linkage2 = computeMI(p00x,p01x,p10x,p11x);
            graph.write(i,j,linkage2);

            for (int k = j + 1; k < ell; ++k) {
                int n1x1, nx11;
                int nikX = countXOR(i, k);
                int njkX = countXOR(j, k);
                n1x1 = (one[i] + one[k] - nikX) / 2;
                nx11 = (one[j] + one[k] - njkX) / 2;

                int n000, n001, n010, n011, n100, n101, n110, n111;
                n111 = countAND3(i, j, k);
                n110 = n11x - n111;
                n101 = n1x1 - n111;
                n011 = nx11 - n111;
                n100 = n10x - n101;
                n010 = n01x - n011;
                n001 = one[k] - n1x1 - n011;
                n000 = n00x - n001;

                if (n000 + n001 + n010 + n011 + n100 + n101 + n110 + n111 != nCurrent)
                    cout << "Total Sum Error!" << endl;

                vector<double> p(8);
                p[0] = 1.0 * n000 / nCurrent;
                p[1] = 1.0 * n001 / nCurrent;
                p[2] = 1.0 * n010 / nCurrent;
                p[3] = 1.0 * n011 / nCurrent;
                p[4] = 1.0 * n100 / nCurrent;
                p[5] = 1.0 * n101 / nCurrent;
                p[6] = 1.0 * n110 / nCurrent;
                p[7] = 1.0 * n111 / nCurrent;
                double linkage3 = compute_three_mi_predict(p);

                sprintf(key, "%03d%03d%03d", i, j, k);
                //cout << key << endl;
                threeMI[key] = linkage3;
            }
        }
    }


    delete []one;
    delete [] key;
}

// from 1 to ell, pick by max edge
void DSMGA2::findClique(int startNode, list<int>& result) {


    result.clear();

    DLLA rest(ell);         // why using DOUBLE link list instead of linked list or just vector
    genOrderELL();
    for (int i=0; i<ell; ++i) {
        if (orderELL[i]==startNode)
            result.push_back(orderELL[i]);
        else
            rest.insert(orderELL[i]);
    }

    double *connection = new double[ell];

    for (DLLA::iterator iter = rest.begin(); iter != rest.end(); ++iter)
        connection[*iter] = graph(startNode, *iter);

    bool second_pushed = false;    // pushed second element into result
    bool third_pushed = false;
    while (!rest.isEmpty()) {
        double max = -INF;
        int index = -1;

        if (second_pushed && !third_pushed) {
            char* key = new char[10];
            int first, second;
            double temp_max = -INF;
            int temp_index = -1;

            if (result.front() < result.back()) {
                first = result.front();
                second = result.back();
            } else {
                first = result.back();
                second = result.front();
            }

            for (DLLA::iterator iter = rest.begin(); iter != rest.end(); ++iter) {
                if (*iter < first)
                    sprintf(key, "%03d%03d%03d", *iter, first, second);
                else if (*iter > first && *iter < second)
                    sprintf(key, "%03d%03d%03d", first, *iter, second);
                else
                    sprintf(key, "%03d%03d%03d", first, second, *iter);
                if (max < threeMI[key]) {
                    max = threeMI[key];
                    index = *iter;
                }
                if (temp_max < connection[*iter]) {
                    temp_max = connection[*iter];
                    temp_index = *iter;
                }
            }

            cout << "ThreeMI Index: " << index << endl;
            cout << "AddMI Index: " << temp_index << endl;

            if (index != temp_index) {
                cout << "ThreeMI!" << endl;
            }

            third_pushed = true;
            delete [] key;
        } else {
            for (DLLA::iterator iter = rest.begin(); iter != rest.end(); ++iter) {
                if (max < connection[*iter]) {
                    max = connection[*iter];
                    index = *iter;
                }
            }
            second_pushed = true;
        }

        rest.erase(index);
        result.push_back(index);

        for (DLLA::iterator iter = rest.begin(); iter != rest.end(); ++iter)
            connection[*iter] += graph(index, *iter);
    }


    delete []connection;

}

double DSMGA2::computeMI(double a00, double a01, double a10, double a11) const {

    double p0 = a00 + a01;
    double q0 = a00 + a10;
    double p1 = 1-p0;
    double q1 = 1-q0;

    double join = 0.0;
    if (a00 > EPSILON)
        join += a00*log(a00);
    if (a01 > EPSILON)
        join += a01*log(a01);
    if (a10 > EPSILON)
        join += a10*log(a10);
    if (a11 > EPSILON)
        join += a11*log(a11);

    double p = 0.0;
    if (p0 > EPSILON)
        p += p0*log(p0);
    if (p1 > EPSILON)
        p += p1*log(p1);


    double q = 0.0;
    if (q0 > EPSILON)
        q += q0*log(q0);
    if (q1 > EPSILON)
        q += q1*log(q1);

    double mi = -p-q+join;

    return mi;
}

void DSMGA2::reduce_p3(vector<double> p, vector<double>& p1, vector<double>& p2) {
    if ((p1.size() != 6) || (p2.size() != 12))
        cout << "length error in reduce_p3!" << endl;
    p2[0]  = p[0] + p[1];
    p2[1]  = p[2] + p[4];
    p2[2]  = p[3] + p[6];
    p2[3]  = p[4] + p[8];
    p2[4]  = p[0] + p[2];
    p2[5]  = p[1] + p[4];
    p2[6]  = p[3] + p[5];
    p2[7]  = p[6] + p[7];
    p2[8]  = p[0] + p[3];
    p2[9]  = p[1] + p[6];
    p2[10] = p[2] + p[5];
    p2[11] = p[4] + p[7];
    p1[0] = p2[0] + p2[1];
    p1[2] = p2[0] + p2[2];
    p1[4] = p2[4] + p2[6];
    for (int i = 1; i < 6; i += 2)
        p1[i] = 1 - p1[i - 1]; 
}

double DSMGA2::compute_h(vector<double> p) {
    double ans = 0.0;
    for (size_t i = 0; i < p.size(); ++i) {
        if (p[i] > EPSILON)
            ans -= p[i] * log(p[i]);
    }
    return ans;
}

double DSMGA2::compute_three_mi(vector<double> p) {
    if (p.size() != 8)
        cout << "The length of input in compute_three_mi is not 8!" << endl;
    vector<double> p2(12);
    vector<double> p1(6);
    reduce_p3(p, p1, p2);
    return compute_h(p1) - compute_h(p2) + compute_h(p);
}

double DSMGA2::compute_three_mi_predict(vector<double> p) {
    vector<double> p2(12);      // probability of 00~11 of bit 12, 13, 23
    vector<double> p1(6);       // probability of bit 1, 2, 3 to be 0 or 1
    vector<double> high(8);     // upper bound
    vector<double> low(8);      // lower bound
    vector<double> predict(8);  // prediction based on upper and lower bound
    vector<double> candi(8);    // candidate

    reduce_p3(p, p1, p2);       // compute p1, p2 given p3

    candi[0] = 1 - p1[1] - p1[3] - p1[5] + p2[3] + p2[7] + p2[11];
    high[0] = min(min(min(p2[0], p2[4]), p2[8]), candi[0]);
    candi[1] = 1 - p1[1] - p1[3] - p1[4] + p2[3] + p2[6] + p2[10];
    high[1] = min(min(min(p2[0], p2[5]), p2[9]), candi[1]);
    candi[2] = 1 - p1[1] - p1[2] - p1[5] + p2[2] + p2[7] + p2[9];
    high[2] = min(min(min(p2[1], p2[4]), p2[10]), candi[2]);
    candi[3] = 1 - p1[1] - p1[2] - p1[4] + p2[2] + p2[6] + p2[8];
    high[3] = min(min(min(p2[1], p2[5]), p2[11]), candi[3]);
    candi[4] = 1 - p1[0] - p1[3] - p1[5] + p2[1] + p2[5] + p2[11];
    high[4] = min(min(min(p2[2], p2[6]), p2[8]), candi[4]);
    candi[5] = 1 - p1[0] - p1[3] - p1[4] + p2[1] + p2[4] + p2[10];
    high[5] = min(min(min(p2[2], p2[7]), p2[9]), candi[5]);
    candi[6] = 1 - p1[0] - p1[2] - p1[5] + p2[0] + p2[5] + p2[9];
    high[6] = min(min(min(p2[3], p2[6]), p2[10]), candi[6]);
    candi[7] = 1 - p1[0] - p1[2] - p1[4] + p2[0] + p2[4] + p2[8];
    high[7] = min(min(min(p2[3], p2[7]), p2[11]), candi[7]);

    low[0] = max(max(p2[0] - high[1] , p2[4] - high[2]) , p2[8] - high[4]);  //000
    low[1] = max(max(p2[0] - high[0] , p2[5] - high[3]) , p2[9] - high[5]);  //001
    low[2] = max(max(p2[1] - high[3] , p2[4] - high[0]) , p2[10] - high[6]); //010
    low[3] = max(max(p2[1] - high[2] , p2[5] - high[1]) , p2[11] - high[7]); //011
    low[4] = max(max(p2[2] - high[5] , p2[6] - high[6]) , p2[8] - high[0]);  //100
    low[5] = max(max(p2[2] - high[4] , p2[7] - high[7]) , p2[9] - high[1]);  //101
    low[6] = max(max(p2[3] - high[7] , p2[6] - high[4]) , p2[10] - high[2]); //110
    low[7] = max(max(p2[3] - high[6] , p2[7] - high[5]) , p2[11] - high[3]); //111

    bool ppp = false;
    for (int i = 0; i < 8; ++i) {
        if (high[i] < p[i]) {
            cout << "high error in index " << i << endl;
            ppp = true;
        }
        if (low[i] > p[i]) {
            cout << "low error in index " << i << endl;
            ppp = true;
        }
        predict[i] = (high[i] + low[i]) / 2;
    }
    if (ppp)
        for (int j = 0; j < 8; ++j)
            cout << p[j] << endl;

    return compute_three_mi(predict);     
}

void DSMGA2::selection () {
    tournamentSelection ();
}

// tournamentSelection without replacement
void DSMGA2::tournamentSelection () {
    int i, j;

    int randArray[selectionPressure * nCurrent];

    for (i = 0; i < selectionPressure; i++)
        myRand.uniformArray (randArray + (i * nCurrent), nCurrent, 0, nCurrent - 1);

    for (i = 0; i < nCurrent; i++) {

        int winner = 0;
        double winnerFitness = -INF;

        for (j = 0; j < selectionPressure; j++) {
            int challenger = randArray[selectionPressure * i + j];
            double challengerFitness = population[challenger].getFitness();

            if (challengerFitness > winnerFitness) {
                winner = challenger;
                winnerFitness = challengerFitness;
            }

        }
        selectionIndex[i] = winner;     // set of winners of each selection
    }
}
