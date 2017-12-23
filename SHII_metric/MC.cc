#include <iomanip>
#include "MC.h"

namespace _MC {
    
    MC::MC(AnyOption* opt1) {
        opt = opt1;
        
        cout << "In testing phase" << endl;
        
        outdir = opt->getValue("outdir");
        probGraphFile = opt->getValue("probGraphFile");
        string seedFile = opt->getValue("seedFile");
        //eta = strToInt(opt->getValue("eta"));
        //    graphType = UNDIRECTED;
        graphType = DIRECTED;
        //flag_in = opt->getValue("");
        
        setModel(opt);
        
        AM = new HashTreeCube(89041);
        revAM = NULL;
        
        countIterations = strToInt(opt->getValue("mcruns"));
        countRandSeedIterations = strToInt(opt->getValue("randSeedIter"));
        seedRatio = strToFloat(opt->getValue("seedRatio"));
        //    edgeWeight = strToFloat(opt->getValue("edgeWeight"));
        //	testingActionsFile = (const char*) opt->getValue("testingActionsFile");
        
        cout << "User specified options: " << endl;
        cout << "model : " << m << " or " << model << endl;
        cout << "outdir : " << outdir << endl;
        cout << "probGraphFile : " << probGraphFile << endl;
        cout << "seedFile : " << seedFile << endl;
        cout << "Number of iterations in MC : " << countIterations << endl;
        
        time (&startTime);
        srand ( startTime );
        time (&stime_mintime);
        
    }
    
    MC::~MC() {
        cout << "total time taken : " << getTime() << endl;
    }
    
    
    void MC::setModel(AnyOption* opt) {
        m = opt->getValue("propModel");
        
        model = LT;
        
        if (m.compare("LT") == 0) {
            model = LT;
            
        } else if (m.compare("IC") == 0) {
            model = IC;
        }
    }
    
    
    
    void MC::doAll() {
        
        readInputData();
        
        cout << "Model " << m << " : " << model << " running" << endl;
        
        openCoverageFiles();
        computeCov();
        outFile.close();
        
    }
    
    void MC::computeCov() {
        // read the seed set from the file
        // for each radius R, generate the table
        // <S_alg, R, cov^R>
        UserList seedSetInCommunity;
        vector<UID> nodesInCommunity;
        //    curSeedSet.clear();
        //    vector<UID> seedsVec;
        //    seedsVec.clear();
        for (unsigned int community_label = 1; community_label<= communitySize; community_label++) {
            nodesInCommunity.clear();
            seedSetInCommunity.clear();
            for (CommunityMap::const_iterator it = nodeLabels.begin(); it != nodeLabels.end(); ++it)
            {
                // if the node belong to the current community and also the seeds in the file
                // put it to the seedset in community @luc 12/07/2015
                if (it->second == community_label)
                {
                    nodesInCommunity.push_back(it->first);
                }
                if (it->second == community_label && curSeedSet.find(it->first) != curSeedSet.end())
                {
                    seedSetInCommunity.insert(it->first);
                }
            }
            unsigned int seedSetSize = int(ceil(float(nodesInCommunity.size()) *seedRatio));
            if(seedSetInCommunity.empty()){
                continue;
            }
            // print out the seeds
//            for (set<UID>::iterator it = seedSetInCommunity.begin(); it != seedSetInCommunity.end(); ++it) {
//                UID v = *it;
//                cout << v << " ";
//            }
//            cout << endl;
        
            // use only one SH spanner and a set of random nodes that are in the same community as the seed set
            // for running information diffuision models @luc 1/3/2016
            
            for (set<UID>::iterator it = seedSetInCommunity.begin(); it != seedSetInCommunity.end(); ++it) {
                UID v = *it;
                UserList oneSeedSet;
                UserList seedNeighborSet; // sample seeds from seedNeighbors
                queue<UID> Q;
                //            float cov = 0;
                // AM is adjacency matrix
                
                // using BFS to add neighbors of the SH spanner to the seedNeighborSet as seed candidates @luc 1/3/2016
                oneSeedSet.insert(v);
                Q.push(v);
                while(Q.empty() == false) {
                    UID q = Q.front();
                    FriendsMap* neighbors = AM->find(q);
                    if (neighbors != NULL) {
                        for (FriendsMap::iterator j = neighbors->begin(); j!=neighbors->end(); ++j) {
                            UID u = j->first;
                            seedNeighborSet.insert(u);
                        }
                        // if the number of seedNeighborSet plus the number of current seed set
                        // is smaller than the predefined setset size,
                        // add all the neighbors to the seed set and push each neighbor to the queue for BFS
                        // @luc 1/3/2016
                        if (Q.size() == 1 && (oneSeedSet.size() + seedNeighborSet.size()  < seedSetSize)){
                            for (set<UID>::iterator j = seedNeighborSet.begin(); j != seedNeighborSet.end(); ++j) {
                                UID u = *j;
                                oneSeedSet.insert(u);
                                Q.push(u);
                            }
                            seedNeighborSet.clear();
                        }
                    }
                    Q.pop();
                }
                vector<UID> seedNeighbors(seedNeighborSet.begin(), seedNeighborSet.end());
                
                // for (vector<UID>::iterator k = seedNeighbors.begin(); k != seedNeighbors.end(); ++k) {
                //     cout << *k <<" ";
                // }
                // cout << endl;
                
                std::pair<double, double> avg_censor_score = std::make_pair(0,0);
                for (int randIter=0; randIter < countRandSeedIterations; randIter++) {
                    
                    UserList randSeedSet(oneSeedSet);
                    vector<UID> seedNeighbors(seedNeighborSet.begin(), seedNeighborSet.end());
                    // sample seeds from seedNeighborSet @luc 1/3/2016
                    while (!seedNeighbors.empty() && randSeedSet.size() < seedSetSize){
                        int r = (unsigned int)((float)seedNeighbors.size() * (((float)(rand() % 1001))/(float)1000));
                        if (r>= seedNeighbors.size()) continue;
                        UID u = seedNeighbors[r];
                        randSeedSet.insert(u);
                        seedNeighbors.erase(seedNeighbors.begin()+r);
                    }
                    
                    std::pair<double, double> censor_score = std::make_pair(0,0);
                    if (model == IC) {
                        if (randIter%100==0) {
                            cout << "IC model community " << community_label << " node "<< v << " for " << randIter <<" iterations" << endl;
                            // print out the rand seeds and their communities @luc 1/3/2016
//                            for (set<UID>::iterator k = randSeedSet.begin(); k != randSeedSet.end(); ++k) {
//                                cout << *k << " " << nodeLabels[*k] << endl;
//                            }
                        }
                        censor_score = ICCov(randSeedSet, community_label);
                    } else if (model == LT) {
                        if (randIter%100==0) {
                            cout << "LT model community " << community_label << " node "<< v << " for " << randIter <<" iterations" << endl;
                             // print out the rand seeds and their communities @luc 1/3/2016
//                            for (set<UID>::iterator k = randSeedSet.begin(); k != randSeedSet.end(); ++k) {
//                                cout << *k << " " << nodeLabels[*k] << endl;
//                            }
                        }
                        censor_score = LTCov(randSeedSet, community_label);
                    }
                    avg_censor_score.first += censor_score.first/countRandSeedIterations;
                    avg_censor_score.second += censor_score.second/countRandSeedIterations;
                    
                }
                outFile << avg_censor_score.first << " " << avg_censor_score.second << endl;
            }
        }
        
        //	outFile.close();
        
    }
    
    
    std::pair<double, double> MC::ICCov(UserList& S, unsigned int community_label) {
        //	double cov = 0;
        
        // initialize ppIn
        //    map<UID, float> ppIn; // what is the prob with which the node is covered
        
        double avg_result_1 = 0;
        double avg_result_2 = 0;
        
        for (int b = 0; b < countIterations; ++b) {
            // Q is the queue in the depth/breadth first search
            queue<UID> Q;
            // activeNodes is the set of nodes that are activated in the current
            // run of Monte Carlo simulation
            UserList activeNodes;
            
            // S is the seed set
            // for each seed node v is S,
            // add it to activeNodes
            // add it to Q as well
            for (UserList::iterator i=S.begin(); i!=S.end(); ++i) {
                UID v = *i;
                //            if (nodeLabels[v] != community_label){
                //                cout << "ERROR1: " << v ;
                //                continue;
                //            }
                Q.push(v);
                activeNodes.insert(v);
            }
            
            while(Q.empty() == false) {
                UID v = Q.front();
//            if (nodeLabels[v] != community_label){
//                cout << "ERROR2: " << v ;
//                continue;
//            }
                // AM is adjacency matrix
                FriendsMap* neighbors = AM->find(v);
                if (neighbors != NULL) {
                    for (FriendsMap::iterator j = neighbors->begin(); j!=neighbors->end(); ++j) {
                        UID u = j->first;
                        
                        // if not actived
                        if (activeNodes.find(u) == activeNodes.end()) {
                            float toss = ((float)(rand() % 1001))/(float)1000;
                            float p = j->second;
                            
                            if (p >= toss) {
                                activeNodes.insert(u);
                                // continue to propagate if the influenced node
                                // belongs to the same community as seed nodes @luc 12/06/2015
//                            if (nodeLabels[u] == community_label) {
//                                Q.push(u);
//                            }
                                Q.push(u);
                            }
                        }
                        
                    }
                }
                
                Q.pop();
                
            }
            // compute the censorship score for this iteration
            double self_cov = 0;
            double total_cov = 0;
            for (UserList::iterator i=activeNodes.begin(); i!=activeNodes.end(); ++i) {
                UID v = *i;
                total_cov += 1;
                
                if (nodeLabels[v] == community_label){
                    self_cov += 1;
                }
            }
            double censor_score_1 = (total_cov - self_cov)/(total_cov);
            double censor_score_2 = (total_cov - self_cov)/(self_cov);
            
            //        double censor_score_1 = (self_cov)/(total_cov);
            //        double censor_score_2 = (self_cov)/(total_cov - self_cov);
            
            avg_result_1 += censor_score_1/countIterations;
            avg_result_2 += censor_score_2/countIterations;
        }
        
        //    double self_cov = 0;
        //    double total_cov = 0;
        //    for (map<UID, float>::iterator i = ppIn.begin(); i!=ppIn.end(); ++i) {
        //        double avg_influenced = (double) i->second/countIterations;
        //        total_cov += avg_influenced;
        //
        //        if (nodeLabels[i->first] == community_label){
        //            self_cov += avg_influenced;
        //        }
        //        //        cout << "ppIN: " << i->first << " " << i->second/countIterations << endl;
        //        //        outFile << i->second/countIterations << " ";
        //    }
        return std::make_pair(avg_result_1, avg_result_2);
        
        
        // compute two things: cov(S) and cov(S+x)
        //	cov = cov/countIterations;
        
        
    }
    
    void MC::clear() {
        covBestNode.clear();
        curSeedSet.clear();
        seedSetNeighbors.clear();
        nodeLabels.clear();
        totalCov = 0;
    }
    
    
    
    std::pair<double, double> MC::LTCov(UserList& S, unsigned int community_label) {
        //    cout << "In LTCov" << endl;
        float tol = 0.00001;
        float cov = 0;
        double avg_result_1 = 0;
        double avg_result_2 = 0;
        // initialize ppIn
        //    map<UID, float> ppIn; // what is the prob with which the node is covered
        
        
        for (int b = 0; b < countIterations; ++b) {
            /* initialize random seed: */
            UserList activeNodes;
            float cov1 = 0;
            // T is the set of nodes that are to be processed
            queue<UID> T;
            // Q is the set of nodes that have been seen until now
            // Thus, Q is a superset of T
            // Q contains the nodes that are processed and to be processed
            map<UID, NodeParams> Q;
            
            //        cov += S.size();
            // cov1 == coverage in one run .. that is, number of nodes reachable
            cov1 += S.size();
            
            for (UserList::iterator i=S.begin(); i!=S.end(); ++i) {
                UID v = *i;
                activeNodes.insert(v);
                
                FriendsMap* neighbors = AM->find(v);
                if (neighbors != NULL) {
                    for (FriendsMap::iterator j = neighbors->begin(); j!=neighbors->end(); ++j) {
                        UID u = j->first;
                        
                        if (S.find(u) == S.end()) {
                            
                            if (Q.find(u) == Q.end()) {
                                // if the node u has not been seen before
                                // create a new NodeParams
                                // create a random threhsold
                                NodeParams& np = Q[u];
                                np.active = false;
                                np.inWeight = j->second;
                                
                                /* generate secret number: */
                                np.threshold = ((float)(rand() % 1001))/(float)1000;
                                T.push(u);
                            } else {
                                NodeParams& np = Q[u];
                                np.inWeight += j->second;
                            } //endif
                        } //endif
                    } //endfor
                } //endif
            } //endfor
            
            while (!T.empty()) {
                UID u = T.front();
                
                //            cout << "T.size " << T.size() << endl;
                
                NodeParams& np = Q.find(u)->second;
                if (np.active == false && np.inWeight >= np.threshold + tol) {
                    activeNodes.insert(u);
                    np.active = true;
                    //                cov++;
                    cov1++;
                    
                    // continue to propagate if the influenced node
                    // belongs to the same community as seed nodes @luc 12/06/2015
                    //                if (nodeLabels[u] != community_label)
                    //                    continue;
                    
                    // add u's neighbors to T
                    FriendsMap* neighbors = AM->find(u);
                    if (neighbors != NULL) {
                        // for each neighbor w of u
                        for (FriendsMap::iterator k = neighbors->begin(); k!=neighbors->end(); ++k) {
                            UID w = k->first;
                            // is w is in S, no need to do anything
                            if (S.find(w) != S.end()) continue;
                            
                            // if w is not in S, locate it in Q
                            map<UID, NodeParams>::iterator it = Q.find(w);
                            
                            if (it == Q.end()) {
                                // if it is not in Q, then
                                NodeParams& np_w = Q[w];
                                np_w.threshold = ((float)(rand() % 1001))/(float)1000;
                                //                            np_w.threshold = (float)rand()/RAND_MAX;
                                np_w.active = false;
                                np_w.inWeight = k->second;
                                T.push(w);
                            } else {
                                // if w is in Q, then
                                NodeParams& np_w = it->second;
                                if (np_w.active == false) {
                                    T.push(w);
                                    np_w.inWeight += k->second;
                                    // ignore the warning now, @luc 10/30/2015
                                    if (np_w.inWeight - 1 > tol) {
                                        cout << "Something wrong, the inweight for a node is > 1. (w, inweight) = " << w << ", " << np_w.inWeight - 1<< endl;
                                    }
                                }
                            }
                        }
                    }
                }
                
                // deletes the first element
                T.pop();
                
            } //endwhile
            //        cout << "Coverage in this iteration: " << cov1 << endl;
            //		cov += cov1/countIterations;
            // compute the censorship score for this iteration
            double self_cov = 0;
            double total_cov = 0;
            for (UserList::iterator i=activeNodes.begin(); i!=activeNodes.end(); ++i) {
                UID v = *i;
                // if active node is in the seed set, continue @luc 12/07/2015
                //            if( S.find(v) != S.end()){
                //                continue;
                //            }
                total_cov += 1;
                
                if (nodeLabels[v] == community_label){
                    self_cov += 1;
                }
            }
            
            //        double alpha = 1 / double(communitySize);
            double censor_score_1 = (total_cov - self_cov)/(total_cov);
            double censor_score_2 = (total_cov - self_cov)/(self_cov);
            //        double censor_score_1 = (self_cov)/(total_cov);
            //        double censor_score_2 = (self_cov)/(total_cov - self_cov);
            
            avg_result_1 += censor_score_1/countIterations;
            avg_result_2 += censor_score_2/countIterations;
        }
        return std::make_pair(avg_result_1, avg_result_2);
        
        // coverage from ppIn
        //    double self_cov = 0;
        //    double total_cov = 0;
        //    for (map<UID, float>::iterator i = ppIn.begin(); i!=ppIn.end(); ++i) {
        //        double avg_influenced = (double) i->second/countIterations;
        //        total_cov += avg_influenced;
        //
        //        if (nodeLabels[i->first] == community_label){
        //            self_cov += avg_influenced;
        //        }
        ////        cout << "ppIN: " << i->first << " " << i->second/countIterations << endl;
        ////        outFile << i->second/countIterations << " ";
        //    }
        //    outFile << endl;
        
        //    double retCov = (double) cov/countIterations;
        
        //    cout << "(cov, retCov, cov1) = " << cov << ", " << retCov << ", " << cov1 << endl;
        
        //    return self_cov/total_cov;
        
        
    }
    
    
    void MC::printVector(vector<UID>& vec, float pp) {
        cout << "AMIT " << pp << " " ;
        for (vector<UID>::iterator i=vec.begin(); i!=vec.end(); ++i) {
            
            cout << *i << " ";
        }
        
        cout << endl;
        
    }
    
    void MC::writeInFile(UID v, float cov, float marginal_gain, int curTimeStep, float actualCov, float actualMG, int countUsers) {
        cout << endl << endl << "Picked a seed node: " << v << ", total: " << curSeedSet.size() << endl;
        outFile << v << " " << cov << " " << marginal_gain << " " << curTimeStep << " " << getCurrentMemoryUsage() << " " << getTime() <<  " " << getTime_cur() << " " << actualCov << " " << actualMG << " " << countUsers << endl;
        cout << v << " " << cov << " " << marginal_gain << " " << curTimeStep << " " << getCurrentMemoryUsage() << " " << getTime() << " " << getTime_cur() <<  " " << actualCov << " " << actualMG << " " << countUsers << endl;
        cout << endl << endl;
    }
    
    void MC::openCoverageFiles() {
        
        if (outFile.is_open()) {
            outFile.close();
        }
        string datasetFile = opt->getValue("datasetFile");
        string m = opt->getValue("propModel");
        //    string filename = opt->getValue("dataName");
        
        string filename = outdir + "/" + datasetFile +"_" + m + ".txt";
        
        outFile.open (filename.c_str());
        
        if (outFile.is_open() == false) {
            cout << "Can't open file " << filename << " for writing" << endl;
            exit(1);
        }
    }
    
    
    
    
    PropModels MC::getModel() {
        return model;
    }
    
    float MC::getTime_cur() const {
        time_t curTime;
        time(&curTime);
        
        float min = ((float)(curTime - stime_mintime))/60;
        return min;
    }
    
    float MC::getTime() const {
        time_t curTime;
        time(&curTime);
        
        float min = ((float)(curTime - startTime))/60;
        return min;
    }
    
    
    void MC::readInputData(float alpha) {
        cout << "in readInputData for model " << model << " with alpha " << alpha << endl;
        
        unsigned int edges = -1;
        unsigned int numUsers = 0;
        string probGraphFile = opt->getValue("probGraphFile");
        cout << "Reading file " << probGraphFile << endl;
        ifstream myfile (probGraphFile.c_str(), ios::in);
        string delim = " \t";	
        if (myfile.is_open()) {
            while (! myfile.eof() )	{
                std::string line;
                getline (myfile,line);
                if (line.empty()) continue;
                edges++;
                std::string::size_type pos = line.find_first_of(delim);
                int	prevpos = 0;
                
                // skip the first three lines in mtx file format
                if (edges <= 2){
                    //                string str = line.substr(prevpos, pos-prevpos);
                    //                numUsers= strToInt(str); // ignore the first line
                    continue;
                }
                // get first user
                string str = line.substr(prevpos, pos-prevpos);
                UID u1 = strToInt(str);
                
                // get the second user
                prevpos = line.find_first_not_of(delim, pos);
                pos = line.find_first_of(delim, prevpos);
                UID u2 = strToInt(line.substr(prevpos, pos-prevpos));
                
                if (u1 == u2) continue;
                
                // get the parameter
                float parameter1 = 0;
                prevpos = line.find_first_not_of(delim, pos);
                pos = line.find_first_of(delim, prevpos);
                if (pos == string::npos) 
                    parameter1 = strToFloat(line.substr(prevpos));
                else
                    parameter1 = strToFloat(line.substr(prevpos, pos-prevpos));
                
                if (parameter1 == 0) continue;
                
                //          CHANGE THE EDGE WEIGHT TO 0.5
                parameter1 = parameter1 * edgeWeight;
                //			users.insert(u1);
                //			users.insert(u2);
                
                
                if (edges % 1000000 == 0) {
                    cout << "(node1, node2, weight,  AM size till now, edges till now, mem) = " << u1 << ", " << u2 << ", " << parameter1 << ", " << AM->size() << ", " << edges << ", " << getCurrentMemoryUsage() << endl;
                }
                
                
                FriendsMap* neighbors = AM->find(u1);
                if (neighbors == NULL) {
                    neighbors = new FriendsMap();
                    neighbors->insert(pair<UID, float>(u2, parameter1));
                    AM->insert(u1, neighbors);
                } else {
                    //FriendsMap::iterator it = neighbors->find(u2);
                    //if (it == neighbors->end()) {
                    neighbors->insert(pair<UID, float>(u2, parameter1));
                    //} else {
                    //	cout << "WARNING: Edge redundant between users " << u1 << " and " << u2 << endl;
                    //}
                }
                
                if (revAM != NULL) {
                    multimap<float, UID> *revNeighbors = revAM->find(u1);
                    if (revNeighbors == NULL) {
                        revNeighbors = new multimap<float, UID>();
                        revNeighbors->insert(std::make_pair(parameter1, u2));
                        revAM->insert(u1, revNeighbors);
                    } else {
                        revNeighbors->insert(std::make_pair(parameter1, u2));
                    }
                }
                /*
                 if (AM_in != NULL) {
                 FriendsMap *inNeighbors = AM_in->find(u2);
                 if (inNeighbors == NULL) {
                 inNeighbors = new FriendsMap();
                 inNeighbors->insert(std::make_pair(u1, parameter1));
                 AM_in->insert(u2, inNeighbors);
                 } else {
                 FriendsMap::iterator it = inNeighbors->find(u1);
                 if (it == inNeighbors->end()) {
                 inNeighbors->insert(std::make_pair(u1, parameter1));
                 }
                 }
                 }*/
                
                // also add the edges u2->u1 but done allocate Edge class to them
                // .. it is just to find friends efficiently
                if (graphType == UNDIRECTED) { 
                    neighbors = AM->find(u2);
                    if (neighbors == NULL) {
                        neighbors = new FriendsMap();
                        neighbors->insert(pair<UID, float>(u1, parameter1 ));
                        AM->insert(u2, neighbors);
                    } else {
                        FriendsMap::iterator it = neighbors->find(u1);
                        if (it == neighbors->end()) {
                            neighbors->insert(pair<UID, float>(u1, parameter1 ));
                        } else {
                            //	cout << "WARNING: Edge redundant between users " << u1 << " and " << u2 << endl;
                        }
                    }
                }
            }
            myfile.close();
        } else {
            cout << "Can't open friendship graph file " << probGraphFile << endl;
        }
        
        // read labels for each UID
        string labelFile = opt->getValue("labelFile");
        cout << "Reading file " << labelFile << endl;
        myfile.open (labelFile.c_str(), ios::in);
        unsigned int community_label = 1;
        if (myfile.is_open()) {
            while (! myfile.eof() )	{
                std::string line;
                getline (myfile,line);
                if (line.empty()) continue;
                std::string::size_type pos = line.find_first_of(delim);
                int	prevpos = 0;
                // get first user
                UID u = strToInt(line.substr(prevpos, pos-prevpos));
                nodeLabels.insert(pair<UID, unsigned int>(u, community_label));
                users.insert(u);
                
                while (pos!=std::string::npos)
                {
                    // get the next user
                    prevpos = line.find_first_not_of(delim, pos);
                    pos = line.find_first_of(delim, prevpos);
                    UID u = strToInt(line.substr(prevpos, pos-prevpos));
                    nodeLabels.insert(pair<UID, unsigned int>(u, community_label));
                    users.insert(u);
                }
                community_label += 1;
            }
            myfile.close();
        }
        
        // read seeds
        // note that seed start index from 1 @luc 12/07/2015
        string seedFile = opt->getValue("seedFile");
        curSeedSet.clear();
        cout << "Reading file " << seedFile << endl;
        
        myfile.open (seedFile.c_str(), ios::in);
        if (myfile.is_open()) {
            while (! myfile.eof() )	{
                std::string line;
                getline (myfile,line);
                if (line.empty()) continue;
                std::string::size_type pos = line.find_first_of(delim);
                int	prevpos = 0;
                // get first user
                UID u = strToInt(line.substr(prevpos, pos-prevpos));
                curSeedSet.insert(u);
                
                while (pos!=std::string::npos)
                {
                    // get the next user
                    prevpos = line.find_first_not_of(delim, pos);
                    pos = line.find_first_of(delim, prevpos);
                    UID u = strToInt(line.substr(prevpos, pos-prevpos));
                    curSeedSet.insert(u);
                }
            }
            myfile.close();
        }
        
        
        
        // generate weight as inverse of indegree for LT model @luc
        if (model==LT) {
            for (UserList::iterator i=  users.begin(); i!= users.end(); ++i) {
                UID u = *i;
                FriendsMap* neighbors = AM->find(u);
                if (neighbors != NULL) {
                    for (FriendsMap::iterator outedge = neighbors->begin(); outedge!=neighbors->end(); ++outedge) {
                        UID v = outedge->first;
                        FriendsMap* v_neighbors = AM->find(v);
                        // set weight
                        outedge->second = 1 / float(v_neighbors->size());
                    }
                }
            }
        }
        communitySize = community_label-1;
        
        //    for (unsigned int community_label = 1; community_label<= communitySize; community_label++) {
        //        for (CommunityMap::const_iterator it = nodeLabels.begin(); it != nodeLabels.end(); ++it)
        //        {
        //            if (it->second == community_label)
        //            {
        //                cout << it->first << " ";
        //            }
        //        }
        //        cout << endl;
        //    }
        
        this->numEdges = edges;
        
        cout << "BuildAdjacencyMatFromFile done" << endl;
        cout << "Size of friendship graph hashtree is : " << AM->size() << endl;
        //cout << "Size of NEW friendship graph hashtree is : " << revAM->size() << endl;
        cout << "Number of users are: " << users.size() << endl;
        cout << "Number of communities are: " << communitySize << endl;
        cout << "Number of edges in the friendship graph are: " << edges << endl;
    }
    
    
    
}
