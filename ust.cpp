// UST
// --- VERSION 1.0 ----
// Author: amatur, Last Edited: Nov 1
#include <cmath>
#include<cstring>
#include <fstream>
#include <iostream>
#include <vector>
#include <assert.h>
#include <stdint.h>
#include <unordered_set>
#include <map>
#include <set>
#include <sstream>
#include <algorithm>
#include <cstdlib>
#include <list>
#include <stack>
#include <unordered_map>
#include <utility>
#include <queue>
#include <deque>
#include <unistd.h>
#include <tuple>
using namespace std;

bool DEBUGMODE = false;
int K = 0;
string UNITIG_FILE = "examples/k11.unitigs.fa";
string OUTPUT_FILENAME;

enum DEBUGFLAG_T { NONE = 0,  UKDEBUG = 6, VERIFYINPUT = 1, INDEGREEPRINT = 2, DFSDEBUGG = 3, PARTICULAR = 4, NODENUMBER_DBG = 5, OLDNEWMAP = 9, PRINTER = 10, SINKSOURCE = 12};

enum ALGOMODE_T { BASIC = 0, INDEGREE_DFS = 1, INDEGREE_DFS_1 = 2, OUTDEGREE_DFS = 3, OUTDEGREE_DFS_1 = 4, INDEGREE_DFS_INVERTED = 5, PLUS_INDEGREE_DFS = 6, RANDOM_DFS = 7, NODEASSIGN = 8, SOURCEFIRST = 9, TWOWAYEXT = 10, PROFILE_ONLY = 11, EPPRIOR=12, GRAPHPRINT = 13, TIGHTUB = 14, BRACKETCOMP = 15};

bool FLG_NEWUB = true;
bool FLG_ABUNDANCE = false;

DEBUGFLAG_T DBGFLAG = NONE; //NODENUMBER_DBG
ALGOMODE_T ALGOMODE = TWOWAYEXT;

string mapmode[] = {"basic", "indegree_dfs", "indegree_dfs_initial_sort_only", "outdegree_dfs", "outdegree_dfs_initial_sort_only", "inverted_indegree_dfs", "plus_indegree_dfs", "random_dfs", "node_assign", "source_first", "twoway", "profile_only", "endpoint_priority", "graph_print", "tight_ub", "tip"
};
string modefilename[] = {"Fwd", "indegree_dfs", "indegree_dfs_initial_sort_only", "outdegree_dfs", "outdegree_dfs_initial_sort_only", "inverted_indegree_dfs", "plus_indegree_dfs", "random_dfs", "node_assign", "source_first", "", "profile_only", "endpoint_priority", "graph_print", "tight_ub", "Tip"
};


typedef tuple<int,int,int, int> fourtuple; // uid, walkid, pos, isTip

bool sort_by_walkId (const fourtuple &lhs, const fourtuple &rhs){
    return get<1>(lhs) < get<1>(rhs);
}
bool sort_by_pos (const fourtuple &lhs, const fourtuple &rhs){
    return get<2>(lhs) < get<2>(rhs);
}
bool sort_by_tipstatus (const fourtuple &lhs, const fourtuple &rhs){
    return get<3>(lhs) < get<3>(rhs);
}

typedef struct {
    int serial = -1;
    int startPosWithKOverlap;
    int endPosWithKOVerlap;
    bool isWalkEnd = false;
    int pos_in_walk = -100;
    int finalWalkId = -1; // renders some walkId as invalid
    int isTip = 0;
} new_node_info_t;

typedef struct {
    int serial;
    int ln;
} unitig_struct_t;

typedef struct {
    //1 means +, 0 means -
    bool left;
    bool right;
    int toNode;
} edge_t;

typedef struct {
    edge_t edge;
    int fromNode;
} edge_both_t;

typedef struct {
    edge_t edge;
    int kmerStartIndex;
    int kmerEndIndex;
} newEdge_t;

int C_ustitch = 0;
int C_twoway_ustitch = 0;
int C_tip_ustitch = 0;
int V_ustitch = 0;
int V_twoway_ustitch = 0;
int V_tip_ustitch = 0;

int isolated_node_count = 0;
int sink_count = 0;
int source_count = 0;
int sharedparent_count = 0;
int sharparentCntRefined = 0;
int onecount = 0;

struct node_sorter {
    int node;
    int sortkey;
    //bool operator() (struct node_sorter  i, struct node_sorter  j) { return (i.sortkey<j.sortkey);}
};
bool sort_by_key (struct node_sorter i, struct node_sorter j) { return (i.sortkey<j.sortkey); }
bool sort_by_key_inverted (struct node_sorter i, struct node_sorter j) { return (i.sortkey>j.sortkey); }

int* global_indegree;
int* global_outdegree;
int* global_plusindegree;
int* global_plusoutdegree;
int* global_issinksource;
//int* global_priority;

map<pair <int, int>, int> inOutCombo;

vector<vector<edge_t> > adjList;
vector<vector<edge_t> > reverseAdjList;

vector<unitig_struct_t> unitigs;
map<int, string> newSequences;
map<int, string> newNewSequences; //int is the unitig id (old id)
set<int> newNewMarker;

vector<list<int> > newToOld;
vector<int> walkFirstNode; //given a walk id, what's the first node of that walk
unordered_map<int, vector<edge_t> > sinkSrcEdges; //int is the unitig id (old id)

string getFileName(const string& s) {

   char sep = '/';

   size_t i = s.rfind(sep, s.length());
   if (i != string::npos) {
      return(s.substr(i+1, s.length() - i));
   }

   return s;
}


inline string plus_strings(const string& a, const string& b, size_t kmersize) {
    if (a == "") return b;
    if (b == "") return a;
    string ret = a + b.substr(kmersize - 1, b.length() - (kmersize - 1));
    return ret;
}

string delSpaces(string &str) {
    str.erase(std::remove(str.begin(), str.end(), ' '), str.end());
    return str;
}

bool charToBool(char c) {
    if (c == '+') {
        return true;
    } else {
        if (c != '-') cout << "Erroneus character in BCALM2 output." << endl;
        return false;
    }
}

string reverseComplement(string base) {
    size_t len = base.length();
    char* out = new char[len + 1];
    out[len] = '\0';
    for (int i = 0; i < len; i++) {
        if (base[i] == 'A') out[len - i - 1] = 'T';
        else if (base[i] == 'C') out[len - i - 1] = 'G';
        else if (base[i] == 'G') out[len - i - 1] = 'C';
        else if (base[i] == 'T') out[len - i - 1] = 'A';
    }
    string outString(out);
    free(out);
    return outString;
}

double readTimer() {
    return clock() / (double) CLOCKS_PER_SEC;
}

inline string currentDateTime() {
    // Get current date/time, format is YYYY-MM-DD HH:mm:ss
    time_t now = time(0);
    struct tm tstruct;
    char buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof (buf), "%Y-%m-%d %X\n", &tstruct);
    return buf;
}


int countOutArcs(int node) {
    return (adjList.at(node)).size();
}


inline char boolToCharSign(bool sign) {
    return (sign == true) ? '+' : '-';
}


// @@ --- ALL PRINTING CODE --- //

void printBCALMGraph(vector<vector<edge_t> > adjList) {
    for (int i = 0; i < adjList.size(); i++) {
        cout << i << "# ";
        for (edge_t edge : adjList.at(i)) {
            cout << boolToCharSign(edge.left) << ":" << edge.toNode << ":" << boolToCharSign(edge.right) << ", ";
        }
        cout << endl;
    }
}

class GroupMerger {
public:
    map<int, bool> fwdVisited;
    map<int, bool> bwdVisited;
    map<int, int> bwdWalkId;
    map<int, int> fwdWalkId;
    GroupMerger() {
    }
    void connectGroups(int from, int to){
        fwdVisited[from] = false;
        bwdVisited[to] = false;
        fwdWalkId[from] = to;
        bwdWalkId[to] = from;
    }
    ~GroupMerger() {
    }
};


class DisjointSet {
    unordered_map<int, int> parent;

public:
    DisjointSet() {
    }
    void make_set(int id) {
        this->parent[id] = -1;
    }

    void Union(int xId, int yId) {
        int xset = find_set(xId);
        int yset = find_set(yId);

        if(xset != yset)
        {
            parent[xset] = yset;
        }
    }

    int find_set(int id) {
        if (parent[id] == -1)
            return id;
        return find_set(parent[id]);
    }
    ~DisjointSet(){
    }

};


class Graph {
public:
    size_t V = adjList.size();
    int countNewNode = 0;
    int time = 0;

    char* color;
    int* p_dfs;
    bool* nodeSign;
    new_node_info_t* oldToNew;
    bool* saturated;
    struct node_sorter * sortStruct;
    bool* countedForLowerBound;
    DisjointSet disSet;
    GroupMerger gmerge;

    Graph() {
        color = new char[V];
        p_dfs = new int[V];
        nodeSign = new bool[V];
        oldToNew = new new_node_info_t[V];
        saturated = new bool[V];
        sortStruct = new struct node_sorter[V];
        global_indegree = new int[V];
        global_outdegree = new int[V];
        global_plusindegree = new int[V];
        global_plusoutdegree = new int[V];
        global_issinksource = new int[V];
        //global_priority = new int[V];
        countedForLowerBound = new bool[V];

        for (int i = 0; i < V; i++) {

            if(ALGOMODE == TWOWAYEXT || ALGOMODE == BRACKETCOMP ){
                disSet.make_set(i);
            }

            oldToNew[i].serial = -1;
            saturated[i] = false;
            sortStruct[i].sortkey = 0;
            sortStruct[i].node = i;
            global_indegree[i] = 0;
            global_outdegree[i] = 0;
            global_plusindegree[i] = 0;
            global_plusoutdegree[i] = 0;
            global_issinksource[i] = 0;
            //global_priority[i] = 0;
            countedForLowerBound[i] = false;
        }
    }

    inline bool sRight(edge_t plusminusedge){
        return !(plusminusedge.right == true);
    }

    inline bool sLeft(edge_t plusminusedge){
        return (plusminusedge.left == true);
    }

    void indegreePopulate(){
        int xc = 0;
        for(vector<edge_t> elist: adjList){
            for(edge_t e: elist){
                global_indegree[e.toNode] += 1;
                sortStruct[e.toNode].sortkey = sortStruct[e.toNode].sortkey + 1;
                if(e.right == true){
                    global_plusindegree[e.toNode] += 1;
                }
                if(e.left == true){
                    global_plusoutdegree[xc] += 1;
                }

            }
            global_outdegree[xc] = elist.size();
            xc++;
        }

        for(int i = 0; i<5; i++){
            for(int j = 0; j<5; j++){
                inOutCombo[make_pair(i,j)] = 0;
            }
        }
        for(int i = 0; i<V; i++){
            pair<int, int> a;
            a = make_pair(global_plusindegree[i], global_plusoutdegree[i]);
            inOutCombo[a] = (inOutCombo.count(a)  ? inOutCombo[a] + 1 : 1  );


            if(DBGFLAG == SINKSOURCE){
                cout<<i<<"is ";
            }

            if(global_plusoutdegree[i] == 0 && global_plusindegree[i] != 0){
                sink_count++;
                global_issinksource[i] = 1;
                //global_priority[i] = 5;
                countedForLowerBound[i] = true;

                if(DBGFLAG == SINKSOURCE){
                    cout<<"sink, ";
                }

            }
            if(global_plusindegree[i] == 0 && global_plusoutdegree[i] != 0){
                source_count++;
                global_issinksource[i] = 1;
                //global_priority[i] = 5;
                countedForLowerBound[i] = true;

                if(DBGFLAG == SINKSOURCE){
                    cout<<"source, ";
                }
            }

            if(global_indegree[i] == 0){
                global_issinksource[i] = 1;
                isolated_node_count++;

                if(DBGFLAG == SINKSOURCE){
                    cout<<"isolated, ";
                }
            }
            if(global_indegree[i] == 1){
                onecount++;
            }

            if(DBGFLAG == SINKSOURCE){
                cout<<endl;
            }

        }

        xc = 0; // current vertex while traversing the adjacency list
        for(vector<edge_t> elist: adjList){
            int neighborCount = 0;
            int spNeighborCount[2];
            spNeighborCount[0]=0;
            spNeighborCount[1]=0;
            stack<int> countedNodes;
            set<pair<int, bool> > countedSides;
            //if(true){
            if(FLG_NEWUB == true){
                //ENDPOINT SIDE UPPER BOUND - improved
                    for(edge_t e_xy: elist){    //central node: all neighbors of x
                        int y = e_xy.toNode;
                        vector<edge_t> adjY = adjList[y];
                        bool eligibleSp = true;

                        //pair<int, bool> pairr;
                        for(edge_t e_from_y : adjY){    // check if this neighbor is speacial
                            //pairr =make_pair(e_from_y.toNode, sRight(e_xy) );
                            if(e_from_y.toNode!=xc){

                                if(sRight(e_xy) == sLeft(e_from_y)){
                                    eligibleSp = false;
                                    break;
                                }
                            }
                        }

                        if(eligibleSp){
                            spNeighborCount[sLeft(e_xy)]++;
                        }
                    }
                    if(spNeighborCount[0]>1){
                        sharedparent_count += spNeighborCount[0] - 1 ;
                    }
                    if(spNeighborCount[1]>1){
                        sharedparent_count += spNeighborCount[1] - 1 ;
                    }
            }
            if(FLG_NEWUB == false){
                //ENDPOINT SIDE UPPER BOUND
                for(edge_t e_xy: elist){
                    int y = e_xy.toNode;
                    vector<edge_t> adjY = adjList[y];
                    bool eligible = true;
                    pair<int, bool> pairr;
                    for(edge_t e_from_y : adjY){
                        pairr =make_pair(e_from_y.toNode, sRight(e_xy) );
                        if(e_from_y.toNode!=xc){

                            if(sRight(e_xy) == sLeft(e_from_y)){
                                eligible = false;
                                break;
                            }

                        }

                    }

                    if(eligible){
                        neighborCount++;
                    }
                }

                if(global_issinksource[xc] == 1){
                    if(neighborCount>1){
                        sharedparent_count += neighborCount - 1 ;
                    }
                }else{
                    if(neighborCount>2){
                        sharedparent_count += neighborCount - 2 ;
                    }
                }
            }
            //sharedparent_count_wrong =sharedparent_count;

            //if(true){
            if(1==0){
                // OLDER UPPER BOUND CALC

                int neighborCount = 0;
                for(edge_t e_xy: elist){
                    int y = e_xy.toNode;

                    if(!countedForLowerBound[y]){
                        vector<edge_t> adjY = adjList[y];
                        bool eligible = true;
                        for(edge_t e_from_y : adjY){
                            if(e_from_y.toNode!=xc){
                                if(sRight(e_xy) == sLeft(e_from_y) ){
                                    eligible = false;
                                    break;
                                }
                            }
                        }
                        if(eligible){
                            countedForLowerBound[y] = true;
                            //global_priority[y] = 4;
                            neighborCount++;
                            countedNodes.push(y);
                        }
                    }
                }

                if(global_issinksource[xc] == 1){
                    if(neighborCount>1){
                        sharedparent_count += neighborCount - 1 ;
                    }else{
                        while(!countedNodes.empty()){
                            countedForLowerBound[countedNodes.top()] = false;
                            countedNodes.pop();
                        }
                    }
                }else{
                    if(neighborCount>2){
                        sharedparent_count += neighborCount - 2 ;
                    }else{
                        while(!countedNodes.empty()){
                            countedForLowerBound[countedNodes.top()] = false;
                            countedNodes.pop();
                        }
                    }
                }
            }

            xc++;
        }

        //check if not ALGOMODE == INDEGREE_DFS_INVERTED
        delete [] global_indegree;
        delete [] global_outdegree;
        delete [] global_plusindegree;
        delete [] global_plusoutdegree;
    }


    void DFS_visit(int u) {
        if(ALGOMODE == BRACKETCOMP){
            if(global_issinksource[u]==1){
                vector<edge_t> adju = adjList.at(u);
                vector<edge_t> myvector;
                for (edge_t e : adju) {
                    myvector.push_back(e);
                }
                sinkSrcEdges[u] = myvector;
                return;
            }
        }

        stack<edge_t> s;
        edge_t uEdge;
        uEdge.toNode = u;
        s.push(uEdge);

        while (!s.empty()) {
            edge_t xEdge = s.top();

            int x = xEdge.toNode;
            s.pop();

            if (color[x] == 'w') {
                //Original DFS code
                time = time + 1;
                color[x] = 'g';
                s.push(xEdge);
                vector<edge_t> adjx = adjList.at(x);
                if(ALGOMODE == RANDOM_DFS){
                    random_shuffle ( adjx.begin(), adjx.end() );
                }

//                if(ALGOMODE == EPPRIOR){
//                    sort( adjx.begin( ), adjx.end( ), [ ]( const edge_t& lhs, const edge_t& rhs )
//                         {
//                             return global_priority[lhs.toNode]   <  global_priority[rhs.toNode]  ;
//                         });
//                }

                if(ALGOMODE == INDEGREE_DFS){
                    sort( adjx.begin( ), adjx.end( ), [ ]( const edge_t& lhs, const edge_t& rhs )
                         {
                             return global_indegree[lhs.toNode] < global_indegree[rhs.toNode];
                         });
                }

                if(ALGOMODE == PLUS_INDEGREE_DFS){
                    sort( adjx.begin( ), adjx.end( ), [ ]( const edge_t& lhs, const edge_t& rhs )
                         {
                             return global_indegree[lhs.toNode]  - global_plusindegree[lhs.toNode] > global_indegree[lhs.toNode]  - global_plusindegree[rhs.toNode];
                         });
                }

                if(ALGOMODE == INDEGREE_DFS_INVERTED){
                    sort( adjx.begin( ), adjx.end( ), [ ]( const edge_t& lhs, const edge_t& rhs )
                         {
                             return global_indegree[lhs.toNode] > global_indegree[rhs.toNode];
                         });
                }
                if (ALGOMODE == OUTDEGREE_DFS){
                    if(p_dfs[x] == -1){
                        sort( adjx.begin( ), adjx.end( ), [ ]( const edge_t& lhs, const edge_t& rhs )
                             {
                                 return global_outdegree[lhs.toNode] < global_outdegree[rhs.toNode];
                             });
                    }else if (nodeSign[x] == false){
                        sort( adjx.begin( ), adjx.end( ), [ ]( const edge_t& lhs, const edge_t& rhs )
                             {
                                 return global_outdegree[lhs.toNode] - global_plusoutdegree[lhs.toNode] < global_outdegree[rhs.toNode] - global_plusoutdegree[rhs.toNode];
                             });
                    }else if (nodeSign[x] == true){
                        sort( adjx.begin( ), adjx.end( ), [ ]( const edge_t& lhs, const edge_t& rhs )
                             {
                                 return global_plusoutdegree[lhs.toNode] < global_plusoutdegree[rhs.toNode];
                             });
                    }
                }


                // Now our branching code ::
                // For a white x
                // Consider 2 case:
                // Case 1. p[x] = -1, it can happen in two way, x is the first one ever in this connected component, or no one wanted to take x
                // either way, if p[x] = -1, i can be representative of a new node in new graph
                // Case 2. p[x] != -1, so x won't be the representative/head of a newHome. x just gets added to its parent's newHome.
                int u = unitigs.at(x).ln; //unitig length

                if (p_dfs[x] == -1) {

                    list<int> xxx;
                    xxx.push_back(x);
                    newToOld.push_back(xxx);
                    oldToNew[x].serial = countNewNode++; // countNewNode starts at 0, then keeps increasing
                    oldToNew[x].finalWalkId = oldToNew[x].serial;


                    //added while doing bracket comp
                    walkFirstNode.push_back(x);

                    oldToNew[x].pos_in_walk = 1;
                    oldToNew[x].startPosWithKOverlap = 1;
                    if (u < K) {
                        oldToNew[x].endPosWithKOVerlap = 1; // do we actually see this? yes
                        if(DBGFLAG == UKDEBUG){
                            cout<< "node: "<< x<<"u< k ***** u = "<<u<<endl;
                        }
                    } else {
                        oldToNew[x].endPosWithKOVerlap = u - K + 1;
                    }

                } else {

                    newToOld[oldToNew[p_dfs[x]].serial].push_back(x);
                    oldToNew[x].serial = oldToNew[p_dfs[x]].serial;
                    oldToNew[x].finalWalkId = oldToNew[x].serial;


                    if(ALGOMODE==TWOWAYEXT || ALGOMODE==BRACKETCOMP ){
                        disSet.Union(x, p_dfs[x]);
                    }


                    oldToNew[x].startPosWithKOverlap = oldToNew[p_dfs[x]].endPosWithKOVerlap + 1;
                    oldToNew[x].pos_in_walk = oldToNew[p_dfs[x]].pos_in_walk + 1;
                    if (u < K) {
                        oldToNew[x].endPosWithKOVerlap = oldToNew[x].startPosWithKOverlap + 1; // do we actually see this? yes
                        if(DBGFLAG == UKDEBUG){
                            cout<< "node: "<< x<<"u< k ***** u = "<<u<<endl;
                        }
                    } else {
                        oldToNew[x].endPosWithKOVerlap = u - K + (oldToNew[x].startPosWithKOverlap); //check correctness
                    }
                }

                // x->y is the edge, x is the parent we are extending
                for (edge_t yEdge : adjx) { //edge_t yEdge = adjx.at(i);
                    int y = yEdge.toNode;

                    if(ALGOMODE == BRACKETCOMP){
                        if(global_issinksource[y] == true){
                            continue;
                        }
                    }


                    //Normal DFS
                    if (color[y] == 'w') {
                        s.push(yEdge);
                    }


                    //handle self-loop, self-loop will always be an extra edge
                    // redundant, because we remove self-loop anyway
                    if (y == x) {
                        edge_both_t e;
                        e.edge = yEdge;
                        e.fromNode = x;
                    } else if (saturated[x]) {
                        // Since x is saturated, we only add resolveLater edges
                        // no need to check for consistency
                        if (y != p_dfs[x]) {
                            edge_both_t e;
                            e.edge = yEdge;
                            e.fromNode = x;
                        }
                    } else {
                        // If x has space to take a child, meaning x is not saturated
                        // hunting for potential child

                        if (color[y] == 'w' && p_dfs[y] == -1) {
                            // y has white color & got no parent => means it's homeless, so let's see if we can take it as a child of x
                            //But just see if it is eligible to be a child, i.e. is it consistent (sign check)?

                            //2 case, Does x's child have grandparent?
                            // If No:
                            if (p_dfs[x] == -1 && ALGOMODE != NODEASSIGN) {
                                // case 1: child has no grandparent
                                // so extend path without checking any sign

                                nodeSign[x] = yEdge.left;
                                nodeSign[y] = yEdge.right;
                                p_dfs[y] = x;
                                saturated[x] = true; //found a child

                            } else if (nodeSign[x] == yEdge.left) {
                                // case 2: child (=y) has grandparent, i.e. x's parent exists
                                nodeSign[y] = yEdge.right;
                                p_dfs[y] = x;
                                saturated[x] = true; //found a child

                            } else {
                                // do we reach this case?
                                edge_both_t e;
                                e.edge = yEdge;
                                e.fromNode = x;
                            }

                        } else {

                            //merger
                            if(ALGOMODE == TWOWAYEXT || ALGOMODE == BRACKETCOMP){
                                // y is not white
                                bool consistentEdge = (nodeSign[y] == yEdge.right && (p_dfs[x]==-1 || (p_dfs[x]!=-1&& nodeSign[x] == yEdge.left)) );
                                if(p_dfs[y]==-1 && consistentEdge && oldToNew[x].serial != oldToNew[y].serial){

                                    //cout<<"x: "<<x<<":" <<disSet.find_set(x)<<" ";
                                    //cout<<"y: "<<y<<":" <<disSet.find_set(y) <<endl;

                                    //not in same group already, prevent cycle
                                    if(disSet.find_set(x)!=disSet.find_set(y)){
                                        nodeSign[x] = yEdge.left;
                                        nodeSign[y] = yEdge.right;
                                        p_dfs[y] = x;
                                        saturated[x] = true; //found a child
                                        // oldToNew[y].serial

                                        disSet.Union(x, y);
                                    gmerge.connectGroups(oldToNew[x].serial,oldToNew[y].serial );

                                    }
                                }
                            }

                            if (y != p_dfs[x]) {
                                edge_both_t e;
                                e.edge = yEdge;
                                e.fromNode = x;
                            }
                        }
                    }
                }
            } else if (color[x] == 'g') {
                time = time + 1;
                color[x] = 'b';
            }
        }
    }



    void DFS() {

        if(ALGOMODE == NODEASSIGN){
            for (int i=0; i<V; i++) {
                nodeSign[i] = true;
                if(global_plusindegree[i]< global_indegree[i] - global_plusindegree[i]){
                    nodeSign[i] = true;
                }
            }
        }

        if(ALGOMODE == SOURCEFIRST){
            for (int i = 0; i < V; i++) {
                sortStruct[i].node = i;
                sortStruct[i].sortkey = global_issinksource[i];
            }
            vector<struct node_sorter> myvector (sortStruct, sortStruct+V);
            sort (myvector.begin(), myvector.end(), sort_by_key_inverted);
            copy(myvector.begin(), myvector.end(), sortStruct);
        }

//        if(ALGOMODE == EPPRIOR){
//            for (int i = 0; i < V; i++) {
//                sortStruct[i].node = i;
//                sortStruct[i].sortkey = global_priority[i];
//            }
//            vector<struct node_sorter> myvector (sortStruct, sortStruct+V);
//            sort (myvector.begin(), myvector.end(), sort_by_key_inverted);
//            copy(myvector.begin(), myvector.end(), sortStruct);
//        }

        if(ALGOMODE == INDEGREE_DFS_INVERTED){
            for (int i = 0; i < V; i++) {
                sortStruct[i].node = i;
                sortStruct[i].sortkey = global_indegree[i];
            }
            vector<struct node_sorter> myvector (sortStruct, sortStruct+V);
            sort (myvector.begin(), myvector.end(), sort_by_key_inverted);
            //random_shuffle ( myvector.begin(), myvector.end() );
            copy(myvector.begin(), myvector.end(), sortStruct);

        }

        if (ALGOMODE == INDEGREE_DFS || ALGOMODE == INDEGREE_DFS_1 ){
            for (int i = 0; i < V; i++) {
                sortStruct[i].node = i;
                sortStruct[i].sortkey = global_indegree[i];
            }
            vector<struct node_sorter> myvector (sortStruct, sortStruct+V);
            sort (myvector.begin(), myvector.end(), sort_by_key);
            copy(myvector.begin(), myvector.end(), sortStruct);

            if(DBGFLAG == INDEGREEPRINT){
                cout<<"print in degrees"<<endl;
                for(int i = 0; i<V; i++){
                    cout<<sortStruct[i].node<<"->"<<sortStruct[i].sortkey<<endl;
                }
            }
        }


        if(ALGOMODE == RANDOM_DFS){
            for (int i = 0; i < V; i++) {
                sortStruct[i].node = i;
                sortStruct[i].sortkey = global_indegree[i];
            }
            vector<struct node_sorter> myvector (sortStruct, sortStruct+V);
            sort (myvector.begin(), myvector.end(), sort_by_key);
            random_shuffle ( myvector.begin(), myvector.end() );
            copy(myvector.begin(), myvector.end(), sortStruct);

        }


        if (ALGOMODE == OUTDEGREE_DFS || ALGOMODE == OUTDEGREE_DFS_1){
            for (int i = 0; i < V; i++) {
                sortStruct[i].node = i;
                sortStruct[i].sortkey = global_outdegree[i];
            }
            vector<struct node_sorter> myvector (sortStruct, sortStruct+V);
            sort (myvector.begin(), myvector.end(), sort_by_key);
            copy(myvector.begin(), myvector.end(), sortStruct);
        }

        double time_a = readTimer();
        for (int i = 0; i < V; i++) {
            color[i] = 'w';
            p_dfs[i] = -1;
        }
        cout<<"Basic V loop time: "<<readTimer() - time_a<<" sec"<<endl;


        time_a = readTimer();
        for (int j = 0; j < V; j++) {
            int i;
            if(ALGOMODE == OUTDEGREE_DFS || ALGOMODE == OUTDEGREE_DFS_1 || ALGOMODE == INDEGREE_DFS || ALGOMODE == INDEGREE_DFS_1 || ALGOMODE == SOURCEFIRST){
                i = sortStruct[j].node;
            }else{
                i = j;
            }



            if (color[i] == 'w') {
                if(DBGFLAG == DFSDEBUGG ){
                    cout<<"visit start of node: "<<i<<endl;
                }
                DFS_visit(i);
            }
        }
        cout<<"DFS time: "<<readTimer() - time_a<<" sec"<<endl;


        /***MERGE START***/
        bool* merged = new bool[countNewNode];
        for (int i = 0; i<countNewNode; i++) {
            merged[i] = false;
        }

        if(ALGOMODE == TWOWAYEXT || ALGOMODE == BRACKETCOMP){
            ofstream uidSequenceFile;
            uidSequenceFile.open("uidSeq"+modefilename[ALGOMODE]+".txt");


            for ( const auto& p: gmerge.fwdWalkId)
            {
                if(gmerge.fwdVisited[p.first] == false){

                    int fromnode =p.first;
                    int tonode = p.second;
                    deque<int> lst;

                    lst.push_back(fromnode);
                    lst.push_back(tonode);

                    gmerge.fwdVisited[fromnode] = true;
                    gmerge.bwdVisited[tonode] = true;


                    if(gmerge.fwdVisited.count(tonode)>0){
                        while(gmerge.fwdVisited[tonode] == false){
                            gmerge.fwdVisited[tonode] = true;
                            tonode = gmerge.fwdWalkId[tonode];
                            gmerge.bwdVisited[tonode] = true;

                            lst.push_back(tonode);
                            if(gmerge.fwdVisited.count(tonode)==0)
                                break;
                        }
                    }
                    if(gmerge.bwdVisited.count(fromnode)>0){
                        while(gmerge.bwdVisited[fromnode] == false){
                            gmerge.bwdVisited[fromnode] = true;
                            fromnode = gmerge.bwdWalkId[fromnode];
                            gmerge.fwdVisited[fromnode] = true;

                            lst.push_front(fromnode);
                            if(gmerge.bwdVisited.count(fromnode)==0)
                                break;
                        }
                    }


                    int headOfThisWalk = walkFirstNode[lst.at(0)]; //CHECK AGAIN
                    assert(!lst.empty());
                    int commonWalkId = lst.at(0);

                    int posOffset = 1;

                    int lastWalk = -1;
                    for(auto i: lst){
                        // i is new walk id before merging
                        merged[i] = true;


                        walkFirstNode[i] = headOfThisWalk;

                        // travesing the walk list of walk ID i
                        for(int uid: newToOld[i]){
                            oldToNew[uid].serial = commonWalkId;
                            oldToNew[uid].finalWalkId = commonWalkId;
                            oldToNew[uid].pos_in_walk = posOffset++;
                        }
                    }
                    oldToNew[newToOld[lst.back()].back()].isWalkEnd = true;

                    V_twoway_ustitch ++;

                }
            }
            for (int newNodeNum = 0; newNodeNum<countNewNode; newNodeNum++){


                if(merged[newNodeNum] == false){
                    oldToNew[newToOld[newNodeNum].back()].isWalkEnd = true;

                    V_twoway_ustitch++;
                }
            }
        }



        //sorter of all walks and printing them
        vector<fourtuple> sorter;
        for(int uid = 0 ; uid< V; uid++){
            new_node_info_t nd = oldToNew[uid];
            sorter.push_back(make_tuple(uid, nd.finalWalkId, nd.pos_in_walk, nd.isTip));
        }
        //stable_sort(sorter.begin(),sorter.end(),sort_by_tipstatus);
        stable_sort(sorter.begin(),sorter.end(),sort_by_pos);
        stable_sort(sorter.begin(),sorter.end(),sort_by_walkId);

        ofstream uidSequence;
        string uidSeqFilename = "uidSeq.usttemp"; //"uidSeq"+ mapmode[ALGOMODE] +".txt"
       uidSequence.open(uidSeqFilename);

        int finalUnitigSerial = 0;
        for(fourtuple n : sorter){
            int uid = get<0>(n);
            int bcalmid = unitigs.at(uid).serial;
//                        int finalWalkId = get<1>(n);
//                        int pos_in_walk = get<2>(n);
//                        int isTip = get<3>(n);
//                        cout<<uid<<" " <<finalWalkId<<" "<<pos_in_walk<<" "<<isTip<<" "<<oldToNew[uid].isWalkEnd<< " was merged: "<< merged[oldToNew[uid].finalWalkId]<< endl;
            uidSequence << finalUnitigSerial <<" "<< bcalmid << endl;
            finalUnitigSerial++;
        }
        uidSequence.close();


        //system();

        //keep the sequences only
        system(("awk '!(NR%2)' "+UNITIG_FILE+" > seq.usttemp").c_str());
        system("sort -n -k 2 -o uidSeq.usttemp uidSeq.usttemp");
        if(FLG_ABUNDANCE){
            system(("awk '(NR%2)' "+UNITIG_FILE+" | cut -f 5 -d ':' | cut -f 1 -d 'L' > count.usttemp").c_str()); // get a separate count file
            system("paste -d' ' uidSeq.usttemp seq.usttemp count.usttemp > merged.usttemp ");
            system("sort -n -k 1 -o merged.usttemp merged.usttemp");
            system(("cat  merged.usttemp  | awk '{for (i=4;i<=NF;i+=1) print $i}' > "+getFileName(UNITIG_FILE)+".ust.counts").c_str());
        }else{
            system("paste -d' ' uidSeq.usttemp seq.usttemp > merged.usttemp ");
            system("sort -n -k 1 -o merged.usttemp merged.usttemp");
        }
        system("cat  merged.usttemp  | cut -d' ' -f3 >  seq.usttemp");


        ifstream sequenceStringFile ("seq.usttemp");
        ofstream ustOutputFile (OUTPUT_FILENAME);
        //ofstream smallKFile("smallK.fa");

        //both string and abundance sort
        //keep string only and output
        //open the string file
        if(100==100){

            int lastWalk = -1;
            string walkString = "";
            string unitigString = "";
            for(fourtuple n : sorter){
                int uid = get<0>(n);
                int finalWalkId = get<1>(n);
                int pos_in_walk = get<2>(n);

                //for each line in file
                string sequenceFromFile = "";//getline
                getline (sequenceStringFile,sequenceFromFile);
                if(nodeSign[uid] == false){
                    unitigString =  reverseComplement(sequenceFromFile);
                }else{
                    unitigString =  sequenceFromFile;
                }

                if(finalWalkId!=lastWalk){
                    if(lastWalk != -1){
                        //print previous walk
                        //ustOutputFile<<">"<<lastWalk << " " << uid<<" " <<finalWalkId<<" "<<pos_in_walk<<" "<<endl;
                        if(walkString.length()>=K){
                            ustOutputFile<<">"<<endl;
                            C_twoway_ustitch+=walkString.length();

                            ustOutputFile<< walkString<<endl;
                        }else{
//                            smallKFile<<">\n"<< walkString+"A" <<endl;
//                            smallKFile<<">\n"<< walkString+"C" <<endl;
//                            smallKFile<<">\n"<< walkString+"G" <<endl;
//                            smallKFile<<">\n"<< walkString+"T" <<endl;
//                            smallKFile<<">\n"<< "A"+walkString <<endl;
//                            smallKFile<<">\n"<< "C"+walkString <<endl;
//                            smallKFile<<">\n"<< "G"+walkString <<endl;
//                            smallKFile<<">\n"<< "T"+walkString <<endl;
                        }
                    }

                    //start a new walk
                    // cout<<"Walk: (" <<finalWalkId<<" ) = ";
                    walkString = "";
                    lastWalk = finalWalkId;
                }
                walkString = plus_strings(walkString, unitigString, K);

                //ustOutputFile<<">"<<uid<<" " <<finalWalkId<<" "<<pos_in_walk<<endl;
            }
             if(walkString.length()>=K){
               ustOutputFile<<">"<<endl;
               C_twoway_ustitch+=walkString.length();

               ustOutputFile<< walkString<<endl;
           }else{
//               smallKFile<<">\n"<< walkString+"A" <<endl;
//               smallKFile<<">\n"<< walkString+"C" <<endl;
//               smallKFile<<">\n"<< walkString+"G" <<endl;
//               smallKFile<<">\n"<< walkString+"T" <<endl;
//               smallKFile<<">\n"<< "A"+walkString <<endl;
//               smallKFile<<">\n"<< "C"+walkString <<endl;
//               smallKFile<<">\n"<< "G"+walkString <<endl;
//               smallKFile<<">\n"<< "T"+walkString <<endl;
           }

            sequenceStringFile.close();
            //smallKFile.close();
            system("rm -rf *.usttemp");
            system("rm -rf uidSeq.txt");
        }

        // clean up
        delete [] merged;

        /// TWOWAYEXT DONE:

        delete []  global_issinksource;
        //delete []  global_priority;
        freeDFS();
    }

    inline void freeDFS(){
        delete [] saturated;
        delete [] p_dfs;
    }

    ~Graph() {
        delete [] color;
        delete [] nodeSign;
        delete [] oldToNew;

        delete [] sortStruct;
        delete [] countedForLowerBound;
    }
};

int maximumUnitigLength(){
    int m = 0;
    for(unitig_struct_t u: unitigs){
        if(u.ln > m){
            m = u.ln;
        }
    }
    return m;
}

void printNewGraph(Graph &G){
    list<int> *newToOld = new list<int>[G.countNewNode];

    for (int i = 0; i < G.V; i++) {
        newToOld[G.oldToNew[i].serial].push_back(i);
        cout << "old " << i << "-> new" << G.oldToNew[i].serial << endl;
    }
    for (int i = 0; i < G.countNewNode; i++) {
        list<int> adj = newToOld[i];

        cout<<"new " << i<<": old (index) ";
        for(int val : adj){
            cout<<val<<" ";
        }
        cout<<endl;
    }
    delete [] newToOld;
}


void formattedOutputForwardExt(Graph &G){
    string plainOutput = "plainOutput"+modefilename[ALGOMODE]+".txt";
    ofstream plainfile;
    plainfile.open(plainOutput);

    string stitchedUnitigs = "stitchedUnitigs"+modefilename[ALGOMODE]+".fa";
    ofstream myfile;
    myfile.open (stitchedUnitigs);
    //>0 LN:i:13 KC:i:12 km:f:1.3  L:-:0:- L:-:2:-  L:+:0:+ L:+:1:-
    for (int newNodeNum = 0; newNodeNum<G.countNewNode; newNodeNum++){
        myfile << '>' << newNodeNum <<" LN:i:"<<newSequences[newNodeNum].length()<<" ";
        myfile<<endl;

        plainfile<<newSequences[newNodeNum];
        myfile<<newSequences[newNodeNum];

        plainfile<<endl;
        myfile<<endl;
    }
    //myfile << '>' << newNodeNum <<">0 LN:i:13 KC:i:12 km:f:1.3  L:-:0:- L:-:2:-  L:+:0:+ L:+:1:- " ;
    myfile.close();
    plainfile.close();
}

int read_unitig_file(const string& unitigFileName, vector<unitig_struct_t>& unitigs) {
    ifstream unitigFile;
    unitigFile.open(unitigFileName);

    string line;

    int nodeNum;
    char lnline[20];
    char kcline[20];
    char kmline[20];
    char edgesline[100000];
    bool doCont = false;

    int smallestK = 9999999;
    getline(unitigFile, line);

    do {
        unitig_struct_t unitig_struct;

        if(FLG_ABUNDANCE){
            //>3 LN:i:24 ab:Z:10 10 10 10 10 7 7 7 7 7 3 3 3 3   L:-:0:+ L:-:1:+  L:+:0:-
            edgesline[0] = '\0';
            sscanf(line.c_str(), "%*c %d %s", &unitig_struct.serial, lnline);

			if(	line.find("ab:Z") == string::npos){
				cout<<"Incorrect input format. Try using flag -a 0."<<endl;
				exit(3);
			}

            sscanf(lnline, "%*5c %d", &unitig_struct.ln);


            int abpos = line.find("ab") + 5;
            int Lpos = line.find("L:");

            if(Lpos < 0){
                Lpos = line.length() ;
            }
            // initialize string stream
            //cout<<line.substr(abpos, Lpos - abpos);
            stringstream ss(line.substr(abpos, Lpos - abpos));
            string abun;

           sscanf(line.substr(Lpos, line.length() - Lpos).c_str(), "%[^\n]s", edgesline);

            if(unitig_struct.ln < smallestK){
                           smallestK = unitig_struct.ln ;
                       }
            if(unitig_struct.ln < K){
                printf("Wrong k! Try again with correct k value. \n");
                exit(2);
            }

        }else{
            edgesline[0] = '\0';
            sscanf(line.c_str(), "%*c %d %s  %s  %s %[^\n]s", &unitig_struct.serial, lnline, kcline, kmline, edgesline);

			if(	line.find("KC") == string::npos){
				cout<<"Incorrect input format. Try using flag -a 1."<<endl;
				exit(3);
			}

            //>0 LN:i:13 KC:i:12 km:f:1.3  L:-:0:- L:-:2:-  L:+:0:+ L:+:1:-
            sscanf(lnline, "%*5c %d", &unitig_struct.ln);


			if(unitig_struct.ln < smallestK){
                           smallestK = unitig_struct.ln ;
                       }
			if(unitig_struct.ln < K){
                printf("Wrong k! Try again with correct k value. \n");
                exit(2);
            }

        }



        char c1, c2;
        stringstream ss(edgesline);

        vector<edge_t> edges;
        while (getline(ss, line, ' ')) {
            if (delSpaces(line).length() != 0) {
                if(DBGFLAG==VERIFYINPUT){
                    cout<<line<<endl;
                }

                sscanf(line.c_str(), "%*2c %c %*c %d  %*c  %c", &c1, &nodeNum, &c2); //L:-:0:-
                edge_t newEdge;

                bool DELSELFLOOP=true;
                if(DELSELFLOOP){
                    if((unitig_struct.serial)!= nodeNum){
                        newEdge.left = charToBool(c1);
                        newEdge.right = charToBool(c2);
                        newEdge.toNode = nodeNum;
                        edges.push_back(newEdge);
                    }
                }else{
                    newEdge.left = charToBool(c1);
                    newEdge.right = charToBool(c2);
                    newEdge.toNode = nodeNum;
                    edges.push_back(newEdge);
                }

            }

        }
        adjList.push_back(edges);



        doCont = false;
        while (getline(unitigFile, line)) {
            if (line.substr(0, 1).compare(">")) {
                //unitig_struct.sequence = unitig_struct.sequence + line;
                unitigs.push_back(unitig_struct);
            } else {
                doCont = true;
                break;
            }
        }
    } while (doCont);


    unitigFile.close();

    if(smallestK > K ){
        cout<<"\n :::: :::: :::: :::: !!!!!!!!! WARNING !!!!!!!!!!! :::: :::: :::: ::::"<<endl;
        cout<<"The length of the smallest string we found was " << smallestK << ". Please make sure you are using the correct value of 'k' to ensure correctness of output."<<endl;
         cout << "------------------------------------------------------"<<endl;
    }

    //cout << "Complete reading input unitig file (bcalm2 file)." << endl;
    return EXIT_SUCCESS;
}

int main(int argc, char** argv) {


    #if defined(__unix__) || (defined(__APPLE__) && defined(__MACH__))
    #else
        cout<<"UST does not support this OS."<<endl;
        exit(2);
    #endif



    //    string debugFileName = "debug.txt";
    //    ofstream debugFile;
    //    debugFile.open(debugFileName);


    const char* nvalue = "" ;

    int c ;

    ///*
    if(DEBUGMODE==false){
        while( ( c = getopt (argc, argv, "i:k:a:") ) != -1 )
        {
            switch(c)
            {
                case 'a':
                    if(optarg) {
						if(strcmp(optarg, "0")==0 || strcmp(optarg, "1")==0){
                            FLG_ABUNDANCE = static_cast<bool>(std::atoi(optarg));
                        }else{
							fprintf(stderr, "Error: use either -a 0 or -a 1 \n");
                            exit(EXIT_FAILURE);
						}

                    }

                    break;
                case 'i':
                    if(optarg) nvalue = optarg;
                    break;
                //case 'f':
                    //if(optarg) {
                        //FLG_NEWUB = static_cast<bool>(std::atoi(optarg));
                    //}
                    //break;
//                case 'm':
//                    if(optarg) {
//                        ALGOMODE = static_cast<ALGOMODE_T>(std::atoi(optarg));
//                    }
//                    break;
                //case 'd':
                    //if(optarg) {
                        //DBGFLAG = static_cast<DEBUGFLAG_T>(std::atoi(optarg));
                    //}
                    //break;
                case 'k':
                    if(optarg) {
                        K = std::atoi(optarg) ;
                        if(K<=0){
                            fprintf(stderr, "Error: Specify a positive k value.\n");
                            exit(EXIT_FAILURE);
                        }
                    }else{
                        fprintf(stderr, "Usage: %s -k <kmer size> -i <input-file-name> [-a <0: for no counts, 1: to output separate count, default = 0>]\n",
                                argv[0]);
                        exit(EXIT_FAILURE);
                    }
                    break;
                default: //
                    fprintf(stderr, "Usage: %s -k <kmer size> -i <input-file-name> [-a <0: for no counts, 1: to output separate count, default = 0>]\n",
                            argv[0]);
                    exit(EXIT_FAILURE);

            }
        }

        if(K==0 || strcmp(nvalue, "")==0){
            fprintf(stderr, "Usage: %s -k <kmer size> -i <input-file-name> [-a <0: for no counts, 1: to output separate count, default = 0>]\n",
                    argv[0]);
            exit(EXIT_FAILURE);
        }



        UNITIG_FILE = string(nvalue);
    }
    //*/

    ifstream infile(UNITIG_FILE);
    if(!infile.good()){
        fprintf(stderr, "Error: File named \"%s\" cannot be opened.\n", UNITIG_FILE.c_str());
        exit(EXIT_FAILURE);
    }
    OUTPUT_FILENAME = getFileName(UNITIG_FILE)+".ust.fa";
    ofstream globalStatFile;
    system(("rm -rf "+getFileName(UNITIG_FILE)+".stats.txt").c_str());
    globalStatFile.open((getFileName(UNITIG_FILE)+".stats.txt").c_str(), std::fstream::out | std::fstream::app);


    double startTime = readTimer();
    cout << "## Please wait while we read input file.... " << UNITIG_FILE << ": k = "<<K<<endl;
    if (EXIT_FAILURE == read_unitig_file(UNITIG_FILE, unitigs)) {
        return EXIT_FAILURE;
    }
    infile.close();
    double TIME_READ_SEC = readTimer() - startTime;
    cout<<"Done. TIME to read file "<<TIME_READ_SEC<<" sec."<<endl;


    Graph G;

    //count total number of edges
    int E_bcalm = 0;
    for (int i = 0; i < G.V; i++) {
        E_bcalm += adjList[i].size();
    }
    int V_bcalm = G.V;
    int numKmers = 0;
    int C_bcalm = 0;

    for (unitig_struct_t unitig : unitigs) {
        C_bcalm += unitig.ln;
        numKmers +=  unitig.ln - K + 1;
    }

    //    if(DBGFLAG == NODENUMBER_DBG){
    //        cout<<"Total Nodes: "<<V<<" Edges: "<<E<<" K-mers: "<<numKmers<<endl;
    //    }

    cout<<"## Please wait while we gather info about lower bound.... "<<endl;
    double time_a = readTimer();
    G.indegreePopulate();

    cout<<"Done. TIME to gather information about unitig graph: "<<readTimer() - time_a<<" sec."<<endl;


    int walkstarting_node_count = ceil((sharedparent_count + sink_count + source_count)/2.0) + isolated_node_count;
    int charLowerbound = C_bcalm-(K-1)*(G.V - walkstarting_node_count*1.0);
    float upperbound = (1-((C_bcalm-(K-1)*(G.V - walkstarting_node_count*1.0))/C_bcalm))*100.0;

//    printf( "%d\t\
//           %d\t\
//           %d\t\
//           %d\t\
//           %d\t\
//           %d\t\
//           %.2f\t\
//           %.2f%%\t\
//           %d\t\
//           %d\t\
//           %d\t\
//           %d\t\
//           %.2f%%\t",
//           K,
//           numKmers,
//           V_bcalm,
//           E_bcalm,
//           C_bcalm,
//           charLowerbound,
//           (charLowerbound*2.0)/numKmers,
//           upperbound,
//           isolated_node_count,
//           sink_count,
//           source_count,
//           sharedparent_count,
//           sharedparent_count*100.0/V_bcalm
//           );


    globalStatFile << "K" <<  "=" << K << endl;
    globalStatFile << "N_KMER" <<  "=" << numKmers << endl;
    globalStatFile << "V_BCALM" <<  "=" << V_bcalm << endl;
    globalStatFile << "E_BCALM" <<  "=" << E_bcalm << endl;
    globalStatFile << "C_BCALM" <<  "=" << C_bcalm << endl;
    globalStatFile << "C_LB" <<  "=" << charLowerbound << endl;
    globalStatFile << "V_LB" <<  "=" << walkstarting_node_count << endl;
    globalStatFile << "N_ISOLATED" <<  "=" << isolated_node_count << endl;
    globalStatFile << "N_SINK" <<  "=" << sink_count << endl;
    globalStatFile << "N_SOURCE" <<  "=" << source_count << endl;

    if(FLG_NEWUB)
        globalStatFile << "N_SPECIAL_NEW" <<  "=" << sharedparent_count << endl;
    if(!FLG_NEWUB)
        globalStatFile << "N_SPECIAL_OLD" <<  "=" << sharedparent_count << endl;

    //OPTIONAL;DERIVABLE
    globalStatFile << "BITSKMER_LB" <<  "=" << (charLowerbound*2.0)/numKmers << endl;
    globalStatFile << "PERCENT_UB" <<  "=" << upperbound <<  "%" << endl;
    globalStatFile << "PERCENT_N_SPECIAL" <<  "=" << sharedparent_count*100.0/V_bcalm <<"%" << endl;


    for (auto i = inOutCombo.begin(); i != inOutCombo.end(); i++) {
        //printf("%.2f%%\t", (i->second)*100.0/V_bcalm);
        globalStatFile << "PERCENT_DEGREE_"<<(i->first).first << "_" << (i->first).second <<  "=" << (i->second)*100.0/V_bcalm <<"%" << endl;
    }
    //printf("%.2f%%\t\
           %.2f%%\t",
     //      isolated_node_count*100.0/V_bcalm,
     //      (sink_count+source_count)*100.0/V_bcalm);
    //OPTIONAL; derivable
    globalStatFile << "PERCENT_N_ISOLATED" <<  "=" << isolated_node_count*100.0/V_bcalm <<"%" << endl;
    globalStatFile << "PERCENT_N_DEADEND" <<  "=" << (sink_count+source_count)*100.0/V_bcalm <<"%" << endl;

    // Iterating the map and printing ordered values
    //    for (auto i = inOutCombo.begin(); i != inOutCombo.end(); i++) {
    //        cout << "(" << i->first.first<< ", "<< i->first.second << ")" << " := " << i->second << '\n';
    //    }

    if(ALGOMODE == PROFILE_ONLY){
        //printf("\n");
        return 0;
    }


    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@//
    //##################################################//
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@//
    //##################################################//

    //STARTING DFS
    cout<<"## Please wait while DFS step is going on... "<<endl;
    G.DFS();


    if(DBGFLAG == PRINTER){
        printBCALMGraph(adjList);
        printNewGraph(G);
        for(int i = 0; i< G.countNewNode; i++){
            cout<<"new ->" <<i<<" ";
            for(int x: newToOld[i]){
                cout<<x<<" ";
            }
            cout<<endl;
        }
    }


    double TIME_TOTAL_SEC = readTimer() - startTime;

    // For collecting stats
    //int U_MAX = maximumUnitigLength();

    time_a = readTimer();


    if(ALGOMODE==BRACKETCOMP){
        C_ustitch = C_tip_ustitch;
        V_ustitch  = V_tip_ustitch;
    }else if(ALGOMODE == TWOWAYEXT){
        V_ustitch = V_twoway_ustitch;
        C_ustitch = C_twoway_ustitch;
    }else if(ALGOMODE == BASIC){
        formattedOutputForwardExt(G);
        V_ustitch = G.countNewNode;
    }else{
        formattedOutputForwardExt(G);
        V_ustitch = G.countNewNode;
    }


    cout<<"Done. TIME to output: "<<readTimer() - time_a<<" sec."<<endl;


    float percent_saved_c = (1-(C_ustitch*1.0/C_bcalm))*100.0;
    float ustitchBitsPerKmer;
    if(ALGOMODE== BRACKETCOMP){
        ustitchBitsPerKmer =C_ustitch*3.0/numKmers;
    }else{
        ustitchBitsPerKmer =C_ustitch*2.0/numKmers;
    }


//    printf("%s\t",  mapmode[ALGOMODE].c_str());
//    printf("%d\t\
//           %.2f%%\t\
//           %.2f%%\t\
//           %d\t\
//           %.2f\t",
//           V_ustitch,
//           percent_saved_c,
//           upperbound - percent_saved_c,
//           C_ustitch,
//           ustitchBitsPerKmer
//           );
//    printf("%.2f\t\
//           %.2f\t",
//           TIME_READ_SEC,
//           TIME_TOTAL_SEC
//           );
//    printf("\n");
//    printf("\n");

    //globalStatFile << "USTITCH_MODE" <<  "=" << mapmode[ALGOMODE].c_str() << endl;
    globalStatFile << "V_USTITCH_" << mapmode[ALGOMODE].c_str() <<  "=" << V_ustitch << endl;
    globalStatFile << "PERCENT_C_SAVED_" << mapmode[ALGOMODE].c_str() <<  "=" << percent_saved_c << "%" << endl;
    globalStatFile << "PERCENT_C_GAP_WITH_UB_" << mapmode[ALGOMODE].c_str() <<  "=" << upperbound - percent_saved_c << "%" << endl;
    globalStatFile << "C_USTITCH_" << mapmode[ALGOMODE].c_str() <<   "=" <<C_ustitch << endl;
    globalStatFile << "BITSKMER_USTITCH_" << mapmode[ALGOMODE].c_str() <<  "=" <<ustitchBitsPerKmer << endl;
    globalStatFile << "TIME_READINPUT_SEC_" << mapmode[ALGOMODE].c_str() <<  "=" <<TIME_READ_SEC << endl;
    globalStatFile << "TIME_TOTAL_SEC_" << mapmode[ALGOMODE].c_str() <<  "=" <<TIME_TOTAL_SEC << endl;

    globalStatFile.close();

     cout << "------------ UST completed successfully. Output is in file "<<OUTPUT_FILENAME << " ------------"<<endl;
    cout << "Total number of unique "<<K<<"-mers " <<  "= " << numKmers << endl;
     cout << "Lower bound on any SPSS representation " <<  "= " << (charLowerbound)*1.0/numKmers << " neucleotide/k-mer"<< endl;
    cout << "Size of UST output" <<  "= " <<ustitchBitsPerKmer/2.0 << " neucleotide/k-mer"<< endl;
     //cout << "Gap with lower bound" << "= " << upperbound - percent_saved_c << "%" << endl;

    cout << "------------------------------------------------------"<<endl;


    return EXIT_SUCCESS;
}
