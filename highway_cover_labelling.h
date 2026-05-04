#ifndef HGHWAY_LABELING_H_
#define HGHWAY_LABELING_H_

#include <map>
#include "libraries.h"

class HighwayLabelling {
 public:
  HighwayLabelling(NetworKit::Graph &g, int l, int ordering_type, int changes, int dyntype);
    void ConstructDirWeighHL();
    void FarhanConstruction();
    void ConstructUnweightedNaiveLabelling();
    void ConstructWeightedNaiveLabelling();
    void ConstructDirWeighNL();
  int GetNumberOfNodes(){return V;};
  long LabellingSize();
  long DirectedLabellingSize();
  long NaiveLabellingSize();
  long DirectedNaiveLabellingSize();
  dist min(dist a, dist b);

  // Returns distance vetween vertices v and w if they are connected.
  dist QueryDistance(vertex s, vertex t);
  bool QueryDistanceBound(vertex s, vertex t, dist upperbound);
  bool DirectedQueryDistanceBound(vertex s, vertex t, dist ub);
    dist DirectedQueryDistance(vertex s, vertex t);
  dist NaiveQueryDistance(vertex s, vertex t);
  dist DirectedNaiveQueryDistance(vertex s, vertex t);
  dist BFSQuery(vertex s, vertex t);
    dist DijkstraQuery(vertex s, vertex t);
    dist DirectedDijkstra(vertex s, vertex t);
    void GetLandmarks(std::vector<vertex>& lndmrks);
    void GetIncrementalLandmarks(std::vector<vertex>& lndmrks);
  void SetLandmarks(std::vector<vertex>& lndmrks);
  void StoreIndex(std::string filename);
  dist SPQuery(vertex s, vertex t);
  // INCREMENTAL
  void AddLandmarkUnweighted(vertex r);
    void AddLandmarkDirected(vertex r);

  // DECREMENTAL
  void RemoveLandmarkUnweighted(vertex r);
    void RemoveLandmarkDirected(vertex r);
  std::vector<vertex> reverse_ordering;
  vertex L;
  std::set<vertex> landmarks;
    std::set<vertex> landmark_pool;


private:
  vertex V;  // total number of vertices
  vertex E; // total number of edges
   // total number of landmarks
  NetworKit::Graph &graph;
    std::vector<std::vector<std::pair<vertex,dist>>> landmarks_distances;
    std::vector<std::vector<std::pair<vertex,dist>>> in_landmarks_distances;
    std::vector<std::vector<std::pair<vertex,dist>>> out_landmarks_distances;
  std::unordered_map<vertex, std::unordered_map<vertex, dist>> highway;
  std::vector<vertex> ordering;
  std::vector<std::vector<dist>> naive_labeling;
  std::vector<std::vector<dist>> in_naive_labeling;
  std::vector<std::vector<dist>> out_naive_labeling;
  std::vector<vertex> landmarks_incremental;
    // temp structures
    std::vector<bool> settled;
    std::vector<dist> dij_distances;
    std::vector<std::pair<dist,bool>> dij_distances_pair;
    std::vector<bool> is_landmark;
};

template<typename Iter, typename RandomGenerator>
Iter select_randomly(Iter start, Iter end, RandomGenerator& g) {
    std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
    std::advance(start, dis(g));
    return start;
}

template<typename Iter>
Iter select_randomly(Iter start, Iter end) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    return select_randomly(start, end, gen);
}

HighwayLabelling::HighwayLabelling(NetworKit::Graph &g, int l, int ordering_type, int changes, int dyntype)
    : graph(g) {
    V = 0;
    E = 0;
    L = l;
    graph = g;

    V = this->graph.numberOfNodes();
    
    E = this->graph.numberOfEdges();
    auto *ordering_rank = new std::pair<double, vertex>[graph.numberOfNodes()];
    this->ordering.resize(graph.numberOfNodes());
    this->reverse_ordering.resize(graph.numberOfNodes());

    this->graph.parallelForNodes([&](vertex i) {
        assert(graph.hasNode(i));
        this->ordering[i] = null_vertex;
        this->reverse_ordering[i] = null_vertex;
        ordering_rank[i] = {0, i};
    });


    double centr_time = 0.0;

    if (ordering_type == 0) {
        INFO("BY DEGREE");
        NetworKit::DegreeCentrality *rank = new NetworKit::DegreeCentrality(graph);
        rank->run();
        this->graph.forNodes([&](vertex i) {
            assert(graph.hasNode(i));
            ordering_rank[i] = std::make_pair(rank->score(i), i);
        });
        delete rank;
    }

    if (ordering_type == 1) {
        INFO("RANDOM ORDERING");
        std::vector<vertex> t(graph.numberOfNodes());
        for (vertex i = 0; i < graph.numberOfNodes(); i++) t[i] = i;
        std::shuffle(t.begin(), t.end(), std::mt19937());
        for (vertex i = 0; i < graph.numberOfNodes(); i++) ordering_rank[i] = std::make_pair(t[i], i);
    }

    if (ordering_type == 2) {
        INFO("BY APX BETW");
        mytimer local_constr_timer;
        double max_time = 30.0;
        double cumulative_time = 0.0;
        double fract = 0.33;
        double n_samples = round(std::pow((double) graph.numberOfNodes(), fract));

        while (cumulative_time < max_time && n_samples < (double) graph.numberOfNodes()) {
            local_constr_timer.restart();

            std::cout << "fract: " << fract << " " << n_samples << " SAMPLES\n";
            NetworKit::EstimateBetweenness *rank = new NetworKit::EstimateBetweenness(graph, n_samples, false, true);

            rank->run();

            this->graph.forNodes([&](vertex i) {
                assert(graph.hasNode(i));
                assert(i<graph.numberOfNodes());
                ordering_rank[i] = std::make_pair(rank->score(i), i);
            });
            delete rank;
            cumulative_time += local_constr_timer.elapsed();
            n_samples *= 2;
        }
    }
    if (ordering_type == 3) {
        srandom(2024);
        INFO("BY DISTANCE-" + std::to_string((vertex) log2(V)) +" BOUNDED DOMINATING SET");
        std::set<vertex> to_cover;
        for (vertex v: graph.nodeRange()) to_cover.insert(v);
        vertex landmark_counter = 0;
        while (!to_cover.empty() && landmark_counter < L) {
            std::queue<std::pair<vertex, dist> > que;
            auto it = to_cover.begin();
            std::advance(it, random() % to_cover.size());
            vertex source = *it;
            landmark_counter++;
            que.push(std::make_pair(source, 0));
            while (!que.empty()) {
                auto p = que.front();
                que.pop();
                if (p.second > ((vertex) log2(V))) break;
                ordering_rank[p.first] = std::make_pair(to_cover.size(), p.first);
                to_cover.erase(p.first);
                for (vertex v: graph.neighborRange(p.first)) {
                    if (to_cover.find(v) != to_cover.end())
                        que.push(std::make_pair(v, p.second + 1));
                }
            }
            //ordering_rank[source] = std::make_pair(graph.numberOfNodes()+1+landmark_counter, source);
        }
    }

    if (ordering_type == 4) {
        INFO("BY APPROX CLOSENESS");
        NetworKit::ApproxCloseness *rank = new NetworKit::ApproxCloseness(graph, 1000);
        rank->run();
        this->graph.forNodes([&](vertex i) {
            assert(graph.hasNode(i));
            ordering_rank[i] = std::make_pair(rank->score(i), i);
        });
        delete rank;
    }

    std::sort(ordering_rank, ordering_rank + graph.numberOfNodes(),
              [](const std::pair<double, vertex> &a, const std::pair<double, vertex> &b) {
                  if (a.first == b.first)
                      return a.second > b.second;
                  else {
                      return a.first > b.first;
                  }
              });
    for (size_t count = 0; count < graph.numberOfNodes(); count++) {
        this->reverse_ordering[count] = ordering_rank[count].second;
        this->ordering[ordering_rank[count].second] = count;
    }
    landmarks.clear();
    is_landmark.resize(V);
    if(dyntype == 1) {
        for (vertex i = 0; i < L ; i++) {
            landmarks.insert(reverse_ordering[i]);
            is_landmark[reverse_ordering[i]] = true;
        }
        for(vertex i = L; i < L+changes; i++){
            landmarks_incremental.push_back(reverse_ordering[i]);
        }
    }
    else if(dyntype == 2) {
        for (vertex i = 0; i < L ; i++) {
            landmarks.insert(reverse_ordering[i]);
            is_landmark[reverse_ordering[i]] = true;
        }
    }
    else if(dyntype == 3) {
         for (vertex i = 0; i < L ; i++) {
             landmarks.insert(reverse_ordering[i]);
             is_landmark[reverse_ordering[i]] = true;
         }

         for(vertex i = L; i < L+changes/2; i++){
             landmarks_incremental.push_back(reverse_ordering[i]);
         }
    }
    else if(dyntype == 4) {
        for (vertex i = 0; i < L ; i++) {
            landmark_pool.insert(reverse_ordering[i]);
        }
        for (vertex i = 0; i < L/2 ; i++) {
            vertex landmark = *select_randomly(landmark_pool.begin(), landmark_pool.end());
            landmarks.insert(landmark);
            is_landmark[reverse_ordering[i]] = true;
            landmark_pool.erase(landmark_pool.find(landmark));
        }
    }
    delete[] ordering_rank;
    //L = L/2;
    for(const auto &v: landmarks) std::cout << v << " ";
    std::cout << "\n";
    settled.resize(V);
    dij_distances.resize(V, null_distance);
    dij_distances_pair.resize(V, std::make_pair(null_distance,true));
    settled.resize(V,false);
}

void HighwayLabelling::GetLandmarks(std::vector<vertex> &lndmrks) {
    for(auto v: landmarks){
        lndmrks.push_back(v);
    }
}

void HighwayLabelling::GetIncrementalLandmarks(std::vector<vertex> &lndmrks) {
    for(auto v: landmarks_incremental){
        lndmrks.push_back(v);
    }
}

void HighwayLabelling::SetLandmarks(std::vector<vertex> &lndmrks) {
    landmarks.clear();
    L = lndmrks.size();
    for(auto v: graph.nodeRange()) {
        is_landmark[v] = false;
    }
    for(auto v: lndmrks){
        landmarks.insert(v);
        is_landmark[v] = true;
    }
}

long HighwayLabelling::LabellingSize() {
  long size = 0;
  for (int i = 0; i < V; i++) {
    for (int j = 0; j < landmarks_distances[i].size(); j++) {
      if(landmarks_distances[i][j].second != null_distance)
        size+=2;
    }
  }

  for (const vertex & v: landmarks) {
    for (const vertex & w: landmarks) {
      if(highway[v][w] != null_distance)
       size++;
    }
  }

  return size;
}

long HighwayLabelling::DirectedLabellingSize() {
  long size = 0;
  for (int i = 0; i < V; i++) {
    for (int j = 0; j < in_landmarks_distances[i].size(); j++) {
      if(in_landmarks_distances[i][j].second != null_distance)
        size+=2;
    }
    for (int j = 0; j < out_landmarks_distances[i].size(); j++) {
      if(out_landmarks_distances[i][j].second != null_distance)
        size+=2;
    }
  }
  auto temp = size;
  for (const vertex & v : landmarks) {
    auto it_v = highway.find(v);
    if (it_v == highway.end()) continue;

    for (const vertex & w : landmarks) {
        auto it_w = it_v->second.find(w);
        if (it_w != it_v->second.end() && it_w->second != null_distance) {
                size++;
            }
        }
    } 
    std::cout << "highway size " << size- temp << "\n"; 
  return size;
}

long HighwayLabelling::NaiveLabellingSize() {
    long size = 0;
    for (int i = 0; i < V; i++) {
        for (int j = 0; j < L; j++) {
            if(naive_labeling[i][j] != null_distance)
                size++;
        }
    }

    return size;
}

long HighwayLabelling::DirectedNaiveLabellingSize() {
    long size = 0;
    for (int i = 0; i < V; i++) {
        for (int j = 0; j < L; j++) {
            if(in_naive_labeling[i][j] != null_distance)
                size++;
            if(out_naive_labeling[i][j] != null_distance)
                size++;
        }
    }

    return size;
}

void HighwayLabelling::FarhanConstruction(){
    landmarks_distances.resize(V);

    for(int i = 0; i < V; i++) {
        landmarks_distances[i].clear();
    }
    highway.clear();
    for(const auto & l1: landmarks) {
        highway[l1] = std::unordered_map<vertex, dist> ();
    }
    ProgressStream hl_bar(L);

    hl_bar.label() << "Unweighted highway labeling construction";
    std::vector<vertex> reached_vertices;
        for(const vertex & b : landmarks) {
            highway[b][b] = 0;
            dij_distances[b] = 0;
            landmarks_distances[b].emplace_back(b,0);
            reached_vertices.push_back(b);
            std::queue<vertex> q[2];
            q[0].push(b); q[0].push(-1);
            int use=0;

            while (!q[0].empty()){
                vertex u = q[use].front();
                q[use].pop();
                if(u==-1){
                    use = 1-use;
                    q[use].push(-1);
                    continue;
                }
                for(auto v: graph.neighborRange(u)){
                    if(dij_distances[v]!= null_distance) continue;
                    reached_vertices.push_back(v);
                    dij_distances[v] = dij_distances[u] +1;
                    if(use == 1 || is_landmark[v]){
                        q[1].push(v);
                    } else{
                        q[0].push(v);
                        std::vector<std::pair<vertex,dist>>::iterator insertion_index;
                        insertion_index = std::upper_bound(landmarks_distances[v].begin(),
                                                                landmarks_distances[v].end(), std::make_pair(b,dij_distances[v]));
                        
                        landmarks_distances[v].insert(insertion_index, std::make_pair(b,dij_distances[v]));
                    }
                }
            }
            while (!q[1].empty()) {
                vertex u = q[1].front();
                q[1].pop();
                if(u == -1){continue;}
                for (auto v: graph.neighborRange(u)) {
                    if(dij_distances[v]!= null_distance) continue;
                    dij_distances[v] = dij_distances[u] +1;
                    reached_vertices.push_back(v);
                    
                    q[1].push(v);
                }
            }
            for(const auto & v: landmarks){
                highway[b][v] = dij_distances[v];
                highway[v][b] = dij_distances[v];
                if (dij_distances[v] == null_distance) {
                    throw new std::runtime_error("experiment fails");
                }
            }
            for(const auto & v: reached_vertices){
                dij_distances[v] = null_distance;
            }
            reached_vertices.clear();

            ++hl_bar;
        }
}

void HighwayLabelling::ConstructDirWeighHL() {
    // Initialization
    in_landmarks_distances.resize(V);
    out_landmarks_distances.resize(V);
    for(int i = 0; i < V; i++) {
        in_landmarks_distances[i].clear();
        out_landmarks_distances[i].clear();
    }
    highway.clear();
    for(const auto & l1: landmarks) {
        highway[l1] = std::unordered_map<vertex, dist> ();
    }

    // Start computing Highway Labelling (HL)
    ProgressStream hl_bar(L);

    hl_bar.label() << "Directed weighted highway labeling construction";
    std::vector<vertex> reached_vertices;
    // FORWARD
    for(const vertex & b : landmarks) {
        highway[b][b] = 0;
        dij_distances_pair[b] = std::make_pair(0,true);
        reached_vertices.push_back(b);
        out_landmarks_distances[b].emplace_back(b,0);
        std::priority_queue<std::pair<std::pair<dist, bool>, vertex>, std::vector<std::pair<std::pair<dist, bool>, vertex>>,
                PQFlagComparator> pq;
        pq.push(std::make_pair(std::make_pair(0,true),b));
        while (!pq.empty()) {
            vertex v = pq.top().second;
            std::pair<dist,bool> f = pq.top().first;
            pq.pop();
            if(settled[v]) continue;
            if (is_landmark[v] && v != b) {
                highway[b][v] = f.first;
                f.second = false;
            }
            else{
                if(f.second && v != b){
                    out_landmarks_distances[v].emplace_back(b,f.first);
                }
            }
            settled[v] = true;
            for (auto w: graph.neighborRange(v)) {
                    std::pair<dist,bool> temp = {dij_distances_pair[v].first + graph.weight(v, w), f.second};
                    if (dij_distances_pair[w] > temp) {
                        dij_distances_pair[w] = temp;
                        reached_vertices.push_back(w);
                        pq.push(std::make_pair(dij_distances_pair[w], w));
                    }
                }
            
        }
        for (const auto &v: reached_vertices) {
            dij_distances_pair[v] = std::make_pair(null_distance,true);
            settled[v] = false;
        }
        reached_vertices.clear();
        // REVERSE
        dij_distances_pair[b] = std::make_pair(0,true);
        reached_vertices.push_back(b);
        in_landmarks_distances[b].emplace_back(b,0);
        while(!pq.empty()) pq.pop();
        pq.push(std::make_pair(std::make_pair(0,true),b));
        while (!pq.empty()) {
            vertex v = pq.top().second;
            std::pair<dist,bool> f = pq.top().first;
            pq.pop();
            if(settled[v]) continue;
            if (is_landmark[v] && v != b) {
                highway[v][b] = f.first;
                f.second = false;
            }
            else{
                if(f.second && v != b){
                    in_landmarks_distances[v].emplace_back(b,f.first);
                }
            }
            settled[v] = true;

            for (auto w: graph.inNeighborRange(v)) {
                    std::pair<dist,bool> temp = {dij_distances_pair[v].first + graph.weight(w,v), f.second};
                    if (dij_distances_pair[w] > temp) {
                        dij_distances_pair[w] = temp;
                        reached_vertices.push_back(w);
                        pq.push(std::make_pair(dij_distances_pair[w], w));
                    }
                }
            
        }
        for (const auto &v: reached_vertices) {
            dij_distances_pair[v] = std::make_pair(null_distance,true);
            settled[v] = false;
        }
        reached_vertices.clear();
        ++hl_bar;
    }
    // for (const auto& [u, inner_map] : highway) {
    //     std::cout << "From vertex " << u << ":\n";

    //     for (const auto& [v, d] : inner_map) {
    //         std::cout << "  -> " << v << " (dist = " << d << ")\n";
    //     }
    // }
    
    //  for (int v = 0; v < V; v++) {
    //     std::cout << "L-In(" << v << "): [";
    //     for (int i = 0; i < in_landmarks_distances[v].size(); i++) {
    //         std::cout << "(" << in_landmarks_distances[v][i] << "," << in_distances[v][i] << "), ";
    //     }
    //     std::cout << "]" << std::endl;
    //     std::cout << "L-Out(" << v << "): [";
    //     for (int i = 0; i < out_landmarks_distances[v].size(); i++) {
    //         std::cout << "(" << out_landmarks_distances[v][i] << "," << out_distances[v][i] << "), ";
    //     }
    //     std::cout << "]" << std::endl;
    // }
}

void HighwayLabelling::ConstructUnweightedNaiveLabelling() {
    // Initialization
    naive_labeling.resize(V);
    for(vertex i = 0; i < V; i++) {
        naive_labeling[i].resize(L);
        for(vertex j = 0; j < L; j++)
            naive_labeling[i][j] = null_distance;
    }


    // Start computing Naive Labelling
    ProgressStream na_bar(L);

    na_bar.label() << "Unweighted naive labeling construction";
    for(vertex s = 0; s < L; s++){
        dist *P = new dist[V];
        for(vertex j = 0; j < V; j++)
            P[j] = null_distance;
        std::queue<vertex> que;
        que.push(reverse_ordering[s]);
        naive_labeling[reverse_ordering[s]][s] = 0; P[reverse_ordering[s]] = 0;

        while(!que.empty()){
            vertex v = que.front();
            que.pop();
            naive_labeling[v][s] = P[v];

            for(vertex w: graph.neighborRange(v)){
                if(P[w]==null_distance){
                    P[w] = P[v] + 1;
                    que.push(w);
                }
            }
        }
    ++na_bar;
    }

}

void HighwayLabelling::ConstructWeightedNaiveLabelling() {
    // Initialization
    naive_labeling.resize(V);
    for(vertex i = 0; i < V; i++) {
        naive_labeling[i].resize(L);
        for(vertex j = 0; j < L; j++)
            naive_labeling[i][j] = null_distance;
    }
    ProgressStream na_bar(L);

    na_bar.label() << "Weighted naive labeling construction";
    for(vertex b = 0; b < L; b++) {
        bool* settled = new bool[V];

        naive_labeling[reverse_ordering[b]][b] = 0;
        std::priority_queue<std::pair<dist,vertex>, std::vector<std::pair<dist,vertex>>,
                PQComparator> pq;
        pq.push(std::make_pair(0,reverse_ordering[b]));
        while(!pq.empty()){
            vertex v = pq.top().second;
            pq.pop();

            for(auto w: graph.neighborRange(v)){
                if(naive_labeling[w][b] > naive_labeling[v][b] + graph.weight(v,w)){
                    naive_labeling[w][b] = naive_labeling[v][b] + (dist)graph.weight(v,w);

                    settled[w] = settled[v];
                    pq.push(std::make_pair(naive_labeling[w][b], w));
                }

            }
        }
        ++na_bar;
        delete [] settled;
    }
    
}

void HighwayLabelling::ConstructDirWeighNL() {
    // Initialization
    in_naive_labeling.resize(V);
    out_naive_labeling.resize(V);
    for(vertex i = 0; i < V; i++) {
        in_naive_labeling[i].resize(L);
        out_naive_labeling[i].resize(L);
        for(vertex j = 0; j < L; j++){
            in_naive_labeling[i][j] = null_distance;
            out_naive_labeling[i][j] = null_distance;
            }
    }
    ProgressStream na_bar(L);

    na_bar.label() << "Directed weighted naive labeling construction";
    for(vertex b = 0; b < L; b++) {
        bool* settled = new bool[V];
        std::vector<vertex> reached_vertices;
        out_naive_labeling[reverse_ordering[b]][b] = 0;
        std::priority_queue<std::pair<dist,vertex>, std::vector<std::pair<dist,vertex>>,
                PQComparator> pq;
        pq.push(std::make_pair(0,reverse_ordering[b]));
        while(!pq.empty()){
            vertex v = pq.top().second;
            pq.pop();
            reached_vertices.push_back(v);
            for(auto w: graph.neighborRange(v)){
                if(out_naive_labeling[w][b] > out_naive_labeling[v][b] + graph.weight(v,w)){
                    out_naive_labeling[w][b] = out_naive_labeling[v][b] + (dist)graph.weight(v,w);

                    settled[w] = settled[v];
                    pq.push(std::make_pair(out_naive_labeling[w][b], w));
                }

            }
        }

        for(const auto & v: reached_vertices){
            settled[v] = false;
        }
        reached_vertices.clear();
        in_naive_labeling[reverse_ordering[b]][b] = 0;
        while(!pq.empty()) pq.pop();
        pq.push(std::make_pair(0,reverse_ordering[b]));
        while(!pq.empty()){
            vertex v = pq.top().second;
            pq.pop();

            for(auto w: graph.inNeighborRange(v)){
                if(in_naive_labeling[w][b] > in_naive_labeling[v][b] + graph.weight(w,v)){
                    in_naive_labeling[w][b] = in_naive_labeling[v][b] + (dist)graph.weight(w,v);

                    settled[w] = settled[v];
                    pq.push(std::make_pair(in_naive_labeling[w][b], w));
                }

            }
        }

        ++na_bar;
        delete [] settled;
    }
    
}

dist HighwayLabelling::min(dist a, dist b) {
  return (a < b) ? a : b;
}

bool HighwayLabelling::QueryDistanceBound(vertex s, vertex t, dist upperbound) {
    dist m = null_distance;

    // Replace VLAs with std::vector
    std::vector<vertex> uni1(landmarks_distances[s].size(), 0);
    std::vector<vertex> uni2(landmarks_distances[t].size(), 0);

    vertex i = 0, j = 0, i1 = 0, j1 = 0;

    while (i < landmarks_distances[s].size() && j < landmarks_distances[t].size()) {
        if (landmarks_distances[s][i].first < landmarks_distances[t][j].first) {
            uni1[i1] = i; i++; i1++;
        } else if (landmarks_distances[t][j].first < landmarks_distances[s][i].first) {
            uni2[j1] = j; j++; j1++;
        } else {
            m = std::min(m, landmarks_distances[s][i].second + landmarks_distances[t][j].second);
            if (m <= upperbound) return true;
            i++; j++;
        }
    }

    while (i < landmarks_distances[s].size()) {
        uni1[i1] = i; i++; i1++;
    }

    while (j < landmarks_distances[t].size()) {
        uni2[j1] = j; j++; j1++;
    }

    i = 0;
    while (i < i1) {
        for (j = 0; j < landmarks_distances[t].size(); j++) {
            m = std::min(m,
                landmarks_distances[s][uni1[i]].second +
                highway[landmarks_distances[s][uni1[i]].first][landmarks_distances[t][j].first] +
                landmarks_distances[t][j].second);
            if (m <= upperbound) return true;
        }
        i++;
    }

    j = 0;
    while (j < j1) {
        for (i = 0; i < landmarks_distances[s].size(); i++) {
            m = std::min(m,
                landmarks_distances[s][i].second +
                highway[landmarks_distances[s][i].first][landmarks_distances[t][uni2[j]].first] +
                landmarks_distances[t][uni2[j]].second);
            if (m <= upperbound) return true;
        }
        j++;
    }

    return m <= upperbound;
}

bool HighwayLabelling::DirectedQueryDistanceBound(vertex s, vertex t, dist ub) {
    dist m = null_distance;
    
    // Replace VLAs with std::vector
    std::vector<vertex> uni1(in_landmarks_distances[s].size(), 0);
    std::vector<vertex> uni2(out_landmarks_distances[t].size(), 0);

    vertex i = 0, j = 0, i1 = 0, j1 = 0;

    while (i < in_landmarks_distances[s].size() && j < out_landmarks_distances[t].size()) {
        if (in_landmarks_distances[s][i].first < out_landmarks_distances[t][j].first) {
            uni1[i1] = i; i++; i1++;
        } else if (out_landmarks_distances[t][j].first < in_landmarks_distances[s][i].first) {
            uni2[j1] = j; j++; j1++;
        } else {
            if(m > in_landmarks_distances[s][i].second + out_landmarks_distances[t][j].second){
                m = in_landmarks_distances[s][i].second + out_landmarks_distances[t][j].second;
                if (m <= ub) return true;
            }
            i++; j++;
        }
    } 
    while (i < in_landmarks_distances[s].size()) {
        uni1[i1] = i; i++; i1++;
    }

    while (j < out_landmarks_distances[t].size()) {
        uni2[j1] = j; j++; j1++;
    }

    i = 0;
    while (i < i1) {
        for (j = 0; j < out_landmarks_distances[t].size(); j++) {
            auto it_outer = highway.find(in_landmarks_distances[s][uni1[i]].first);
            if (it_outer != highway.end()) {
                auto it_inner = it_outer->second.find(out_landmarks_distances[t][j].first);
                if (it_inner != it_outer->second.end()) {
                    if(m > in_landmarks_distances[s][uni1[i]].second +
                        it_inner->second +
                        out_landmarks_distances[t][j].second){
                            m = in_landmarks_distances[s][uni1[i]].second +
                        it_inner->second +
                        out_landmarks_distances[t][j].second;
                            if (m <= ub) return true;
                        }
                }
            }
        }
        i++;
    }
    j = 0;
    while (j < j1) {
        for (i = 0; i < in_landmarks_distances[s].size(); i++) {

            auto it_outer = highway.find(in_landmarks_distances[s][i].first);
            if (it_outer != highway.end()) {

                auto it_inner = it_outer->second.find(out_landmarks_distances[t][uni2[j]].first);
                if (it_inner != it_outer->second.end()) {
                    if(m > in_landmarks_distances[s][i].second +
                                    it_inner->second +
                                    out_landmarks_distances[t][uni2[j]].second){
                                        m = in_landmarks_distances[s][i].second +
                                    it_inner->second +
                                    out_landmarks_distances[t][uni2[j]].second;
                            if (m <= ub) return true;
                                    }
                }
            }
        }
        j++;
    }

    return m <= ub;
}

dist HighwayLabelling::QueryDistance(vertex s, vertex t) {
    dist m = null_distance;

    std::vector<vertex> uni1(landmarks_distances[s].size(), 0);
    std::vector<vertex> uni2(landmarks_distances[t].size(), 0);

    vertex i = 0, j = 0, i1 = 0, j1 = 0;

    while (i < landmarks_distances[s].size() && j < landmarks_distances[t].size()) {
        if (landmarks_distances[s][i].first < landmarks_distances[t][j].first) {
            uni1[i1] = i; i++; i1++;
        } else if (landmarks_distances[t][j].first < landmarks_distances[s][i].first) {
            uni2[j1] = j; j++; j1++;
        } else {
            m = std::min(m, landmarks_distances[s][i].second + landmarks_distances[t][j].second);
            i++; j++;
        }
    }

    while (i < landmarks_distances[s].size()) {
        uni1[i1] = i; i++; i1++;
    }

    while (j < landmarks_distances[t].size()) {
        uni2[j1] = j; j++; j1++;
    }

    i = 0;
    while (i < i1) {
        for (j = 0; j < landmarks_distances[t].size(); j++) {
            m = std::min(m,
                landmarks_distances[s][uni1[i]].second +
                highway[landmarks_distances[s][uni1[i]].first][landmarks_distances[t][j].first] +
                landmarks_distances[t][j].second);
        }
        i++;
    }

    j = 0;
    while (j < j1) {
        for (i = 0; i < landmarks_distances[s].size(); i++) {
            m = std::min(m,
                landmarks_distances[s][i].second +
                highway[landmarks_distances[s][i].first][landmarks_distances[t][uni2[j]].first] +
                landmarks_distances[t][uni2[j]].second);
        }
        j++;
    }

    return m;
}


dist HighwayLabelling::DirectedQueryDistance(vertex s, vertex t) {
    dist m = null_distance;
    
    // Replace VLAs with std::vector
    std::vector<vertex> uni1(in_landmarks_distances[s].size(), 0);
    std::vector<vertex> uni2(out_landmarks_distances[t].size(), 0);

    vertex i = 0, j = 0, i1 = 0, j1 = 0;

    while (i < in_landmarks_distances[s].size() && j < out_landmarks_distances[t].size()) {
        if (in_landmarks_distances[s][i].first < out_landmarks_distances[t][j].first) {
            uni1[i1] = i; i++; i1++;
        } else if (out_landmarks_distances[t][j].first < in_landmarks_distances[s][i].first) {
            uni2[j1] = j; j++; j1++;
        } else {
            m = std::min(m, in_landmarks_distances[s][i].second + out_landmarks_distances[t][j].second);
            i++; j++;
        }
    } 
    while (i < in_landmarks_distances[s].size()) {
        uni1[i1] = i; i++; i1++;
    }

    while (j < out_landmarks_distances[t].size()) {
        uni2[j1] = j; j++; j1++;
    }

    i = 0;
    while (i < i1) {
        for (j = 0; j < out_landmarks_distances[t].size(); j++) {
            auto it_outer = highway.find(in_landmarks_distances[s][uni1[i]].first);
            if (it_outer != highway.end()) {
                auto it_inner = it_outer->second.find(out_landmarks_distances[t][j].first);
                if (it_inner != it_outer->second.end()) {
                    m = std::min(m,
                        in_landmarks_distances[s][uni1[i]].second +
                        it_inner->second +
                        out_landmarks_distances[t][j].second);
                }
            }
        }
        i++;
    }
    j = 0;
    while (j < j1) {
        for (i = 0; i < in_landmarks_distances[s].size(); i++) {

            auto it_outer = highway.find(in_landmarks_distances[s][i].first);
            if (it_outer != highway.end()) {

                auto it_inner = it_outer->second.find(out_landmarks_distances[t][uni2[j]].first);
                if (it_inner != it_outer->second.end()) {

                    m = std::min(m,
                                    in_landmarks_distances[s][i].second +
                                    it_inner->second +
                                    out_landmarks_distances[t][uni2[j]].second);
                }
            }
        }
        j++;
    }

    return m;
}


dist HighwayLabelling::NaiveQueryDistance(vertex s, vertex t) {
    dist m = null_distance;
    for(vertex i = 0; i < L; i++){
            m = min(m, naive_labeling[s][i]+naive_labeling[t][i]);
    }
    return m;
}

dist HighwayLabelling::DirectedNaiveQueryDistance(vertex s, vertex t) {
    dist m = null_distance;
    for(vertex i = 0; i < L; i++){
            m = min(m, in_naive_labeling[s][i]+out_naive_labeling[t][i]);
    }
    return m;
}

dist HighwayLabelling::BFSQuery(vertex s, vertex t) {
    // Initialization
    dist *s_to_vertices = new dist[V];
    dist *t_to_vertices = new dist[V];
    for(vertex j = 0; j < V; j++) {
        s_to_vertices[j] = null_distance;
        t_to_vertices[j] = null_distance;
    }
    std::queue<vertex> que;
    que.push(s);
    s_to_vertices[s] = 0;

    while(!que.empty()){
        vertex v = que.front();
        que.pop();

        for(vertex w: graph.neighborRange(v)){
            if(s_to_vertices[w]==null_distance){
                s_to_vertices[w] = s_to_vertices[v] + 1;
                que.push(w);
            }
        }
    }

    que.push(t);
    t_to_vertices[t] = 0;

    while(!que.empty()){
        vertex v = que.front();
        que.pop();

        for(vertex w: graph.neighborRange(v)){
            if(t_to_vertices[w]==null_distance){
                t_to_vertices[w] = t_to_vertices[v] + 1;
                que.push(w);
            }
        }
    }

    dist m = null_distance;
    for(auto l: landmarks){
        m = min(s_to_vertices[l] + t_to_vertices[l], m);
    }

    delete [] s_to_vertices;
    delete [] t_to_vertices;
    return m;
}

dist HighwayLabelling::DijkstraQuery(vertex s, vertex t) {
    dist* s_distances = new dist[V];
    for(vertex i = 0; i < V; i++){
        s_distances[i] = null_distance;
    }
    s_distances[s] = 0;
    std::priority_queue<std::pair<dist,vertex>, std::vector<std::pair<dist,vertex>>,
            PQComparator> pq;
    pq.push(std::make_pair(0,s));
    while(!pq.empty()){
        vertex v = pq.top().second;
        pq.pop();
        for(auto w: graph.neighborRange(v)){
            if(s_distances[w] > s_distances[v] + graph.weight(v,w)){
                s_distances[w] = s_distances[v] + (dist)graph.weight(v,w);

                pq.push(std::make_pair(s_distances[w], w));
            }

        }
    }

    dist* t_distances = new dist[V];
    for(vertex i = 0; i < V; i++){
        t_distances[i] = null_distance;
    }
    t_distances[t] = 0;
    std::priority_queue<std::pair<dist,vertex>, std::vector<std::pair<dist,vertex>>,
            PQComparator> pqt;
    pqt.push(std::make_pair(0,t));
    while(!pqt.empty()){
        vertex v = pqt.top().second;
        pqt.pop();
        for(auto w: graph.neighborRange(v)){
            if(t_distances[w] > t_distances[v] + graph.weight(v,w)){
                t_distances[w] = t_distances[v] + (dist)graph.weight(v,w);

                pqt.push(std::make_pair(t_distances[w], w));
            }

        }
    }

//    auto s_to_vertices = NetworKit::Dijkstra(graph, s, false, false);
//    s_to_vertices.run();
//    auto t_to_vertices = NetworKit::Dijkstra(graph, t, false, false);
//    t_to_vertices.run();
//    dist m = null_distance;
//    for(vertex l = 0; l < L; l++){
//        m = min(s_to_vertices.distance(reverse_ordering[l]) + t_to_vertices.distance(reverse_ordering[l]), m);
//    }
    dist m = null_distance;
    vertex minhub = 0;
    for(const vertex & l: landmarks){
        if(m > s_distances[l] + t_distances[l]) minhub = l;
        m = min(m, s_distances[l] + t_distances[l]);
    }
    //std::cout << "minhub " << minhub << std::endl;
    delete [] s_distances;
    delete [] t_distances;
    return m;
}

dist HighwayLabelling::DirectedDijkstra(vertex s, vertex t) {
    dist* s_distances = new dist[V];
    for(vertex i = 0; i < V; i++){
        s_distances[i] = null_distance;
    }
    s_distances[s] = 0;
    std::set<vertex> rlndm;
    std::priority_queue<std::pair<dist,vertex>, std::vector<std::pair<dist,vertex>>,
            PQComparator> pq;
    pq.push(std::make_pair(0,s));
    while(!pq.empty()){
        vertex v = pq.top().second;
        dist delta = pq.top().first;
        pq.pop();
        if(delta > s_distances[v]) continue;
        if(is_landmark[v]) {
            rlndm.insert(v);
        }
        if(rlndm.size() == landmarks.size()) break;
        for(auto w: graph.neighborRange(v)){
            if(s_distances[w] > s_distances[v] + graph.weight(v,w)){
                s_distances[w] = s_distances[v] + (dist)graph.weight(v,w);
                pq.push(std::make_pair(s_distances[w], w));
            }

        }
    }

    dist* t_distances = new dist[V];
    for(vertex i = 0; i < V; i++){
        t_distances[i] = null_distance;
    }
    rlndm.clear();
    t_distances[t] = 0;
    std::priority_queue<std::pair<dist,vertex>, std::vector<std::pair<dist,vertex>>,
            PQComparator> pqt;
    pqt.push(std::make_pair(0,t));
    while(!pqt.empty()){
        vertex v = pqt.top().second;
        dist delta = pqt.top().first;
        pqt.pop();
        if(delta > t_distances[v]) continue;

        if(is_landmark[v]) {
            rlndm.insert(v);
        }
        if(rlndm.size() == landmarks.size()) break;
        for(auto w: graph.inNeighborRange(v)){
            if(t_distances[w] > t_distances[v] + graph.weight(w,v)){
                t_distances[w] = t_distances[v] + (dist)graph.weight(w,v);

                pqt.push(std::make_pair(t_distances[w], w));
            }

        }
    }

    dist m = null_distance;
    vertex v = null_vertex;
    dist from_s = null_distance;
    dist to_t = null_distance;
    for(const vertex & l: landmarks){
        if(m > s_distances[l] + t_distances[l]) {
            v = l;
            m = s_distances[l] + t_distances[l];
            from_s = s_distances[l];
            to_t = t_distances[l];
        }
        //m = min(m, s_distances[l] + t_distances[l]);
    }
    delete [] s_distances;
    delete [] t_distances;
    //std::cout << "hub " << v << " from s " << from_s << " to " << to_t << "\n";
    return m;
}

void HighwayLabelling::AddLandmarkUnweighted(vertex r) {
    // compute Highway distances from previous landmarks to r
    // add distances in Highway between r and landmarks directly connected to r (i.e. in L(r))

    landmarks.insert(r);
    is_landmark[r] = true;
    highway[r] = std::unordered_map<vertex, dist> ();
    for(size_t i = 0; i < landmarks_distances[r].size(); i++){
        highway[r][landmarks_distances[r][i].first] = landmarks_distances[r][i].second;
        highway[landmarks_distances[r][i].first][r] = landmarks_distances[r][i].second;
    }
    highway[r][r] = 0;
    // fill distances in Highway between r and landmarks not covering r (i.e. not in L(r))
    for(const vertex & l: landmarks){
        if(highway[r].find(l) == highway[r].end()){
            highway[r][l] = null_distance;
            for(size_t j = 0; j < landmarks_distances[r].size(); j++){
                highway[r][l] = std::min(highway[r][l],highway[r][landmarks_distances[r][j].first] + highway[landmarks_distances[r][j].first][l]);
            }
            highway[l][r] = highway[r][l];
        }
    }
    landmarks_distances[r].clear();
    landmarks_distances[r].emplace_back(r,0);
    // search from r to fill L a la Akiba (prune when reaching landmarks different from r, compute query to vertices not in L)
    std::vector<vertex> reached_vertices;

    dij_distances[r] = 0;
    reached_vertices.push_back(r);
    std::vector<vertex> reached_landmarks;
    std::queue<vertex> q;
    std::vector<vertex> obsolete_positions;
    for(auto w: graph.neighborRange(r)){
            if (is_landmark[w]) {
                 settled[w] = true;

                reached_landmarks.push_back(w);
                continue;
            }
            dij_distances[w] = dij_distances[r] + 1;
            q.push(w);
            reached_vertices.push_back(w);

        }
    while(!q.empty()){
        vertex v = q.front();
        q.pop();
        if(settled[v]) continue;
        //dist query_dist = QueryDistance(r,v);
        if(is_landmark[v] && v != r) {
            settled[v] = true;
            reached_landmarks.push_back(v);
            continue;
        }
        if(QueryDistanceBound(r,v,dij_distances[v])){
            settled[v] = true;
            continue;
        }

        //auto insertion_index = std::upper_bound(landmarks_distances[v].begin(), landmarks_distances[v].end(), r);
        vertex i = 0;
        for(; i < landmarks_distances[v].size(); i++) {
            if(landmarks_distances[v][i].first > r) {
                break;
            }
            if(landmarks_distances[v][i].second>=dij_distances[v]+highway[r][landmarks_distances[v][i].first]) {
                //lndm_rch_vrtx[landmarks_distances[v][i]].emplace_back(std::make_pair(distances[v][i],v),i);
                obsolete_positions.push_back(i);
            }
        }
        landmarks_distances[v].insert(landmarks_distances[v].begin() + i, std::make_pair(r,dij_distances[v]));
        i += 1;
        for(; i < landmarks_distances[v].size(); i++) {
            if(landmarks_distances[v][i].second>=dij_distances[v]+highway[r][landmarks_distances[v][i].first]){
                //lndm_rch_vrtx[landmarks_distances[v][i]].emplace_back(std::make_pair(distances[v][i],v),i);
                obsolete_positions.push_back(i);
            }
        }
            for (auto it = obsolete_positions.rbegin(); it != obsolete_positions.rend(); ++it) {
                const auto& j = *it;

                landmarks_distances[v].erase(landmarks_distances[v].begin() + j);
            }
        obsolete_positions.clear();

        for(auto w: graph.neighborRange(v)){
            if(dij_distances[w] == null_distance){
                dij_distances[w] = dij_distances[v] + 1;
                q.push(w);
                reached_vertices.push_back(w);
            }
        }
    }
    for(const auto & v: reached_landmarks) settled[v] = false;
    for(const auto & v: reached_vertices){
        dij_distances[v] = null_distance;
        settled[v] = false;
    }

    
    L++;
    
}

void HighwayLabelling::AddLandmarkDirected(vertex r) {
    // compute Highway distances from previous landmarks to r
    // add distances in Highway between r and landmarks directly connected to r (i.e. in L(r))
    landmarks.insert(r);
    is_landmark[r] = true;
    highway[r] = std::unordered_map<vertex, dist> ();
    for(size_t i = 0; i < out_landmarks_distances[r].size(); i++){
        highway[out_landmarks_distances[r][i].first][r] = out_landmarks_distances[r][i].second;
    }
    for(size_t i = 0; i < in_landmarks_distances[r].size(); i++){
        highway[r][in_landmarks_distances[r][i].first] = in_landmarks_distances[r][i].second;
    }

    highway[r][r] = 0;
    // fill distances in Highway between r and landmarks not covering r (i.e. not in L(r))
    for (const vertex & l : landmarks) {

        // ---------- r -> l ----------
        {
            auto it_rl = highway[r].find(l);
            dist best = null_distance;

            if (it_rl != highway[r].end()) {
                best = it_rl->second;
            }

            for (size_t j = 0; j < in_landmarks_distances[r].size(); j++) {
                vertex mid = in_landmarks_distances[r][j].first;

                auto it_rm = highway[r].find(mid);
                auto it_ml = highway[mid].find(l);

                if (it_rm != highway[r].end() &&
                    it_ml != highway[mid].end()) {

                    best = std::min(best, it_rm->second + it_ml->second);
                }
            }

            if (best != null_distance) {
                highway[r][l] = best;
            }
        }

        // ---------- l -> r ----------
        {
            auto it_lr = highway[l].find(r);
            dist best = null_distance;

            if (it_lr != highway[l].end()) {
                best = it_lr->second;
            }

            for (size_t j = 0; j < out_landmarks_distances[r].size(); j++) {
                vertex mid = out_landmarks_distances[r][j].first;

                auto it_lm = highway[l].find(mid);
                auto it_mr = highway[mid].find(r);

                if (it_lm != highway[l].end() &&
                    it_mr != highway[mid].end()) {

                    best = std::min(best, it_lm->second + it_mr->second);
                }
            }

            if (best != null_distance) {
                highway[l][r] = best;
            }
        }
    }

    // FORWARD

    out_landmarks_distances[r].clear();
    out_landmarks_distances[r].emplace_back(r,0);
    // search from r to fill L a la Akiba (prune when reaching landmarks different from r, compute query to vertices not in L)
    std::vector<vertex> reached_vertices;
    dij_distances_pair[r] = std::make_pair(0,true);
    reached_vertices.push_back(r);
    std::priority_queue<std::pair<std::pair<dist, bool>, vertex>, std::vector<std::pair<std::pair<dist, bool>, vertex>>,
                PQFlagComparator> pq;
    std::vector<vertex> obsolete_positions;

    for(auto w: graph.neighborRange(r)){
        if (dij_distances_pair[w].first > dij_distances_pair[r].first + graph.weight(r, w)) {
                        dij_distances_pair[w] = std::make_pair(dij_distances_pair[r].first + static_cast<dist>(graph.weight(r, w)), true);
                        reached_vertices.push_back(w);
                        pq.push(std::make_pair(dij_distances_pair[w], w));
                    }
    }
    settled[r] = true;
    while(!pq.empty()){
        vertex v = pq.top().second;
        std::pair<dist,bool> f = pq.top().first;
        pq.pop();
        if(settled[v]) continue;
        if(is_landmark[v] && v != r){
            settled[v] = true;
            reached_vertices.push_back(v);
            if(highway[r].find(v)!=highway[r].end()){
                highway[r][v] = min(highway[r][v], f.first);
            } else{
                highway[r][v] = f.first;
            }
            continue;
        }
        if(DirectedQueryDistanceBound(r,v, dij_distances_pair[v].first)){
            settled[v] = true;
            continue;
        }
        vertex i = 0;
        for(; i < out_landmarks_distances[v].size(); i++){
            if(out_landmarks_distances[v][i].first > r) break;
            auto position = highway[out_landmarks_distances[v][i].first].find(r);
            if(position != highway[out_landmarks_distances[v][i].first].end() && out_landmarks_distances[v][i].second >= dij_distances_pair[v].first + position->second){
                obsolete_positions.push_back(i);
            }
        }
        out_landmarks_distances[v].insert(out_landmarks_distances[v].begin() + i, std::make_pair(r,dij_distances_pair[v].first));
        i ++;
        for(; i < out_landmarks_distances[v].size(); i++){
            auto position = highway[out_landmarks_distances[v][i].first].find(r);
            if(position != highway[out_landmarks_distances[v][i].first].end() && out_landmarks_distances[v][i].second >= dij_distances_pair[v].first + position->second){
                obsolete_positions.push_back(i);
            }
        }
        for (auto it = obsolete_positions.rbegin(); it != obsolete_positions.rend(); ++it) {
                const auto& j = *it;
                out_landmarks_distances[v].erase(out_landmarks_distances[v].begin() + j);
            }
        obsolete_positions.clear();
        settled[v] = true;
        for (auto w: graph.neighborRange(v)) {
                    std::pair<dist,bool> temp = {dij_distances_pair[v].first + graph.weight(v, w), f.second};
                    if (dij_distances_pair[w] > temp) {
                        dij_distances_pair[w] = temp;
                        reached_vertices.push_back(w);
                        pq.push(std::make_pair(dij_distances_pair[w], w));
                    }
                }
    }

    for(const auto & v: reached_vertices){
        settled[v] = false;
        dij_distances_pair[v] = std::make_pair(null_distance,true);
    }


    // REVERSE
    in_landmarks_distances[r].clear();
    in_landmarks_distances[r].emplace_back(r,0);
    // search from r to fill L a la Akiba (prune when reaching landmarks different from r, compute query to vertices not in L)
    dij_distances_pair[r] = std::make_pair(0,true);
    reached_vertices.clear();
    reached_vertices.push_back(r);
    while(!pq.empty()) pq.pop();
    obsolete_positions.clear();

    for(auto w: graph.inNeighborRange(r)){
       if (dij_distances_pair[w].first > dij_distances_pair[r].first + graph.weight(w, r)) {
                        dij_distances_pair[w] = std::make_pair(dij_distances_pair[r].first + static_cast<dist>(graph.weight(w,r)), true);
                        reached_vertices.push_back(w);
                        pq.push(std::make_pair(dij_distances_pair[w], w));
                    }
    }
    settled[r] = true;
    while(!pq.empty()){
        vertex v = pq.top().second;
        std::pair<dist,bool> f = pq.top().first;
        pq.pop();
        if(settled[v]) continue;
        if(is_landmark[v] && v != r){
            settled[v] = true;
            if(highway[v].find(r)!=highway[v].end()){
                highway[v][r] = min(highway[v][r], f.first);
            } else{
                highway[v][r] = f.first;
            }
            reached_vertices.push_back(v);
            continue;
        }
        if(DirectedQueryDistanceBound(v,r, dij_distances_pair[v].first)){
            settled[v] = true;
            continue;
        }
        vertex i = 0;
        for(; i < in_landmarks_distances[v].size(); i++){
            if(in_landmarks_distances[v][i].first > r) break;
            auto position = highway[r].find(in_landmarks_distances[v][i].first);
            if(position != highway[r].end() && in_landmarks_distances[v][i].second >= dij_distances_pair[v].first + position->second){
                obsolete_positions.push_back(i);
            }
        }
        
        in_landmarks_distances[v].insert(in_landmarks_distances[v].begin() + i, std::make_pair(r,dij_distances_pair[v].first));
        i ++;
        for(; i < in_landmarks_distances[v].size(); i++){
            auto position = highway[r].find(in_landmarks_distances[v][i].first);
            if(position != highway[r].end() && in_landmarks_distances[v][i].second >= dij_distances_pair[v].first + position->second){
                obsolete_positions.push_back(i);
            }
        }
        for (auto it = obsolete_positions.rbegin(); it != obsolete_positions.rend(); ++it) {
                const auto& j = *it;
                in_landmarks_distances[v].erase(in_landmarks_distances[v].begin() + j);
            }
        obsolete_positions.clear();
        settled[v] = true;
        for (auto w: graph.inNeighborRange(v)) {
                    std::pair<dist,bool> temp = {dij_distances_pair[v].first + graph.weight(w, v), f.second};
                    if (dij_distances_pair[w] > temp) {
                        dij_distances_pair[w] = temp;
                        reached_vertices.push_back(w);
                        pq.push(std::make_pair(dij_distances_pair[w], w));
                    }
                }
    }

    for(const auto & v: reached_vertices){
        settled[v] = false;
        dij_distances_pair[v] = std::make_pair(null_distance,true);
    }

    reached_vertices.clear();
    L++;
    
    // std::cout << "Current landmarks ";
    // for (const auto& v : landmarks) std::cout << v << " ";
    // std::cout << "\n";

    // for (const auto& [u, inner_map] : highway) {
    //     std::cout << "From vertex " << u << ":\n";

    //     for (const auto& [v, d] : inner_map) {
    //         std::cout << "  -> " << v << " (dist = " << d << ")\n";
    //     }
    // }

    // for (int v = 0; v < V; v++) {

    //     std::cout << "L-In(" << v << "): [";
    //     for (int i = 0; i < in_landmarks_distances[v].size(); i++) {
    //         std::cout << "("
    //                 << in_landmarks_distances[v][i].first << ","
    //                 << in_landmarks_distances[v][i].second << "), ";
    //     }
    //     std::cout << "]" << std::endl;

    //     std::cout << "L-Out(" << v << "): [";
    //     for (int i = 0; i < out_landmarks_distances[v].size(); i++) {
    //         std::cout << "("
    //                 << out_landmarks_distances[v][i].first << ","
    //                 << out_landmarks_distances[v][i].second << "), ";
    //     }
    //     std::cout << "]" << std::endl;
    // }
}

void HighwayLabelling::RemoveLandmarkUnweighted(vertex r) {

    L--;
    landmarks.erase(r);
    is_landmark[r] = false;
    std::vector<std::pair<vertex,dist>> affected_landmarks;

    // clear label of r and treat it as affected
    landmarks_distances[r].clear();
    // search from r for affected landmarks

    std::vector<vertex> reached_vertices;
    dij_distances[r] = 0;
    reached_vertices.push_back(r);
    std::queue<vertex> q;
    q.push(r);
    while (!q.empty()) {
        vertex v = q.front();
        q.pop();
        if (settled[v]) continue;
        if (is_landmark[v]) {
            settled[v] = true;
            if(dij_distances[v] > highway[r][v]) continue;
            affected_landmarks.emplace_back(v,dij_distances[v]);
            continue;
        }
        for(size_t i = 0; i < landmarks_distances[v].size(); i++){
            if(landmarks_distances[v][i].first == r){
                landmarks_distances[v].erase(landmarks_distances[v].begin()+i);
                break;
            }
            if(landmarks_distances[v][i].first > r){
                break;
            }

        }
        for (auto w: graph.neighborRange(v)) {
            if (dij_distances[w] == null_distance) {
                dij_distances[w] = dij_distances[v] + 1;
                settled[w] = settled[v];
                q.push(w);
                reached_vertices.push_back(w);
            }

        }
    }

    for(const auto & v: reached_vertices){
        dij_distances[v] = null_distance;
        settled[v] = false;
    }
    reached_vertices.clear();

    // remove entries of r in highway
    for(const vertex & l: landmarks){
        highway[l].erase(r);
    }
    highway[r].clear();
    highway.erase(r);
    // search from r rooted in each affected landmark to potentially cover vertices
    for(const auto & ld: affected_landmarks){
        const vertex &l = ld.first;
        const dist &d = ld.second;
        dij_distances[l] = 0;
        dij_distances[r] = d;
        reached_vertices.push_back(r);
        std::queue<vertex> q;
        q.push(r);
        settled[l] = true;
        reached_vertices.push_back(l);

        while (!q.empty()) {
            vertex v = q.front();
            q.pop();
            if (settled[v]) continue;
            if (is_landmark[v]) {
                settled[v] = true;
                continue;
            }
            if(QueryDistanceBound(v,l, dij_distances[v])){
                settled[v] = true;
                continue;
            }
            auto it = landmarks_distances[v].begin();
            bool inserted = false;

            for (; it != landmarks_distances[v].end(); ++it) {
                if (it->first == l) {
                    inserted = true;
                    break;
                }
                if (it->first > l) {
                    landmarks_distances[v].insert(it, std::make_pair(l,dij_distances[v]));
                    inserted = true;
                    break;
                }
            }

            if (!inserted) {
                landmarks_distances[v].emplace_back(l,dij_distances[v]);
            }

            for (auto w: graph.neighborRange(v)) {
                if (dij_distances[w] == null_distance) {
                    dij_distances[w] = dij_distances[v] + 1;
                    q.push(w);
                    reached_vertices.push_back(w);
                }

            }
        }

        for(const auto & v: reached_vertices){
            dij_distances[v] = null_distance;
            settled[v] = false;
        }
        reached_vertices.clear();

    }
}



void HighwayLabelling::RemoveLandmarkDirected(vertex r) {

    L--;
    landmarks.erase(r);
    is_landmark[r] = false;

    std::vector<std::pair<vertex,dist>> in_affected_landmarks;

    // clear label of r and treat it as affected
    in_landmarks_distances[r].clear();
    // search from r for affected landmarks

    std::vector<vertex> in_reached_vertices;
    dij_distances[r] = 0;
    in_reached_vertices.push_back(r);
    std::priority_queue<std::pair<dist, vertex>, std::vector<std::pair<dist, vertex>>,
            PQComparator> pq;
    pq.push(std::make_pair(0, r));
    while (!pq.empty()) {
        vertex v = pq.top().second;
        pq.pop();
        if (settled[v]) continue;
        if (is_landmark[v]) {
            settled[v] = true;
            if(dij_distances[v] > highway[r][v]) continue;
            in_affected_landmarks.emplace_back(v,dij_distances[v]);
            continue;
        }
        for(size_t i = 0; i < out_landmarks_distances[v].size(); i++){
            if(out_landmarks_distances[v][i].first == r){
                out_landmarks_distances[v].erase(out_landmarks_distances[v].begin()+i);
                break;
            }
            if(out_landmarks_distances[v][i].first > r){
                break;
            }

        }
        settled[v] = true;
        for (auto w: graph.neighborRange(v)) {
            if (dij_distances[w] > dij_distances[v] + graph.weight(v, w)) {
                dij_distances[w] = dij_distances[v] + (dist) graph.weight(v, w);
                pq.push(std::make_pair(dij_distances[w], w));
                in_reached_vertices.push_back(w);
            }

        }
    }

    for(const auto & v: in_reached_vertices){
        dij_distances[v] = null_distance;
        settled[v] = false;
    }
    in_reached_vertices.clear();
    std::vector<std::pair<vertex,dist>> out_affected_landmarks;
    std::vector<vertex> out_reached_vertices;

    // clear label of r and treat it as affected
    out_landmarks_distances[r].clear();
    // search from r for affected landmarks

    dij_distances[r] = 0;
    out_reached_vertices.push_back(r);
    while(!pq.empty()) pq.pop();
    pq.push(std::make_pair(0, r));
    while (!pq.empty()) {
        vertex v = pq.top().second;
        pq.pop();
        if (settled[v]) continue;
        if (is_landmark[v]) {
            settled[v] = true;
            if(dij_distances[v] > highway[v][r]) continue;
            out_affected_landmarks.emplace_back(v,dij_distances[v]);
            continue;
        }
        for(size_t i = 0; i < in_landmarks_distances[v].size(); i++){
            if(in_landmarks_distances[v][i].first == r){
                in_landmarks_distances[v].erase(in_landmarks_distances[v].begin()+i);
                break;
            }
            if(in_landmarks_distances[v][i].first > r){
                break;
            }
        }
        settled[v] = true;
        for (auto w: graph.inNeighborRange(v)) {
            if (dij_distances[w] > dij_distances[v] + graph.weight(w, v)) {
                dij_distances[w] = dij_distances[v] + (dist) graph.weight(w, v);
                pq.push(std::make_pair(dij_distances[w], w));
                out_reached_vertices.push_back(w);
            }

        }
    }

    for(const auto & v: out_reached_vertices){
        dij_distances[v] = null_distance;
        settled[v] = false;
    }
    out_reached_vertices.clear();
    // remove entries of r in highway
    for(const vertex & l: landmarks){
        highway[l].erase(r);
    }
    highway[r].clear();
    highway.erase(r);
    for(const auto & ld: in_affected_landmarks){
        const vertex &l = ld.first;
        const dist &d = ld.second;
        if(DirectedQueryDistanceBound(r,l, d)) continue;
        auto insertion_index = std::upper_bound(in_landmarks_distances[r].begin(),
                                                in_landmarks_distances[r].end(), std::make_pair(l,d));
        in_landmarks_distances[r].insert(insertion_index, std::make_pair(l,d));
    }
    for(const auto & ld: out_affected_landmarks){
        const vertex &l = ld.first;
        const dist &d = ld.second;
        if(DirectedQueryDistanceBound(l,r, d)) continue;
        auto insertion_index = std::upper_bound(out_landmarks_distances[r].begin(),
                                                    out_landmarks_distances[r].end(), std::make_pair(l,d));
        out_landmarks_distances[r].insert(insertion_index, std::make_pair(l,d));
    }

    // search from r rooted in each affected landmark to potentially cover vertices
    for(const auto & ld: in_affected_landmarks){
        const vertex &l = ld.first;
        const dist &d = ld.second;
        dij_distances_pair[l] = std::make_pair(0,true);
        dij_distances_pair[r] = std::make_pair(d,true);
        in_reached_vertices.push_back(r);
        std::priority_queue<std::pair<std::pair<dist, bool>, vertex>, std::vector<std::pair<std::pair<dist, bool>, vertex>>,
                PQFlagComparator> pq;
        for (auto w: graph.inNeighborRange(r)) {
            if (w == l) continue;
            dij_distances_pair[w] = std::make_pair(dij_distances_pair[r].first + (dist) graph.weight(w, r), true);
            in_reached_vertices.push_back(w);
            pq.push(std::make_pair(dij_distances_pair[w], w));

        }
        settled[l] = true;
        settled[r] = true;
        in_reached_vertices.push_back(l);

        while (!pq.empty()) {
            vertex v = pq.top().second;
            std::pair<dist,bool> f = pq.top().first;
            pq.pop();
            if (settled[v]) continue;
            if(is_landmark[v] && v != l){
                settled[v] = true;
                if(highway[v].find(l)!=highway[v].end()){
                    highway[v][l] = min(highway[v][l], f.first);
                } else{
                    highway[v][l] = f.first;
                }
                in_reached_vertices.push_back(v);
                continue;
            }
            if(DirectedQueryDistanceBound(v,l, dij_distances_pair[v].first)){
                settled[v] = true;
                continue;
            }
            if(in_landmarks_distances[v].rbegin()->first < l || in_landmarks_distances[v].empty()) {
                in_landmarks_distances[v].emplace_back(l,dij_distances_pair[v].first);
            }
            else {
                for(size_t i = 0; i < in_landmarks_distances[v].size(); i++) {
                    if(in_landmarks_distances[v][i].first < l) {
                        continue;
                    }
                    if(in_landmarks_distances[v][i].first > l) {
                        in_landmarks_distances[v].insert(in_landmarks_distances[v].begin() + i, std::make_pair(l,dij_distances_pair[v].first));
                        break;
                    }
                    if(in_landmarks_distances[v][i].first == l) {
                        in_landmarks_distances[v][i].second = dij_distances_pair[v].first;
                        break;
                    }
                }
            }
            settled[v] = true;

            for (auto w: graph.inNeighborRange(v)) {
                    std::pair<dist,bool> temp = {dij_distances_pair[v].first + graph.weight(w, v), f.second};
                    if (dij_distances_pair[w] > temp) {
                        dij_distances_pair[w] = temp;
                        in_reached_vertices.push_back(w);
                        pq.push(std::make_pair(dij_distances_pair[w], w));
                    }
                }
        }

        for(const auto & v: in_reached_vertices){
            dij_distances_pair[v] = std::make_pair(null_distance,true);
            settled[v] = false;
        }
        in_reached_vertices.clear();

    }

    // search from r rooted in each affected landmark to potentially cover vertices
    for(const auto & ld: out_affected_landmarks){
        const vertex &l = ld.first;
        const dist &d = ld.second;
        dij_distances_pair[l] = std::make_pair(0,true);
        dij_distances_pair[r] = std::make_pair(d,true);
        out_reached_vertices.push_back(r);
        std::priority_queue<std::pair<std::pair<dist, bool>, vertex>, std::vector<std::pair<std::pair<dist, bool>, vertex>>,
                PQFlagComparator> pq;
        for (auto w: graph.neighborRange(r)) {
            if (w == l) continue;
            dij_distances_pair[w] = std::make_pair(dij_distances_pair[r].first + (dist) graph.weight(r,w), true);
            out_reached_vertices.push_back(w);
            pq.push(std::make_pair(dij_distances_pair[w], w));
        }
        settled[l] = true;
        settled[r] = true;
        out_reached_vertices.push_back(l);

        while (!pq.empty()) {
            vertex v = pq.top().second;
            std::pair<dist,bool> f = pq.top().first;
            pq.pop();
            if(settled[v]) continue;
            if(is_landmark[v] && v != l){
                settled[v] = true;
                out_reached_vertices.push_back(v);
                if(highway[l].find(v)!=highway[l].end()){
                    highway[l][v] = min(highway[l][v], f.first);
                } else{
                    highway[l][v] = f.first;
                }
                continue;
            }
            if(DirectedQueryDistanceBound(l,v, dij_distances_pair[v].first)){
                settled[v] = true;
                continue;
            }

            if(out_landmarks_distances[v].empty() || out_landmarks_distances[v].rbegin()->first < l) {
                    out_landmarks_distances[v].emplace_back(l,dij_distances_pair[v].first);
            }
            else{
                for(size_t i = 0; i < out_landmarks_distances[v].size(); i++) {
                                if(out_landmarks_distances[v][i].first < l) {
                                    continue;
                                }
                                if(out_landmarks_distances[v][i].first > l) {
                                    out_landmarks_distances[v].insert(out_landmarks_distances[v].begin() + i, std::make_pair(l,dij_distances_pair[v].first));
                                    break;
                                }
                                if(out_landmarks_distances[v][i].first == l) {
                                    out_landmarks_distances[v][i].second = dij_distances_pair[v].first;
                                    break;
                                }
                            }
            }
            
            settled[v] = true;
            for (auto w: graph.neighborRange(v)) {
                    std::pair<dist,bool> temp = {dij_distances_pair[v].first + graph.weight(v, w), f.second};
                    if (dij_distances_pair[w] > temp) {
                        dij_distances_pair[w] = temp;
                        out_reached_vertices.push_back(w);
                        pq.push(std::make_pair(dij_distances_pair[w], w));
                    }
                }
        }

        for(const auto & v: out_reached_vertices){
            dij_distances_pair[v] = std::make_pair(null_distance, true);
            settled[v] = false;
        }
        out_reached_vertices.clear();

    }
    // std::cout << "Current landmarks ";
    // for(const auto& v: landmarks) std::cout << v << " ";
    // std::cout << "\n";
    // for (const auto& [u, inner_map] : highway) {
    //     std::cout << "From vertex " << u << ":\n";

    //     for (const auto& [v, d] : inner_map) {
    //         std::cout << "  -> " << v << " (dist = " << d << ")\n";
    //     }
    // }
    
    //  for (int v = 0; v < V; v++) {
    //     std::cout << "L-In(" << v << "): [";
    //     for (int i = 0; i < in_landmarks_distances[v].size(); i++) {
    //         std::cout << "(" << in_landmarks_distances[v][i] << "," << in_distances[v][i] << "), ";
    //     }
    //     std::cout << "]" << std::endl;
    //     std::cout << "L-Out(" << v << "): [";
    //     for (int i = 0; i < out_landmarks_distances[v].size(); i++) {
    //         std::cout << "(" << out_landmarks_distances[v][i] << "," << out_distances[v][i] << "), ";
    //     }
    //     std::cout << "]" << std::endl;
    // }
}

// inline void HighwayLabelling::StoreIndex(std::string filename) {
//         std::ofstream ofs(std::string("index/")+std::string(filename) + std::string("index"));
//     for (int i = 0; i < V; i++) {
//         vertex C = 0;
//         for (int j = 0; j < distances[i].size(); j++) {
//             if(distances[i][j] != null_distance)
//                 C++;
//         }
//         ofs.write((char*)&C, sizeof(C));
//         for (int j = 0; j < distances[i].size(); j++) {
//             if(distances[i][j] != null_distance) {
//                 ofs.write((char*)&j, sizeof(j));
//                 ofs.write((char*)&distances[i][j], sizeof(distances[i][j]));
//             }
//         }
//     }

//     for (const vertex & v: landmarks) {
//         for (const vertex & w: landmarks) {
//             if(highway[v][w] != null_distance)
//                 ofs.write((char*)&highway[v][w], sizeof(highway[v][w]));
//         }
//     }
//     ofs.close();
// }

#endif  // HGHWAY_LABELING_H_
