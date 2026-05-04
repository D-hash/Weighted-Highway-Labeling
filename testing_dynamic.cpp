#include <iostream>
#include "highway_cover_labelling.h"
#include <vector>
#include <boost/program_options.hpp>
#include "networkit/components/ConnectedComponents.hpp"
#include <networkit/io/EdgeListReader.hpp>
#include <networkit/graph/GraphTools.hpp>
#include "mytimer.h"
#include <networkit/distance/Diameter.hpp>
#include  <random>
#include  <iterator>
#include <routingkit/contraction_hierarchy.h>
#include <fstream>
#include <stdexcept>
#include <cassert>
#include <utility>
#include <networkit/graph/Graph.hpp>
#include <sstream>
#include <string>
#include <algorithm> // std::remove_if

using namespace std;
double median(std::vector<double>& arr) { //SORTS
    size_t n = arr.size() / 2;
    if (n % 2 == 0) {
        std::nth_element(arr.begin(),arr.begin() + n/2,arr.end());
        std::nth_element(arr.begin(),arr.begin() + (n - 1) / 2,arr.end());
        return (double) (arr[(n-1)/2]+ arr[n/2])/2.0;
    }

    else{
        std::nth_element(arr.begin(),arr.begin() + n / 2,arr.end());
        return (double) arr[n/2];
    }
    assert(false);
}


namespace {
inline std::string trim(std::string s) {
    auto notSpace = [](unsigned char c) { return !std::isspace(c); };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), notSpace));
    s.erase(std::find_if(s.rbegin(), s.rend(), notSpace).base(), s.end());
    return s;
}

// Split by comma, keep empty fields if any (we'll validate)
inline std::vector<std::string> splitCSVLine(const std::string& line) {
    std::vector<std::string> out;
    std::string field;
    std::istringstream iss(line);
    while (std::getline(iss, field, ',')) {
        out.push_back(trim(field));
    }
    return out;
}
} // namespace

void read_csv(const std::string& source, NetworKit::Graph& g,
              bool directed = false) {
    /* --------------------------------------------------
     * First pass: determine
     *   - number of data lines
     *   - max node id
     *   - weighted vs unweighted (2 vs 3 columns)
     * -------------------------------------------------- */
    std::ifstream countStream(source);
    if (!countStream)
        throw std::runtime_error("Error opening file: " + source);

    std::size_t totalLines = 0;
    NetworKit::node maxNode = 0;
    bool seenUnweighted = false;
    bool seenWeighted = false;

    std::string line;
    while (std::getline(countStream, line)) {
        line = trim(line);
        if (line.empty()) continue;
        if (!line.empty() && line[0] == '#') continue; // comment/header

        auto cols = splitCSVLine(line);
        if (!(cols.size() == 2 || cols.size() == 3)) {
            throw std::runtime_error("Invalid CSV line (expected 2 or 3 columns): " + line);
        }

        // Parse u,v
        NetworKit::node u, v;
        try {
            u = static_cast<NetworKit::node>(std::stoull(cols[0]));
            v = static_cast<NetworKit::node>(std::stoull(cols[1]));
        } catch (...) {
            throw std::runtime_error("Invalid node id in CSV line: " + line);
        }

        if (u > maxNode) maxNode = u;
        if (v > maxNode) maxNode = v;

        if (cols.size() == 3) seenWeighted = true;
        else seenUnweighted = true;

        ++totalLines;
    }

    if (totalLines == 0)
        throw std::runtime_error("Empty CSV graph file (no data lines): " + source);

    if (seenWeighted && seenUnweighted) {
        throw std::runtime_error("Mixed weighted/unweighted rows in CSV file: " + source);
    }

    const bool isWeighted = seenWeighted;
    const NetworKit::count n = static_cast<NetworKit::count>(maxNode) + 1;

    std::cout << "[CSV] Total data lines in file: " << totalLines << std::endl;
    std::cout << "[CSV] Detected " << (isWeighted ? "weighted" : "unweighted")
              << " graph" << std::endl;

    countStream.close();

    /* --------------------------------------------------
     * Second pass: build graph
     * -------------------------------------------------- */
    std::ifstream ifs(source);
    if (!ifs)
        throw std::runtime_error("Error reopening file: " + source);

    NetworKit::Graph tmp(n, isWeighted, directed);

    ProgressStream reader(totalLines);
    reader.label() << "Reading "
                   << (isWeighted ? "weighted " : "unweighted ")
                   << (directed ? "directed " : "undirected ")
                   << "graph in " << source << " (CSV format)";

    while (std::getline(ifs, line)) {
        line = trim(line);
        if (line.empty()) continue;
        if (!line.empty() && line[0] == '#') continue;

        auto cols = splitCSVLine(line);
        if (!(cols.size() == 2 || cols.size() == 3)) continue; // or throw

        ++reader;

        NetworKit::node u, v;
        NetworKit::edgeweight w = 1.0;

        try {
            u = static_cast<NetworKit::node>(std::stoull(cols[0]));
            v = static_cast<NetworKit::node>(std::stoull(cols[1]));
            if (isWeighted && cols.size() == 3) {
                w = static_cast<NetworKit::edgeweight>(std::stod(cols[2]));
            }
        } catch (...) {
            // skip malformed line (or throw if you prefer strictness)
            continue;
        }

        if (u == v) continue;

        // First edge wins for duplicates
        tmp.addEdge(u, v, w, /*checkMultiEdge=*/true);
    }

    tmp.removeSelfLoops();
    tmp.removeMultiEdges();

    g = std::move(tmp);
}


void read_nde(const std::string& source, NetworKit::Graph& g,
              bool directed = false) {
    /* --------------------------------------------------
     * First pass: determine
     *   - number of lines
     *   - max node id
     *   - weighted vs unweighted
     * -------------------------------------------------- */
    std::ifstream countStream(source);
    if (!countStream)
        throw std::runtime_error("Error opening file: " + source);

    std::size_t totalLines = 0;
    NetworKit::node maxNode = 0;
    bool isWeighted = false;
    bool seenUnweighted = false;
    bool seenWeighted = false;

    std::string line;
    while (std::getline(countStream, line)) {
        if (line.empty()) continue;

        std::istringstream iss(line);
        NetworKit::node u, v;
        NetworKit::edgeweight w;

        if (!(iss >> u >> v))
            throw std::runtime_error("Invalid NDE line: " + line);

        if (u > maxNode) maxNode = u;
        if (v > maxNode) maxNode = v;

        if (iss >> w) {
            seenWeighted = true;
        } else {
            seenUnweighted = true;
        }

        ++totalLines;
    }

    if (totalLines == 0)
        throw std::runtime_error("Empty NDE file: " + source);

    if (seenWeighted && seenUnweighted) {
        throw std::runtime_error(
            "Mixed weighted/unweighted edges in NDE file: " + source);
    }

    isWeighted = seenWeighted;

    std::cout << "[NDE] Total lines in file: " << totalLines << std::endl;
    std::cout << "[NDE] Detected "
              << (isWeighted ? "weighted" : "unweighted")
              << " graph" << std::endl;

    const NetworKit::count n =
        static_cast<NetworKit::count>(maxNode) + 1;

    countStream.close();

    /* --------------------------------------------------
     * Second pass: build graph
     * -------------------------------------------------- */
    std::ifstream ifs(source);
    if (!ifs)
        throw std::runtime_error("Error reopening file: " + source);

    NetworKit::Graph tmp(n, isWeighted, directed);

    ProgressStream reader(totalLines);
    reader.label() << "Reading "
                   << (isWeighted ? "weighted " : "unweighted ")
                   << (directed ? "directed " : "undirected ")
                   << "graph in " << source << " (NDE format)";

    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        ++reader;

        std::istringstream iss(line);
        NetworKit::node u, v;
        NetworKit::edgeweight w = 1.0;

        iss >> u >> v;
        if (isWeighted)
            iss >> w;

        if (u == v) continue;

        // First edge wins
        tmp.addEdge(u, v, w, /*checkMultiEdge=*/true);
    }

    tmp.removeSelfLoops();
    tmp.removeMultiEdges();

    g = std::move(tmp);
}


void read_hist(const std::string& source, NetworKit::Graph& g) {
    /* --------------------------------------------------
     * First pass: count lines
     * -------------------------------------------------- */
    std::ifstream countStream(source);
    if (!countStream)
        throw std::runtime_error("Error opening file: " + source);

    std::size_t totalLines = 0;
    std::string line;
    while (std::getline(countStream, line))
        ++totalLines;

    if (totalLines == 0)
        throw std::runtime_error("Empty HIST file: " + source);

    const std::size_t dataLines = totalLines - 1; // minus header
    countStream.close();
    // Print number of lines read
    std::cout << "[HIST] Total lines in file: " << totalLines << std::endl;
    /* --------------------------------------------------
     * Second pass: actual reading
     * -------------------------------------------------- */
    std::ifstream ifs(source);
    if (!ifs)
        throw std::runtime_error("Error reopening file: " + source);

    NetworKit::count n = 0, m = 0;
    int weightedFlag = 0, directedFlag = 0;

    if (!(ifs >> n >> m >> weightedFlag >> directedFlag))
        throw std::runtime_error("Invalid HIST header in: " + source);

    assert((weightedFlag == 0 || weightedFlag == 1) &&
           (directedFlag == 0 || directedFlag == 1));

    const bool weighted = (weightedFlag != 0);
    const bool directed = (directedFlag != 0);

    NetworKit::Graph tmp(n, weighted, directed);

    ProgressStream reader(dataLines);
    reader.label() << "Reading "
                   << (weighted ? "weighted " : "unweighted ")
                   << (directed ? "directed " : "undirected ")
                   << "graph in " << source << " (HIST format)";

    long long time;
    NetworKit::node u, v;
    NetworKit::edgeweight w;

    while (ifs >> time >> u >> v >> w) {
        ++reader;

        if (u >= n || v >= n) continue;
        if (u == v) continue;

        // First weight wins
        tmp.addEdge(u, v, weighted ? w : 1.0, /*checkMultiEdge=*/true);
    }

    tmp.removeSelfLoops();
    tmp.removeMultiEdges();

    g = std::move(tmp);
}


double average(std::vector<double> & arr) {

    auto const count = static_cast<double>(arr.size());
    double sum = 0;
    for(double value: arr) sum += value;
    return sum / count;
}

int main(int argc, char **argv) {
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");

    desc.add_options()
            ("graph_location,g", po::value<std::string>(), "Input Graph File Location")
            ("landmarks,l", po::value<int>(), "Number of Landmarks")
            ("num_queries,q", po::value<int>(), "Number of Queries to Be Performed")
            ("ordering,o",po::value<int>(), "Type of Node Ordering [DEGREE(0) RANDOM(1)]")
            ("index_file,i",po::value<string>(), "Name of the file that will store the index")
            ("graph_changes,c",po::value<int>(), "Number of landmarks' changes")
            ("vertex_id,v",po::value<int>(), "Starting ID of the vertices; either 0 or 1")
            ("dynamic_exp,d",po::value<int>(), "Type of dynamic experiment [INC(1) DEC(2) FUL(3)]")
            ("contraction_hierarchy,k",po::value<int>(), "Use contraction hierarchy (0 = NO | 1 = YES)")
            ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    if(vm.empty()){
        std::cout << desc << "\n";
        throw std::runtime_error("empty argument array");
    }
    std::string graph_location;

    if(vm.count("graph_location")){
        graph_location = vm["graph_location"].as<std::string>();
    }


    int L = -1;

    if (vm.count("landmarks")){
        L = vm["landmarks"].as<int>();
    }

    if (L < 1){
        std::cout << desc << "\n";
        throw std::runtime_error("L must be at least 1");
    }

    int ordering = -1;

    if (vm.count("ordering")){
        ordering = vm["ordering"].as<int>();
    }
    if(ordering != 0 && ordering != 1 && ordering != 2 && ordering != 3){
        std::cout << desc << "\n";
        throw std::runtime_error("Wrong ordering selection (0 or 1 or 2)");
    }

    int num_queries = -1;
    if (vm.count("num_queries")){
        num_queries = vm["num_queries"].as<int>();
    }
    if(num_queries < 2){
        std::cout << desc << "\n";
        throw std::runtime_error("Wrong num queries");
    }

   std::string index_file;
    if(vm.count("index_file")){
        index_file = vm["index_file"].as<std::string>();
    }
    int graph_changes = 0;
    if(vm.count("graph_changes")){
        graph_changes = vm["graph_changes"].as<int>();
    }

    int starting_vertex_id = 0;
    if(vm.count("vertex_id")){
        starting_vertex_id = vm["vertex_id"].as<int>();
    }

    int dyn_type = 0;
    if(vm.count("dynamic_exp")){
        dyn_type = vm["dynamic_exp"].as<int>();
    }

    int contraction_hierarchy = 0;
    if(vm.count("contraction_hierarchy")){
        contraction_hierarchy = vm["contraction_hierarchy"].as<int>();
    }
    std::cout << "Reading " << graph_location << " building with L = " << L << " num_queries = "<<num_queries
     << " graph_changes " << graph_changes << " ordering "<<ordering<<" index_file " << index_file
     << " dynamic experiment type " << dyn_type << " CH " << contraction_hierarchy ? "YES" : "NO";
    std::cout<<std::endl;

    NetworKit::Graph graph;
    
    if(graph_location.find(".hist") != std::string::npos){
        read_hist(graph_location, graph);
    } else {
        NetworKit::EdgeListReader edl(' ', starting_vertex_id, "#", true);
        graph = edl.read(graph_location);
    }
    // graph.removeMultiEdges();
    // graph.removeSelfLoops();
    // if(!graph.isDirected()) {
    //     const NetworKit::Graph &graph_handle = graph;
    //     NetworKit::ConnectedComponents *cc = new NetworKit::ConnectedComponents(graph_handle);
    //     graph = cc->extractLargestConnectedComponent(graph_handle, true);
    //     graph.shrinkToFit();
    //     graph.indexEdges();
    // }
    graph.shrinkToFit();
    //graph.indexEdges();
    std::cout << "Graph has " << graph.numberOfNodes() << " vertices and " << graph.numberOfEdges()
              << " edges \n";

    vertex diameter = 0;
    HighwayLabelling *hl = new HighwayLabelling(graph, L, ordering, graph_changes, dyn_type);

    mytimer timer;

    // construct labelling
    double hl_constr_time = 0.0;
    timer.restart();
    if(graph.isDirected()){
        hl->ConstructDirWeighHL();
    }else
        hl->FarhanConstruction();
    hl_constr_time = timer.elapsed();
    long hl_size = graph.isDirected() ? hl->DirectedLabellingSize() : hl->LabellingSize();
    cout << "HL constr time " << hl_constr_time << " | HL size " << hl_size << "\n";

    int modifications = graph_changes;
    std::vector<double> incremental_time;
    ProgressStream updates(modifications);

    updates.label() << "Updates";
    if(dyn_type == 1) {
        std::vector<vertex> incland;
        hl->GetIncrementalLandmarks(incland);
        for(const vertex& newlandmark: incland) {
            if(hl->landmarks.find(newlandmark) != hl->landmarks.end())
                throw new std::runtime_error("old landmark");
            timer.restart();
            if(graph.isDirected())
                hl->AddLandmarkDirected(newlandmark);
            else hl->AddLandmarkUnweighted(newlandmark);
            incremental_time.push_back(timer.elapsed());
            ++updates;
        }
    }
    else if(dyn_type == 2) {
        while (graph_changes--) {
            int removelandmark = *select_randomly(hl->landmarks.begin(), hl->landmarks.end());

            timer.restart();
            if(graph.isDirected())
                hl->RemoveLandmarkDirected(removelandmark);
            else hl->RemoveLandmarkUnweighted(removelandmark);
            incremental_time.push_back(timer.elapsed());
            ++updates;
        }
    }

    else if(dyn_type == 3){
        int ins_changes = graph_changes / 2;
        int del_changes = graph_changes / 2;
        std::vector<vertex> incland;
        hl->GetIncrementalLandmarks(incland);
        while(graph_changes--) {
            if ((rand() % 2 == 0 || del_changes == 0) && ins_changes) {
                if(hl->landmarks.find(incland[incland.size() - ins_changes]) != hl->landmarks.end()) {
                    std::cout << "Repeated landmark\n";
                    throw new std::runtime_error("old landmark");   
                }
                timer.restart();
                if(graph.isDirected())
                    hl->AddLandmarkDirected(incland[incland.size() - ins_changes]);
                else hl->AddLandmarkUnweighted(incland[incland.size() - ins_changes]);
                ins_changes--;
                incremental_time.push_back(timer.elapsed());
                ++updates;

            } else {
                int removelandmark = *select_randomly(hl->landmarks.begin(), hl->landmarks.end());
                timer.restart();
                if(graph.isDirected())
                    hl->RemoveLandmarkDirected(removelandmark);
                else hl->RemoveLandmarkUnweighted(removelandmark);
                del_changes--;
                incremental_time.push_back(timer.elapsed());
                ++updates;
            }
        }
    }
    std::cout << "Executed total of " << modifications << " landmark modifications | Mean time " << average(incremental_time) << "\n";
    long dhl_size = graph.isDirected() ? hl->DirectedLabellingSize() : hl->LabellingSize();
    std::cout << "Current DHL size " << dhl_size << "\n";
    std::cout << "Current number of landmarks " << hl->L << "\n";

    //hl->StoreIndex(graph_location+"_"+std::to_string(L)+"_"+std::to_string(dyn_type)+"_"+std::to_string(modifications)+"dynamic");
    vector<double> hl_query_time;
    vector<double> naive_query_time;
    vector<double> sssp_query_time;
    vector<pair<int,int>> queried_nodes;
    vector<double> scratch_query_time;
    vector<double> ch_query_time(num_queries);
    vector<dist> queried_distances;
    long q = num_queries;
    ProgressStream dyn_query_bar(q);
    dyn_query_bar.label() << "Query computation";
    while(num_queries--){
        int s = NetworKit::GraphTools::randomNode(graph);
        int t = NetworKit::GraphTools::randomNode(graph);
        queried_nodes.emplace_back(s,t);
        double hl_q_time = 0.0;
        timer.restart();
        dist  hld= graph.isDirected() ? hl->DirectedQueryDistance(s,t) : hl->QueryDistance(s,t);
        hl_q_time += timer.elapsed();
        queried_distances.push_back(hld);
        hl_query_time.push_back(hl_q_time);
        dist ssspd;
        double sssp_time = 0.0;
        timer.restart();
        if(q - num_queries < 10) {
            if (graph.isDirected())
                ssspd = hl->DirectedDijkstra(s, t);
            else ssspd = hl->BFSQuery(s, t);
            sssp_time += timer.elapsed();
            sssp_query_time.push_back(sssp_time);
            if (ssspd != hld) {
                cout << "Error between " << s << " and " << t << "\n";
                cout << "HL distance " << (int) hld << " | SSSP distance " << (int) ssspd << "\n";
                throw new std::runtime_error("experiment fails");
            }
        }
        ++dyn_query_bar;
    }
    cout << "HL median time " << median(hl_query_time) << " | mean time " << average(hl_query_time) << "\n";
    cout << "SSSP median time " << median(sssp_query_time) << " | mean time " << average(sssp_query_time) << "\n";
    vector<vertex> new_landmarks;
    hl->GetLandmarks(new_landmarks);
    delete hl;

    HighwayLabelling *scratch_hl = new HighwayLabelling(graph, L, ordering, 0, dyn_type);

    scratch_hl->SetLandmarks(new_landmarks);

    double scratch_hl_constr_time = 0.0;
    timer.restart();
    if(graph.isDirected()){
        scratch_hl->ConstructDirWeighHL();
    }
    else 
        scratch_hl->FarhanConstruction();
    scratch_hl_constr_time = timer.elapsed();
    long scratch_hl_size = graph.isDirected() ? scratch_hl->DirectedLabellingSize() : scratch_hl->LabellingSize();
    cout << "Scratch HL constr time " << scratch_hl_constr_time << " | HL size " << scratch_hl_size << "\n";
    //scratch_hl->StoreIndex(graph_location+"_"+std::to_string(L)+"_"+std::to_string(dyn_type)+"_"+std::to_string(modifications)+"_scratch");
    ProgressStream scratch_query_bar(q);
    scratch_query_bar.label() << "Scratch query computation";
    int i = 0;
    for(const auto & p: queried_nodes){
        double scratch_hl_q_time = 0.0;
        timer.restart();
        dist  shld=  graph.isDirected() ? scratch_hl->DirectedQueryDistance(p.first,p.second) : scratch_hl->QueryDistance(p.first,p.second);
        scratch_hl_q_time += timer.elapsed();
        scratch_query_time.push_back(scratch_hl_q_time);
        if(queried_distances[i] != shld) {
            cout << "Error between " << p.first << " and " << p.second << "\n";
            cout << "Updated HL distance " << (int) queried_distances[i] << " | Scratch HL distance " << (int) shld << "\n";
            throw new std::runtime_error("experiment fails");
        }
        ++i;
        ++scratch_query_bar;
    }
    delete scratch_hl;


    timer.restart();
    double ch_cnst_time = 0.0;

    if(contraction_hierarchy) {
        unsigned node_count = graph.numberOfNodes();
        std::vector<unsigned> tail;
        std::vector<unsigned> head;
        std::vector<unsigned> weight;
        graph.forNodes([&](vertex from) {
            graph.forNeighborsOf(from, [&] (vertex to) {
                tail.push_back(from);
                head.push_back(to);
                weight.push_back(graph.weight(from, to));
                tail.push_back(to);
                head.push_back(to);
                weight.push_back(graph.weight(from, to));
            });
        });
        RoutingKit::ContractionHierarchy ch = RoutingKit::ContractionHierarchy::build(node_count, tail, head, weight);

        ch_cnst_time = timer.elapsed();
        ch.save_file("index/"+graph_location+"_"+std::to_string(L)+"_ch");
        std::cout << "CH in time " << ch_cnst_time << "\n";

        tail.clear();
        head.clear();
        weight.clear();

        ProgressStream ch_query_bar(q);
        ch_query_bar.label() << "CH query computation";
        i = 0;
       RoutingKit::ContractionHierarchyQuery chquery(ch);
        for(const auto & p: queried_nodes){
            double ch_hl_q_time = 0.0;
            dist ch_d = null_distance;
            dist sl = null_distance;
            dist lt = null_distance;
            timer.restart();

            for(vertex l: new_landmarks) {
                chquery.reset().add_source(p.first).add_target(l).run();
                sl = chquery.get_distance();
                chquery.reset().add_source(l).add_target(p.second).run();
                lt = chquery.get_distance();

                if(ch_d > sl+lt) {
                    ch_d = sl+lt;
                }
            }
            //scratch_hl->DirectedDijkstra(p.first, p.second);
            ch_hl_q_time += timer.elapsed();
            ch_query_time.push_back(ch_hl_q_time);
            if(queried_distances[i] != ch_d) {
                cout << "Error between " << p.first << " and " << p.second << "\n";
                cout << "Updated HL distance " << (int) queried_distances[i] << " | Scratch HL distance " << (int) ch_d << "\n";
                throw new std::runtime_error("experiment fails");
            }
            ++i;
            ++ch_query_bar;
        }
    }
    std::string order_string;

    switch(ordering){
        case (0):
            order_string = "DEG";
            break;
        case (1):
            order_string = "RND";
            break;
        case (2):
            order_string = "BET";
            break;
        case (3):
            order_string = "DTS";
            break;
        default:
            throw new std::runtime_error("problem on order string");


    }

    std::time_t rawtime;
    std::tm* timeinfo;
    char buffer [80];

    std::time(&rawtime);
    timeinfo = std::localtime(&rawtime);
    std::strftime(buffer,80,"%d-%m-%Y-%H-%M-%S",timeinfo);
    std::string tmp_time = buffer;
    std::string shortenedName = graph_location.substr(0, 16);
    stringstream string_object_name;
    string_object_name<<L;
    std::string l_string;
    string_object_name>>l_string;
    std::string prefix = graph.isDirected() ? "dir_" : "undir_";
    std::string timestampstring = prefix +"dynamic_landmarks_with_ch_"+std::to_string(dyn_type)+"_"+shortenedName+"_"+l_string+"_"+"_"+order_string+"_"+"_"+tmp_time;

    std::string logFile = timestampstring +".csv";


    std::ofstream ofs(logFile);


    ofs << "G,V,E,L,ORD,Q,DIAM,C,HLCT,HLIS,DHLUT,DHLIS,DHLQT,SHLCT,SHLIS,SHLQT,CHCT,CHQT\n";
    ofs << graph_location << "," << graph.numberOfNodes() << "," << graph.numberOfEdges() << "," << L << "," <<
        order_string << "," << q << "," << diameter << "," << modifications  << "," << hl_constr_time << "," << hl_size << ","
        << average(incremental_time) << "," << dhl_size << "," << average(hl_query_time) << "," << scratch_hl_constr_time
        << "," << scratch_hl_size << "," << average(scratch_query_time) << "," << ch_cnst_time << "," << average(ch_query_time) <<
        "\n";
    ofs.close();
    exit(EXIT_FAILURE);
}
