#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <limits>
#include <queue>
#include <math.h>

using namespace std;

#define CLUSTER_RADIUS 100  // in km
#define EARTH_RADIUS 6371  
#define EARTH_MASS 5.972e24
#define G 6.674e-11

enum class NodeType {
    GROUND,
    LEO,
    HEO
};

struct Position{
    double x;
    double y;
    double z;
};

struct Normal{
    double a;
    double b;
    double c;
};

struct Node {
    int id;
    int cpuCapacity;
    int memCapacity;
    int cpuAvailable;
    int memAvailable;
    double velocity; // 0 for ground nodes
    double communicationRange; // this value is for satellites only
    int height; // this value is for satellites only
    NodeType type;
    Position position;
    Normal normalVector;
};

struct Cluster{
  vector<Node *> nodes;
  Node *sourceNode;
  Position center;   
};

struct Link {
    Node *node1;
    Node *node2;
    int bandwidthCapacity;
    int bandwidthAvailable;
    double delay;
    bool active;
    Link * revLink;
};

struct VNF{
    int id;
    int cpuRequirement;
    int memRequirement;
    int delay;
};

struct SET{
    int id;
    vector<VNF *> vnfs;
    vector<Node *> deployedNodes;
};

struct SFCR{ 
    int id;
    Node *sourceNode;
    vector<SET *> sets;
    int bandwidthRequirement; 
    int maxDelay;
    int lifetime;
};

struct Path{
    vector<Node *> nodes;
    int cost;
};

class Graph{
    private:
        unordered_map<int, vector<pair<int, double>>> graph;
    public:
        void add_edge(int u, int v, double w){
            graph[u].push_back({v, w});
        }

        vector<pair<double, vector<int>>> dijkstra(int src, int length, int k = 3) {
            typedef pair<double, pair<int, vector<int>>> pq_entry;
            priority_queue<pq_entry, vector<pq_entry>, greater<pq_entry>> pq;

            vector<pair<double, vector<int>>> shortest_paths;
            pq.push({0.0, {src, {src}}});

            unordered_map<int, bool> visited_nodes;

            while (!pq.empty() && shortest_paths.size() < k) {
                auto pq_entry = pq.top();
                double cost = pq_entry.first;
                int node = pq_entry.second.first;
                vector<int> path = pq_entry.second.second;
                pq.pop();

                if (visited_nodes[node]) continue; // Skip if node is already visited
                visited_nodes[node] = true;

                if (path.size() == length + 1) {
                    shortest_paths.push_back({cost, path});
                } else if (path.size() < length + 1) {
                    for (const auto& neighbor : graph[node]) {
                        int neighbor_node = neighbor.first;
                        double weight = neighbor.second;
                        vector<int> new_path = path;
                        if (!visited_nodes.count(neighbor_node)) { // Only consider unvisited nodes
                            new_path.push_back(neighbor_node);
                            pq.push({cost + weight, {neighbor_node, new_path}});
                        }
                    }
                }
            }

            return shortest_paths;
        }
};

Node *find_node(vector<Node *> &nodes, int id)
{
    for (Node *node : nodes)
    {
        if (node->id == id)
            return node;
    }
    return nullptr;
}

Link *find_link(vector<Link *> links, Node *node1, Node *node2)
{
    for (Link *link : links)
    {
        if ((link->node1->id == node1->id && link->node2->id == node2->id) || (link->node1->id == node2->id && link->node2->id == node1->id))
            return link;
    }
    return nullptr;
}

VNF *find_vnf(vector<VNF *> &vnfs, int id)
{
    for (VNF *vnf : vnfs)
    {
        if (vnf->id == id)
            return vnf;
    }
    return nullptr;
}

SET *find_set(vector<SET *> &sets, int id)
{
    for (SET *set : sets)
    {
        if (set->id == id)
            return set;
    }
    return nullptr;
}

Link *create_link(Node *node1, Node *node2, int bandwidth, double delay, vector<Link *> &links)
{
    Link *link = new Link();
    link->node1 = node1;
    link->node2 = node2;
    link->bandwidthCapacity = bandwidth;
    link->bandwidthAvailable = bandwidth;
    link->delay = delay;
    link->active = true;
    links.push_back(link);

    Link *revLink = new Link();
    revLink->node1 = node2;
    revLink->node2 = node1;
    revLink->bandwidthCapacity = bandwidth;
    revLink->bandwidthAvailable = bandwidth;
    revLink->delay = delay;
    revLink->active = true;

    link->revLink = revLink;
    revLink->revLink = link;

    links.push_back(revLink);

    return link;
}

void printNodes(vector<Node *> &nodes)
{
    cout<<"Nodes : \n";
    for (Node *node : nodes)
    {
        cout << "Node " << node->id << " : \n" << 
        " Type: " << (node->type == NodeType::GROUND ? "Ground" : (node->type == NodeType::LEO ? "LEO" : "HEO")) << "\n" <<
        " Height: " << node->height <<
        " cpu Capacity: " << node->cpuCapacity << 
        " Memory Capacity: " << node->memCapacity << 
        " Cpu Available: " << node->cpuAvailable << 
        " Memory Available: " << node->memAvailable << 
        " Communication Range: " << node->communicationRange << 
        " Velocity: " << node->velocity <<
        " Normal Vector: " << " A: " << node->normalVector.a << " B: " << node->normalVector.b << " C: " << node->normalVector.c <<
        " \nPosition:"
        " X: " << node->position.x << " Y: " << node->position.y << " Z: " << node->position.z << endl;
    }
}

void printLinks(vector<Link *> &links)
{
    cout<<"Links : \n";
    for (Link *link : links)
    {
        cout << "Link between " << link->node1->id << " and " << link->node2->id << " : \n" << 
        "Bandwidth Capacity: " << link->bandwidthCapacity << 
        " Bandwidth Available: " << link->bandwidthAvailable << 
        " Delay: " << link->delay << 
        " Active: " << link->active << endl;
    }
}

void printVNFs(vector<VNF *> &vnfs)
{
    cout<<"VNFs: \n";
    for (VNF *vnf : vnfs)
    {
        cout << "VNF " << vnf->id << " : \n" << 
        "Cpu Requirement: " << vnf->cpuRequirement << 
        " Memory Requirement: " << vnf->memRequirement << 
        " Delay: " << vnf->delay << endl;
    }
}

void printSFCR(vector<SFCR *> &sfcrs)
{
    cout<<"SFCRS: \n";
    cout<<"Number of SFCRs: "<<sfcrs.size()<<endl;
    for (SFCR *sfc : sfcrs)
    {
        cout << "SFCR " << sfc->id << " : \n" << 
        "Source Node: " << sfc->sourceNode->id << 
        " Bandwidth Requirement: " << sfc->bandwidthRequirement << 
        " Max Delay: " << sfc->maxDelay << 
        " Lifetime: " << sfc->lifetime << endl;
    }
}

void printCluster(Cluster *cluster){
    cout<<"Cluster: \n";
    cout<<"Source Node: "<<cluster->sourceNode->id<<endl;
    cout<<"Center: X: "<<cluster->center.x<<" Y: "<<cluster->center.y<<" Z: "<<cluster->center.z<<endl;
    cout<<"Nodes: ";
    for(Node *node : cluster->nodes){
        cout<<node->id<<" ";
    }
    cout<<endl;
}

void initializeNodesAndLinks(vector<Node *> &nodes, vector<Link *> &links)
{
    int numNodes, type;
    ifstream thefile("groundnodes.txt");

    if (!thefile.is_open())
    {
        cerr << "Error : Unable to open file" << endl;
        exit(EXIT_FAILURE);
    }

    //ground nodes input
    thefile >> numNodes;
    for(int i = 0; i < numNodes; i++){
        Node *node = new Node();
        node->type = NodeType::GROUND;
        thefile >> node->id >> node->cpuCapacity >> node->memCapacity >> type >> node->position.x >> node->position.y >> node->position.z;

        node->cpuAvailable = node->cpuCapacity;
        node->memAvailable = node->memCapacity;
        node->velocity = 0;
        node->communicationRange = 0;
        node->normalVector.a = 0;
        node->normalVector.b = 0;
        node->normalVector.c = 0;
        nodes.push_back(node);
    }

    //satellites input
    thefile.close();
    thefile.open("satellites.txt");
    thefile >> numNodes;
    for(int i = 0; i < numNodes; i++){
        Node *node = new Node();
        thefile >> node->id >> node->cpuCapacity >> node->memCapacity >> type >> node->position.x >> node->position.y >> node->position.z >> node->height >> node->communicationRange >> node->normalVector.a >> node->normalVector.b >> node->normalVector.c;
        node->cpuAvailable = node->cpuCapacity;
        node->memAvailable = node->memCapacity;
        node->velocity = sqrt(G*EARTH_MASS/(EARTH_RADIUS + node->height));
        if(type == 1){
            node->type = NodeType::LEO;
        }
        else if(type == 2){
            node->type = NodeType::HEO;
        }
        nodes.push_back(node);
    }

    thefile.close();

    thefile.open("links.txt");
    int numLinks;
    thefile >> numLinks;

    for (int i = 0; i < numLinks; ++i) 
    {
        int id1, id2, bandwidthCapacity;
        double delay;
        thefile >> id1 >> id2 >> delay >> bandwidthCapacity;
        Node *node1 = find_node(nodes, id1);
        Node *node2 = find_node(nodes, id2);
        //if a satellite is present in either end of the link, the delay is calculated based on radio wave propagation speed
        if(delay == 0){
            double distance = sqrt(pow(node1->position.x - node2->position.x, 2) + pow(node1->position.y - node2->position.y, 2) + pow(node1->position.z - node2->position.z, 2));
            delay = distance*1000/200000000; // delay in milliseconds
            // cout<<"delay: "<<delay<<endl;

        }
        // inactive links between ground nodes and satellites will be provided in input
        Link *link = create_link(node1, node2, bandwidthCapacity, delay, links); // create an active two way link and adds it to the links vector

        //configure active and inactive links based on the type of nodes
        if(link->node1->type == NodeType::GROUND && link->node2->type == NodeType::GROUND){
            link->active = true;
            link->revLink->active = true;
        }
        else{
            link->active = false;
            link->revLink->active = false;
        }
    }
    thefile.close();
}

void initializeVNFs(vector<VNF *> &vnfs)
{
    ifstream thefile("vnfs.txt");

    if (!thefile.is_open())
    {
        cerr << "Error: Unable to open file " << endl;
        exit(EXIT_FAILURE);
    }

    int numVNFs;
    thefile >> numVNFs;

    for (int i = 1; i <= numVNFs; ++i)
    {
        VNF *vnf = new VNF();
        thefile >> vnf->id >> vnf->cpuRequirement >> vnf->memRequirement >> vnf->delay;
        vnfs.push_back(vnf);
    }

    thefile.close();
}

void initializeSFCR(vector<SFCR *> &sfcrs, vector<VNF *> &allVNFs, vector<Node *> &nodes)
{
    ifstream thefile("sfcs.txt");

    if (!thefile.is_open())
    {
        cerr << "Error: Unable to open file " << endl;
        exit(EXIT_FAILURE);
    }

    int numSFCR;
    thefile >> numSFCR;

    for (int i = 1; i <= numSFCR; ++i)
    {
        SFCR *sfc = new SFCR();
        sfc->id = i;

        int nVNF = 0, sourceId;
        thefile >> nVNF;
        thefile >> sourceId >> sfc->lifetime >> sfc->bandwidthRequirement >> sfc->maxDelay;

        sfc->sourceNode = find_node(nodes, sourceId);

        for (int j = 1; j <= nVNF; ++j)
        {
            SET *set;
            int vnfID, setNumber;
            thefile >> vnfID >> setNumber;
            if (setNumber > sfc->sets.size())
            {
                set = new SET();
                set->id = setNumber;
                sfc->sets.push_back(set);
            }
            else
            {
                set = sfc->sets[setNumber-1];
            }

            VNF *it = find_vnf(allVNFs, vnfID);
            if (it)
            {
                set->vnfs.push_back(it);
            }
            else
            {
                cerr << "Error: VNF with ID " << vnfID << " not found." << endl;
                exit(EXIT_FAILURE);
            }
        }   
    
        sfcrs.push_back(sfc);
    }

    thefile.close();
}

bool compareDistance(Node* node1,  Node* node2, Position center){
    double d1 = sqrt(pow(node1->position.x - center.x, 2) + pow(node1->position.y - center.y, 2) + pow(node1->position.z - center.z, 2));
    double d2 = sqrt(pow(node2->position.x - center.x, 2) + pow(node2->position.y - center.y, 2) + pow(node2->position.z - center.z, 2));
    return d1 < d2; // ascending order of distance
}

Cluster *create_cluster(Node *sourceNode, vector<Node *> &nodes){
    //find the cluster where the source node resides
    Cluster *cluster = new Cluster();
    cluster->sourceNode = sourceNode;
    for(Node *node : nodes){
        if(node->type == NodeType::GROUND){
            double distance = sqrt(pow(node->position.x - sourceNode->position.x, 2) + pow(node->position.y - sourceNode->position.y, 2));
            if( distance <= CLUSTER_RADIUS){
                cluster->nodes.push_back(node);
            }
        }
    }

    // find the center of the cluster by averaging all the coordinates of cluster ground nodes
    int x = 0, y = 0, z = 0, n = cluster->nodes.size();
    for(Node *node : cluster->nodes){
        x += node->position.x;
        y += node->position.y;
        z += node->position.z;
    }
    cluster->center.x = x/n;
    cluster->center.y = y/n;
    cluster->center.z = z/n;
}

vector<Node *> algorithm1(SFCR *sfc, Cluster *cluster, vector<Node *> &nodes, vector<Link *> &links, vector<VNF *> &VNFS){
    vector<Node *> available_satellites; 
    //finding available satellites for the sfcr
    for(Node *node : nodes){
        if(node->type == NodeType::HEO || node->type == NodeType::LEO){
            double c = (node->normalVector.a * cluster->center.x + node->normalVector.b * cluster->center.y + node->normalVector.c * cluster->center.z) / sqrt(pow(node->normalVector.a, 2) + pow(node->normalVector.b, 2) + pow(node->normalVector.c, 2));
            double ratio = ((2*pow(EARTH_RADIUS,2)) + pow(node->height,2) - pow(node->communicationRange, 2) + (2*EARTH_RADIUS*node->height))/((2*sqrt(EARTH_RADIUS*EARTH_RADIUS - c*c))*(EARTH_RADIUS+node->height));
            if(ratio > 1 || ratio < -1) 
                continue;
            // cout<<"ratio: "<<ratio<<endl;
            double visibility_window = 2*(EARTH_RADIUS+node->height)*acos(ratio)*1000/node->velocity;
            // cout<<"visibility window: "<<visibility_window<<endl; 
            if (visibility_window > sfc->lifetime)
                available_satellites.push_back(node);
        }
    }
    sort(available_satellites.begin(), available_satellites.end(), [&](Node *node1, Node *node2){
        return compareDistance(node1, node2, cluster->center);
    });

    //make links from each satellite to ground nodes in cluster active
    for(Node *satellite : available_satellites){
        for(Node *node : cluster->nodes){
            Link *link = find_link(links, satellite, node);
            if(link){
                link->active = true;
                link->revLink->active = true;
            }
            else{
                cout<<"link between satellite "<<satellite->id<<" and ground node "<<node->id<<"not found\n";
            }
        }
    }

    return available_satellites;
}

vector<pair<double,vector<int>>> kShortestPath(Node* source, vector<Node*>& nodes, vector<Link*>& links, int maxNodes, int K = 3) {
    Graph g;
    for(Link* link : links){
        g.add_edge(link->node1->id, link->node2->id, link->delay);
    }
    return g.dijkstra(source->id, maxNodes-1, K);
}

bool deployVNFsOnPath(vector<Node *> path, vector<SET *> sets)
{
    vector<SET *> parallelSets;
    vector<SET *> seriesSets;

    for (SET *set : sets)
    {
        if (set->vnfs.size() > 1)
        {
            parallelSets.push_back(set);
        }
        else
        {
            seriesSets.push_back(set);
        }
    }

    // Deploy
    for (SET *set : parallelSets)
    {
        bool found = false;
        for (VNF *vnf : set->vnfs)
        {
            found = false;
            for (Node *node : path)
            {
                if (node->cpuAvailable >= vnf->cpuRequirement && node->memAvailable >= vnf->memRequirement)
                {
                    node->cpuAvailable -= vnf->cpuRequirement;
                    node->memAvailable -= vnf->memRequirement;
                    set->deployedNodes.push_back(node);
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                // if VNF placement failed
                for (SET *set_ : parallelSets)
                {
                    for (Node *node : set_->deployedNodes)
                    {
                        node->cpuAvailable += vnf->cpuRequirement;
                        node->memAvailable += vnf->memRequirement;
                    }
                }
                break;
            }
        }
        if (!found)
        {
            break;
        }
    }

    // series VNF sets
    for (SET *set : seriesSets)
    {
        bool found = false;
        for (VNF *vnf : set->vnfs)
        {
            found = false;
            for (Node *node : path)
            {
                if (node->cpuAvailable >= vnf->cpuRequirement && node->memAvailable >= vnf->memRequirement)
                {
                    node->cpuAvailable -= vnf->cpuRequirement;
                    node->memAvailable -= vnf->memRequirement;
                    set->deployedNodes.push_back(node);
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                // if VNF placement failed
                for (SET *set_ : sets)
                {
                    for (Node *node : set_->deployedNodes)
                    {
                        node->cpuAvailable += vnf->cpuRequirement;
                        node->memAvailable += vnf->memRequirement;
                    }
                }
                return false; // rejected
            }
        }
    }

    return true; // successfully deployed
}

void algorithm2(SFCR *sfcr, Cluster *cluster, vector<Node *> &available_satellites, vector<Link *> &links){
    vector<Node *> totalAvailableNodes;
    vector<Link *> totalAvailableLinks;
    for(Node *node : cluster->nodes){ // cluster->nodes should contain the source node as well
        totalAvailableNodes.push_back(node);
    }

    for(Node *node : available_satellites){
        totalAvailableNodes.push_back(node);
    }

    // sorting available nodes based on available resources
    sort(totalAvailableNodes.begin(), totalAvailableNodes.end(), [&](Node *node1, Node *node2){
        return node1->cpuAvailable + node1->memAvailable > node2->cpuAvailable + node2->memAvailable;
    });

    //find the available links for this set of clusters and satellite
    //Step 1: find the available links between ground nodes in the cluster
    for(Node *node1 : cluster->nodes){
        for(Node *node2 : cluster->nodes){
            Link *link = find_link(links, node1, node2);
            if(link){
                totalAvailableLinks.push_back(link);
                totalAvailableLinks.push_back(link->revLink);
            }
        }
    } 

    //Step 2: find the available links between satellite and ground nodes in the cluster
    for(Node *node1 : cluster->nodes){
        for(Node *node2 : available_satellites){
            Link *link = find_link(links, node1, node2);
            if(link){
                totalAvailableLinks.push_back(link);
                totalAvailableLinks.push_back(link->revLink);
            }
        }
    }

    //Find the k shortest paths with from the source node to other nodes where the number of nodes is number of vnfs and weight is the link delay
    int numberOfVNFs = 0;
    for(SET *set : sfcr->sets){
        for(VNF *vnf : set->vnfs){
            numberOfVNFs++;
        } 
    } 
    
    cout<<"k shortest paths: \n";
    auto kShortestPaths = kShortestPath(cluster->sourceNode, totalAvailableNodes, totalAvailableLinks, numberOfVNFs);
    cout<<"k shortest paths found\n";
    cout<<"Number of sets: "<<sfcr->sets.size()<<"\n";
    cout<<"Number of VNfs: "<<numberOfVNFs<<endl;
    //print the k shortest paths
    for (int idx = 0; idx < kShortestPaths.size(); ++idx) {
        cout << "Path " << idx + 1 << ": Cost=" << kShortestPaths[idx].first << ", Path=";
        for (int node : kShortestPaths[idx].second) {
            cout << node << " ";
        }
        cout << endl;
    }

    vector<pair<double, vector<Node *>> > Kpaths;
    for(int i = 0; i < kShortestPaths.size(); i++){
        vector<Node *> path;
        for(int node : kShortestPaths[i].second){
            path.push_back(find_node(totalAvailableNodes, node));
        }
        Kpaths.push_back({kShortestPaths[i].first, path});
    }

    for(int i = 0; i < Kpaths.size(); i++){
        vector<Node *> path = Kpaths[i].second;
        if(deployVNFsOnPath(path, sfcr->sets)){
            cout<<"Succesfully deployed onto Path " << i+1 << endl;
            break;
        }
        cout<<"Unsuccesfull deployment on Path " << i+1 << endl;
    }
}

int main(){
    vector <Node *> nodes;
    vector <Link *> links;
    vector <VNF *> VNFS;
    vector<SFCR *> SFCRS;

    initializeNodesAndLinks(nodes, links);
    initializeVNFs(VNFS);
    initializeSFCR(SFCRS, VNFS, nodes);
    // printNodes(nodes);
    // printLinks(links);
    // printVNFs(VNFS);
    // printSFCR(SFCRS);
    for(SFCR *sfc: SFCRS){
        Cluster *cluster = create_cluster(sfc->sourceNode, nodes);
        vector<Node *> available_satellites = algorithm1(sfc, cluster, nodes, links, VNFS);
        // cout<<"number of satellites available: "<<available_satellites.size()<<endl;
        // printCluster(cluster);
        cout<<"satellites found"<<endl;
        // printNodes(available_satellites);
        algorithm2(sfc, cluster, available_satellites, links);
        break;
    }

    return 0;
}