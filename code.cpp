#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <tuple>
#include <limits>
#include <queue>
#include <math.h>

using namespace std;

#define CLUSTER_RADIUS 100 // in km
#define EARTH_RADIUS 6371
#define EARTH_MASS 5.972e24
#define G 6.674e-11

enum class NodeType
{
    GROUND,
    LEO,
    HEO
};

struct Position
{
    double x;
    double y;
    double z;
};

struct Normal
{
    double a;
    double b;
    double c;
};

struct Node
{
    int id;
    int cpuCapacity;
    int memCapacity;
    int cpuAvailable;
    int memAvailable;
    double velocity;           // 0 for ground nodes
    double communicationRange; // this value is for satellites only
    int height;                // this value is for satellites only
    NodeType type;
    Position position;
    Normal normalVector;
};

struct Cluster
{
    vector<Node *> nodes;
    Node *sourceNode;
    Position center;
};

struct Link
{
    Node *node1;
    Node *node2;
    int bandwidthCapacity;
    int bandwidthAvailable;
    double delay;
    bool active;
    Link *revLink;
};

struct VNF
{
    int id;
    int cpuRequirement;
    int memRequirement;
    int delay;
};

struct SET
{
    int id;
    vector<VNF *> vnfs;
    vector<Node *> deployedNodes;
};

struct SFCR
{
    int id;
    Node *sourceNode;
    vector<SET *> sets;
    int bandwidthRequirement;
    int maxDelay;
    int lifetime;
};

struct Path
{
    vector<Node *> nodes;
    int cost;
};

struct Compare {
    bool operator()(const std::pair<double, std::vector<Node *>>& a, const std::pair<double, std::vector<Node *>>& b) const {
        return a.first > b.first; // Compare based on the first element (double value)
    }
};

class Graph
{
private:
    unordered_map<Node *, vector<pair<Node *, double>>> graph;

public:
    void add_edge(Node *u, Node *v, double w)
    {
        graph[u].push_back({v, w});
        // graph[v].push_back({u, w});
    }

    void consolidate_k_shortest(priority_queue<pair<double, vector<Node *>>> &k_shortest, vector<Node *> &new_path, double new_cost, int k)
    {
        if (k_shortest.size() < k)
        {
            k_shortest.push({new_cost, new_path});
        }
        else if (new_cost < k_shortest.top().first)
        {
            k_shortest.pop();
            k_shortest.push({new_cost, new_path});
        }
    }

    // priority_queue<pair<double, vector<Node *>>> BFS(Node *src, int min_nodes, int max_nodes, int k = 10)
    std::priority_queue<std::pair<double, std::vector<Node *>>, std::vector<std::pair<double, std::vector<Node *>>>, Compare> BFS(Node *src, int min_nodes, int max_nodes, int k)
    {
        vector<tuple<Node *, vector<Node *>, unordered_set<Node *>, double>> curr = {make_tuple(src, vector<Node *>({src}), unordered_set<Node *>{src}, 0.0)};
        // priority_queue<pair<double, vector<Node *>>> k_shortest;

        std::priority_queue<std::pair<double, std::vector<Node *>>, std::vector<std::pair<double, std::vector<Node *>>>, Compare> k_shortest;

        while (!curr.empty())
        {
            vector<tuple<Node *, vector<Node *>, unordered_set<Node *>, double>> next_nodes;
            for (auto &current_node : curr)
            {
                Node *start_node = get<0>(current_node);
                vector<Node *> path_so_far = get<1>(current_node);
                unordered_set<Node *> visited = get<2>(current_node);
                double path_cost = get<3>(current_node);
                for (auto &neighbour : graph[start_node])
                {
                    Node *neighbor = neighbour.first;
                    // cout<<start_node->id<<" "<<neighbor->id<<endl;
                    double LinkWeight = neighbour.second;
                    if (visited.find(neighbor) != visited.end())
                    {
                        continue;
                    }
                    vector<Node *> new_path(path_so_far.begin(), path_so_far.end());
                    new_path.push_back(neighbor);
                    unordered_set<Node *> new_visited(visited.begin(), visited.end());
                    new_visited.insert(neighbor);
                    double new_cost = path_cost + LinkWeight;
                    if (new_path.size() >= min_nodes)
                    {
                        // consolidate_k_shortest(k_shortest, new_path, new_cost, k);
                        k_shortest.push({new_cost, new_path});
                        // for (int i = 0; i < new_path.size(); i++)
                        // {
                        //     cout<<new_path[i]->id<<" ";
                        // }
                        // cout<<endl;
                    }
                    if (new_path.size() <= max_nodes)
                    {
                        next_nodes.push_back(make_tuple(neighbor, new_path, new_visited, new_cost));

                        
                    }
                }
            }
            curr = next_nodes;
        }
        return k_shortest;
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
        if ((link->node1->id == node1->id && link->node2->id == node2->id) /*|| (link->node1->id == node2->id && link->node2->id == node1->id)*/)
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
    cout << "Nodes : \n";
    for (Node *node : nodes)
    {
        cout << "Node " << node->id << " : \n"
             << " Type: " << (node->type == NodeType::GROUND ? "Ground" : (node->type == NodeType::LEO ? "LEO" : "HEO")) << "\n"
             << " Height: " << node->height << " cpu Capacity: " << node->cpuCapacity << " Memory Capacity: " << node->memCapacity << " Cpu Available: " << node->cpuAvailable << " Memory Available: " << node->memAvailable << " Communication Range: " << node->communicationRange << " Velocity: " << node->velocity << " Normal Vector: "
             << " A: " << node->normalVector.a << " B: " << node->normalVector.b << " C: " << node->normalVector.c << " \nPosition:"
                                                                                                                      " X: "
             << node->position.x << " Y: " << node->position.y << " Z: " << node->position.z << endl;
    }
}

void printLinks(vector<Link *> &links)
{
    cout << "Links : \n";
    for (Link *link : links)
    {
        cout << "Link between " << link->node1->id << " and " << link->node2->id << " : \n"
             << "Bandwidth Capacity: " << link->bandwidthCapacity << " Bandwidth Available: " << link->bandwidthAvailable << " Delay: " << link->delay << " Active: " << link->active << endl;
    }
}

void printVNFs(vector<VNF *> &vnfs)
{
    cout << "VNFs: \n";
    for (VNF *vnf : vnfs)
    {
        cout << "VNF " << vnf->id << " : \n"
             << "Cpu Requirement: " << vnf->cpuRequirement << " Memory Requirement: " << vnf->memRequirement << " Delay: " << vnf->delay << endl;
    }
}

void printSFCR(vector<SFCR *> &sfcrs)
{
    cout << "SFCRS: \n";
    cout << "Number of SFCRs: " << sfcrs.size() << endl;
    for (SFCR *sfc : sfcrs)
    {
        cout << "SFCR " << sfc->id << " : \n\t"
             << "Source Node: " << sfc->sourceNode->id << " Bandwidth Requirement: " << sfc->bandwidthRequirement << " Max Delay: " << sfc->maxDelay << " Lifetime: " << sfc->lifetime << endl;
        for (SET *set : sfc->sets)
        {
            cout << "\tSet " << set->id << " : \tVNFs:";
            for (VNF *vnf : set->vnfs)
            {
                cout << vnf->id << " ";
            }
            cout<<endl;
        }
    }
}

void printCluster(Cluster *cluster)
{
    cout << "Cluster: \n";
    cout << "\tSource Node: " << cluster->sourceNode->id << endl;
    cout << "\tCenter: X: " << cluster->center.x << " Y: " << cluster->center.y << " Z: " << cluster->center.z << endl;
    cout << "\tCluster Nodes: ";
    for (Node *node : cluster->nodes)
    {
        cout << node->id << " ";
    }
    cout << endl;
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

    // ground nodes input
    thefile >> numNodes;
    for (int i = 0; i < numNodes; i++)
    {
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

    // satellites input
    thefile.close();
    thefile.open("satellites.txt");
    thefile >> numNodes;
    for (int i = 0; i < numNodes; i++)
    {
        Node *node = new Node();
        thefile >> node->id >> node->cpuCapacity >> node->memCapacity >> type >> node->position.x >> node->position.y >> node->position.z >> node->height >> node->communicationRange >> node->normalVector.a >> node->normalVector.b >> node->normalVector.c;
        node->cpuAvailable = node->cpuCapacity;
        node->memAvailable = node->memCapacity;
        // cout<<"height: "<<node->height;
        node->velocity = EARTH_RADIUS  * sqrt(9.8 / (EARTH_RADIUS + node->height)*1000);
        // cout<<" velocity: "<<node->velocity<<endl;
        if (type == 1)
        {
            node->type = NodeType::LEO;
        }
        else if (type == 2)
        {
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
        // if a satellite is present in either end of the link, the delay is calculated based on radio wave propagation speed
        if (delay == 0)
        {
            double distance = sqrt(pow(node1->position.x - node2->position.x, 2) + pow(node1->position.y - node2->position.y, 2) + pow(node1->position.z - node2->position.z, 2));
            delay = distance * 1000 / 200000000; // delay in milliseconds
            // cout<<"delay: "<<delay<<endl;
        }
        // inactive links between ground nodes and satellites will be provided in input
        Link *link = create_link(node1, node2, bandwidthCapacity, delay, links); // create an active two way link and adds it to the links vector

        // configure active and inactive links based on the type of nodes
        if (link->node1->type == NodeType::GROUND && link->node2->type == NodeType::GROUND)
        {
            link->active = true;
            link->revLink->active = true;
        }
        else
        {
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
                set = sfc->sets[setNumber - 1];
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

bool compareDistance(Node *node1, Node *node2, Position center)
{
    double d1 = sqrt(pow(node1->position.x - center.x, 2) + pow(node1->position.y - center.y, 2) + pow(node1->position.z - center.z, 2));
    double d2 = sqrt(pow(node2->position.x - center.x, 2) + pow(node2->position.y - center.y, 2) + pow(node2->position.z - center.z, 2));
    return d1 < d2; // ascending order of distance
}

Cluster *create_cluster(Node *sourceNode, vector<Node *> &nodes)
{
    // find the cluster where the source node resides
    Cluster *cluster = new Cluster();
    cluster->sourceNode = sourceNode;
    for (Node *node : nodes)
    {
        if (node->type == NodeType::GROUND)
        {
            double distance = sqrt(pow(node->position.x - sourceNode->position.x, 2) + pow(node->position.y - sourceNode->position.y, 2));
            if (distance <= CLUSTER_RADIUS)
            {
                cluster->nodes.push_back(node);
            }
        }
    }

    // find the center of the cluster by averaging all the coordinates of cluster ground nodes
    int x = 0, y = 0, z = 0, n = cluster->nodes.size();
    for (Node *node : cluster->nodes)
    {
        x += node->position.x;
        y += node->position.y;
        z += node->position.z;
    }
    cluster->center.x = x / n;
    cluster->center.y = y / n;
    cluster->center.z = z / n;
}

vector<Node *> algorithm1(SFCR *sfc, Cluster *cluster, vector<Node *> &nodes, vector<Link *> &links, vector<VNF *> &VNFS)
{
    vector<Node *> available_satellites;
    // finding available satellites for the sfcr
    for (Node *node : nodes)
    {
        if (node->type == NodeType::HEO || node->type == NodeType::LEO)
        {
            double c = (node->normalVector.a * cluster->center.x + node->normalVector.b * cluster->center.y + node->normalVector.c * cluster->center.z) / sqrt(pow(node->normalVector.a, 2) + pow(node->normalVector.b, 2) + pow(node->normalVector.c, 2));
            double ratio = ((2 * pow(EARTH_RADIUS, 2)) + pow(node->height, 2) - pow(node->communicationRange, 2) + (2 * EARTH_RADIUS * node->height)) / ((2 * sqrt(EARTH_RADIUS * EARTH_RADIUS - c * c)) * (EARTH_RADIUS + node->height));
            if (ratio > 1 || ratio < -1)
                continue;
            // cout<<"ratio: "<<ratio<<endl;
            double visibility_window = 2 * (EARTH_RADIUS + node->height) * acos(ratio) / node->velocity;
            // cout<<"visibility window: "<<visibility_window<< " sfc lifetime "<<sfc->lifetime<<endl;
            if (visibility_window > sfc->lifetime)
                available_satellites.push_back(node);
        }
    }
    sort(available_satellites.begin(), available_satellites.end(), [&](Node *node1, Node *node2)
         { return compareDistance(node1, node2, cluster->center); });

    // make links from each satellite to ground nodes in cluster active
    for (Node *satellite : available_satellites)
    {
        for (Node *node : cluster->nodes)
        {
            Link *link = find_link(links, satellite, node);
            if (link)
            {
                link->active = true;
                link->revLink->active = true;
            }
            else
            {
                cout << "link between satellite " << satellite->id << " and ground node " << node->id << "not found\n";
            }
        }
    }

    return available_satellites;
}

std::priority_queue<std::pair<double, std::vector<Node *>>, std::vector<std::pair<double, std::vector<Node *>>>, Compare> kShortestPath(Node *source, vector<Link *> &links, int minNodes, int maxNodes, vector<Graph> &graphs, int K = 30)
{
    Graph g;
    for (Link *link : links)
    {
        g.add_edge(link->node1, link->node2, link->delay);
        // cout<<"adding link between: "<< link->node1->id<<" "<<link->node2->id<<endl;
    }
    graphs.push_back(g);
    return g.BFS(source, minNodes, maxNodes, K);
}

bool deployVNFsOnPath(vector<Node *> path, vector<SET *> sets)
{
    vector<SET *> parallelSets;
    vector<SET *> seriesSets;

    unordered_map<int, bool> deployedNodes;

    for (Node *node : path)
    {
        deployedNodes[node->id] = false;
    }

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
                if (deployedNodes[node->id])
                    continue;
                if (node->cpuAvailable >= vnf->cpuRequirement && node->memAvailable >= vnf->memRequirement)
                {
                    node->cpuAvailable -= vnf->cpuRequirement;
                    node->memAvailable -= vnf->memRequirement;
                    set->deployedNodes.push_back(node);
                    found = true;
                    deployedNodes[node->id] = true;
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
                return false; // path rejected
            }
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
                if (deployedNodes[node->id])
                    continue;
                if (node->cpuAvailable >= vnf->cpuRequirement && node->memAvailable >= vnf->memRequirement)
                {
                    node->cpuAvailable -= vnf->cpuRequirement;
                    node->memAvailable -= vnf->memRequirement;
                    set->deployedNodes.push_back(node);
                    found = true;
                    deployedNodes[node->id] = true;
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

void algorithm2(SFCR *sfcr, Cluster *cluster, vector<Node *> &available_satellites, vector<Link *> &links, vector<Graph> &graphs, int K)
{
    vector<Node *> totalAvailableNodes;
    vector<Link *> totalAvailableLinks;
    for (Node *node : cluster->nodes)
    { // cluster->nodes should contain the source node as well
        totalAvailableNodes.push_back(node);
    }

    for (Node *node : available_satellites)
    {
        totalAvailableNodes.push_back(node);
    }

    // sorting available nodes based on available resources
    sort(totalAvailableNodes.begin(), totalAvailableNodes.end(), [&](Node *node1, Node *node2)
         { return node1->cpuAvailable + node1->memAvailable > node2->cpuAvailable + node2->memAvailable; });

    // find the available links for this set of clusters and satellite
    // Step 1: find the available links between ground nodes in the cluster
    for (Node *node1 : cluster->nodes)
    {
        for (Node *node2 : cluster->nodes)
        {
            Link *link = find_link(links, node1, node2);
            if (link)
            {
                totalAvailableLinks.push_back(link);
                // totalAvailableLinks.push_back(link->revLink);
            }
        }
    }

    // Step 2: find the available links between satellite and ground nodes in the cluster
    for (Node *node1 : cluster->nodes)
    {
        for (Node *node2 : available_satellites)
        {
            Link *link = find_link(links, node1, node2);
            if (link)
            {
                totalAvailableLinks.push_back(link);
                // totalAvailableLinks.push_back(link->revLink);
            }
        }
    }

    // Find the k shortest paths with from the source node to other nodes where the number of nodes is number of vnfs and weight is the link delay
    int numberOfVNFs = 0;
    for (SET *set : sfcr->sets)
    {
        for (VNF *vnf : set->vnfs)
        {
            numberOfVNFs++;
        }
    }
    cout << "Number of sets: " << sfcr->sets.size() << endl;
    cout << "Number of VNfs: " << numberOfVNFs << endl;
    cout << "Number of Nodes: " << totalAvailableNodes.size() << endl;
    cout << "Number of Links: " << totalAvailableLinks.size() << endl;
    cout << "\nFinding k shortest paths: \n";

    std::priority_queue<std::pair<double, std::vector<Node *>>, std::vector<std::pair<double, std::vector<Node *>>>, Compare> kShortestPaths = kShortestPath(cluster->sourceNode, totalAvailableLinks, numberOfVNFs-1, numberOfVNFs, graphs, K);

    // print the k shortest paths
    cout<<"K shortest Paths found:"<< endl;

    for (int idx = 0; idx < K; ++idx)
    {
        cout<<"Deploying nodes in ";
        cout << "Path " << idx + 1 << ": Cost=" << kShortestPaths.top().first << ", Path=";
        vector<Node *> path = kShortestPaths.top().second;
        for (Node *node : path)
        {
            cout << node->id << " ";
        }
        if(deployVNFsOnPath(path, sfcr->sets))
        {
            cout << " VNFs deployed successfully" << endl;
            break;
        }
        else
        {
            cout << " VNFs deployment failed" << endl;
        }
        cout << endl;
        kShortestPaths.pop();

    }
}

int main()
{
    vector<Node *> nodes;
    vector<Link *> links;
    vector<VNF *> VNFS;
    vector<SFCR *> SFCRS;
    vector<Graph> graphs;
    initializeNodesAndLinks(nodes, links);
    initializeVNFs(VNFS);
    initializeSFCR(SFCRS, VNFS, nodes);
    // printNodes(nodes);
    // printLinks(links);
    // printVNFs(VNFS);
    // printSFCR(SFCRS);
    cout << "Total number of SFCRS: " << SFCRS.size() << endl;
    for (SFCR *sfc : SFCRS)
    {
        cout << "\n\nSFCR: " << sfc->id << " Source Node id: " << sfc->sourceNode->id << endl;
        // Create the Cluster
        Cluster *cluster = create_cluster(sfc->sourceNode, nodes);
        // printCluster(cluster);
        // Find the Available Satellites for the Cluster
        vector<Node *> available_satellites = algorithm1(sfc, cluster, nodes, links, VNFS);
        cout << "number of satellites available: " << available_satellites.size() << endl;
        // printNodes(available_satellites);
        // Find KSP and Deploy the Nodes
        int K = 3;
        algorithm2(sfc, cluster, available_satellites, links, graphs, K);
    }

    cout << "size of graphs: " << graphs.size() << endl;
    return 0;
}
