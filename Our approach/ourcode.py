import math

G = 6.6743 * math.pow(10, -11 - 9)  # m^3/kg/s^2
M = 5.972 * math.pow(10, 24)  # kg
R = 6371  # km
K = 10
CLUSTER_RADIUS = 100  # km


class Position:
    def __init__(self, x: float, y: float, z: float):
        self.x = x
        self.y = y
        self.z = z

    def __repr__(self) -> str:
        return f"({self.x}, {self.y}, {self.z})"


class Normal:
    def __init__(self, a, b, c):
        self.a = a
        self.b = b
        self.c = c

    def __repr__(self) -> str:
        return f"({self.a}, {self.b}, {self.c})"


class Node:
    def __init__(self, id: int, data: list):
        self.id = id
        self.cpuCapacity = data[0]
        self.memoryCapacity = data[1]
        self.cpuAvailable = data[0]
        self.memoryAvailable = data[1]
        self.coordinates = data[2]
        self.edges = dict()


class GroundNode(Node):
    def __init__(self, id: int, data: list):
        super().__init__(id, data)
        self.nodeType = "Ground"

    def __repr__(self) -> str:
        return (
            "Node ID: "
            + str(self.id)
            + " | Type: "
            + self.nodeType
            + " | CPU Capacity: "
            + str(self.cpuAvailable)
            + " | Memory Capacity: "
            + str(self.memoryCapacity)
            + " | Memory Available: "
            + str(self.memoryAvailable)
            + "| Coordinates: "
            + str(self.coordinates)
            + "\n"
        )


class Satellite(Node):
    def __init__(self, id: int, data: list):
        super().__init__(id, data)
        self.nodeType = "Satellite"
        self.height = data[3]
        self.velocity = math.sqrt(G * M / (R + self.height))
        self.communicationRange = data[4]
        if data[5] == 1:
            self.nodeType = "LEO"
        elif data[5] == 2:
            self.nodeType = "HEO"
        self.normal = data[6]

    def __repr__(self) -> str:
        return (
            "Node ID: "
            + str(self.id)
            + " | Type: "
            + self.nodeType
            + " | CPU: "
            + str(self.cpuCapacity)
            + " | Memory: "
            + str(self.memoryCapacity)
            + " | Coordinates: "
            + str(self.coordinates)
            + " | Height: "
            + str(self.height)
            + " | Velocity: "
            + str(self.velocity)
            + " | Communication Range: "
            + str(self.communicationRange)
            + " | Normal: "
            + str(self.normal)
        )


class Link:
    def __init__(self, id: int, data: list[Node, Node, float, int]):
        self.id = id
        self.endpoints = (data[0], data[1])
        self.latency = data[2]
        self.bandwidth = data[3]

    def __repr__(self) -> str:
        return (
            "Link ID: "
            + str(self.id)
            + " | Start Node: "
            + str(self.endpoints[0].id)
            + " | End Node: "
            + str(self.endpoints[1].id)
            + " | Bandwidth: "
            + str(self.bandwidth)
            + " | Latency: "
            + str(self.latency)
        )


class Cluster:
    def __init__(self, nodes: dict, source: GroundNode) -> None:
        self.nodes = nodes
        self.sourceNode = source
        n = len(nodes)
        self.center = Position(
            sum([nodes[node_id].coordinates.x for node_id in nodes]) / n,
            sum([nodes[node_id].coordinates.y for node_id in nodes]) / n,
            sum([nodes[node_id].coordinates.z for node_id in nodes]) / n,
        )

    def __repr__(self) -> str:
        return (
            "Cluster Center: "
            + str(self.center)
            + " | Source Node: "
            + str(self.sourceNode.id)
            + " | Nodes: "
            + str([node_id for node_id in self.nodes])
            + "\n"
        )


class VNF:
    def __init__(self, id: int, data: list[int, int, int]):
        self.id = id
        self.cpuRequirement = data[0]
        self.memoryRequirement = data[1]
        self.delay = data[2]

    def __repr__(self) -> str:
        return (
            "VNF ID: "
            + str(self.id)
            + " | CPU: "
            + str(self.cpuRequirement)
            + " | Memory: "
            + str(self.memoryRequirement)
            + " | Delay: "
            + str(self.delay)
        )


class SET:
    def __init__(self, id: int, vnfs: list[VNF]):
        self.id = id
        self.vnfs = vnfs
        self.highestDelay = 0
        self.deployedNodes = dict()

    def addVNF(self, vnf: VNF):
        self.vnfs.append(vnf)
        self.highestDelay = max(self.vnfs, key=lambda x: x.delay).delay

    def __repr__(self) -> str:
        return (
            "SET ID: "
            + str(self.id)
            + " | VNFs: "
            + str([vnf.id for vnf in self.vnfs])
            + " | Highest Delay: "
            + str(self.highestDelay)
            + "\n"
        )


class SFCR:
    def __init__(self, id: int, source: GroundNode, data: list[dict, int, int, int]):
        self.id = id
        self.sourceNode = source
        self.sets = dict()
        self.bandwidthRequirement = data[1]
        self.MaxDelay = data[2]
        self.lifetime = data[3]

        for set_id in data[0]:
            self.sets[set_id] = data[0][set_id]

    def __repr__(self) -> str:
        return (
            "SFCR ID: "
            + str(self.id)
            + " | Source Node: "
            + str(self.sourceNode.id)
            + " | Bandwidth: "
            + str(self.bandwidthRequirement)
            + " | Max Delay: "
            + str(self.MaxDelay)
            + " | Lifetime: "
            + str(self.lifetime)
            + " | \nSETs: "
            + str([self.sets[setid] for setid in self.sets])
            + "\n"
        )


def get_distance(p1: Position, p2: Position) -> float:
    return math.sqrt(
        math.pow(p1.x - p2.x, 2) + math.pow(p1.y - p2.y, 2) + math.pow(p1.z - p2.z, 2)
    )


def initializeNodesAndLinks(nodes: dict, links: list):
    # initialize ground nodes
    f = open("./groundnodes.txt", "r")
    lines = [line.strip().split() for line in f.readlines()]

    number_of_ground_nodes = int(lines[0][0])
    for i in range(1, number_of_ground_nodes + 1):
        node_id = int(lines[i][0])
        cpuCapacity = int(lines[i][1])
        memoryCapacity = int(lines[i][2])
        coordinates = Position(
            float(lines[i][4]), float(lines[i][5]), float(lines[i][6])
        )
        node = GroundNode(node_id, [cpuCapacity, memoryCapacity, coordinates])
        nodes[node.id] = node
        # print(node)
    f.close()

    # initialize satellite nodes
    f = open("./satellites.txt", "r")
    lines = [line.strip().split() for line in f.readlines()]
    number_of_satellites = int(lines[0][0])
    for i in range(1, number_of_satellites + 1):
        node_id = int(lines[i][0])
        cpuCapacity = int(lines[i][1])
        memoryCapacity = int(lines[i][2])
        nodeType = int(lines[i][3])
        coordinates = Position(
            float(lines[i][4]), float(lines[i][5]), float(lines[i][6])
        )
        height = float(lines[i][7])
        communicationRange = int(lines[i][8])
        normal_vector = Normal(float(lines[i][9]), float(lines[i][10]), float(lines[i][11]))
        node = Satellite(
            node_id,
            [
                cpuCapacity,
                memoryCapacity,
                coordinates,
                height,
                communicationRange,
                nodeType,
                normal_vector,
            ],
        )
        nodes[node.id] = node
        # print(node)
    f.close()

    # initialize links
    f = open("links.txt", "r")
    lines = [line.strip().split() for line in f.readlines()]
    number_of_links = int(lines[0][0])
    for i in range(1, number_of_links + 1):
        node1 = nodes[int(lines[i][0])]
        node2 = nodes[int(lines[i][1])]
        # print(node1, node2)
        latency = float(lines[i][2])
        bandwidth = int(lines[i][3])
        link = Link(i, [node1, node2, latency, bandwidth])
        nodes[node1.id].edges[node2.id] = link
        nodes[node2.id].edges[node1.id] = link
        links.append(link)
        # print(link)
    f.close()


def initializeVNFs(vnfs: dict):
    f = open("vnfs.txt", "r")
    lines = [line.strip().split() for line in f.readlines()]
    number_of_vnfs = int(lines[0][0])
    for i in range(1, number_of_vnfs + 1):
        vnf_id = int(lines[i][0])
        cpuRequirement = int(lines[i][1])
        memoryRequirement = int(lines[i][2])
        delay = int(lines[i][3])
        vnf = VNF(vnf_id, [cpuRequirement, memoryRequirement, delay])
        vnfs[vnf_id] = vnf
        # print(vnf)
    f.close()


def initializeSFCRs(sfcrs: dict, vnfs: dict, nodes: dict):
    f = open("sfcs.txt", "r")
    number_of_sfcs = int(f.readline().strip())
    for i in range(1, number_of_sfcs + 1):

        sfcr_id = i
        number_of_vnfrs = int(f.readline().strip())
        sfc_info = f.readline().strip().split()
        source_node = nodes[int(sfc_info[0])]
        lifetime = int(sfc_info[1])
        bandwidth_requirement = int(sfc_info[2])
        max_delay = int(sfc_info[3])

        sets = dict()
        for j in range(number_of_vnfrs):
            line = f.readline().strip().split()
            vnf = vnfs[int(line[0])]
            set_id = int(line[1])
            if set_id not in sets:
                sets[set_id] = SET(set_id, [])
                sets[set_id].addVNF(vnf)
            else:
                sets[set_id].addVNF(vnf)
            # print(sets[set_id])

        sfcr = SFCR(
            sfcr_id, source_node, [sets, bandwidth_requirement, max_delay, lifetime]
        )
        sfcrs[sfcr_id] = sfcr
        # print(sfcr)


def find_cluster_nodes(nodes: dict, source: GroundNode):
    cluster_nodes = dict()
    for node_id in nodes:
        node = nodes[node_id]
        if node.nodeType != "Ground":
            continue
        distance = get_distance(node.coordinates, source.coordinates)
        if distance <= CLUSTER_RADIUS:
            cluster_nodes[node_id] = node
    return cluster_nodes


def find_available_satellites(sfc: SFCR, cluster: Cluster, nodes: dict):
    available_satellites = dict()
    for node_id in nodes:
        node = nodes[node_id]
        if node.nodeType == "Ground":
            continue
        distance = get_distance(node.coordinates, cluster.center)
        shortest_distance = (
            (node.normal.a * cluster.center.x)
            + (node.normal.b * cluster.center.y)
            + (node.normal.c * cluster.center.z)
        ) / math.sqrt(
            math.pow(node.normal.a, 2)
            + math.pow(node.normal.b, 2)
            + math.pow(node.normal.c, 2)
        )
        ratio = (
            (2 * R * R)
            + math.pow(node.height, 2)
            - math.pow(node.communicationRange, 2)
            + (2 * R * node.height)
        ) / (
            2
            * math.sqrt(R * R + shortest_distance * shortest_distance)
            * (R + node.height)
        ) 
        visibility_window = 2 * (R + node.height) * math.acos(ratio) * 1000/ node.velocity

        # print(distance, node.communicationRange, visibility_window, sfc.lifetime)
        if distance - node.communicationRange <= 4000:
            if visibility_window >= sfc.lifetime:
                available_satellites[node_id] = node
            else: print("satellite was in range but visibility window was less")
    return available_satellites


def revise_k_shortest(k_shortest: set, new_path, new_cost, k):
    if len(k_shortest) < k:
        k_shortest.add((new_path, new_cost))
        return
    max_existing = max(k_shortest, key=lambda x: x[1])
    if max_existing[1] > new_cost:
        k_shortest.remove(max_existing)
        k_shortest.add((new_path, new_cost))


def BFS(
    min_nodes: int, max_nodes: int, source: GroundNode, totalNodes: dict, K: int
) -> set:
    curr = [(source, tuple([source]), {source}, 0)]
    k_shortest = set()
    i = 0
    print("Entered BFS")
    for node in totalNodes.values():
        if node.nodeType == "LEO":
            print("leo exists in totalnodes")
            break
    # return k_shortest
    while curr:
        curr.sort(key=lambda x: x[3])
        print(f"{i}th iteration")
        i += 1
        next_nodes = []
        for start_node, path_so_far, visited, path_cost in curr:
            if k_shortest and max(k_shortest, key=lambda x: x[1])[1] < path_cost:
                continue
            for neighbor_id in start_node.edges.keys():
                if neighbor_id not in totalNodes:
                    continue
                link = start_node.edges[neighbor_id]
                neighbor = totalNodes[neighbor_id]
                if neighbor in visited:
                    continue
                new_path = path_so_far + (neighbor,)
                new_visited = visited.union({neighbor})
                new_cost = path_cost + link.latency

                if len(new_path) >= min_nodes:
                    revise_k_shortest(k_shortest, new_path, new_cost, K)
                if len(new_path) <= max_nodes:
                    next_nodes.append((neighbor, new_path, new_visited, new_cost))
        curr = next_nodes
    return k_shortest

def KSP(K: int, sfc: SFCR, totalNodes: dict) -> set:
    min_number_of_nodes = len(sfc.sets)
    max_number_of_nodes = sum([len(sfc.sets[set_id].vnfs) for set_id in sfc.sets])
    source = totalNodes[sfc.sourceNode.id]
    return BFS(min_number_of_nodes, max_number_of_nodes, source, totalNodes, K)


def deallocate_resources(sets: dict):
    for set_id in sets:
        s = sets[set_id]
        for node_id in s.deployedNodes:
            node = s.deployedNodes[node_id][0]
            for i in range(len(s.deployedNodes[node_id][1])):
                node.cpuAvailable += s.deployedNodes[node_id][1][i].cpuRequirement
                node.memoryAvailable += s.deployedNodes[node_id][1][i].memoryRequirement

def deployVNFsOnPath(
    path: list[Node], sets: dict, max_delay: int, linkDelay: int
) -> bool:
    parallel_sets = dict()
    series_sets = dict()
    unavailable_nodes = set()
    totalProcessingDelay = 0

    for set_id in sets:
        s = sets[set_id]
        totalProcessingDelay += s.highestDelay
        if len(s.vnfs) == 1:
            series_sets[set_id] = s
        else:
            parallel_sets[set_id] = s

    totalDelay = totalProcessingDelay + linkDelay
    print("Total Delay: ", totalDelay, "Max Delay: ", max_delay)
    if totalDelay > max_delay:
        print("Total Delay Exceeds Maximum Delay")
        return False

    for set_id in parallel_sets:
        s = parallel_sets[set_id]
        for vnf in s.vnfs:
            found = False
            print("Deploying VNF ", vnf.id, " on Path")
            if s.deployedNodes:
                for node_id in s.deployedNodes:
                    node = s.deployedNodes[node_id][0]
                    if (
                        node.cpuAvailable >= vnf.cpuRequirement
                        and node.memoryAvailable >= vnf.memoryRequirement
                    ):
                        node.cpuAvailable -= vnf.cpuRequirement
                        node.memoryAvailable -= vnf.memoryRequirement
                        s.deployedNodes[node_id][1].append(vnf)
                        found = True
                        break
            if found:
                continue

            for node in path:
                if len(unavailable_nodes) == len(path):
                    print("No more nodes available, Deployement Fails")
                    break
                if node.id in unavailable_nodes:
                    continue
                if (
                    node.cpuAvailable >= vnf.cpuRequirement
                    and node.memoryAvailable >= vnf.memoryRequirement
                ):
                    node.cpuAvailable -= vnf.cpuRequirement
                    node.memoryAvailable -= vnf.memoryRequirement
                    s.deployedNodes[node.id] = [node, [vnf]]
                    unavailable_nodes.add(node.id)
                    found = True
                    break

            if not found:
                print("VNF Deployment Failed due to Resource Constraints")
                deallocate_resources(parallel_sets)
                return False

    for set_id in series_sets:
        s = series_sets[set_id]
        for vnf in s.vnfs:
            found = False
            print("Deploying VNF ", vnf.id, " on Path")
            for node in path:
                if len(unavailable_nodes) == len(path):
                    print("No more nodes available, Deployement Fails")
                    break
                if node.id in unavailable_nodes:
                    continue
                if (
                    node.cpuAvailable >= vnf.cpuRequirement
                    and node.memoryAvailable >= vnf.memoryRequirement
                ):
                    node.cpuAvailable -= vnf.cpuRequirement
                    node.memoryAvailable -= vnf.memoryRequirement
                    s.deployedNodes[node.id] = [node, [vnf]]
                    unavailable_nodes.add(node.id)
                    found = True
                    break

            if not found:
                print("VNF Deployment Failed due to Resource Constraints")
                deallocate_resources(parallel_sets)
                deallocate_resources(series_sets)
                return False
    return True


def algorithm2(sfc: SFCR, totalNodes: dict, K: int) -> bool:
    kshortestpaths = KSP(K, sfc, totalNodes)
    kshortestpaths = sorted(kshortestpaths, key=lambda x: x[1])
    [
        print(
            "path = ",
            [(node.id, node.nodeType) for node in path],
            "Link Delay = ",
            cost,
        )
        for path, cost in kshortestpaths
    ]

    i = 0
    for path, cost in kshortestpaths:
        i += 1
        if cost < sfc.MaxDelay:
            if(deployVNFsOnPath(path, sfc.sets, sfc.MaxDelay, cost)):
                print("SFC Deployed on Path", i, ":", [node.id for node in path] , "Successfully")
                return True
        else:
            print("SFC Deployment Failed Due Link Delay > Max Delay")
            return False
    return False


def main():
    nodes = dict()
    links = list()
    vnfs = dict()
    sfcrs = dict()

    initializeNodesAndLinks(nodes, links)
    initializeVNFs(vnfs)
    initializeSFCRs(sfcrs, vnfs, nodes)
    print("read input 👍")
    heos = 0
    leos = 0
    for id in nodes:
        if nodes[id].nodeType == "LEO":
            leos += 1
        elif nodes[id].nodeType == "HEO":
            heos += 1
    
    print("LEOs: ", leos, "HEOs: ", heos)

    number_of_SFCS_accepted = 0
    for sfc_id in sfcrs:
        print()
        print(sfcrs[sfc_id])
        # print("SFC ID: ", sfc_id)
        sfc = sfcrs[sfc_id]
        source_node = sfcrs[sfc_id].sourceNode
        cluster_nodes = find_cluster_nodes(nodes, source_node)
        cluster = Cluster(cluster_nodes, source_node)
        # print("Cluster : ", cluster)
        available_satellites = find_available_satellites(sfc, cluster, nodes)
        
        # print("satellites available: ", str([(Satellite.id, Satellite.nodeType) for Satellite in available_satellites.values()]))
        # continue
        totalNodes = dict(cluster_nodes.items() | available_satellites.items())
        # print([node_id for node_id in totalNodes])
        isSFCDeployed = algorithm2(sfc, totalNodes, K)

        if isSFCDeployed:
            print("SFC Deployed Successfully")
            number_of_SFCS_accepted += 1
    print("Acceptance Ratio: ", number_of_SFCS_accepted / len(sfcrs))
    # break

main()