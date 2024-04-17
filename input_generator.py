import random
import math

MIN_NUMBER_OF_CLUSTERS = 10
MAX_NUMBER_OF_CLUSTERS = 30
CLUSTER_RADIUS = 100
MIN_NODES_PER_CLUSTER = 8
MAX_NODES_PER_CLUSTER = 30
MIN_EDGE_WEIGHT = 10
MAX_EDGE_WEIGHT = 50
MIN_EDGE_BANDWIDTH = 5
MAX_EDGE_BANDWIDTH = 10
MIN_LEO_HEIGHT = 500
MAX_LEO_HEIGHT = 2000
MIN_HEO_HEIGHT = 5000
MAX_HEO_HEIGHT = 35000
EARTH_RADIUS = 6371
LEO_HEO_RATIO = 25


class Link:
    def __init__(self, node1, node2, weight) -> None:
        self.endpoints = (node1, node2)
        self.bandwidth = random.randint(MIN_EDGE_BANDWIDTH, MAX_EDGE_BANDWIDTH)
        self.weight = weight

    def __repr__(self) -> str:
        return (
            f"{self.endpoints[0]} {self.endpoints[1]} {self.weight} {self.bandwidth}\n"
        )


class Node:
    def __init__(self, node_id) -> None:
        self.node_id = node_id
        self.cpu_capacity = random.randint(20, 100)
        self.memory_capacity = random.randint(20, 100)
        # The key contains the node_id of the neighbour.
        # Value contains a Link object for the edge between them
        self.edges = dict()


class GroundNode(Node):
    def __init__(self, node_id, coordinates) -> None:
        super().__init__(node_id)
        self.coordinates = coordinates
        self.node_type = 0

    def __repr__(self) -> str:
        return f"{self.node_id} {self.cpu_capacity} {self.memory_capacity} {self.node_type} {self.coordinates[0]} {self.coordinates[1]} {self.coordinates[2]}\n"


class Satellite(Node):
    def __init__(self, node_id, cluster_coordinates=None, node_type=None):
        super().__init__(node_id)
        self.node_type = random.choice([1]*int(LEO_HEO_RATIO) + [2]) if node_type is None else node_type
        self.height = None
        self.communication_range = 0
        if self.node_type == 1:
            self.height = random.randint(MIN_LEO_HEIGHT, MAX_LEO_HEIGHT)
            self.communication_range = int(math.sqrt((self.height * self.height) + (4*CLUSTER_RADIUS * CLUSTER_RADIUS))) + random.randint(1000, 1500)
        else:
            self.height = random.randint(MIN_HEO_HEIGHT, MAX_HEO_HEIGHT)
            self.communication_range = int(math.sqrt((self.height * self.height) + (CLUSTER_RADIUS * CLUSTER_RADIUS)))

        if cluster_coordinates is None:
            x, y, z = get_random_point_at_height(self.height)
        else:
            x, y, z = cluster_coordinates
            height_ratio = (self.height+EARTH_RADIUS)/EARTH_RADIUS
            x = x*height_ratio
            y = y*height_ratio
            z = z*height_ratio

        self.starting_coordinates = (x, y, z)
        a = random.choice([random.randint(-5, -1), random.randint(1, 5)])
        b = random.choice([random.randint(-10, -1), random.randint(1, 10)])
        c = (-x*a - y*b)/z
        self.normal_vector = (a, b, c )

        

    def __repr__(self) -> str:
        return f"{self.node_id} {self.cpu_capacity} {self.memory_capacity} {self.node_type} {self.starting_coordinates[0]} {self.starting_coordinates[1]} {self.starting_coordinates[2]} {self.height} {self.communication_range} {self.normal_vector[0]} {self.normal_vector[1]} {self.normal_vector[2]}\n"


class VNF:
    def __init__(self, vnf_id):
        self.vnf_id = vnf_id
        self.cpu_requirement = random.randint(1, 5)
        self.memory_requirement = random.randint(1, 5)
        self.delay = random.randint(5, 10)

    def __repr__(self) -> str:
        return f"{self.vnf_id} {self.cpu_requirement} {self.memory_requirement} {self.delay}\n"


class SFC:
    def __init__(self, sfc_id, node_id):
        self.sfc_id = sfc_id
        self.source_node = node_id
        self.vnf_list = list()
        self.vnf_count = 0
        self.bandwidth_requirement = random.randint(1, 5)
        self.max_tolerated_delay = random.randint(50, 100)
        self.lifetime = random.randint(5000, 10000)

    def __repr__(self) -> str:
        repr_str = f"{self.vnf_count}\n{self.source_node} {self.lifetime} {self.bandwidth_requirement} {self.max_tolerated_delay}\n"
        for index, vnf_list in enumerate(self.vnf_list):
            for vnf in vnf_list:
                repr_str += f"{vnf.vnf_id} {index + 1}\n"
        return repr_str


def get_random_point_at_height(height):
    axis_coordinates = []
    x = random.randint(height//3, 2*height//3)
    axis_coordinates.append(x)
    y = random.randint(0, 2*height//3)
    axis_coordinates.append(y)
    axis_coordinates.append(math.sqrt(height**2 - x**2 - y**2))
    random.shuffle(axis_coordinates)
    return tuple(axis_coordinates)


def add_satellite(new_satellite, satellites, links, node_lookup):
    
    for groundnode in node_lookup.values():
        if get_distance(groundnode.coordinates, new_satellite.starting_coordinates) > new_satellite.communication_range + 1000:
            continue
        edge = Link(
            groundnode.node_id,
            new_satellite.node_id,
            get_node_satellite_distance(groundnode, new_satellite)
            / 200000,
        )
        groundnode.edges[new_satellite.node_id] = edge
        new_satellite.edges[groundnode.node_id] = edge
        links.append(edge)
    # for old_satellite in satellites:
    #     edge = Link(
    #         old_satellite.node_id,
    #         new_satellite.node_id,
    #         get_satellite_distance(old_satellite, new_satellite) / 200000,
    #     )
    #     old_satellite.edges[new_satellite.node_id] = edge
    #     new_satellite.edges[old_satellite.node_id] = edge
    #     links.append(edge)
    satellites.append(new_satellite)
    

def get_distance(p1, p2):
    x1, y1, z1 = p1
    x2, y2, z2 = p2
    return math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2) * 1000


def get_node_satellite_distance(node, satellite):
    return get_distance(node.coordinates, satellite.starting_coordinates)


def get_satellite_distance(s1, s2):
    return get_distance(s1.starting_coordinates, s2.starting_coordinates)


def revise_k_shortest(k_shortest, new_path, new_cost, k):
    if len(k_shortest) < k:
        k_shortest.add((new_path, new_cost))
        return
    max_existing = max(k_shortest, key=lambda x: x[1])
    if max_existing[1] > new_cost:
        k_shortest.remove(max_existing)
        k_shortest.add((new_path, new_cost))


def rec_bfs(min_nodes, max_nodes, src: GroundNode, node_lookup, k):
    curr = [(src, tuple(), {src}, 0)]
    k_shortest = set()

    while curr:
        next_nodes = []
        for start_node, path_so_far, visited, path_cost in curr:
            if k_shortest and max(k_shortest, key=lambda x: x[1])[1] < path_cost:
                continue
            for neighbor_id in start_node.edges.keys():
                if neighbor_id not in node_lookup:
                    continue
                link = start_node.edges[neighbor_id]
                neighbor = node_lookup[neighbor_id]
                if neighbor in visited:
                    continue
                new_path = path_so_far + (neighbor,)
                new_visited = visited.union({neighbor})
                new_cost = path_cost + link.weight
                if len(new_path) >= min_nodes:
                    revise_k_shortest(k_shortest, new_path, new_cost, k)
                if len(new_path) <= max_nodes:
                    next_nodes.append((neighbor, new_path, new_visited, new_cost))
        curr = next_nodes
    return k_shortest


def find_k_shortest_paths(k, sfc: SFC, node_lookup):
    min_number_of_nodes = len(sfc.vnf_list)
    max_number_of_nodes = sfc.vnf_count
    src = node_lookup[sfc.source_node]
    return rec_bfs(min_number_of_nodes, max_number_of_nodes, src, node_lookup, k)


def main():
    number_of_clusters = random.randint(MIN_NUMBER_OF_CLUSTERS, MAX_NUMBER_OF_CLUSTERS)
    cluster_list = dict()
    link_list = []
    node_id = 0
    node_lookup = dict()
    for i in range(number_of_clusters):
        x, y, z = get_random_point_at_height(EARTH_RADIUS)
        if (x, y, z) in cluster_list:
            continue
        curr_cluster = []

        for j in range(random.randint(MIN_NUMBER_OF_CLUSTERS, MAX_NUMBER_OF_CLUSTERS)):
            x1 = x + random.randint(-CLUSTER_RADIUS // 3, 0)
            y1 = y + random.randint(-CLUSTER_RADIUS // 3, 0)
            z1 = math.sqrt((EARTH_RADIUS**2 - x1**2 - y1**2))
            new_node = GroundNode(node_id, (x1, y1, z1))
            node_lookup[node_id] = new_node
            node_id += 1
            curr_cluster.append(new_node)
            if len(curr_cluster) > 1:
                connections = random.sample(
                    curr_cluster[:-1],
                    random.choice([1, 1, 1, 1, 2, 3]) % len(curr_cluster),
                )
                for neighbor in connections:
                    edge = Link(
                        neighbor.node_id,
                        new_node.node_id,
                        random.randint(MIN_EDGE_WEIGHT, MAX_EDGE_WEIGHT),
                    )
                    new_node.edges[neighbor.node_id] = edge
                    neighbor.edges[new_node.node_id] = edge
                    link_list.append(edge)
        cluster_list[(x, y, z)] = curr_cluster

    number_of_clusters = len(cluster_list.keys())
    number_of_nodes = node_id

    # print(node_lookup.values())
    file = open("./groundnodes.txt", "w")
    # print(file)
    file.write(f"{number_of_nodes}\n")
    for node in node_lookup.values():
        file.write(f"{node}")

    for _ in range(number_of_nodes // 2, 3 * number_of_nodes // 2):
        node_id1, node_id2 = random.sample(range(number_of_nodes), 2)
        edge = Link(
            node_id1, node_id2, random.randint(MIN_EDGE_WEIGHT, MAX_EDGE_WEIGHT)
        )
        node_lookup[node_id1].edges[node_lookup[node_id2].node_id] = edge
        node_lookup[node_id2].edges[node_lookup[node_id1].node_id] = edge
        link_list.append(edge)

    satellites = []

    number_of_satellites = random.randint(200, 300)
    for _ in range(number_of_satellites):
        new_satellite = Satellite(node_id)
        node_id += 1
        add_satellite(new_satellite, satellites, link_list, node_lookup)

    
    for cluster in cluster_list.keys():
        new_satellite = Satellite(node_id, cluster, 2)
        node_id += 1
        add_satellite(new_satellite, satellites, link_list, node_lookup)
        print(cluster)
        
    number_of_satellites += len(cluster_list)
        
        
    # number_of_satellite_links = random.randint(number_of_satellites, number_of_satellites*2)
    # for _ in range(number_of_satellite_links):
    #     a,b = random.sample(satellites, 2)
    #     edge = Link(a.node_id, b.node_id, random.randint(MIN_EDGE_WEIGHT, MAX_EDGE_WEIGHT))
    #     a.edges[b.node_id] = edge
    #     b.edges[a.node_id] = edge
    #     link_list.append(edge)

    # print(satellites)
    file = open("./satellites.txt", "w")
    # print(file)
    file.write(f"{number_of_satellites}\n")
    for satellite in satellites:
        file.write(f"{satellite}")

    file = open("./links.txt", "w")
    # print(file)
    file.write(f"{len(link_list)}\n")
    for link in link_list:
        file.write(f"{link}")

    number_of_sfcs = random.randint(6, 10)
    number_of_vfns = random.randint(15, 20)

    sfc_list = [
        SFC(i, random.randint(0, number_of_nodes - 1)) for i in range(number_of_sfcs)
    ]
    vnf_list = [VNF(i) for i in range(number_of_vfns)]

    for vnf in vnf_list:
        sfcs_to_insert = random.sample(sfc_list, random.randint(0, 6))
        for sfc in sfcs_to_insert:
            if not sfc.vnf_list or random.random() < 0.8:
                sfc.vnf_list.append([vnf])
            else:
                index = random.randint(0, len(sfc.vnf_list) - 1)
                for i in range(len(sfc.vnf_list)):
                    if len(sfc.vnf_list[(i + index) % len(sfc.vnf_list)]) > 2:
                        continue
                    else:
                        sfc.vnf_list[(i + index) % len(sfc.vnf_list)].append(vnf)
                        break
            sfc.vnf_count += 1

    for sfc in sfc_list:
        if sfc.vnf_count == 0:
            new_vnf = VNF(number_of_vfns)
            vnf_list.append(new_vnf)
            number_of_vfns += 1
            sfc.vnf_list.append([new_vnf])
            sfc.vnf_count += 1

    # print(sfc_list)
    file = open("./sfcs.txt", "w")
    # print(file)
    file.write(f"{len(sfc_list)}\n")
    for sfc in sfc_list:
        file.write(f"{sfc}")

    file = open("./vnfs.txt", "w")
    # print(file)
    file.write(f"{len(vnf_list)}\n")
    for vnf in vnf_list:
        file.write(f"{vnf}")

    # print(find_k_shortest_paths(10, sfc_list[0], node_lookup))


main()


# groundnodes.txt
# satellites.txt


# vnf1, vnf2, vnf3, vnf4


# {vnf1, vnf2} {vnf3} {vnf4, vnf5}


# ({vnf1, vnf2} {vnf3} {vnf4, vnf5})

# 10 sfcs, max 4 sets, max 3 vnfs per set


#
# number of nodes has bounds
# k such shortest paths
