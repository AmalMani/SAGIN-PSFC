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


class Link:
    def __init__(
        self, node1, node2, weight=random.randint(MIN_EDGE_WEIGHT, MAX_EDGE_WEIGHT)
    ) -> None:
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
    def __init__(self, node_id):
        super().__init__(node_id)

        self.node_type = random.choice([1, 2])
        self.height = None
        self.communication_range = 0
        if self.node_type == 1:
            self.height = random.randint(MIN_LEO_HEIGHT, MAX_LEO_HEIGHT)
            self.communication_range = self.height + random.randint(100, 500)
        else:
            self.height = random.randint(MIN_HEO_HEIGHT, MAX_HEO_HEIGHT)
            self.communication_range = self.height + random.randint(1000, 10000)

        a, b, c = random.randint(1, 10), random.randint(1, 10), random.randint(1, 10)
        self.normal_vector = (
            a,
            b,
            c,
        )  # Normal vector to the plane, plane passes through origin

        while True:
            x = random.randint(-self.height // 4, self.height // 4)
            k = -1 * a * x
            m = (EARTH_RADIUS + self.height) ** 2 - x**2
            check = 4 * (b**2) * (k**2) - 4 * (b**2 + c**2) * (k**2 - (c**2) * m)
            if check > 0:
                y = (2 * b * k + math.sqrt(check)) / (2 * (b**2 + c**2))
                z = (k - b * y) / c
                self.starting_coordinates = (x, y, z)
                break

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
        self.lifetime = random.randint(500, 1000)

    def __repr__(self) -> str:
        repr_str = f"{self.vnf_count}\n{self.source_node} {self.lifetime} {self.bandwidth_requirement} {self.max_tolerated_delay}\n"
        for index, vnf_list in enumerate(self.vnf_list):
            for vnf in vnf_list:
                repr_str += f"{vnf.vnf_id} {index + 1}\n"
        return repr_str


def main():
    number_of_clusters = random.randint(MIN_NUMBER_OF_CLUSTERS, MAX_NUMBER_OF_CLUSTERS)
    cluster_list = dict()
    link_list = []
    node_id = 0
    node_lookup = dict()
    for i in range(number_of_clusters):
        x = random.randint(
            -EARTH_RADIUS + CLUSTER_RADIUS, EARTH_RADIUS - CLUSTER_RADIUS
        )
        y = random.randint(
            -EARTH_RADIUS + CLUSTER_RADIUS + abs(x),
            EARTH_RADIUS - CLUSTER_RADIUS - abs(x),
        )
        z = math.sqrt((EARTH_RADIUS**2 - x**2 - y**2))
        if (x, y, z) in cluster_list:
            continue
        curr_cluster = []

        for j in range(random.randint(MIN_NUMBER_OF_CLUSTERS, MAX_NUMBER_OF_CLUSTERS)):
            x1 = x + random.randint(-CLUSTER_RADIUS // 2, CLUSTER_RADIUS // 2)
            y1 = y + random.randint(-CLUSTER_RADIUS // 3, CLUSTER_RADIUS // 3)
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
                    edge = Link(neighbor.node_id, new_node.node_id)
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
        edge = Link(node_id1, node_id2)
        node_lookup[node_id1].edges[node_lookup[node_id2].node_id] = edge
        node_lookup[node_id2].edges[node_lookup[node_id1].node_id] = edge
        link_list.append(edge)

    satellites = []

    number_of_satellites = random.randint(200, 300)
    for _ in range(number_of_satellites):
        new_satellite = Satellite(node_id)
        node_id += 1
        for groundnode_id in range(number_of_nodes):
            edge = Link(groundnode_id, new_satellite.node_id, 0)
            node_lookup[groundnode_id].edges[new_satellite.node_id] = edge
            new_satellite.edges[groundnode_id] = edge
            link_list.append(edge)

        satellites.append(new_satellite)

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


main()


# groundnodes.txt
# satellites.txt


# vnf1, vnf2, vnf3, vnf4


# {vnf1, vnf2} {vnf3} {vnf4, vnf5}


# ({vnf1, vnf2} {vnf3} {vnf4, vnf5})

# 10 sfcs, max 4 sets, max 3 vnfs per set
