# %%
import numpy as np
from typing import List
from numpy.typing import NDArray
import os
import matplotlib.pyplot as plt
 
class Mesh:
    nlocal: int
    nnodes: int
    nelem: int
    nedges: int
    nodes: NDArray[np.float64]
    elem: NDArray[np.int32]
    edges: NDArray[np.int32]
    domains: List[str]
 
    def __init__(self, path):
        self.path = path
        with open(path, "r") as f:
            self.nnodes = int(f.readline().split(" ")[3])
            self.nodes = np.zeros((self.nnodes, 2))
            for i in range(self.nnodes):
                line = f.readline()
                parts = [s.strip() for s in line.split(':')]
                self.nodes[i] = [float(val.strip())
                                 for val in parts[1].split()]
 
            self.nedges = int(f.readline().split(" ")[3])
            self.edges = np.zeros((self.nedges, 2), dtype=int)
            for i in range(self.nedges):
                line = f.readline()
                parts = [s.strip() for s in line.split(':')]
                self.edges[i] = [int(val.strip()) for val in parts[1].split()]
 
            line = f.readline()
            self.nlocal = 3 if "triangles" in line else 4
            self.nelem = int(line.split(" ")[3])
            self.elem = np.zeros((self.nelem, self.nlocal), dtype=np.int32)
            for i in range(self.nelem):
                line = f.readline()
                parts = [s.strip() for s in line.split(':')]
                self.elem[i] = [int(val.strip()) for val in parts[1].split()]
 
            self.domains = f.readlines()
 
    def __str__(self):
        info = [
            f"Mesh: {self.path}\n",
            f"├─nlocal: {self.nlocal}\n",
            f"├─nnodes: {self.nnodes}\n",
            f"├─nelem : {self.nelem}\n",
            f"├─nedges: {self.nedges}\n",
            f"├─nodes: array(shape={self.nodes.shape}, dtype=np.float64)\n",
            f"├─elem : array(shape={self.elem.shape}, dtype=np.int32)\n",
            f"└─edges: array(shape={self.edges.shape}, dtype=np.int32)"
        ]
        return "".join(info)
 
    def unfuck(self):
        visited = np.zeros(self.nnodes, dtype=bool)
        visited[self.elem] = True
        if(np.all(visited)):
            print("mesh was already good")
            return
        renumber = np.cumsum(visited) - 1
        self.elem = renumber[self.elem]
        self.edges = renumber[self.edges]
        self.nodes = self.nodes[visited,:]
        print(f"removed {self.nnodes - visited.sum()} disconnected node(s)")
        self.nnodes = len(self.nodes)        
 
    def write(self, path):
        s = []
        s.append(f"Number of nodes {self.nnodes:d} \n")
        for i in range(self.nnodes):
            s.append(f"{i:6d} : {self.nodes[i,0]:14.7e} {self.nodes[i,1]:14.7e} \n")
 
        s.append(f"Number of edges {self.nedges:d} \n")
        for i in range(self.nedges):
            s.append(f"{i:6d} : {self.edges[i,0]:6d} {self.edges[i,1]:6d} \n")
 
        if(self.nlocal == 3):
            s.append(f"Number of triangles {self.nelem:d} \n")
            for i in range(self.nelem):
                s.append(f"{i:6d} : {self.elem[i,0]:6d} {self.elem[i,1]:6d} {self.elem[i,2]:6d}\n")
        elif(self.nlocal == 4):
            s.append(f"Number of quads {self.nelem:d} \n")
            for i in range(self.nelem):
                s.append(f"{i:6d} : {self.elem[i,0]:6d} {self.elem[i,1]:6d} {self.elem[i,2]:6d} {self.elem[i,3]:6d}\n")
        else:
            print("nlocal must be 3 or 4")
            return
        
        s += self.domains
        with open(path, "w") as f:
            f.write("".join(s))
 
 
    def __repr__(self) -> str:
        return self.__str__()
 
if __name__ == "__main__":
    mesh = Mesh("../genmesh/data/mesh.txt")
    print(mesh)
    mesh.unfuck()
    print("Mesh generated and fixed, merci Michel !")
    deleted_nodes = mesh.nnodes - mesh.nodes.shape[0]
    print(f"Number of nodes deleted: {deleted_nodes}")
    mesh.write("mesh.txt")