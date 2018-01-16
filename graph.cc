#include "graph.hh"

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <string>

static void fail(std::string msg) {
    std::cerr << msg << std::endl;
    exit(1);
}

Graph::Graph(unsigned int n) {
    this->n = n;
    adjmat = {n, std::vector<unsigned int>(n, false)};
}

Graph induced_subgraph(struct Graph& g, std::vector<int> vv) {
    Graph subg(vv.size());
    for (int i=0; i<subg.n; i++)
        for (int j=0; j<subg.n; j++)
            subg.adjmat[i][j] = g.adjmat[vv[i]][vv[j]];
    return subg;
}

void add_edge(Graph& g, int v, int w) {
    if (v != w) {
        g.adjmat[v][w] = 1;
        g.adjmat[w][v] = 1;
    } else {
        fail("Loop detected!");
    }
}

struct Graph readDimacsGraph(char* filename) {
    struct Graph g(0);

    FILE* f;
    
    if ((f=fopen(filename, "r"))==NULL)
        fail("Cannot open file");

    char* line = NULL;
    size_t nchar = 0;

    int nvertices = 0;
    int medges = 0;
    int v, w;
    int edges_read = 0;

    while (getline(&line, &nchar, f) != -1) {
        if (nchar > 0) {
            switch (line[0]) {
            case 'p':
                if (sscanf(line, "p edge %d %d", &nvertices, &medges)!=2)
                    fail("Error reading a line beginning with p.\n");
                g = Graph(nvertices);
                break;
            case 'e':
                if (sscanf(line, "e %d %d", &v, &w)!=2)
                    fail("Error reading a line beginning with e.\n");
                add_edge(g, v-1, w-1);
                edges_read++;
                break;
            }
        }
    }

    if (medges>0 && edges_read != medges) fail("Unexpected number of edges.");

    fclose(f);
    return g;
}

struct Graph readLadGraph(char* filename) {
    struct Graph g(0);
    FILE* f;
    
    if ((f=fopen(filename, "r"))==NULL)
        fail("Cannot open file");

    int nvertices = 0;
    int w;

    if (fscanf(f, "%d", &nvertices) != 1)
        fail("Number of vertices not read correctly.\n");
    g = Graph(nvertices);

    for (int i=0; i<nvertices; i++) {
        int edge_count;
        if (fscanf(f, "%d", &edge_count) != 1)
            fail("Number of edges not read correctly.\n");
        for (int j=0; j<edge_count; j++) {
            if (fscanf(f, "%d", &w) != 1)
                fail("An edge was not read correctly.\n");
            add_edge(g, i, w);
        }
    }

    fclose(f);
    return g;
}

int read_word(FILE *fp) {
    unsigned char a[2];
    if (fread(a, 1, 2, fp) != 2)
        fail("Error reading file.\n");
    return (int)a[0] | (((int)a[1]) << 8);
}

struct Graph readBinaryGraph(char* filename)
{
    struct Graph g(0);
    FILE* f;
    
    if ((f=fopen(filename, "rb"))==NULL)
        fail("Cannot open file");

    int nvertices = read_word(f);
    g = Graph(nvertices);

    for (int i=0; i<nvertices; i++) {
        read_word(f);  // ignore label
    }

    for (int i=0; i<nvertices; i++) {
        int len = read_word(f);
        for (int j=0; j<len; j++) {
            int target = read_word(f);
            read_word(f);  // ignore label
            add_edge(g, i, target);
        }
    }
    fclose(f);
    return g;
}

struct Graph readGraph(char* filename, char format) {
    struct Graph g(0);
    if (format=='D') g = readDimacsGraph(filename);
    else if (format=='L') g = readLadGraph(filename);
    else if (format=='B') g = readBinaryGraph(filename);
    else fail("Unknown graph format\n");
    return g;
}

auto calculate_degrees(const Graph & g) -> std::vector<int>
{
    std::vector<int> degree(g.n, 0);
    for (int v=0; v<g.n; v++)
        for (int w=0; w<g.n; w++)
            if (g.adjmat[v][w]) degree[v]++;
    return degree;
}

