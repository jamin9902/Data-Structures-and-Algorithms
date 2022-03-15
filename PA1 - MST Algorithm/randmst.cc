#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <algorithm>
#include <math.h>  
#include <chrono>

using namespace std;

struct Vertex {

    float x;
    float y;
    float z;
    float a;
    Vertex * parent;
    int rank;

    Vertex(float x = 0.0, float y = 0.0, float z = 0.0, float a = 0.0)
    {
        this->x = x;
        this->y = y;
        this->z = z;
        this->a = a;
    }
};

struct Edge {

    float weight;
    std::pair<Vertex *, Vertex *> vertices;

    Edge(float w, Vertex * u, Vertex * v) {
        weight = w;
        vertices = std::make_pair(u, v);
    }
};


struct Graph
{

    int dim, n;
    std::vector<Edge> edges;
    std::vector<Vertex> vertices;
    long totalEdges;
    long edgesPruned;

    void addEdge(Edge e) { edges.push_back(e); }

    void addVertex(Vertex v) { vertices.push_back(v); }

    float calcEuclidean(Vertex * v1, Vertex * v2)
    {
        return sqrt(pow(v1->x - v2->x,2) + pow(v1->y - v2->y,2) + pow(v1->z - v2->z,2) + pow(v1->a - v2->a,2));
    }

    Graph(int dim, int n)
    {   
        this->dim = dim;
        this->n = n;
        totalEdges = 0;
        edgesPruned = 0;
        for (int i = 0; i < n; i++)
        {  
            float x = (float) rand()/RAND_MAX;
            float y = (float) rand()/RAND_MAX;
            float z = (float) rand()/RAND_MAX;
            float a = (float) rand()/RAND_MAX;
            if(dim == 0){
                addVertex(Vertex());
            }
            if(dim == 2){
                addVertex(Vertex(x,y));
            }
            if(dim == 3){
                addVertex(Vertex(x,y,z));
            }
            if(dim == 4){
                addVertex(Vertex(x,y,z,a));
            }
        }
        for (int i = 0; i < n; i++)
        {
            for (int j = i + 1; j < n; j++)
            {   
                totalEdges += 1;
                if(dim == 0){
                    float w = (float) rand()/RAND_MAX;
                    if(w < 4.3*pow(n, -0.82)){
                        addEdge(Edge(w, &vertices[i], &vertices[j]));
                    } else edgesPruned += 1;
                }
                if(dim == 2){
                    float w = calcEuclidean(&vertices[i], &vertices[j]);
                    if(w < 2.8*pow(n, -0.49)){
                        addEdge(Edge(w, &vertices[i], &vertices[j]));
                    } else edgesPruned += 1;
                }
                if(dim == 3){
                    float w = calcEuclidean(&vertices[i], &vertices[j]);
                    if(w < 1.65*pow(n, -0.3)){  
                        addEdge(Edge(w, &vertices[i], &vertices[j]));
                    } else edgesPruned += 1;
                }
                if(dim == 4){
                    float w = calcEuclidean(&vertices[i], &vertices[j]);
                    if(w < 1.62*pow(n, -0.23)){
                        addEdge(Edge(w, &vertices[i], &vertices[j]));
                    } else edgesPruned += 1;
                }
            }
        }
    }

    void makeset(Vertex * x){
        x->parent = x;
        x->rank = 0;
    }

    Vertex * find(Vertex * x){
        if(x != x->parent){
            x->parent = find(x->parent);
        }
        return x->parent;
    }

    Vertex * link(Vertex * root1, Vertex * root2){
        if(root1->rank == root2->rank) {
            root1->parent = root2;
            root2->rank += 1;
            return root2;
        } 
        else if(root1->rank > root2->rank){
            root2->parent = root1;
            return root1;
        }
        else {
            root1->parent = root2;
            return root2;
        }
    }

    void union_procedure(Vertex * x, Vertex * y){
        link(find(x), find(y));
    }


    vector<Edge> MSTkruskal();
    void rerun_kruskal();
};

bool edges_sorter(Edge u, Edge v) { return (u.weight < v.weight); }

void Graph::rerun_kruskal(){
    totalEdges = 0;
    edgesPruned = 0;
    edges.clear(); 
    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {   
            totalEdges += 1;
            if(dim == 0){
                float w = (float) rand()/RAND_MAX;
                addEdge(Edge(w, &vertices[i], &vertices[j]));
            }
            if(dim == 2){
                float w = calcEuclidean(&vertices[i], &vertices[j]);
                addEdge(Edge(w, &vertices[i], &vertices[j]));
            }
            if(dim == 3){
                float w = calcEuclidean(&vertices[i], &vertices[j]);
                addEdge(Edge(w, &vertices[i], &vertices[j]));
            }
            if(dim == 4){
                float w = calcEuclidean(&vertices[i], &vertices[j]);
                addEdge(Edge(w, &vertices[i], &vertices[j]));
            }
        }
    }
    MSTkruskal();
    return;
}

vector<Edge> Graph::MSTkruskal(){
    vector<Edge> X;
    std::sort(edges.begin(), edges.end(), edges_sorter);   


    for (int v = 0; v != vertices.size(); v++){
        makeset(&vertices[v]);
    }

    for (int e = 0; e != edges.size(); e++){
        if(find(edges[e].vertices.first) != find(edges[e].vertices.second)){
            X.push_back(edges[e]);
            union_procedure(edges[e].vertices.first, edges[e].vertices.second);
        }
    }
    
    if (X.size() != (vertices.size() - 1)){
        rerun_kruskal();
    }
    return X;
}


std::vector<float> find_max_edge_weights(long dimension, long numtrials){
    std::vector<float> max_weights;
    for (int i=5; i<2045; i+= 5){

        float max_edge_weight = 0;
        for(int j=0; j<numtrials; j++){
            Graph G(dimension, i);
            vector<Edge> X = G.MSTkruskal();
            
            if(max_edge_weight < X[X.size()-1].weight){
                max_edge_weight = X[X.size()-1].weight;
            }
        }
        max_weights.push_back(max_edge_weight);
        // cout << max_edge_weight << "," << endl;
    }
    return max_weights;
}

float standard_run(long dimension, long numtrials, long numpoints){
    
    
    float sum_results = 0;
    float sum_edges_pruned = 0;
    float sum_edges = 0;

    for(int i=0; i<numtrials; i++){

        Graph G(dimension, numpoints);  
        vector<Edge> X = G.MSTkruskal();
        float weight_sum = 0;

        for (int e = 0; e != X.size(); e++){
            weight_sum += X[e].weight;
        }
        sum_results += weight_sum;
        sum_edges_pruned += G.edgesPruned;
        sum_edges += G.totalEdges;
    }
    // cout << "EDGES total: " << sum_edges / numtrials  << " EDGES pruned: " << sum_edges_pruned / numtrials << endl;
    return sum_results / numtrials;
}

void run_all_n(long numtrials){
    vector<int> dimensions_arr = {2};
    vector<int> test_n = {128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144};


    for (int d=0; d != dimensions_arr.size(); d++){
        for (int i=0; i != test_n.size(); i++){
            auto start = std::chrono::high_resolution_clock::now();
            float avg_results = standard_run(dimensions_arr[d], numtrials, test_n[i]);
            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
            cout << "dimensions: " << dimensions_arr[d] << " numtrials: " << numtrials <<" numpoints: " << test_n[i] << " results: " << avg_results << " microseconds total: " << duration.count() << " average microsecond per trial: " << duration.count()/numtrials << endl;
        }
    }
}


int main(int argc, char * argv[]){

    srand((unsigned)time(NULL));
    long flag = strtol(argv[1], NULL, 10);
    long numpoints = strtol(argv[2], NULL, 10);
    long numtrials = strtol(argv[3], NULL, 10);
    long dimension = strtol(argv[4], NULL, 10);


    // std::vector<float> weights = find_max_edge_weights(dimension, numtrials);
    // run_all_n(numtrials);

    float result = standard_run(dimension, numtrials, numpoints);

    cout << result << " " << numpoints << " " << numtrials << " " << dimension << endl;
    
    return 0;
}

