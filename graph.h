// graph.h <Starter Code>
// < Wilbert Limson >
//
// Basic graph class using adjacency matrix representation.  Currently
// limited to a graph with at most 100 vertices.
//
// University of Illinois at Chicago
// CS 251: Fall 2021
// Project #7 - Openstreet Maps
//
// is filled with 100 vertices  which the limitation is due to
// the underlying implementation based on an adjacency matrix

#pragma once

#include <algorithm>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace std;

template <typename VertexT, typename WeightT>
class graph {
 private:
  // vector of pair
  typedef pair<VertexT, WeightT> WeigPair;
  typedef map<VertexT, vector<WeigPair>> VecMap;
  VecMap pos;
  //     int total;

 public:
  //
  // constructor:
  //
  // Constructs an empty graph where n is the max # of vertices
  // you expect the graph to contain.
  //
  // NOTE: the graph is implemented using an adjacency matrix.
  // If n exceeds the dimensions of this matrix, an exception
  // will be thrown to let you know that this implementation
  // will not suffice.
  //
  graph() {
    // default constructor
    //     total = 0;
    return;
  }

  //
  // NumVertices
  //
  // Returns the # of vertices currently in the graph.
  //
  int NumVertices() const {
    int num = pos.size();
    return num;
  }

  //
  // NumEdges
  //
  // Returns the # of edges currently in the graph.
  //
  int NumEdges() const {
    // return the number of total when it's found
    int total_1 = 0;
    auto temp = pos.begin();
    while (temp != pos.end()) {
      vector<WeigPair> S = temp->second;
      total_1 += S.size();
      temp++;
    }
    return total_1;
  }

  //
  // addVertex
  //
  // Adds the vertex v to the graph if there's room, and if so
  // returns true.  If the graph is full, or the vertex already
  // exists in the graph, then false is returned.
  //

  // to check whether the vertex has reached the last position of the  map
  int LookupVertex(VertexT v) {
    auto a = pos.find(v);
    if (a == pos.end()) {
      return 1;
    } else if (a != pos.end()) {
      return 0;
    } else {
      return 0;
    }
  }

  bool addVertex(VertexT v) {
    //
    // is the vertex already in the graph?  If so, we do not
    // insert again otherwise Vertices may fill with duplicates:
    //
    auto temp = pos.find(v);
    if (temp == pos.end()) {
      // if it doesn't match at the end it will return false
      vector<WeigPair> temp;
      pos.emplace(v, temp);
      return true;
    } else if (temp != pos.end()) {
      // make a temporary position of the vector to store it
      return false;
    }

    return false;
  }

  graph(const graph& other) {
    this->pos = other.pos;
    //     this->total = other.total;
  }

  graph operator=(const graph& other) {
    this->pos = other.pos;
    //     this->total = other.total;
    return *this;
  }

  //
  // addEdge
  //
  // Adds the edge (from, to, weight) to the graph, and returns
  // true.  If the vertices do not exist or for some reason the
  // graph is full, false is returned.
  //
  // NOTE: if the edge already exists, the existing edge weight
  // is overwritten with the new edge weight.
  //
  bool addEdge(VertexT from, VertexT to, WeightT weight) {
    //
    // we need to search the Vertices and find the position
    // of each vertex; this will denote the row and col to
    // access in the adjacency matrix:
    //
    WeigPair temp_pair;
    int row = LookupVertex(from);

    if (row == 1) {  // not found:
      return false;
    }

    int col = LookupVertex(to);

    if (col == 1) {  // not found:
      return false;
    }

    auto current_edge = pos.find(from);
    int temp = 0;
    temp = current_edge->second.size();
    int i = 0;

    while (i < temp) {
      if (current_edge->second.at(i).first == to) {
        current_edge->second.at(i).second = weight;
        //         this->total++;
        return true;
      }
      i++;
    }

    temp_pair = make_pair(to, weight);
    current_edge->second.push_back(temp_pair);

    //     this->total++;
    return true;
  }

  //
  // getWeight
  //
  // Returns the weight associated with a given edge.  If
  // the edge exists, the weight is returned via the reference
  // parameter and true is returned.  If the edge does not
  // exist, the weight parameter is unchanged and false is
  // returned.
  //
  bool getWeight(VertexT from, VertexT to, WeightT& weight) const {
    //
    // we need to search the Vertices and find the position
    // of each vertex; this will denote the row and col to
    // access in the adjacency matrix:
    //

    WeigPair temp_pair;

    auto row = pos.find(from);

    if (row == pos.end()) {  // not found:
      return false;
    }

    auto col = pos.find(to);

    if (col == pos.end()) {  // not found:
      return false;
    }

    auto current_edge = pos.find(from);
    int temp = 0;
    temp = current_edge->second.size();
    int i = 0;

    while (i < temp) {
      if (current_edge->second.at(i).first == to) {
        weight = current_edge->second.at(i).second;
        return true;
      }
      i++;
    }

    return false;
  }

  //
  // neighbors
  //
  // Returns a set containing the neighbors of v, i.e. all
  // vertices that can be reached from v along one edge.
  // Since a set is returned, the neighbors are returned in
  // sorted order; use foreach to iterate through the set.
  //
  set<VertexT> neighbors(VertexT v) const {
    set<VertexT> S;

    //
    // we need to search the Vertices and find the position
    // of v, that will be the row we need in the adjacency
    // matrix:
    //
    auto row = pos.find(v);

    if (row == pos.end()) {  // not found:
      return S;
    }

    //
    // we found the row, so loop along the row and for every
    // edge that exists, add the column vertex to V:
    //
    // NOTE: how many columns are there?  The # of vertices.
    //
    auto temp = pos.find(v);
    int tot = temp->second.size();
    for (int col = 0; col < tot; ++col) {
      VertexT dest = temp->second.at(col).first;  // dest vertex is here:
      S.insert(dest);
    }

    return S;
  }

  //
  // getVertices
  //
  // Returns a vector containing all the vertices currently in
  // the graph.
  //
  vector<VertexT> getVertices() const {
    vector<VertexT> S;

    auto temp = pos.begin();
    while (temp != pos.end()) {
      S.push_back(temp->first);
      temp++;
    }
    return S;
  }

  //
  // dump
  //
  // Dumps the internal state of the graph for debugging purposes.
  //
  // Example:
  //    graph<string,int>  G(26);
  //    ...
  //    G.dump(cout);  // dump to console
  //
  void dump(ostream& output) const {
    //     output << "***************************************************" <<
    //     endl; output << "********************* GRAPH ***********************"
    //     << endl;

    //     output << "**Num vertices: " << this->NumVertices() << endl;
    //     output << "**Num edges: " << this->NumEdges() << endl;

    //     output << endl;
    //     output << "**Vertices:" << endl;

    //     for (int i = 0; i < this->NumVertices(); ++i) {
    //       output << " " << i << ". " << this->Vertices[i] << endl;
    //     }

    //     output << endl;
    //     output << "**Edges:" << endl;
    //     for (int row = 0; row < this->NumVertices(); ++row) {
    //       output << " row " << row << ": ";
    //       output << endl;
    //     }
    //     output << "**************************************************" <<
    //     endl;
  }
};
