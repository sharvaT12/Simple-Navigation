//
// References:
// TinyXML: https://github.com/leethomason/tinyxml2
// OpenStreetMap: https://www.openstreetmap.org
// OpenStreetMap docs:
//   https://wiki.openstreetmap.org/wiki/Main_Page
//   https://wiki.openstreetmap.org/wiki/Map_Features
//   https://wiki.openstreetmap.org/wiki/Node
//   https://wiki.openstreetmap.org/wiki/Way
//   https://wiki.openstreetmap.org/wiki/Relation
//
// application is where we used to navigate two people at two different building
// this is the main screen

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iomanip> /*setprecision*/
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <stack>
#include <string>
#include <vector>
#include "dist.h"
#include "graph.h"
#include "osm.h"
#include "tinyxml2.h"

using namespace std;
using namespace tinyxml2;

// to infinity and beyond redefine INF as a global variable
const double INF = numeric_limits<double>::max();

class prioritize {
 public:
  bool operator()(pair<long long, double> s1,
                  pair<long long, double> s2) const {
    // huffmanNode* will be different
    return s1.second > s2.second;
  }
};

// MS 7 function
BuildingInfo searchBuilding(string query, vector<BuildingInfo>& Buildings) {
  BuildingInfo ret;

  for (auto S : Buildings) {
    if ((S.Abbrev.find(query) != string::npos) ||
        (S.Fullname.find(query) != string::npos)) {
      ret = S;
      return ret;
    }
  }
  return ret;
}

// MS 8 function
BuildingInfo MiddleDistance(BuildingInfo& building_1, BuildingInfo& building_2,
                            vector<BuildingInfo>& Buildings) {
  BuildingInfo temp;

  // standard min algorithm
  double min = INF;
  double latitude = 0.0;
  double longtitude = 0.0;

  // nearest building MS8 midpoint coordinates
  Coordinates midpoint;

  midpoint = centerBetween2Points(building_1.Coords.Lat, building_1.Coords.Lon,
                                  building_2.Coords.Lat, building_2.Coords.Lon);

  latitude = midpoint.Lat;
  longtitude = midpoint.Lon;
  // for every distance you calculate
  // if distance < min
  // min = distance
  // return min
  for (auto S : Buildings) {
    double distance_1 = 0.0;
    double latitude_2 = S.Coords.Lat;
    double longtitude_2 = S.Coords.Lon;
    distance_1 =
        distBetween2Points(latitude, longtitude, latitude_2, longtitude_2);
    if (distance_1 < min) {  // coba tambah nanti
      min = distance_1;
      temp = S;
    }
  }
  return temp;
}

// mildstone 9
long long nearestNode(BuildingInfo& building_1, vector<FootwayInfo>& Footways,
                      map<long long, Coordinates>& Nodes) {
  double min = INF;
  long long temp;

  // for every distance you calculate
  for (auto S : Footways) {
    double length = 0.0;
    //     int size = 0;
    //     size = S.Nodes.size();
    for (unsigned i = 0; i < S.Nodes.size(); i++) {
      length =
          distBetween2Points(Nodes.at(S.Nodes[i]).Lat, Nodes.at(S.Nodes[i]).Lon,
                             building_1.Coords.Lat, building_1.Coords.Lon);
      if (length < min) {
        min = length;
        temp = Nodes.at(S.Nodes[i]).ID;
      }
    }
  }

  return temp;
}

// mildstone 10
// from lab 11
void Dijkstra(graph<long long, double>& G, long long startV, long long endV,
              map<long long, double>& distances, map<long long, long long>& M) {
  vector<long long> visited;
  set<long long> verticeVis;

  priority_queue<pair<long long, double>, vector<pair<long long, double>>,
                 prioritize>
      tmp;  // tmp
  vector<long long> vertices = G.getVertices();

  for (auto& e : vertices) {
    pair<long long, double> p1;
    p1.first = e;
    p1.second = INF;
    tmp.push(p1);
    distances.emplace(p1);
    distances[e] = INF;
    M[e] = 0;
  }

  pair<long long, double> initializeV;
  initializeV.first = startV;
  initializeV.second = 0;
  tmp.push(initializeV);
  distances[startV] = 0;

  // lab
  visited.push_back(startV);
  verticeVis.insert(startV);

  while (!tmp.empty()) {
    pair<long long, double> currentV;
    currentV.first = tmp.top().first;
    currentV.second = tmp.top().second;
    tmp.pop();

    if (currentV.second == INF) {
      break;
    }

    else if (!verticeVis.count(currentV.first)) {
      visited.push_back(currentV.first);
      verticeVis.insert(currentV.first);
    }

    else if (currentV.first == endV) {
      visited.push_back(endV);
      break;
    }

    else if (verticeVis.count(currentV.first) == 1) {
      continue;
    }

    set<long long> neighbor = G.neighbors(currentV.first);

    for (auto& e : neighbor) {
      double edgeweight = INF;

      if (G.getWeight(currentV.first, e, edgeweight)) {
        ;
      }

      double altPathDist = currentV.second + edgeweight;

      if (altPathDist < distances.at(e)) {
        distances[e] = altPathDist;
        pair<long long, double> e1;
        M[e] = currentV.first;
        e1.first = e;
        e1.second = distances[e];
        tmp.push(e1);
      }
    }
  }
}

// MS11
// to get the path from A to B, repeadtedly pop off from the stack until its
// empty
vector<long long> getPath(graph<long long, double>& G, long long startV,
                          long long endVertex,
                          map<long long, double>& distances,
                          map<long long, long long>& M) {
  // create vector

  // create the stack to push late
  stack<long long> visited;

  // vector to push back the value when is not empty
  vector<long long> verticeVis;

  // currV = endVertex
  long long currV = endVertex;
  long long temp = visited.top();

  // while currV is not null{
  // push currV to stack
  // curV = currV->predV
  // }
  // while stack is not empty{
  // currV = pop from stack
  // push back currV to vector
  // }
  while (currV != 0) {
    visited.push(currV);
    currV = M[currV];
  }

  while (!visited.empty()) {
    verticeVis.push_back(temp);
    visited.pop();
  }

  return verticeVis;
}

//
// Implement your creative component application here
// TO DO: add arguments
//
void creative(map<long long, Coordinates>& Nodes, vector<FootwayInfo>& Footways,
              vector<BuildingInfo>& Buildings, graph<long long, double>& G) {
  string person1Building, person2Building;

  cout << endl;
  cout << "Enter person 1's building (partial name or abbreviation), or #> ";
  getline(cin, person1Building);

  while (person1Building != "#") {
    cout << "Enter person 2's building (partial name or abbreviation)> ";
    getline(cin, person2Building);

    // MILDSTONE 7
    BuildingInfo ms7 = searchBuilding(person1Building, Buildings);
    BuildingInfo ms7_1 = searchBuilding(person2Building, Buildings);

    // MILDSTONE 8
    BuildingInfo ms8 = MiddleDistance(ms7, ms7_1, Buildings);

    // mildstone 9
    // FOR first person
    long long ms9 = nearestNode(ms7, Footways, Nodes);
    // for second person
    long long ms9_1 = nearestNode(ms7, Footways, Nodes);
    // for the total destination
    long long ms9_2 = nearestNode(ms7, Footways, Nodes);

    //
    // TO DO: lookup buildings, find nearest start and dest nodes, find center
    // run Dijkstra's alg from each start, output distances and paths to
    // destination:
    //

    vector<long long> ms11;
    vector<long long> ms11_1;
    vector<long long> ms11_2;

    map<long long, long long> M;
    map<long long, long long> M_1;
    map<long long, long long> M_2;

    map<long long, double> distances;
    map<long long, double> distances_1;
    map<long long, double> distances_2;

    // MS 10
    // for the first person
    Dijkstra(G, ms9, ms9_2, distances, M);
    ms11 = getPath(G, ms9, ms9_2, distances, M);

    // for the second person
    Dijkstra(G, ms9_1, ms9_2, distances_1, M_1);
    ms11 = getPath(G, ms9_1, ms9_2, distances_1, M_1);

    // for the total distance
    Dijkstra(G, ms9, ms9_1, distances_2, M_2);
    ms11 = getPath(G, ms9, ms9_1, distances_2, M_2);

    if (ms7.Coords.Lon == 0.0 && ms7.Coords.Lat == 0.0) {
      cout << "Person 1's building not found" << endl;
    }

    else if (ms7_1.Coords.Lon == 0.0 && ms7_1.Coords.Lat == 0.0) {
      cout << "Person 2's building not found" << endl;
    }

    else {
      if (distances_2[ms9_1] >= INF) {
        cout << "Person 1's point:" << endl;
        cout << " " << ms7.Fullname << endl;
        cout << " (" << ms7.Coords.Lat << ", " << ms7.Coords.Lon << ")" << endl;

        cout << "Person 2's point:" << endl;
        cout << " " << ms7_1.Fullname << endl;
        cout << " (" << ms7_1.Coords.Lat << ", " << ms7_1.Coords.Lon << ")"
             << endl;

        cout << "Destination Building:" << endl;
        cout << " " << ms8.Fullname << endl;
        cout << " (" << ms8.Coords.Lat << ", " << ms8.Coords.Lon << ")" << endl;
        cout << endl;

        cout << "Nearest P1 node:" << endl;
        cout << " " << ms9 << endl;
        cout << " (" << Nodes[ms9].Lat << ", " << Nodes[ms9].Lon << ")" << endl;

        cout << "Nearest P2 node:" << endl;
        cout << " " << ms9_1 << endl;
        cout << " (" << Nodes[ms9_1].Lat << ", " << Nodes[ms9_1].Lon << ")"
             << endl;

        cout << "Nearest destination node:" << endl;
        cout << ms9_2 << endl;
        cout << " (" << Nodes[ms9_2].Lat << ", " << Nodes[ms9_2].Lon << ")"
             << endl;
        cout << endl;

        cout << endl;
        cout << "At least one person was unable to reach the destination "
                "building. Finding next closest building..."
             << endl;
        cout << endl;

        cout << "Sorry, destination unreachable." << endl;
        cout << endl;
      }
    }

    cout << endl;
    cout << "Enter person 1's building (partial name or abbreviation), or #> ";
    getline(cin, person1Building);
  }
}

//
// Implement your standard application here
// TO DO: add a parameter for the graph you make.
//
void application(map<long long, Coordinates>& Nodes,
                 vector<FootwayInfo>& Footways, vector<BuildingInfo>& Buildings,
                 graph<long long, double>& G) {
  string person1Building, person2Building;

  cout << endl;
  cout << "Enter person 1's building (partial name or abbreviation), or #> ";
  getline(cin, person1Building);

  while (person1Building != "#") {
    cout << "Enter person 2's building (partial name or abbreviation)> ";
    getline(cin, person2Building);

    // MILDSTONE 7
    BuildingInfo ms7 = searchBuilding(person1Building, Buildings);
    BuildingInfo ms7_1 = searchBuilding(person2Building, Buildings);

    // MILDSTONE 8
    BuildingInfo ms8 = MiddleDistance(ms7, ms7_1, Buildings);

    // mildstone 9
    // FOR first person
    long long ms9 = nearestNode(ms7, Footways, Nodes);
    // for second person
    long long ms9_1 = nearestNode(ms7, Footways, Nodes);
    // for the total destination
    long long ms9_2 = nearestNode(ms7, Footways, Nodes);

    //
    // TO DO: lookup buildings, find nearest start and dest nodes, find center
    // run Dijkstra's alg from each start, output distances and paths to
    // destination:
    //

    vector<long long> ms11;
    vector<long long> ms11_1;
    vector<long long> ms11_2;

    map<long long, long long> M;
    map<long long, long long> M_1;
    map<long long, long long> M_2;

    map<long long, double> distances;
    map<long long, double> distances_1;
    map<long long, double> distances_2;

    // MS 10
    // for the first person
    Dijkstra(G, ms9, ms9_2, distances, M);
    ms11 = getPath(G, ms9, ms9_2, distances, M);

    // for the second person
    Dijkstra(G, ms9_1, ms9_2, distances_1, M_1);
    ms11 = getPath(G, ms9_1, ms9_2, distances_1, M_1);

    // for the total distance
    Dijkstra(G, ms9, ms9_1, distances_2, M_2);
    ms11 = getPath(G, ms9, ms9_1, distances_2, M_2);

    if (ms7.Coords.Lon == 0.0 && ms7.Coords.Lat == 0.0) {
      cout << "Person 1's building not found" << endl;
    }

    else if (ms7_1.Coords.Lon == 0.0 && ms7_1.Coords.Lat == 0.0) {
      cout << "Person 2's building not found" << endl;
    }

    else {
      if (distances_2[ms9_1] >= INF) {
        cout << "Person 1's point:" << endl;
        cout << " " << ms7.Fullname << endl;
        cout << " (" << ms7.Coords.Lat << ", " << ms7.Coords.Lon << ")" << endl;

        cout << "Person 2's point:" << endl;
        cout << " " << ms7_1.Fullname << endl;
        cout << " (" << ms7_1.Coords.Lat << ", " << ms7_1.Coords.Lon << ")"
             << endl;

        cout << "Destination Building:" << endl;
        cout << " " << ms8.Fullname << endl;
        cout << " (" << ms8.Coords.Lat << ", " << ms8.Coords.Lon << ")" << endl;
        cout << endl;

        cout << "Nearest P1 node:" << endl;
        cout << " " << ms9 << endl;
        cout << " (" << Nodes[ms9].Lat << ", " << Nodes[ms9].Lon << ")" << endl;

        cout << "Nearest P2 node:" << endl;
        cout << " " << ms9_1 << endl;
        cout << " (" << Nodes[ms9_1].Lat << ", " << Nodes[ms9_1].Lon << ")"
             << endl;

        cout << "Nearest destination node:" << endl;
        cout << ms9_2 << endl;
        cout << " (" << Nodes[ms9_2].Lat << ", " << Nodes[ms9_2].Lon << ")"
             << endl;
        cout << endl;

        cout << endl;
        cout << "At least one person was unable to reach the destination "
                "building. Finding next closest building..."
             << endl;
        cout << endl;

        cout << "Sorry, destination unreachable." << endl;
        cout << endl;
      }
      //       else {
      //         cout << "Person 1's point:" << endl;
      //         cout << " " << ms7.Fullname << endl;
      //         cout << " (" << ms7.Coords.Lat << ", " << ms7.Coords.Lon << ")"
      //         << endl;

      //         cout << "Person 2's point:" << endl;
      //         cout << " " << ms7_1.Fullname << endl;
      //         cout << " (" << ms7_1.Coords.Lat << ", " << ms7_1.Coords.Lon <<
      //         ")"
      //              << endl;

      //         cout << "Destination Building:" << endl;
      //         cout << " " << ms8.Fullname << endl;
      //         cout << " (" << ms8.Coords.Lat << ", " << ms8.Coords.Lon << ")"
      //         << endl; cout << endl;

      //         cout << "Nearest P1 node:" << endl;
      //         cout << " " << ms9 << endl;
      //         cout << " (" << Nodes[ms9].Lat << ", " << Nodes[ms9].Lon << ")"
      //         << endl;

      //         cout << "Nearest P2 node:" << endl;
      //         cout << " " << ms9_1 << endl;
      //         cout << " (" << Nodes[ms9_1].Lat << ", " << Nodes[ms9_1].Lon <<
      //         ")"
      //              << endl;

      //         cout << "Nearest destination node:" << endl;
      //         cout << ms9_2 << endl;
      //         cout << " (" << Nodes[ms9_2].Lat << ", " << Nodes[ms9_2].Lon <<
      //         ")"
      //              << endl;
      //         cout << endl;

      //         cout << endl;
      //         cout << "At least one person was unable to reach the
      //         destination "
      //                 "building. Finding next closest building..."
      //              << endl;
      //         cout << endl;

      //         cout << endl;
      //         cout << "Person 1's distance to dest: " << distances.at(ms9_2)
      //              << " miles" << endl;

      //         int k = 0;
      //         int size_1 = ms11.size();

      //         cout << "Path: ";

      //         for (auto a : ms11) {
      //           if (k == size_1 - 1) {
      //             break;
      //           }
      //           cout << a << "->";
      //           k++;
      //         }

      //         cout << ms9_2 << endl;

      //         cout << "Person 2's distance to dest: " <<
      //         distances_1.at(ms9_2)
      //              << " miles" << endl;

      //         int L = 0;
      //         int size_2 = ms11_1.size();

      //         cout << "Path: ";

      //         for (auto b : ms11_1) {
      //           if (L == size_2 - 1) {
      //             break;
      //           }
      //           cout << b << "->";
      //           L++;
      //         }

      //         cout << ms9_2 << endl;
      //       }
    }

    //
    // another navigation?
    //
    cout << endl;
    cout << "Enter person 1's building (partial name or abbreviation), or #> ";
    getline(cin, person1Building);
  }
}

int main() {
  // maps a Node ID to it's coordinates (lat, lon)
  map<long long, Coordinates> Nodes;
  // info about each footway, in no particular order
  vector<FootwayInfo> Footways;
  // info about each building, in no particular order
  vector<BuildingInfo> Buildings;
  // declare graph variable
  graph<long long, double> G;
  XMLDocument xmldoc;

  cout << "** Navigating UIC open street map **" << endl;
  cout << endl;
  cout << std::setprecision(8);

  string def_filename = "map.osm";
  string filename;

  cout << "Enter map filename> ";
  getline(cin, filename);

  if (filename == "") {
    filename = def_filename;
  }

  //
  // Load XML-based map file
  //
  if (!LoadOpenStreetMap(filename, xmldoc)) {
    cout << "**Error: unable to load open street map." << endl;
    cout << endl;
    return 0;
  }

  //
  // Read the nodes, which are the various known positions on the map:
  //
  int nodeCount = ReadMapNodes(xmldoc, Nodes);

  //
  // Read the footways, which are the walking paths:
  //
  int footwayCount = ReadFootways(xmldoc, Footways);

  //
  // Read the university buildings:
  //
  int buildingCount = ReadUniversityBuildings(xmldoc, Nodes, Buildings);

  //
  // Stats
  //
  assert(nodeCount == (int)Nodes.size());
  assert(footwayCount == (int)Footways.size());
  assert(buildingCount == (int)Buildings.size());

  cout << endl;
  cout << "# of nodes: " << Nodes.size() << endl;
  cout << "# of footways: " << Footways.size() << endl;
  cout << "# of buildings: " << Buildings.size() << endl;

  // Todo MS 5 : add verticies
  // loop through this map (for each)
  for (auto a : Nodes) {
    G.addVertex(a.first);
  }

  // Todo MS 6: add edges
  for (auto a : Footways) {
    for (int i = 0; i < a.Nodes.size() - 1; i++) {
      double distance;
      distance = distBetween2Points(
          Nodes.at(a.Nodes[i]).Lat, Nodes.at(a.Nodes[i]).Lon,
          Nodes.at(a.Nodes[i + 1]).Lat, Nodes.at(a.Nodes[i + 1]).Lon);
      G.addEdge(a.Nodes[i], a.Nodes[i + 1], distance);
      G.addEdge(a.Nodes[i + 1], a.Nodes[i], distance);
    }
  }

  //
  // TO DO: build the graph, output stats:
  //

  // todo uncomment bellow after ms 6
  cout << "# of vertices: " << G.NumVertices() << endl;
  cout << "# of edges: " << G.NumEdges() << endl;
  cout << endl;

  //
  // Menu
  //
  string userInput;
  cout << "Enter \"a\" for the standard application or "
       << "\"c\" for the creative component application> ";
  getline(cin, userInput);
  if (userInput == "a") {
    // TO DO: add argument for the graph you make.
    application(Nodes, Footways, Buildings, G);
  } else if (userInput == "c") {
    // TO DO: add arguments
    creative(Nodes, Footways, Buildings, G);
  }
  //
  // done:
  //
  cout << "** Done **" << endl;
  return 0;
}
