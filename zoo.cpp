// IDENTIFIER:3E33912F8BAA7542FC4A1585D2DB6FE0312725B9
#include "zoo.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <sys/types.h>
#include <vector>

void Zoo::getMode(int argc, char *argv[]) {
  // These are used with getopt_long()
  opterr = false; // Let us handle all error output for command line options
  int choice;
  int index = 0;
  option long_options[] = {
      // TODO: Fill in two lines, for the "mode" ('m') and
      // the "help" ('h') options.
      {"help", no_argument, nullptr, 'h'},
      {"mode", required_argument, nullptr, 'm'},
      {nullptr, 0, nullptr, '\0'},
  }; // long_options[]
  while ((choice = getopt_long(argc, argv, "hm:", long_options, &index)) !=
         -1) {
    switch (choice) {
    case 'h':
      printHelp(argv);
      exit(0);

    case 'm': { // Need a block here to declare a variable inside a case
      std::string arg{optarg};
      // Checks if it is a valid mode
      if (arg != "MST" && arg != "FASTTSP" && arg != "OPTTSP") {
        std::cerr << "Not a correct option\n";
        exit(1);
      }
      mode = arg;
      break;
    } // case 'm'
    default:
      std::cerr << "Error: invalid option" << std::endl;
      exit(1);
    } // switch ..choice
  }   // while
  readCoordinates();
  chooseAlgorithm();
} // getMode()

void Zoo::printHelp(char *argv[]) {
  std::cout << "Usage: " << argv[0] << " [-o map|list | -q  |  -s  |  -h\n";
  std::cout << "This program is to help you learn command-line processing,\n";
  std::cout
      << "reading data into a vector, the difference between resize and\n";
  std::cout << "reserve and how to properly read until end-of-file."
            << std::endl;
} // printHelp()

void Zoo::readCoordinates() {
  std::cin >> NodeCount;
  int x;
  int y;
  vertices.reserve(NodeCount);
  while (std::cin >> x) {
    std::cin >> y;
    // Border if x or y is 0
    //
    if ((y == 0 && x <= 0) || (x == 0 && y <= 0)) {
      coordinate cord(x, y, WALL);
      vertices.push_back(cord);
    } else if (x < 0 && y < 0) {
      coordinate cord(x, y, WILD);
      vertices.push_back(cord);
    } else {
      coordinate cord(x, y, SAFER);
      vertices.push_back(cord);
    }
  }
}

void Zoo::chooseAlgorithm() {
  if (mode == "MST") {
    primTable.resize(NodeCount, PrimData());
    solveMST();
    printMST();
  } else if (mode == "FASTTSP") {
    greedyNearest.resize(NodeCount, PrimData());
    partB();
    printFASTTSP();
  } else
    greedyNearest.resize(NodeCount, PrimData());
  partC();
}

double Zoo::Euclerian(size_t vertice1, size_t vertice2) {
  return std::sqrt(((static_cast<double>(vertices[vertice1].x) -
                     static_cast<double>(vertices[vertice2].x)) *
                    (static_cast<double>(vertices[vertice1].x) -
                     static_cast<double>(vertices[vertice2].x))) +
                   ((static_cast<double>(vertices[vertice1].y) -
                     static_cast<double>(vertices[vertice2].y)) *
                    (static_cast<double>(vertices[vertice1].y) -
                     static_cast<double>(vertices[vertice2].y))));
}

double Zoo::checkDistance(size_t vertice1, size_t vertice2) {
  // check which zone it is in.
  // If one is in the wild, the other point has to be a wild or the border
  // if (vertice1 == 19 && vertice2 == 8)
  //   std::cout << vertices[vertice1].category << " "
  //             << vertices[vertice2].category << "\n";

  if ((vertices[vertice1].category == WILD &&
       vertices[vertice2].category == SAFER) ||
      (vertices[vertice1].category == SAFER &&
       vertices[vertice2].category == WILD)) {
    return -1;
  }
  // Use Euclerian Distance
  return std::sqrt(((static_cast<double>(vertices[vertice1].x) -
                     static_cast<double>(vertices[vertice2].x)) *
                    (static_cast<double>(vertices[vertice1].x) -
                     static_cast<double>(vertices[vertice2].x))) +
                   ((static_cast<double>(vertices[vertice1].y) -
                     static_cast<double>(vertices[vertice2].y)) *
                    (static_cast<double>(vertices[vertice1].y) -
                     static_cast<double>(vertices[vertice2].y))));
  // If
}

void Zoo::solveMST() {
  // // 1. Set the d value of the starting vertex to 0
  primTable[0].distance = 0;
  // 5. Reapeat steps 2-4 until j is true for every vertex (i.e., ever vertex
  // has been added to the tree)
  while (true) {
    // 2. From the set of verices for which k is false, select the vertex v that
    // has teh smlalest value of d
    int smallestIndex = -1;
    double smallestDistance = -1;
    for (size_t j = 0; j < NodeCount; j++) {
      // If smallestDistance is unitalized
      if (smallestDistance == -1 && !primTable[j].k) {
        smallestIndex = static_cast<int>(j);
        smallestDistance = primTable[j].distance;
      } else if (!primTable[j].k && primTable[j].distance < smallestDistance) {
        smallestIndex = static_cast<int>(j);
        smallestDistance = primTable[j].distance;
      }
    }
    if (smallestIndex == -1)
      break;
    // 3. Set the value of k for this vertex to true.
    primTable[static_cast<size_t>(smallestIndex)].k = true;
    //! add to total weight, This might not work.
    weightTotal += primTable[static_cast<size_t>(smallestIndex)].distance;
    // 4. For each vertex w adjact to v for which kw is false, check whether dw
    // is greater than the edge weight that connects v and w. If it is, set d to
    // the weight of thee edge that connects v and w, and set p to vertex v.
    // int count = 0; FORGOT WHAT THESE ARE FOR
    for (size_t a = 0; a < NodeCount; a++) {
      // check that the vertice is false, and if the distance is less than
      if (!primTable[a].k) {
        double distance = checkDistance(static_cast<size_t>(smallestIndex), a);
        // Checks that is valid and if so it is less than current distance. If
        // so change distance and change index.
        if (distance != -1 && distance < primTable[a].distance) {
          primTable[a].distance = distance;
          primTable[a].index = smallestIndex;
        } // else if (distance != -1)
          // count++; FORGOT WHAT THESE ARE FOR
      }
    }
    // returns error if there was no adjacent .
  }
}

void Zoo::printMST() {
  std::cout << std::fixed << std::setprecision(2) << weightTotal << "\n";
  for (int32_t i = 1; i < static_cast<int32_t>(NodeCount); i++) {
    if (i < primTable[static_cast<size_t>(i)].index)
      std::cout << i << " " << primTable[static_cast<size_t>(i)].index << "\n";
    else
      std::cout << primTable[static_cast<size_t>(i)].index << " " << i << "\n";
  }
}

// void Zoo::printFASTTSP() {
//   std::cout << std::fixed << std::setprecision(2) << weightTSP << "\n";
//   int curr = 0;
//   std::cout << 0;
//   curr = greedyNearest[0].index;
//   while (curr != 0) {
//     std::cout << " " << curr;
//     curr = greedyNearest[static_cast<size_t>(curr)].index;
//   }
//   std::cout << "\n";
// }

void Zoo::printFASTTSP() {
  std::cout << std::fixed << std::setprecision(2) << weightTSP << "\n";
  std::cout << bestPathSeen[0];
  for (size_t i = 1; i < bestPathSeen.size(); i++) {
    std::cout << " " << bestPathSeen[i];
  }
  std::cout << "\n";
}

void Zoo::partB() {
  // 1. Initialize a partial tour with a vertexx i, chosen arbitrarily
  greedyNearest[0].k = true;
  greedyNearest[0].index = 1;
  greedyNearest[0].distance = Euclerian(0, 1);
  // 2. Choose another arbitrary vertex j and set the initial partial tour to i
  // to j to i
  greedyNearest[1].k = true;
  greedyNearest[1].index = 0;
  greedyNearest[1].distance = Euclerian(0, 1);
  // 3. Arbitrarily select a vertex k that is currently not in the partial tour.
  // 4. Find the best place to insert vertexx k into the partial tour to
  // minimize cost
  int size = 2;
  double minDistance = -1;
  int32_t minIndex = -1;
  for (int j = 2; j < static_cast<int>(NodeCount); j++) {
    // Loops through all the possible options
    for (int i = 0; i < size; i++) {
      // Calculate the minDistance it would be if you were to reconnect. First
      // connection between where your checking and new vertex. Second: vertex
      // distance from where the where you checking is point to and new vertex.
      double totalDistance =
          Euclerian(static_cast<size_t>(j), static_cast<size_t>(i)) +
          Euclerian(static_cast<size_t>(j),
                    static_cast<size_t>(
                        greedyNearest[static_cast<size_t>(i)].index)) -
          greedyNearest[static_cast<size_t>(i)].distance;
      // Euclerian(static_cast<size_t>(i),
      // static_cast<size_t>(greedyNearest[static_cast<size_t>(i)].index));
      // Make sure a vertex cant create a path to self
      if (i != j) {
        // Initalize the first value
        if (i == 0) {
          minDistance = totalDistance;
          minIndex = i;
        }
        // Find the smallest Distance.
        else if (minDistance > totalDistance) {
          minDistance = totalDistance;
          minIndex = i;
        }
      }
    }
    size++;
    // Now connect everything after finding the smallest Distance.
    // Start by setting the values of the vertex
    greedyNearest[static_cast<size_t>(j)].index =
        greedyNearest[static_cast<size_t>(minIndex)].index;
    greedyNearest[static_cast<size_t>(j)].distance = Euclerian(
        static_cast<size_t>(j),
        static_cast<size_t>(greedyNearest[static_cast<size_t>(j)].index));
    greedyNearest[static_cast<size_t>(j)].k = true;
    // Reatch the old point
    greedyNearest[static_cast<size_t>(minIndex)].index = j;
    greedyNearest[static_cast<size_t>(minIndex)].distance =
        Euclerian(static_cast<size_t>(minIndex), static_cast<size_t>(j));
  }
  for (int i = 0; i < static_cast<int>(NodeCount); i++) {
    weightTSP += greedyNearest[static_cast<size_t>(i)].distance;
  }
  int curr = 0;
  bestPathSeen.push_back(0);
  curr = greedyNearest[0].index;
  while (curr != 0) {
    bestPathSeen.push_back(curr);
    curr = greedyNearest[static_cast<size_t>(curr)].index;
  }
}

double Zoo::calculateCurrCost(size_t permLength) {
  double totalcost = 0;
  for (size_t i = 0; i < permLength - 1; i++) {
    totalcost += Euclerian(i, i + 1);
  }
  return totalcost;
}

double Zoo::connectingArms(size_t permLength) {
  // beginning vertex = 0
  // final vertex - permLength - 1
  double minValue = -1;
  for (size_t i = permLength; i < pathBuilding.size(); i++) {
    // Initalizes minValue, length of first arm (first vertex) + length of
    // second arm (final vertex)
    double value = Euclerian(pathBuilding[i], 0) +
                   Euclerian(pathBuilding[i], permLength - 1);
    if (i == permLength)
      minValue = value;
    else if (minValue > value)
      minValue = value;
  }
  return minValue;
}

double Zoo::unvistedMST(size_t permLength) {
  size_t count = bestPathSeen.size() - permLength;
  // resize so it fits.
  std::vector<PrimData> primTable;
  primTable.resize(count, PrimData());
  // primTable.resize(count, PrimData());
  primTable[0].distance = 0;
  double weightMST = 0;
  while (true) {
    int smallestIndex = -1;
    double smallestDistance = -1;
    for (size_t j = 0; j < count; j++) {
      // If smallestDistance is unitalized
      if (smallestDistance == -1 && !primTable[j].k) {
        smallestIndex = static_cast<int>(j);
        smallestDistance = primTable[j].distance;
      } else if (!primTable[j].k && primTable[j].distance < smallestDistance) {
        smallestIndex = static_cast<int>(j);
        smallestDistance = primTable[j].distance;
      }
    }
    if (smallestIndex == -1)
      break;
    primTable[static_cast<size_t>(smallestIndex)].k = true;
    //! This is wrong because im not changing the smallestDistance
    weightMST += primTable[static_cast<size_t>(smallestIndex)].distance;
    for (size_t a = 0; a < count; a++) {
      // check that the vertice is false, and if the distance is less than
      if (!primTable[a].k) {
        double distance = Euclerian(permLength + smallestIndex, permLength + a);
        // Checks that is valid and if so it is less than current distance. If
        // so change distance and change index.
        if (distance < primTable[a].distance) {
          primTable[a].distance = distance;
          primTable[a].index = smallestIndex;
        } // else if (distance != -1)
          // count++; FORGOT WHAT THESE ARE FOR
      }
    }
  }
  return weightMST;
}

bool Zoo::promising(size_t permLength) {
  if (permLength <= 4)
    return true;
  double unvistedMSTs = unvistedMST(permLength);
  double arms = connectingArms(permLength);
  double test = unvistedMSTs + arms + currPathLength;
  std::cout << "k: " << permLength << " currcost: " << currPathLength
            << " MST: " << unvistedMSTs << " arms: " << arms
            << " Total: " << test << " Best: " << lengthOfBestPath << "\n";
  return ((unvistedMSTs + arms + currPathLength) < lengthOfBestPath);
}

void Zoo::genPerms(size_t permLength) {
  // We have reached leaf node everything has been fixed
  if (permLength == pathBuilding.size()) {
    // Do something with the path
    // get total path length totalpath length < bestSoFar.
    // currCost += (closing edge)
    currPathLength += Euclerian(0, pathBuilding.size() - 1);
    // Check if it's better
    if (currPathLength < lengthOfBestPath) {
      lengthOfBestPath = currPathLength;
      // copy over the values to best path seen
      for (size_t i = 0; i < bestPathSeen.size(); i++) {
        std::cout << "size best Path:" << bestPathSeen.size();
        std::cout << "size best Path:" << bestPathSeen.size();
        bestPathSeen[i] = pathBuilding[i];
      }
    }
    //  update things (best path and best cost are changed)
    // currCost -= (closing edge) (should be close to 0 or 0 at the end)
    currPathLength -= Euclerian(0, pathBuilding.size() - 1);
    return;
  } // if ..complete path
  // If not promising go back and make another choice
  if (!promising(permLength)) {
    return;
  } // if ..not promising

  for (size_t i = permLength; i < lengthOfBestPath; ++i) {
    // Pick someone new to be part of the path
    std::swap(pathBuilding[permLength], pathBuilding[i]);
    // curCost =+ (cost of the new edge) (a,b,c) (a+ b) + (b + c) Use perm
    currPathLength += calculateCurrCost(permLength);
    // length to help with index.
    genPerms(permLength + 1);
    // currCost -= (cost of that same edge)
    currPathLength -= calculateCurrCost(permLength);
    std::swap(pathBuilding[permLength], pathBuilding[i]);
  } // for ..unpermuted elements
}

void Zoo::partC() {
  // builds bestPathSeen and best length.
  partB();
  lengthOfBestPath = weightTSP;
  // Add zero to the start
  for (auto x : bestPathSeen) {
    pathBuilding.push_back(x);
  }
  currPathLength = 0;
  // start the genPerms
  genPerms(1);
  std::cout << lengthOfBestPath << "\n";
  std::cout << bestPathSeen[0];
  for (size_t i = 1; i < bestPathSeen.size(); i++) {
    std::cout << " " << bestPathSeen[i];
  }
  std::cout << "\n";
  return;
}
