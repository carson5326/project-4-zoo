// IDENTIFIER:3E33912F8BAA7542FC4A1585D2DB6FE0312725B9
#include <cstdint>
#include <getopt.h>
#include <limits>
#include <string>
#include <vector>

class Zoo {
public:
  enum Category { SAFER, WILD, WALL };

  struct coordinate {
    int x;
    int y;
    Category category;
    coordinate(int x, int y, Zoo::Category category)
        : x(x), y(y), category(category) {}
  };

  struct PrimData {

    double distance;
    int32_t index;
    bool k;
    PrimData()
        : distance{std::numeric_limits<double>::infinity()}, index{-1},
          k{false} {}
  };

  Zoo() {
    weightTotal = 0;
    weightTSP = 0;
  };
  void getMode(int argc, char *argv[]);

  void printHelp(char *argv[]);

  void readCoordinates();

  void chooseAlgorithm();

  double Euclerian(size_t vertice1, size_t vertice2);

  double checkDistance(size_t vertice1, size_t vertice2);

  void solveMST();

  void printMST();

  void printFASTTSP();

  void partB();

  // PartC: functions
  double connectingArms(size_t permLength);

  double unvistedMST(size_t permLength);

  bool promising(size_t permLength);

  double calculateCurrCost(size_t permLength);

  void genPerms(size_t permLength);

  void partC();

private:
  // The mode the prgram will run in MST = 1, FASTTSP = 2, and OPTTSP == 3
  std::string mode;
  uint32_t NodeCount;
  // Stores coordinates
  std::vector<coordinate> vertices;
  // Part A:
  std::vector<PrimData> primTable;
  double weightTotal;

  // Part B:
  std::vector<PrimData> greedyNearest;
  double weightTSP;

  // Part C: Variables
  std::vector<size_t> pathBuilding;
  // Running total
  double currPathLength;
  std::vector<size_t> bestPathSeen;
  double lengthOfBestPath;
};
