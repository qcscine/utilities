/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Geometry/ElementInfo.h>
#include <Utils/IO/ChemicalFileFormats/PdbStreamHandler.h>
#include <gmock/gmock.h>
#include <iostream>
#include <sstream>
#include <vector>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

class PdbStreamHandlerTest : public Test {
 public:
  PdbStreamHandler handler;
};

std::string tyrosine = "HETATM    1  N   UNK     1       1.171  -2.180   0.228  1.00  0.00           N \n"
                       "HETATM    2  C   UNK     2       1.491  -0.984  -0.582  1.00  0.00           C\n"
                       "HETATM    3  C   UNK     3       2.946  -0.962  -1.046  1.00  0.00           C\n"
                       "HETATM    4  O   UNK     4       3.915  -1.497  -0.536  1.00  0.00           O\n"
                       "HETATM    5  O   UNK     5       3.185  -0.264  -2.178  1.00  0.00           O\n"
                       "HETATM    6  C   UNK     6       1.205   0.310   0.201  1.00  0.00           C\n"
                       "HETATM    7  H   UNK     7       1.674   0.272   1.206  1.00  0.00           H\n"
                       "HETATM    8  H   UNK     8       1.681   1.167  -0.320  1.00  0.00           H\n"
                       "HETATM    9  C   UNK     9      -0.262   0.568   0.317  1.00  0.00           C\n"
                       "HETATM   10  C   UNK    10      -0.964   1.087  -0.775  1.00  0.00           C\n"
                       "HETATM   11  H   UNK    11      -0.432   1.295  -1.711  1.00  0.00           H\n"
                       "HETATM   12  C   UNK    12      -2.324   1.341  -0.687  1.00  0.00           C\n"
                       "HETATM   13  H   UNK    13      -2.871   1.749  -1.544  1.00  0.00           H\n"
                       "HETATM   14  C   UNK    14      -2.990   1.069   0.517  1.00  0.00           C\n"
                       "HETATM   15  C   UNK    15      -2.299   0.548   1.618  1.00  0.00           C\n"
                       "HETATM   16  H   UNK    16      -2.818   0.332   2.558  1.00  0.00           H\n"
                       "HETATM   17  C   UNK    17      -0.936   0.301   1.508  1.00  0.00           C\n"
                       "HETATM   18  H   UNK    18      -0.391  -0.112   2.364  1.00  0.00           H\n"
                       "HETATM   19  O   UNK    19      -4.331   1.343   0.538  1.00  0.00           O\n"
                       "HETATM   20  H   UNK    20      -4.663   1.103   1.394  1.00  0.00           H\n"
                       "HETATM   21  H   UNK    21       0.824  -1.021  -1.482  1.00  0.00           H\n"
                       "HETATM   22  H   UNK    22       4.113  -0.293  -2.390  1.00  0.00           H\n"
                       "HETATM   23  H   UNK    23       1.681  -2.167   1.087  1.00  0.00           H\n"
                       "HETATM   24  H   UNK    24       1.396  -3.005  -0.284  1.00  0.00           H\n"
                       "HETATM    1  O   HOH     1      -0.711   3.428   0.621  0.00  0.00           O\n"
                       "HETATM    2  O   HOH     1      -1.659  -1.373  -0.032  0.00  0.00           O\n"
                       "HETATM    3  H   HOH     1      -0.743   2.872   1.402  0.00  0.00           H\n"
                       "HETATM    4  H   HOH     1      -1.403   4.092   0.670  0.00  0.00           H\n"
                       "HETATM    5  H   HOH     1      -0.764  -1.432  -0.374  0.00  0.00           H\n"
                       "HETATM    6  H   HOH     1      -1.909  -0.450   0.053  0.00  0.00           H\n"
                       "HETATM    7  N   UNK     1       1.171  -2.180   0.228  1.00  0.00           N\n"
                       "CONECT    1    2   23   24 \n"
                       "CONECT    2    1    3    6   21\n"
                       "CONECT    3    2    4    5\n"
                       "CONECT    4    3\n"
                       "CONECT    5    3   22\n"
                       "CONECT    6    2    7    8    9\n"
                       "CONECT    7    6\n"
                       "CONECT    8    6\n"
                       "CONECT    9    6   10   17\n"
                       "CONECT   10    9   11   12\n"
                       "CONECT   11   10\n"
                       "CONECT   12   10   13   14\n"
                       "CONECT   13   12\n"
                       "CONECT   14   12   15   19\n"
                       "CONECT   15   14   16   17\n"
                       "CONECT   16   15\n"
                       "CONECT   17    9   15   18\n"
                       "CONECT   18   17\n"
                       "CONECT   19   14   20\n"
                       "CONECT   20   19\n"
                       "CONECT   21    2\n"
                       "CONECT   22    5\n"
                       "CONECT   23    1\n"
                       "CONECT   24    1\n"
                       "END";

std::string structureWithOverlay = "ATOM      1  CB ATYR A  14     -13.771   2.911 -17.239  0.50 27.81           C\n"
                                   "ATOM      2  CB BTYR A  14     -13.786   2.928 -17.233  0.50 27.85           C\n"
                                   "ATOM      3  CG ATYR A  14     -13.608   4.388 -17.494  0.50 28.31           C\n"
                                   "ATOM      4  CG BTYR A  14     -13.829   4.427 -17.427  0.50 28.20           C\n"
                                   "ATOM      5  CD1ATYR A  14     -12.442   4.893 -18.059  0.50 23.25           C\n"
                                   "ATOM      6  CD1BTYR A  14     -12.714   5.121 -17.879  0.50 25.75           C\n"
                                   "ATOM      7  CD2ATYR A  14     -14.614   5.282 -17.153  0.50 28.36           C\n"
                                   "ATOM      8  CD2BTYR A  14     -14.980   5.152 -17.138  0.50 27.68           C\n"
                                   "ATOM      9  CE1ATYR A  14     -12.289   6.252 -18.283  0.50 27.32           C\n"
                                   "ATOM     10  CE1BTYR A  14     -12.746   6.496 -18.047  0.50 27.49           C\n"
                                   "ATOM     11  CE2ATYR A  14     -14.472   6.636 -17.372  0.50 29.32           C\n"
                                   "ATOM     12  CE2BTYR A  14     -15.021   6.527 -17.300  0.50 28.26           C\n"
                                   "ATOM     13  CZ ATYR A  14     -13.311   7.117 -17.938  0.50 29.06           C\n"
                                   "ATOM     14  CZ BTYR A  14     -13.902   7.193 -17.757  0.50 29.57           C\n"
                                   "ATOM     15  OH ATYR A  14     -13.179   8.469 -18.153  0.50 33.01           O\n"
                                   "ATOM     16  OH BTYR A  14     -13.938   8.561 -17.922  0.50 32.79           O\n"
                                   "TER\n"
                                   "END";

TEST_F(PdbStreamHandlerTest, Selfconsistent) {
  for (const auto& PDBInput : {tyrosine, structureWithOverlay}) {
    std::stringstream in(PDBInput), out;
    auto structure = handler.read(in);
    handler.write(out, structure[0], "Test PDB file");
    auto secondTime = handler.read(out);
    // Make sure that Hydrogens are parsed
    ASSERT_THAT(structure[0].size(), Eq(secondTime[0].size()));

    // Compare AtomCollections
    PositionCollection p1 = secondTime[0].getPositions();
    PositionCollection p2 = structure[0].getPositions();
    for (int i = 0; i < p1.rows(); ++i) {
      ASSERT_TRUE(ElementInfo::Z(secondTime[0].getElement(i)) == ElementInfo::Z(structure[0].getElement(i)));
      for (int j = 0; j < p1.cols(); ++j) {
        ASSERT_THAT(p1(i, j), DoubleNear(p2(i, j), 1e-4));
      }
    }
  }
}

TEST_F(PdbStreamHandlerTest, HydrogensAreSkipped) {
  std::stringstream in(tyrosine), out;

  // Change the setting
  handler.setReadH(false);
  auto structureWithoutHydrogens = handler.read(in);
  ASSERT_THAT(structureWithoutHydrogens[0].size(), Eq(16));
}

TEST_F(PdbStreamHandlerTest, SolventIsSkipped) {
  std::stringstream in(tyrosine), out;

  // Change the setting
  handler.setReadHOH(false);
  auto structureWithoutSolvent = handler.read(in);
  ASSERT_THAT(structureWithoutSolvent[0].size(), Eq(25));
}

TEST_F(PdbStreamHandlerTest, CommentIsSetCorrectly) {
  std::stringstream in(tyrosine), out;
  auto structure = handler.read(in);
  // Pick same comment as in input
  PdbStreamHandler::write(out, structure[0], "comment");
  std::string line;
  std::istringstream f(out.str());
  std::getline(f, line);
  ASSERT_THAT(line, Eq("comment"));
}

TEST_F(PdbStreamHandlerTest, TwoSubstructuresAreReturnedCorrectly) {
  std::stringstream in(structureWithOverlay), out;

  auto structure = handler.read(in);
  ASSERT_THAT(structure[0].size(), Eq(structure[1].size()));
}

} // namespace Tests
} // namespace Utils
} // namespace Scine
