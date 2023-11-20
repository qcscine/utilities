/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Geometry/ElementInfo.h>
#include <Utils/IO/ChemicalFileFormats/PdbStreamHandler.h>
#include <gmock/gmock.h>
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

std::string tyrosine = "HETATM    1  N   UNK     1       1.171  -2.180   0.228  1.00  0.00           N\n"
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

std::string structureWithModels = "MODEL        1\n"
                                  "ATOM      1  N   ILE A   1      -7.158   5.359   0.606  1.00  0.00           N\n"
                                  "ATOM      2  CA  ILE A   1      -5.843   5.515  -0.080  1.00  0.00           C\n"
                                  "ATOM      3  C   ILE A   1      -4.773   4.681   0.627  1.00  0.00           C\n"
                                  "ATOM      4  O   ILE A   1      -4.775   4.541   1.834  1.00  0.00           O\n"
                                  "ATOM      5  CB  ILE A   1      -5.512   7.005   0.020  1.00  0.00           C\n"
                                  "ATOM      6  CG1 ILE A   1      -4.205   7.285  -0.725  1.00  0.00           C\n"
                                  "ATOM      7  CG2 ILE A   1      -5.350   7.396   1.490  1.00  0.00           C\n"
                                  "ATOM      8  CD1 ILE A   1      -3.902   8.784  -0.683  1.00  0.00           C\n"
                                  "ATOM     22  N   CYS A   2      -3.859   4.125  -0.119  1.00  0.00           N\n"
                                  "ATOM     23  CA  CYS A   2      -2.785   3.297   0.503  1.00  0.00           C\n"
                                  "ATOM     24  C   CYS A   2      -1.717   2.955  -0.538  1.00  0.00           C\n"
                                  "ATOM     25  O   CYS A   2      -1.962   2.227  -1.480  1.00  0.00           O\n"
                                  "ATOM     26  CB  CYS A   2      -3.469   2.018   1.013  1.00  0.00           C\n"
                                  "ATOM     27  SG  CYS A   2      -4.843   1.547  -0.078  1.00  0.00           S\n"
                                  "ATOM     32  N   VAL A   3      -0.533   3.477  -0.374  1.00  0.00           N\n"
                                  "ATOM     33  CA  VAL A   3       0.556   3.186  -1.353  1.00  0.00           C\n"
                                  "ATOM     34  C   VAL A   3       1.888   3.026  -0.633  1.00  0.00           C\n"
                                  "ATOM     35  O   VAL A   3       2.940   2.958  -1.238  1.00  0.00           O\n"
                                  "ATOM     36  CB  VAL A   3       0.589   4.401  -2.271  1.00  0.00           C\n"
                                  "ATOM     37  CG1 VAL A   3       1.897   4.417  -3.067  1.00  0.00           C\n"
                                  "ATOM     38  CG2 VAL A   3      -0.597   4.339  -3.236  1.00  0.00           C\n"
                                  "ENDMDL                                                                        \n"
                                  "MODEL        2                                                                \n"
                                  "ATOM      1  N   ILE A   1      -6.479   5.364  -0.431  1.00  0.00           N\n"
                                  "ATOM      2  CA  ILE A   1      -5.250   6.198  -0.291  1.00  0.00           C\n"
                                  "ATOM      3  C   ILE A   1      -4.192   5.440   0.512  1.00  0.00           C\n"
                                  "ATOM      4  O   ILE A   1      -3.960   5.716   1.672  1.00  0.00           O\n"
                                  "ATOM      5  CB  ILE A   1      -5.705   7.452   0.455  1.00  0.00           C\n"
                                  "ATOM      6  CG1 ILE A   1      -6.924   8.049  -0.253  1.00  0.00           C\n"
                                  "ATOM      7  CG2 ILE A   1      -4.571   8.480   0.468  1.00  0.00           C\n"
                                  "ATOM      8  CD1 ILE A   1      -7.583   9.095   0.648  1.00  0.00           C\n"
                                  "ATOM     22  N   CYS A   2      -3.551   4.484  -0.101  1.00  0.00           N\n"
                                  "ATOM     23  CA  CYS A   2      -2.505   3.701   0.620  1.00  0.00           C\n"
                                  "ATOM     24  C   CYS A   2      -1.476   3.162  -0.373  1.00  0.00           C\n"
                                  "ATOM     25  O   CYS A   2      -1.754   2.278  -1.159  1.00  0.00           O\n"
                                  "ATOM     26  CB  CYS A   2      -3.223   2.535   1.326  1.00  0.00           C\n"
                                  "ATOM     27  SG  CYS A   2      -4.785   2.120   0.496  1.00  0.00           S\n"
                                  "ATOM     32  N   VAL A   3      -0.287   3.690  -0.335  1.00  0.00           N\n"
                                  "ATOM     33  CA  VAL A   3       0.777   3.216  -1.265  1.00  0.00           C\n"
                                  "ATOM     34  C   VAL A   3       2.063   2.965  -0.493  1.00  0.00           C\n"
                                  "ATOM     35  O   VAL A   3       3.115   2.740  -1.058  1.00  0.00           O\n"
                                  "ATOM     36  CB  VAL A   3       0.962   4.353  -2.261  1.00  0.00           C\n"
                                  "ATOM     37  CG1 VAL A   3       2.221   4.104  -3.094  1.00  0.00           C\n"
                                  "ATOM     38  CG2 VAL A   3      -0.258   4.425  -3.181  1.00  0.00           C\n"
                                  "ENDMDL                                                                        \n"
                                  "MODEL        3                                                                \n"
                                  "ATOM      1  N   ILE A   1      -5.132   5.503  -2.479  1.00  0.00           N\n"
                                  "ATOM      2  CA  ILE A   1      -4.984   6.235  -1.190  1.00  0.00           C\n"
                                  "ATOM      3  C   ILE A   1      -3.959   5.516  -0.305  1.00  0.00           C\n"
                                  "ATOM      4  O   ILE A   1      -3.180   6.139   0.388  1.00  0.00           O\n"
                                  "ATOM      5  CB  ILE A   1      -6.395   6.234  -0.572  1.00  0.00           C\n"
                                  "ATOM      6  CG1 ILE A   1      -7.027   7.611  -0.783  1.00  0.00           C\n"
                                  "ATOM      7  CG2 ILE A   1      -6.341   5.931   0.933  1.00  0.00           C\n"
                                  "ATOM      8  CD1 ILE A   1      -8.544   7.507  -0.626  1.00  0.00           C\n"
                                  "ATOM     22  N   CYS A   2      -3.948   4.215  -0.330  1.00  0.00           N\n"
                                  "ATOM     23  CA  CYS A   2      -2.962   3.469   0.505  1.00  0.00           C\n"
                                  "ATOM     24  C   CYS A   2      -1.846   2.907  -0.377  1.00  0.00           C\n"
                                  "ATOM     25  O   CYS A   2      -2.069   2.067  -1.226  1.00  0.00           O\n"
                                  "ATOM     26  CB  CYS A   2      -3.755   2.343   1.168  1.00  0.00           C\n"
                                  "ATOM     27  SG  CYS A   2      -4.514   1.301  -0.099  1.00  0.00           S\n"
                                  "ATOM     32  N   VAL A   3      -0.645   3.376  -0.183  1.00  0.00           N\n"
                                  "ATOM     33  CA  VAL A   3       0.500   2.890  -1.006  1.00  0.00           C\n"
                                  "ATOM     34  C   VAL A   3       1.607   2.356  -0.104  1.00  0.00           C\n"
                                  "ATOM     35  O   VAL A   3       2.708   2.086  -0.542  1.00  0.00           O\n"
                                  "ATOM     36  CB  VAL A   3       0.939   4.139  -1.777  1.00  0.00           C\n"
                                  "ATOM     37  CG1 VAL A   3       2.448   4.401  -1.626  1.00  0.00           C\n"
                                  "ATOM     38  CG2 VAL A   3       0.594   3.963  -3.257  1.00  0.00           C\n"
                                  "ENDMDL                                                                        \n"
                                  "END";

TEST_F(PdbStreamHandlerTest, Selfconsistent) {
  for (const auto& PDBInput : {tyrosine, structureWithOverlay}) {
    std::stringstream in(PDBInput), out;
    auto structure = handler.read(in);
    handler.write(out, structure[0], BondOrderCollection(), "Test PDB file");
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

TEST_F(PdbStreamHandlerTest, HydrogensAndSolventAreSkipped) {
  std::stringstream in(tyrosine), out;

  // Change the setting
  handler.setReadH(false);
  handler.parseOnlySolvent(false);
  auto structureWithoutHydrogens = handler.read(in);

  ASSERT_THAT(structureWithoutHydrogens[0].size(), Eq(14));

  // change the setting
  handler.setReadH(true);
  in.clear();
  in.seekg(0);
  auto structureWithHydrogens = handler.read(in);
  ASSERT_THAT(structureWithHydrogens[0].size(), Eq(25));

  // Change the setting
  handler.parseOnlySolvent(true);
  in.clear();
  in.seekg(0);
  auto solvent = handler.read(in);
  ASSERT_THAT(solvent[0].size(), Eq(6));
}

TEST_F(PdbStreamHandlerTest, CommentIsSetCorrectly) {
  std::stringstream in(tyrosine), out;
  auto structure = handler.read(in);
  // Pick same comment as in input
  PdbStreamHandler::write(out, structure[0], BondOrderCollection(), "comment");
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

TEST_F(PdbStreamHandlerTest, ThreeModelsAreReturnedCorrectly) {
  std::stringstream in(structureWithModels), out;

  auto structure = handler.read(in);
  ASSERT_THAT(structure.size(), Eq(3));
  ASSERT_THAT(structure[0].size(), Eq(structure[1].size()));
  ASSERT_THAT(structure[1].size(), Eq(structure[2].size()));
}

} // namespace Tests
} // namespace Utils
} // namespace Scine
