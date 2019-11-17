// $Id: tsp-subtour.cpp 2469 2019-01-06 23:17:46Z unxusr $

/*! \file tsp-subtour.cpp
  \brief Solves the TSP using dinamically generated subtour elimination constraints

  Solver for the traveling salesman problem using lazy constraints

*/

#include <cstdio>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <string>
#include <OsiClpSolverInterface.hpp>
#include <CbcModel.hpp>
#include <OsiAuxInfo.hpp>
#include <CglCutGenerator.hpp>
#include <cmath>
#include <algorithm>
#include <CoinTime.hpp>

// when separating fractional cuts
#define MIN_VIOLATION 1e-5

int verbose = 1;

int time_limit = 3600;

// leave empty if no log log should be generated
char results_file[] = "results.csv";

using namespace std;

typedef struct 
{
  int n;          // number of cities
  int **c = NULL; // distance matrix
  int **x; // indexes of x variables
} TSPInfo;

// receives info of tsp instance and an integer solution, returns size of
// subtour and fills elements in els
int find_subtour( const TSPInfo *tspi, const double *s, int *els, int start );


// receives info of tsp instance and a fractional solution, returns the
// size of the subtour found and fills elements in els
int find_subtour_frac( const TSPInfo *tspi, const double *s, int *els, int start );

// for debugging
void print_sol(int n, const int **x, const double *s);

// check in a solution which is the output arc of node i
int out_arc(int n, const double *s, const int **x, int i);

class CglSubTour : public CglCutGenerator {
  public:
    CglSubTour(const TSPInfo *info);

    CglSubTour(const CglSubTour& rhs);

    virtual CglCutGenerator * clone() const;

    virtual void
      generateCuts(const OsiSolverInterface& si, OsiCuts & cs,
          const CglTreeInfo info = CglTreeInfo());

    const TSPInfo *info;
};

int main(int argc, char **argv)
{
  if (argc < 2) {
    fprintf(stderr, "Enter instance name (.dist file)");
    exit(EXIT_FAILURE);
  }

  double start = CoinCpuTime();

  int n;          // number of cities
  int **c = NULL; // distance matrix
  int **x; // indexes of x variables

  // reading distance matrix
  FILE *f = fopen( argv[1], "r" );
  fscanf(f, "%d", &n);
  c = new int*[n];
  c[0] = new int[n*n];
  for ( int i=1 ; (i<n) ; ++i )
    c[i] = c[i-1]+n;
  fill(c[0], c[0]+n*n, 0);

  for ( int i=0 ; (i<n) ; ++i )
    for ( int j=0 ; (j<n) ; ++j )
      fscanf(f, "%d", &(c[i][j]));

  x = new int*[n];
  x[0] = new int[n*n];
  for ( int i=1 ; (i<n) ; ++i )
    x[i] = x[i-1]+n;  

  fclose(f);

  int *idx = new int[n];
  double *coef = new double[n];
  fill(coef, coef+n, 1.0);
  double *lb  = new double[n*n];
  fill(lb, lb+n*n, 0.0);
  double *ub = new double[n*n];
  fill(ub, ub+n*n, 1.0);
  double *obj = new double[n*n];
  char name[256];
  vector< string > cnames;

  // creating x variables
  OsiClpSolverInterface *mip = new OsiClpSolverInterface();
  mip->messageHandler()->setLogLevel(0);
  int k = 0;
  for ( int i=0 ; (i<n) ; ++i ) {
    for ( int j=0 ; (j<n) ; ++j ) {
      sprintf(name, "x(%d,%d)", i, j);
      cnames.push_back(name);
      obj[k] = c[i][j];
      x[i][j] = k;    
      ++k;
    }
  }
  CoinBigIndex *starts = new CoinBigIndex[n*n+1];
  fill(starts, starts+n*n+1, 0);
  mip->addCols(k, starts, NULL, NULL, lb, ub, obj);

  // out degree constraints
  for ( int i=0 ; (i<n) ; ++i ) {
    k = 0;
    for ( int j=0 ; (j<n) ; ++j ) {
      if (i==j)
        continue;

      idx[k] = x[i][j];
      ++k;
    }
    mip->addRow(k, idx, coef, 1.0, 1.0);
  }
  // in degree constraints
  for ( int i=0 ; (i<n) ; ++i ) {
    k = 0;
    for ( int j=0 ; (j<n) ; ++j ) {
      if (i==j)
        continue;

      idx[k] = x[j][i];
      ++k;      
    }
    mip->addRow(k, idx, coef, 1.0, 1.0);
  }

  mip->setColNames(cnames, 0, cnames.size(), 0);
  for ( int i=0 ; (i<mip->getNumCols()) ; ++i )
    mip->setInteger(i);

  // 2-subtour elimination
  /*
  coef[0] = 1.0; coef[1] = 1.0;
  for ( int i=0 ; (i<n-1) ; ++i ) {
    for ( int j=i+1 ; (j<n) ; ++j ) {
      idx[0] = x[i][j];
      idx[1] = x[j][i];
      mip->addRow(2, idx, coef, -COIN_DBL_MAX, 1.0);
    }
  }*/

  CbcModel model(*mip);

  // initial formulation is incomplete
  OsiBabSolver defaultC;
  defaultC.setSolverType(4);
  model.solver()->setAuxiliaryInfo(&defaultC);
  model.passInSolverCharacteristics(&defaultC);
  model.setKeepNamesPreproc(true);

  TSPInfo info;
  info.n = n;
  info.c = c;
  info.x = x;

  CglSubTour cglst(&info);

  model.addCutGenerator( &cglst, 1, "subtour", true, true );

  model.setMaximumSeconds(time_limit);
  model.branchAndBound();

  double total_time = CoinCpuTime() - start;

  printf("Branch & Bound Finished with status: %d\n", model.status());
  printf("isProvenOptimal: %d\n", model.isProvenOptimal());

  char has_subtour = 0;
  if (model.getColSolution()) 
  {
    int *els = new int[n];
    int sts = find_subtour(&info, model.getColSolution(), els, 0);
    if (sts) 
    {
      printf("subtour with %d elements found: ", sts);
      for ( int i=0 ; (i<sts) ; ++i )
        printf(" %d", els[i]);
      printf("\n");
      has_subtour = 1;
    }
    delete[] els;

    const double *s= model.getColSolution();

    if (has_subtour == 0) {
      int next = 0;
      printf("Optimal route with length %g : 0 ", model.getObjValue() );
      do {
        next = out_arc(n, s, (const int **)x, next);
        printf("-> %d ", next);
      } while (next != 0);
      printf("\n");
      fflush(stdout);
    }
  }

  if (strlen(results_file)){
    char first_line = 1;
    f = fopen(results_file, "r");
    if (f) {
      first_line = 0;
      fclose(f);
    }

    char iname[256]; 
    strcpy(iname, argv[1]);
    char *s = strstr(iname, ".dist");
    if (s)
      *s = 0;

    f = fopen(results_file, "a");
    if (first_line)
      fprintf(f, "instance,time,status,isProvenOptimal,objValue,hasSubTour\n");

    fprintf(f, "%s,%.4f,%d,%d,%g,%d\n", iname, total_time, model.status(), model.isProvenOptimal(), model.getObjValue(), has_subtour );
    fclose(f);
  }

  delete[] idx;
  delete[] coef;

  delete mip;
  // freeing memoruy
  if (c) {
    delete[] c[0];
    delete[] c;
  }
  if (x) {
    delete[] x[0];
    delete[] x;
  }

  exit(0);
}

CglSubTour::CglSubTour(const TSPInfo *tspi)
{
  this->info = tspi;
}


CglSubTour::CglSubTour(const CglSubTour& rhs)
{
  this->info = rhs.info;
}

CglCutGenerator *CglSubTour::clone() const
{
  return new CglSubTour(this->info);
}

void CglSubTour::generateCuts(const OsiSolverInterface& si, OsiCuts & cs,
    const CglTreeInfo tinfo)
{
  const TSPInfo *tspi = this->info;
  const double *s = si.getColSolution();
  int n = tspi->n;
  const int **x = (const int **)tspi->x;

  double maxidist = 0.0;
  for ( int i=0 ; (i<n) ; ++i )
  {
    for ( int j=0 ; (j<n) ; ++j )
    {
      if (i==j)
        continue;
      double v = s[x[i][j]];
      double idist = fabs(v - floor(v+0.5));
      maxidist = max(idist, maxidist);
    }
  }

  static int cut_it = 0;
  ++cut_it;
  
  if (verbose == 2) {
    printf("\nmax frac value: %g, cut generation call: %d\n", maxidist, cut_it);
    fflush(stdout);
  }

  char separate_integer_sol = 0;
  if (maxidist <= si.getIntegerTolerance())
    separate_integer_sol = 1;

  int *idx = new int[n*n];
  double *coef = new double[n*n];
  fill(coef, coef+n*n, 1.0);

  int *els = new int[n];
  char *visited = new char[n];

  int newCuts = 0;

  for ( int startn = 0 ; (startn<n) ; startn++ )
  {
    int nnodes = 0;
    if (separate_integer_sol)
      nnodes = find_subtour(tspi, s, els, startn);
    else {
      //printf("FRAC SOL:\n");
      //print_sol(n, x, s);
      nnodes = find_subtour_frac(tspi, s, els, startn);
      //if (nnodes)
      //  printf("found cut in frac sol\n");
    }

    if (!nnodes)
      break;

    fill(visited, visited+n, 0);
    for ( int i=0 ; (i<nnodes) ; ++i )
      visited[els[i]] = 1;

    OsiRowCut cut;
    int nz = 0;
    for ( int i=0 ; (i<n) ; ++i )
      for (int j=0 ; (j<n) ; ++j )
        if (i!=j && visited[i] && visited[j])
          idx[nz++] = x[i][j];

    cut.setRow(nz, idx, coef);
    cut.setLb(-COIN_DBL_MAX);
    cut.setUb(nnodes-1);

    double violation = cut.violated(s);
    if (violation < MIN_VIOLATION)
      continue;

    int ncuts = cs.sizeRowCuts();
    cs.insertIfNotDuplicate(cut);

    newCuts += (cs.sizeRowCuts() - ncuts);

    if (verbose >= 2) {
      if (cs.sizeRowCuts() > ncuts)
      {
        printf("subtour found with nodes (violation %.5f):", violation);
        for ( int i=0 ; (i<nnodes) ; ++i )
          printf(" %d", els[i]);
        printf("\n");

        for ( int i=0 ; (i<nz) ; ++i )
          printf("%+g %s", coef[i], si.getColName(idx[i]).c_str());
        printf("<= %d\n", nnodes-1); fflush(stdout);
        printf("");
      }
    } // verbose == 2
  }

  if (verbose == 1) {
    printf("%d cut generation call with maxfrac: %g generated %d new sub-tour elimination constraints.\n", cut_it, maxidist, newCuts);
    fflush(stdout);
  }

  delete[] idx;
  delete[] visited;
  delete[] els;
  delete[] coef;
}

int find_subtour( const TSPInfo *tspi, const double *s, int *els, int start )
{
  int n = tspi->n;
  const int **x = (const int **)tspi->x;
  char *visited = new char[n];
  fill(visited, visited+n, 0);
  visited[start] = 1;

  char has_sub = 0;
  int nnodes = 1;
  int node = start;
  while (nnodes<n)
  {
    for (int j=0 ; j<n ; ++j )
    {
      if (j==node)
        continue;

      //printf("(%d, %d) ", node, j); fflush(stdout);
      if (s[x[node][j]] >= 0.99)
      {
        if (visited[j])
        {
          if (nnodes<n)
            has_sub = 1;
          else
            has_sub = 0;
          break;          
        }
        else
        {
          node = j;
          visited[node] = 1;
          nnodes++;
        }
      }
      if (has_sub)
        break;
    }
    if (has_sub)
      break;
  }

  if (has_sub) {
    int nn = 0;
    for ( int i=0 ; (i<n) ; ++i ) 
      if (visited[i])
        els[nn++] = i;

    delete[] visited;
    return nnodes;
  }
  delete[] visited;

  return 0;
}


int find_subtour_frac( const TSPInfo *tspi, const double *s, int *els, int start )
{
  int result = 0;
  // heuristically generates a series of sub-sets, inserting always the most connected 
  // elemented to the already included ones and returns the most disconnected sub-set (if any)

  int n = tspi->n;
  const int **x = (const int **)tspi->x;
  char *visited = new char[n];
  fill(visited, visited+n, 0);
  visited[start] = 1;

  int nnodes = 1;
  els[0] = start;

  // to compute how much each node is connected to the already inserted ones
  double *link = new double[n];

  // to store which of the subsets (selected from el[0] ... el[n-1])
  // is the most violated one
  double most_violated = -DBL_MAX;

  double links_inside = 0.0;
  while (nnodes<n-1)
  {
    // checks for the most connected node to the node already inserted ones
    // that still outside s
    fill(link, link+n, 0.0);

    // how much the already inserted nodes are connnected
    for ( int i=0 ; (i<n) ; ++i ) { // tail
      for ( int j=0 ; (j<n) ; ++j ) { // head
        char tail_in_s = visited[i];
        char head_in_s = visited[j];

        switch (tail_in_s + head_in_s) {
          case 0:
            {
              continue;
            }
          case 1:
            {
              if (head_in_s)
                link[i] += s[x[i][j]];
              else
                link[j] += s[x[i][j]];
              break;
            }
          case 2:
            {
              continue;
            }
        } // switch
      } // head
    } // tail

    int idx_most_connected = -1;
    double link_most_connected = -DBL_MAX;
    for ( int i=0 ; (i<n) ; ++i ) { 
      if (link[i] > link_most_connected) {
        link_most_connected = link[i];
        idx_most_connected = i;
      }
    }
    links_inside += link[idx_most_connected];

    // adding new node
    els[nnodes++] = idx_most_connected;
    visited[idx_most_connected] = 1;

    double violation  = links_inside - (nnodes-1);
    //printf("violation: %g\n", violation);
    if (violation >= MIN_VIOLATION && violation > most_violated) {
      most_violated = violation;
      result = nnodes;
    }
  }

  delete[] visited;
  delete[] link;

  return result;
}

void print_sol(int n, const int **x, const double *s) {
  char col[256], str[256];
  printf("    ");
  for ( int i=0 ; (i<n) ; ++i ) {
    sprintf(str, "%d", i);
    printf("%7s ", str);
  }
  printf("\n");
  for ( int i=0 ; (i<n) ; ++i ) {
    sprintf(str, "%d", i);
    printf("%3s ", str);

    for ( int j=0 ; (j<n) ; ++j ) 
      printf("%.5f ", s[x[i][j]]);
    printf("\n");
  }
  printf("\n");
}


// check in a solution which is the output arc of node i
int out_arc(int n, const double *s, const int **x, int i) {
  for ( int j=0 ; j<n; ++j ) {
    if (s[x[i][j]] >= 0.99)
      return j;
  }
  return -1;
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
