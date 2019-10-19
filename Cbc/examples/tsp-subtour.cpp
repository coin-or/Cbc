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

using namespace std;

typedef struct 
{
  int n;          // number of cities
  int **c = NULL; // distance matrix
  int **x; // indexes of x variables
} TSPInfo;

// receves info of tsp instance and solution, returns size of
// subtour and fills elements in els
int find_subtour( const TSPInfo *tspi, const double *s, int *els );

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

  model.branchAndBound();

  printf("Branch & Bound Finished with status: %d\n", model.status());
  printf("isProvenOptimal: %d\n", model.isProvenOptimal());

  if (model.getColSolution()) 
  {
    int *els = new int[n];
    int sts = find_subtour(&info, model.getColSolution(), els);
    if (sts) 
    {
      printf("subtour with %d elements found: ", sts);
      for ( int i=0 ; (i<sts) ; ++i )
        printf(" %d", els[i]);
      printf("\n");
    }
    delete[] els;

    const double *s= model.getColSolution();

    for ( int i=0 ; i<mip->getNumCols() ; ++i )
      if (fabs(s[i]) >= 1e-5)
        printf("  %s: %g\n", mip->getColName(i).c_str(), s[i]);
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

  printf("max frac value: %g\n", maxidist);

  // processing only integer solutions to include
  // really lazy constraints
  if (maxidist >= si.getIntegerTolerance()+1e-7)
    return;

  int *idx = new int[n*n];
  double *coef = new double[n*n];
  fill(coef, coef+n*n, 1.0);

  int *els = new int[n];

  int nnodes = find_subtour(tspi, s, els);
  char *visited = new char[n];
  fill(visited, visited+n, 0);
  for ( int i=0 ; (i<nnodes) ; ++i )
    visited[els[i]] = 1;

  if (nnodes)
  {
    printf("subtour found with nodes:");
    for ( int i=0 ; (i<nnodes) ; ++i )
      printf(" %d", els[i]);
    printf("\n");

    OsiRowCut cut;
    int nz = 0;
    for ( int i=0 ; (i<n) ; ++i )
      for (int j=0 ; (j<n) ; ++j )
        if (i!=j && visited[i] && visited[j])
          idx[nz++] = x[i][j];
  
    cut.setRow(nz, idx, coef);
    cut.setLb(-COIN_DBL_MAX);
    cut.setUb(nnodes-1);
    cs.insert(cut);

    for ( int i=0 ; (i<nz) ; ++i )
      printf("%+g %s", coef[i], si.getColName(idx[i]).c_str());
    printf("<= %d\n", nnodes-1); fflush(stdout);
    printf("");
  }

  delete[] visited;
  delete[] els;
  delete[] idx;
  delete[] coef;
}

int find_subtour( const TSPInfo *tspi, const double *s, int *els )
{
  int n = tspi->n;
  const int **x = (const int **)tspi->x;
  char *visited = new char[n];
  fill(visited+1, visited+n, 0);
  visited[0] = 1;

  char has_sub = 0;
  int nnodes = 1;
  int node = 0;
  while (true)
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
