#include <iostream>
#include <list>
#include <vector>

#include "fusion.h"

using namespace mosek::fusion;
using namespace monty;

template< typename MT>
void tsp(int n, MT A,  MT C, bool graph_topology, bool remove_1_hop_loops, bool remove_2_hop_loops)
{
  //TAG:begin-tsp-assignment
  Model::t M = new Model(); 
  auto M_ = finally([&]() { M->dispose(); });

  auto x = M->variable( new_array_ptr<int,1>({ n, n}), Domain::binary());

  M->constraint( Expr::sum(x,0), Domain::equalsTo(1.0));
  M->constraint( Expr::sum(x,1), Domain::equalsTo(1.0));

  M->objective(ObjectiveSense::Minimize, Expr::dot(C ,x) );

  //TAG:end-tsp-assignment

  if( graph_topology)
    //TAG:begin-tsp-topology
    M->constraint( x, Domain::lessThan( A ) );
  //TAG:end-tsp-topology
  
  if( remove_1_hop_loops)
    //TAG:begin-tsp-selfloops
    M->constraint( x->diag(), Domain::equalsTo(0.) );
  //TAG:end-tsp-selfloops
  
  if( remove_2_hop_loops)
    //TAG:begin-tsp-2loops
    M->constraint( Expr::add( x, x->transpose()), Domain::lessThan(1.0));
  //TAG:end-tsp-2loops

  M->setLogHandler([=](const std::string & msg) { std::cout << msg << std::flush; } );
  //TAG:begin-tsp-subtour-search
  while(true)
  {

    M->solve();

    typedef std::vector< std::tuple<int,int> > cycle_t;

    std::list< cycle_t > cycles;

    for(int i=0; i<n;i++)
      for(int j=0; j<n;j++)
      {

        if( (*(x->level()))[i*n+j] <= 0.5 ) 
          continue;

        bool found = false;
        for(auto&& c : cycles)
        {
          for( auto&& cc : c )
          {
            if ( i == std::get<0>(cc) || i == std::get<1>(cc) ||
                 j == std::get<0>(cc) || j == std::get<1>(cc) )
            {
              c.push_back( std::make_tuple(i,j) );
              found = true;
              break;
            }
          }
          if (found) break;
        }
        
        if(!found)
          cycles.push_back( cycle_t(1, std::make_tuple(i,j) ) );
      
      }
    //TAG:end-tsp-subtour-search            
    if (cycles.size()==1) break;
    //TAG:begin-tsp-subtours
    for(auto c : cycles)
    {
      int csize = c.size();

      auto tmp = std::shared_ptr<monty::ndarray<int,2> >(new  ndarray<int,2>( shape(csize,2)) );
      for (auto i = 0; i < csize; ++i) 
      { 
          (*tmp)(i,0) = std::get<0>(c[i]);
          (*tmp)(i,1) = std::get<1>(c[i]);
      }

      M->constraint( Expr::sum(x->pick( tmp )), Domain::lessThan( 1.0 * csize - 1 ) );
    }     
    //TAG:end-tsp-subtours
  }

}

int main()
{
  auto A_i = new_array_ptr<int,1>(   {0 ,1 ,2 ,3,1,0,2,0});
  auto A_j = new_array_ptr<int,1>(   {1 ,2 ,3 ,0,0,2,1,3});

  auto C_v = new_array_ptr<double,1>({1.,1.,1.,1.,0.1,0.1,0.1,0.1});
  
  int n= 4;

  tsp(n, Matrix::sparse(n,n,A_i,A_j, 1.), Matrix::sparse(n,n,A_i,A_j,C_v), true,true,false);
  tsp(n, Matrix::sparse(n,n,A_i,A_j, 1.), Matrix::sparse(n,n,A_i,A_j,C_v), true,true,true);

  return 0;
}