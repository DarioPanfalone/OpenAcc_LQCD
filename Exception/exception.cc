
#ifndef EXCEPTION_CC_
#define EXCEPTION_CC_

#include <exception>
#include <cstdlib>

using namespace std;

class myexception0: public exception
{
  virtual const char* what() const throw()
  {
    return "\033[31mERROR: configuration file '" QUOTEME(CONF_FILE) "' does not exists\033[0m\n";
  }
} no_stored_conf;

class myexception1: public exception
{
  virtual const char* what() const throw()
  {
    return "\033[31mERROR: stored configuration has wrong dimensions\033[0m\n";
  }
} stored_conf_not_fit;

class myexception2: public exception
{
  virtual const char* what() const throw()
  {
    return "\033[31mERROR: multistep integrator cannot be used in PURE_GAUGE simulations\033[0m\n";
  }
} no_multistep;

class myexception3: public exception
{
  virtual const char* what() const throw()
  {
    return "\033[31mERROR: nx, ny, nz, nt have to be even to use odd/even preconditioning\033[0m\n";
  }
} odd_sides;

class myexception4: public exception
{
  virtual const char* what() const throw()
  {
  return "\033[31mERROR: volume must be divisible by 2*NUM_THREADS\033[0m\n";
  }
} wrong_num_threads;


void test_param(void)
  {
  #ifdef DEBUG_MODE
  cout << "DEBUG: inside test_param..."<<endl;
  #endif

  if(nx%2 + ny%2 +nz%2 +nt%2 >0)
    {
    throw odd_sides;
    }

  #ifdef PURE_GAUGE
  if(use_multistep!=0) throw no_multistep;
  #endif     

  #ifdef USE_GPU
  if(size % (2*NUM_THREADS)!=0) throw wrong_num_threads;
  #endif

  #ifdef DEBUG_MODE
  cout << "\tterminated test_param..."<<endl;
  #endif
  }

#endif
