#include <RcppArmadillo.h>
#include <boost/math/tools/minima.hpp>
#include <boost/math/tools/toms748_solve.hpp>

#include <functional> // functional
// using namespace std::placeholders;

// [[Rcpp::export]]
double myFun (const double& x, const arma::vec& singular_values, const double& degrees_of_freedom)
{
  // return exp(x) - std::pow(x,2);
  // return exp(x) - std::pow(x, 3) / 3;
  return 2 * arma::accu(1 / (1 + x * singular_values)) - arma::accu(1 / arma::pow(1 + x * singular_values, 2)) - degrees_of_freedom;
  // return 2 * arma::accu(arma::log(1 + x * singular_values) / singular_values) + arma::accu(1 / (singular_values + x * arma::pow(singular_values, 2))) - x * degrees_of_freedom;
}

struct myFunObj
{
    myFunObj (arma::vec singular_values, double degrees_of_freedom) 
    : singular_values (singular_values), 
      degrees_of_freedom (degrees_of_freedom) { }
    // you may need a copy constructor and operator= here too ...

    double operator()(double x) 
    {  
        return myFun(x, singular_values, degrees_of_freedom);
    }
    arma::vec singular_values; 
    double degrees_of_freedom;
};

// [[Rcpp::export]]
double optimBrent (const arma::vec& singular_values, const double& degrees_of_freedom, const double& lower_bound = 0.,
  const double& upper_bound = 1e15) 
{ 
  // int bits = std::numeric_limits<double>::digits;
  // std::pair<double, double> r = boost::math::tools::brent_find_minima(myFunObj(singular_values, degrees_of_freedom), lower_bound, upper_bound, bits);

  boost::uintmax_t max_iter=500;
  boost::math::tools::eps_tolerance<double> tol(30);
 
  std::pair<double, double> r = boost::math::tools::toms748_solve(std::bind(myFun, std::placeholders::_1, singular_values, degrees_of_freedom), lower_bound, upper_bound, tol, max_iter);

  std::cout << "lower: " << r.first << "; upper: " << r.second << std::endl;

  return r.first;
}