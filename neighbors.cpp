//scaling.cpp
#include <Rcpp.h>
#include <algorithm>    // std::set_intersection, std::sort
#include <vector>       // std::vector
using namespace Rcpp;



// // [[Rcpp::export]]
// NumericVector rcpp_intersect_old(NumericVector vec1, NumericVector vec2){
//   int i, j =0;
//   NumericVector result =NumericVector::create();
//   while (i!=vec1.size() && j!=vec2.size())
//   {
//     if (vec1[i]<vec2[j]) ++i;
//     else if (vec2[j]<vec1[i]) ++j;
//     else {
//       result.push_back(vec1[i]);
//       ++i; ++j;
//     }
//   }
//   return result;
// }

// [[Rcpp::export]]
std::vector<int> rcpp_intersect(std::vector<int> x, std::vector<int> y) {
  std::vector<int> res(x.size()+y.size());
  std::vector<int>::iterator it;
  
  
  it=std::set_intersection (x.begin(), x.end(), y.begin(), y.end(), res.begin());
  res.resize(it-res.begin()); 
  return res;
}

// [[Rcpp::export]]
double rcpp_jaccard_similarity(const std::vector<int>& vec1, const std::vector<int>& vec2){
  double inter = rcpp_intersect(vec1,vec2).size();
  double res=inter/double(vec1.size()+vec2.size()-inter); 
  return(res);
}

// [[Rcpp::export]]
double rcpp_generalized_jaccard_similarity(std::vector<double> countvec1, std::vector<double> countvec2, const std::vector<int>& vec1, const std::vector<int>& vec2){
  
  
  // std::cout<<vec1.size()<<'\t';
  std::vector <int> inter = rcpp_intersect(vec1,vec2);
  // for(int i=0; i<inter.size();i++){
  //   std::cout<<inter[i]<<'\n';
  // }
  std::vector<double> tmp1;
  std::vector<double> tmp2;
  tmp1.reserve(inter.size());
  tmp2.reserve(inter.size());
  for (int i=0; i<inter.size(); ++i) {
    tmp1.push_back(countvec1[inter[i]-1]);
    tmp2.push_back(countvec2[inter[i]-1]);
    
  }
  
  
  
  
  
  // std::cout<<std::accumulate(tmp1.begin(), tmp1.end(), 0.0)<<'\t'<<std::accumulate(tmp2.begin(), tmp2.end(), 0.0)<<'\n';
  double res= (std::accumulate(tmp1.begin(), tmp1.end(), 0.0)+std::accumulate(tmp2.begin(), tmp2.end(), 0.0))/(std::accumulate(countvec1.begin(), countvec1.end(), 0.0)+std::accumulate(countvec2.begin(), countvec2.end(), 0.0));
  return(res);
}




// // [[Rcpp::export]]
// double rcpp_jaccard_similarity_old(NumericVector vec1, NumericVector vec2){
//   
//   double inter = intersect(vec1,vec2).size();
//   // double inter = rcpp_intersect_old(vec1,vec2).size();
//   double res=inter/double(vec1.size()+vec2.size()-inter); 
//   return(res);
// }

// [[Rcpp::export]]
List rcpp_generalized_jaccard_neighbors(NumericMatrix counts,List non_zero, int coleng){
  List jaccard_neighbors = List::create();
  for (int i = 0; i < coleng; i++){
    NumericVector temp_jaccard_neighbor=NumericVector::create();
    for(int j = 0; j<coleng; j++){
      NumericVector countsi = counts(_,i);
      NumericVector countsj = counts(_,j);
      temp_jaccard_neighbor.push_back(rcpp_generalized_jaccard_similarity(as<std::vector<double>> (countsi), as<std::vector<double>> (countsj) ,non_zero[i],non_zero[j]));
    }
    jaccard_neighbors.push_back(temp_jaccard_neighbor);
  }
  return(jaccard_neighbors);
}


// [[Rcpp::export]]
List rcpp_jaccard_neighbors(List non_zero, int coleng){
  List jaccard_neighbors = List::create();
  for (int i = 0; i < coleng; i++){
    NumericVector temp_jaccard_neighbor=NumericVector::create();
    for(int j = 0; j<coleng; j++){
      temp_jaccard_neighbor.push_back(rcpp_jaccard_similarity(non_zero[i],non_zero[j]));
    }
    jaccard_neighbors.push_back(temp_jaccard_neighbor);
  }
  return(jaccard_neighbors);
}

// // [[Rcpp::export]]
// List rcpp_scf_cal(List neighbors, int coleng){
//   List scf = List::create();
//   for (int i = 0; i < coleng; i++){
//     NumericVector nodescf=NumericVector::create();
//     for(int j = 0; j<neighbors[i].length(); j++){
//       temp_jaccard_neighbor.push_back(rcpp_jaccard_similarity(non_zero[i],non_zero[j]));
//     }
//     scf.push_back()<-nodescf;
//   }
//   return(scf);
// }
// 
// 
// 

