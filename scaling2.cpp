//scaling.cpp
#include <Rcpp.h>
#include <algorithm>    // std::set_intersection, std::sort
#include <math.h>
#include <cmath>
#include <vector>       // std::vector
#include <valarray>     
#include <numeric>
using namespace Rcpp;




struct logt { double operator() (double d) const { return std::log10(d); } };

// [[Rcpp::export]]
List rcpp_scale(NumericMatrix x,List neighbors, int coleng){
  List scf = List::create();
  for (int i = 0; i < coleng; i++){
    NumericVector nodescf=NumericVector::create();
    NumericVector tmp=neighbors[i];
    for(int j = 0; j<tmp.length(); j++){
      NumericVector obs=x( _ , (tmp[j]-1) );
      NumericVector ref=x(_, i);
      //std::valarray<double> obs(ref,ref.size());
      std::vector<double> obs_v=as<std::vector<double> >(obs);
      std::vector<double> ref_v=as<std::vector<double> >(ref);
      double nO = sum(obs);
      double nR = sum(ref);
      

      std::vector<double>::iterator itr_o;
      std::vector<double>::iterator itr_r;

      std::vector<int> o_zero;
      itr_o=std::find(obs_v.begin(),obs_v.end(),0);
      while (itr_o != obs_v.cend())
      {

        // Do something with iter
        o_zero.push_back(std::distance(obs_v.begin(), itr_o));
        itr_o++;
        itr_o = std::find(itr_o, obs_v.end(), 0);
      }

      std::vector<int> r_zero;
      itr_r=std::find(ref_v.begin(),ref_v.end(),0);
      while (itr_r != ref_v.cend())
      {
        // Do something with iter
        r_zero.push_back(std::distance(ref_v.begin(), itr_r));
        itr_r++;
        itr_r = std::find(itr_r, ref_v.end(), 0);
      }

      std::vector<int> all_zero(o_zero.size()+r_zero.size());
      std::vector<int>::iterator it;
      it=std::set_union (o_zero.begin(), o_zero.end(), r_zero.begin(), r_zero.end(), all_zero.begin());
      all_zero.resize(it-all_zero.begin());

      // 0부터 n까지 sequential vector를 만들어요.
      std::vector<int> all(obs_v.size());
      std::iota(all.begin(), all.end(), 0);
      
      // sequential vector에서 0 index를 제거한 vector를 만들어요. test해봤는데 이렇게하는게 zeroindex가지고 for문 돌리는것보다 빠르거나 별차이안나더라고요.
      std::vector<int> not_all_zero(obs_v.size()-all_zero.size());
      std::set_difference(all.begin(), all.end(), all_zero.begin(), all_zero.end(), not_all_zero.begin());

      //nonzeroindex에 해당하는 vector element를 tmp vector에 하나씩 집어넣고 swap해서 원래 vector에 넣어줘요.
      std::vector<double> tmp1;
      std::vector<double> tmp2;
      tmp1.reserve(not_all_zero.size());
      tmp2.reserve(not_all_zero.size());
      for (int i=0; i<not_all_zero.size(); ++i) {
        tmp1.push_back(obs_v[not_all_zero[i]]);
        tmp2.push_back(ref_v[not_all_zero[i]]);

      }
      obs_v.swap(tmp1);
      ref_v.swap(tmp2);
      
      
      
      
      // 얘도 느려요.
      // std::vector<double> tmp1;
      // std::vector<double> tmp2;
      // tmp1.reserve(obs_v.size() - all_zero.size());
      // tmp2.reserve(ref_v.size() - all_zero.size());
      // for (int i=0; i<obs_v.size(); ++i) {
      //   if (std::find(all_zero.begin(),all_zero.end(),i)!=all_zero.cend()) {
      //     tmp1.push_back(obs_v[i]);
      //     tmp2.push_back(ref_v[i]);
      //   }
      // }
      // obs_v.swap(tmp1);
      // ref_v.swap(tmp2);
      //
      
      
      // vector에서 zero index를 제거해줘요. 더느려요.
      // for(int i=0; i<all_zero.size();i++){
      //   obs_v.erase(obs_v.begin()+all_zero[i]-i);
      //   ref_v.erase(ref_v.begin()+all_zero[i]-i);
      // }


      // vector에서 zero index를 제거하는 while문인데 더 오래걸려요. 여기서는 nonzero index나 zero index vector를 따로만들어주지는않아요.
      // while(check==0){
      //   itr_o=std::find(obs_v.begin(),obs_v.end(),0);
      //   itr_r=std::find(ref_v.begin(),ref_v.end(),0);
      //
      //   if(itr_o!=obs_v.cend()){
      //     obs_v.erase(obs_v.begin()+std::distance(obs_v.begin(), itr_o));
      //     ref_v.erase(ref_v.begin()+std::distance(obs_v.begin(), itr_o));
      //   }else if(itr_r!=ref_v.cend()){
      //     obs_v.erase(obs_v.begin()+std::distance(ref_v.begin(), itr_r));
      //     ref_v.erase(ref_v.begin()+std::distance(ref_v.begin(), itr_r));
      //   }else{
      //     check++;
      //   }
      // };
      
      
      
      std::vector<double> trsf(obs_v.size());
      trsf.reserve(obs_v.size());
      // std::cout<<obs_v[0]<<'\t'<<ref_v[0]<<'\t';
      
      
      
      // 두 벡터를 서로 나누어줘요.
      std::transform(obs_v.begin(), obs_v.end(), ref_v.begin(), trsf.begin(), std::divides<double>());
      // std::cout<<trsf[0]<<'\t';
      
      // library size로 나누어줘요.
      std::transform(trsf.begin(), trsf.end(), trsf.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, nR/nO));
      // std::cout<<trsf[0]<<'\t';
      
      // logtransform을 해줘요.
      std::transform(trsf.begin(), trsf.end(), trsf.begin(), logt());
      // std::cout<<trsf[0]<<'\n';

      //mean을 r구하고 pow로 다시 조화평균값을 내줘요.
      double res_mean = std::accumulate(trsf.begin(), trsf.end(), 0.0)/trsf.size();
      res_mean = pow(10,(res_mean/2));
      
      
      nodescf.push_back(res_mean);
      // nodescf.push_back(nO);
    }
    scf.push_back(nodescf);
  }
  return(scf);
}


// 
// // [[Rcpp::export]]
// List test2(NumericMatrix x,List neighbors, int coleng){
//   List scf = List::create();
//   for (int i = 0; i < coleng; i++){
//     NumericVector nodescf=NumericVector::create();
//     NumericVector tmp=neighbors[i];
//     for(int j = 0; j<tmp.length(); j++){
//       NumericVector obs=x( _ , (tmp[j]-1) );
//       NumericVector ref=x(_, i);
//       //std::valarray<double> obs(ref,ref.size());
//       std::vector<double> obs_v=as<std::vector<double> >(obs);
//       std::vector<double> ref_v=as<std::vector<double> >(ref);
//       double nO = sum(obs);
//       double nR = sum(ref);
// 
// 
//       std::vector<double>::iterator itr_o;
//       std::vector<double>::iterator itr_r;
// 
//       std::vector<int> o_zero;
//       itr_o=std::find(obs_v.begin(),obs_v.end(),0);
//       while (itr_o != obs_v.cend())
//       {
// 
//         // Do something with iter
//         o_zero.push_back(std::distance(obs_v.begin(), itr_o));
//         itr_o++;
//         itr_o = std::find(itr_o, obs_v.end(), 0);
//       }
// 
//       std::vector<int> r_zero;
//       itr_r=std::find(ref_v.begin(),ref_v.end(),0);
//       while (itr_r != ref_v.cend())
//       {
//         // Do something with iter
//         r_zero.push_back(std::distance(ref_v.begin(), itr_r));
//         itr_r++;
//         itr_r = std::find(itr_r, ref_v.end(), 0);
//       }
// 
//       std::vector<int> all_zero(o_zero.size()+r_zero.size());
//       std::vector<int>::iterator it;
//       it=std::set_union (o_zero.begin(), o_zero.end(), r_zero.begin(), r_zero.end(), all_zero.begin());
//       all_zero.resize(it-all_zero.begin());
// 
// 
//       std::vector<int> all(obs_v.size());
//       std::iota(all.begin(), all.end(), 0);
// 
//       std::vector<int> not_all_zero(obs_v.size()-all_zero.size());
//       std::set_difference(all.begin(), all.end(), all_zero.begin(), all_zero.end(), not_all_zero.begin());
// 
// 
//       std::vector<double> tmp1;
//       std::vector<double> tmp2;
//       tmp1.reserve(not_all_zero.size());
//       tmp2.reserve(not_all_zero.size());
//       for (int i=0; i<not_all_zero.size(); ++i) {
//         tmp1.push_back(obs_v[not_all_zero[i]]);
//         tmp2.push_back(ref_v[not_all_zero[i]]);
// 
//       }
//       obs_v.swap(tmp1);
//       ref_v.swap(tmp2);
// 
// 
// 
// 
// 
//       std::valarray<double> trsf(x.size());
//       std::valarray<double> obs_a(obs_v.size());
//       std::copy(begin(obs_v), end(obs_v), begin(obs_a));
//       std::valarray<double> ref_a(ref_v.size());
//       std::copy(begin(ref_v), end(ref_v), begin(ref_a));
// 
//       trsf=std::log10(obs_a/ref_a*(nR/nO));
//       // trsf=trsf*(nR/nO);
//       // trsf=std::log10(trsf);
//       double res_mean=trsf.sum()/obs_v.size();
//       res_mean = pow(10,(res_mean/2));
//       nodescf.push_back(res_mean);
//       
//       
//       // std::vector<double> trsf;
//       // trsf.reserve(obs_v.size());
//       // // std::cout<<obs_v[0]<<'\t'<<ref_v[0]<<'\t';
//       // std::transform(obs_v.begin(), obs_v.end(), ref_v.begin(), std::back_inserter( trsf ), std::divides<double>());
//       // // std::cout<<trsf[0]<<'\t';
//       // std::transform(trsf.begin(), trsf.end(), trsf.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, nR/nO));
//       // // std::cout<<trsf[0]<<'\t';
//       // std::transform(trsf.begin(), trsf.end(), trsf.begin(), logt());
//       // // std::cout<<trsf[0]<<'\n';
//       // //
//       // double res_mean = std::accumulate(trsf.begin(), trsf.end(), 0.0)/trsf.size();
//       
//       
//       nodescf.push_back(nO);
//     }
//     scf.push_back(nodescf);
//   }
//   return(scf);
// }
// 
// 
// // // [[Rcpp::export]]
// // std::vector<int> test2(std::vector<int> x) {
// //   std::vector<int> res;
// //   std::vector<int>::iterator iterO = std::find(x.begin(), x.end(), 0);
// //   while (iterO != x.cend())
// //   {
// //     
// //     // Do something with iter
// //     res.push_back(std::distance(x.begin(), iterO));
// //     iterO++;
// //     iterO = std::find(iterO, x.end(), 0);
// //   }
// //   return res;
// // }
// 
// // [[Rcpp::export]]
// void test3(std::vector<double> x) {
//   std::valarray<double> v(x.size());
//   std::copy(begin(x), end(x), begin(v));
//   std::valarray<double> y(x.size());
//   y=v+1.0;
//   std::valarray<double> z(x.size());
//   z=y/v;
//   for(int i=0; i<x.size();i++){
//     std::cout<<z[i]<<'\t';
//   }
//   std::cout<<'\n';
//   z=z*(3.0/2);
//   for(int i=0; i<x.size();i++){
//     std::cout<<z[i]<<'\t';
//   }
//   std::cout<<'\n';
//   z=std::log10(z);
//   for(int i=0; i<x.size();i++){
//     std::cout<<z[i]<<'\t';
//   }
//   std::cout<<'\n';
//   z=std::log10(y/v*(3.0/2));
//   // z=std::log10(z);
//   for(int i=0; i<x.size();i++){
//     std::cout<<z[i]<<'\t';
//   }
//   double res_mean=z.sum()/x.size();
//   std::cout<<'\n'<<res_mean<<'\n';
// }


// // [[Rcpp::export]]
// double test2(std::vector<double> x, std::vector<double> y){
//   int check=0;
//   std::vector<double>::iterator itr_x;
//   std::vector<double>::iterator itr_y;
//   while(check==0){
//     itr_x=std::find(x.begin(),x.end(),0);
//     itr_y=std::find(y.begin(),y.end(),0);
//     
//     if(itr_x!=x.cend()){
//       x.erase(x.begin()+std::distance(x.begin(), itr_x));
//       y.erase(y.begin()+std::distance(x.begin(), itr_x));
//     }else if(itr_y!=y.cend()){
//       x.erase(x.begin()+std::distance(y.begin(), itr_y));
//       y.erase(y.begin()+std::distance(y.begin(), itr_y));
//     }else{
//       check++;
//     }
//   };
// 
//   std::vector<double> dv;
//   dv.reserve(x.size());
//   std::transform(x.begin(), x.end(), y.begin(), std::back_inserter(dv), std::divides<double>());
//   std::transform(dv.begin(), dv.end(), dv.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, 3));
//   std::transform(dv.begin(), dv.end(), dv.begin(), logt());
//   double res_mean = std::accumulate(dv.begin(), dv.end(), 0.0)/dv.size();
//   return(sqrt(pow(10, res_mean)));
// }






// // [[Rcpp::export]]
// List test2(NumericMatrix x,List neighbors, int coleng){
//   List scf = List::create();
//   for (int i = 0; i < coleng; i++){
//     NumericVector nodescf=NumericVector::create();
//     NumericVector tmp=neighbors[i];
//     for(int j = 0; j<tmp.length(); j++){
//       NumericVector obs=x( _ , (tmp[j]-1) );
//       NumericVector ref=x(_, i);
//       //std::valarray<double> obs(ref,ref.size());
//       std::vector<double> obs_v=as<std::vector<double> >(obs);
//       std::vector<double> ref_v=as<std::vector<double> >(ref);
//       // double nO = sum(obs);
//       // double nR = sum(ref);
//       
//       double nOv = std::accumulate(obs_v.begin(), obs_v.end(), 0.0);
//       double nRv = std::accumulate(ref_v.begin(), ref_v.end(), 0.0);
//       
//       //double logR = log2((obs)/(ref));
//       
//       
//       
//       nodescf.push_back(nOv); 
//     }
//     scf.push_back(nodescf);
//   }
//   return(scf);
// }
// 
// 
// // [[Rcpp::export]]
// List test3(NumericMatrix x,List neighbors, int coleng){
//   List scf = List::create();
//   for (int i = 0; i < coleng; i++){
//     NumericVector nodescf=NumericVector::create();
//     NumericVector tmp=neighbors[i];
//     for(int j = 0; j<tmp.length(); j++){
//       NumericVector obs=x( _ , (tmp[j]-1) );
//       NumericVector ref=x(_, i);
//       //std::valarray<double> obs(ref,ref.size());
//       std::vector<double> obs_v=as<std::vector<double> >(obs);
//       std::vector<double> ref_v=as<std::vector<double> >(ref);
//       double nO = sum(obs);
//       double nR = sum(ref);
//       
//       // double nOv = std::accumulate(obs_v.begin(), obs_v.end(), 0.0);
//       // double nRv = std::accumulate(ref_v.begin(), ref_v.end(), 0.0);
//       
//       //double logR = log2((obs)/(ref));
//       
//       
//       
//       nodescf.push_back(nO); 
//     }
//     scf.push_back(nodescf);
//   }
//   return(scf);
// }