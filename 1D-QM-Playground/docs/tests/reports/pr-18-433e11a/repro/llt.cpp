#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
int main(int argc,char**argv){
  std::ifstream f(argv[1]); std::string line; std::vector<std::vector<double>> rows;
  while(std::getline(f,line)){ if(line.empty()||line[0]=='#') continue; std::istringstream is(line); std::vector<double> r; double v; while(is>>v) r.push_back(v); rows.push_back(r);}
  int order=rows.size(), n=rows[0].size(), ka=order-1; Eigen::MatrixXd B=Eigen::MatrixXd::Zero(n,n);
  for(int j=0;j<n;++j) for(int i=std::max(0,j-ka);i<=j;++i){ B(i,j)=rows[ka+i-j][j]; B(j,i)=B(i,j);}
  Eigen::LLT<Eigen::MatrixXd> llt(B); Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(B);
  std::cout<<argv[1]<<"  Eigen LLT info==Success? "<<(llt.info()==Eigen::Success)<<"  min eig "<<es.eigenvalues()(0)<<"\n";
}
