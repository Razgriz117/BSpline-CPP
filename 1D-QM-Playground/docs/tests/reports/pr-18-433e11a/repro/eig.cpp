#include <Eigen/Eigenvalues>
#include <Eigen/Cholesky>
#include <iostream>
extern "C" { void dsbgv_(char*,char*,int*,int*,int*,double*,int*,double*,int*,double*,double*,int*,double*,int*); }
void test(const std::string &label, Eigen::MatrixXd A, Eigen::MatrixXd B){
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> s(A,B,Eigen::ComputeEigenvectors|Eigen::Ax_lBx);
  Eigen::LLT<Eigen::MatrixXd> llt(B);
  // LAPACK dsbgv on the same dense-as-banded (full band kd=n-1) matrices
  int n=A.rows(), ka=n-1, kb=n-1, ldab=n, ldbb=n, info=0; char jobz='N', uplo='U';
  std::vector<double> AB(n*n,0), BB(n*n,0), W(n), Z(1), work(3*n);
  for(int j=0;j<n;++j) for(int i=std::max(0,j-ka);i<=j;++i){ AB[j*ldab+(ka+i-j)]=A(i,j); BB[j*ldbb+(kb+i-j)]=B(i,j);}
  int ldz=1; dsbgv_(&jobz,&uplo,&n,&ka,&kb,AB.data(),&ldab,BB.data(),&ldbb,W.data(),Z.data(),&ldz,work.data(),&info);
  std::cout<<label<<": GSAES info==Success? "<<(s.info()==Eigen::Success)<<"  LLT(B).info==Success? "<<(llt.info()==Eigen::Success)
           <<"  dsbgv info="<<info<<"  eigenvalues: "<<s.eigenvalues().transpose()<<"\n";
}
int main(){ using Eigen::MatrixXd; MatrixXd I3=MatrixXd::Identity(3,3);
 { MatrixXd B(3,3); B.setZero(); B(0,0)=1; B(1,1)=-1;    B(2,2)=1; test("1 indefinite      ",I3,B);}
 { MatrixXd B(3,3); B.setZero(); B(0,0)=1; B(1,1)=0;     B(2,2)=1; test("2 exactly singular",I3,B);}
 { MatrixXd B(3,3); B.setZero(); B(0,0)=1; B(1,1)=1e-14; B(2,2)=1; test("3 near-singular PD",I3,B);}
 { MatrixXd B(3,3); B.setZero(); B(0,0)=1; B(1,1)=-1e-10;B(2,2)=1; test("4 tiny negative   ",I3,B);}
 { MatrixXd B(3,3); B<<1,1,0, 1,1,0, 0,0,1;                       test("5 rank-deficient  ",I3,B);}
 { MatrixXd B(3,3); B<<1,1,0, 1,1,0, 0,0,1; B(1,1)+=1e-16;        test("6 rank-def + noise",I3,B);}
 { MatrixXd B(3,3); B<<1,1,0, 1,1,0, 0,0,1; B(1,1)+=1e-10;        test("7 rank-def + 1e-10",I3,B);}
}
