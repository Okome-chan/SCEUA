/********************************************/
// SCE-UA Method
// Developed by Ichiro Yoneda
// version 0.0.1
// June 21st 2020
/********************************************/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <math.h>
#include <time.h>

using namespace std;

#define p 10   // 集団の個数
#define m 21    // 各集団における個体の数(>=n+1)
#define n 10    // 探索するパラメータの数
#define q 11    // 親個体の数
#define alpha 1 // 反復回数
#define beta 21 // 反復回数

  void generate_sample();       // パラメータサンプルのランダム生成
  void ObjectiveFunction(double *fx, double **x, int s);     // 目的関数値の計算
  void SortOrderIncrease(double *fx, double **x, int s);     // 目的関数値の昇順による各個体の並べ替え
  void GroupDivision(double ***A, double **x);         // 集団の分割
  void Probably(double *prob);              // 選択確率の計算
  void SelectParent(double ***A, double **U, double *prob, int *selectq, int k);     // 選択確率に従ってq個の親個体を非復元抽出する． 
  void CalcGravPoint(double *G, double **U);         //目的関数値の重心の計算
  void CalObjectiveFunction(double *fx, double **x, int s);  // 目的関数値の計算
  double ObjectiveFunction2(double *R);     // 目的関数
  void CalcRandMituation(double *Ux, double **U, double *G, double *paramin, double *paramax);     // 子個体R=2G-Uqの計算と突然変異
  void ReChild(double ***A, double **U, int *selectq, int k);           // 子個体により置き換わったq個の個体をAkに戻す
  void SortAk(double ***A, int k);                // Akに含まれるm個の個体の並び替え
void MixAk(double **x, double ***A);          // k個すべての集団の混合


// 目的関数値の計算
void CalObjectiveFunction(double *fx, double **x, int s){
  double *calf;
  int i,j;
  calf = new double[n+1];
  
  for(i=1;i<=s;i++){
    for(j=1;j<=n;j++){
      calf[j] = x[i][j];
    }
    fx[i] = ObjectiveFunction2(calf);
  }
  delete[] calf;
}

double ObjectiveFunction2(double *R){
  double fR=0.;
  int i;
  //F1 Sphere function
  /*for(i=1;i<=n;i++){
    fR+=R[i]*R[i];
  }*/
  //F2 Ridge function
  /*double temp;
  for(i=1;i<=n;i++){
    temp=0.;
    for(int j=1;j<=i;j++){
      temp+=R[j];
    }
    fR+=temp*temp;
  }*/
  //F3 Rosenbrock function
  /*for(i=2;i<=n;i++){
    fR+=100.*(R[1]-R[i]*R[i])*(R[1]-R[i]*R[i])+(R[i]-1.)*(R[i]-1.);
  }*/
  //F4 Bohachevsky function
  /*for(i=1;i<=n-1;i++){
    fR+=R[i]*R[i]+2.*R[i+1]*R[i+1]-0.3*cos(3.*M_PI*R[i])-0.4*cos(4.*M_PI*R[i])+0.7;
  }*/
  //F5 Rastrigin function
  /*for(i=1;i<=n;i++){
    fR+=R[i]*R[i]-10.*cos(2.*M_PI*R[i]);
  }
  fR+=10.*n;*/
  //F6 Schwefel function
  for(i=1;i<=n;i++){
    fR-=R[i]*sin(sqrt(abs(R[i])));
    //cout<<std::scientific <<R[i]<<" "<<abs(R[i])<<endl;
  }
  fR+=418.98288728*n;
  //F7 Rastrigin function
  /*double temp=1.;
  for(i=1;i<=n;i++){
      fR+=R[i]*R[i]/4000.;
  }
  for(i=1;i<=n;i++){
    temp*=cos(R[i]/sqrt(i));
  }
  fR=fR-temp+1.;*/
  //F8 Griewank-d function
  /*double temp=1.;
  for(i=1;i<=n;i++){
      fR+=(R[i]-100.)*(R[i]-100.)/4000.;
  }
  for(i=1;i<=n;i++){
    temp*=cos((R[i]-100.)/sqrt(i));
  }
  fR=fR-temp+1.;*/


  //exit(0);

  return(fR);
}


// 値の昇順並べ替え
void SortOrderIncrease(double *fx, double **x, int s){
  double tmp;
  int i,j,k;
  for(i=1;i<=s;i++){
    for(j=i+1;j<=s;j++){
      if(fx[i]>fx[j]){
        tmp = fx[i];
        fx[i] = fx[j];
        fx[j] = tmp;
        for(k=1;k<=n;k++){
          tmp = x[i][k];
          x[i][k] = x[j][k];
          x[j][k] =tmp;
      }}
  }}
}

// 集団の分割
void GroupDivision(double ***A, double **x){
  int k,i,j;
  for(k=1;k<=p;k++){
    for(i=1;i<=m;i++){
      for(j=1;j<=n;j++){
        A[k][i][j] = x[k+p*(i-1)][j];
  }}} 
}

// 選択確率の計算
void Probably(double *prob){
  int i;
  prob[0] = 0.;
  for(i=1;i<=m;i++){
    prob[i] = prob[i-1] + 2.*(m+1.-i)/(m*(m+1.));
  }
}
  
// 選択確率に従ってq個の個体を非復元抽出する． 
void SelectParent(double ***A, double **U, double *prob, int *selectq, int k){
  
  double randomN;
  int count=1;
  int check=0;
  int i,j,ii;
  for(i=0;i<=q;i++){
    selectq[i] = 0;
  }

  while(count<=q){
    randomN = (double)rand()/RAND_MAX;
    for(i=1;i<=m;i++){
      if(randomN<=prob[i] && randomN>prob[i-1]){
        check = 0;
        for(ii=1;ii<count;ii++){
          if(selectq[ii]==i) check = 1;
        }
        if(check == 0){
          for(j=1;j<=n;j++){
            U[count][j] = A[k][i][j];
            ////cout<<U[count][j]<<" ";
          }
          selectq[count] = i;
          count++;
        }
      }
    }
  }
}

// 重心Gの計算
void CalcGravPoint(double *G, double **U){
  int i,j;
  for(j=1;j<=n;j++){
    G[j] = 0.;
  }
  for(i=1;i<q;i++){
    for(j=1;j<=n;j++){
      G[j] += U[i][j];
  }}
  for(j=1;j<=n;j++){
    G[j] = G[j]/((double)q-1.);
  }
}

// 子個体R=2G-Uqの計算と突然変異
void CalcRandMituation(double *Ux, double **U, double *G, double *paramin, double *paramax){
  int check = 0;
  double *R;    // 子個体の配列
  double *z;    // 新しい子個体の配列
  double fR;
  double fUq;
  int i,j;

  R = new double[n+1];            // 子個体のパラメータの配列
  z = new double[n+1];            // 突然変異子個体のパラメータの配列

  for(j=1;j<=n;j++){
    R[j] = 2.*G[j]-U[q][j];
    if(R[j]<paramin[j] || paramax[j]<R[j]){check=1;}
  }
  if(check==0){
    fR = ObjectiveFunction2(R);
  }else{
    for(j=1;j<=n;j++){
      z [j] = (double)rand()/(double)RAND_MAX*(paramax[j]-paramin[j]) + paramin[j];
      R[j] = z[j];
    }
    fR = ObjectiveFunction2(z);
  }
  if(fR < Ux[q]){
    Ux[q] = fR;
    for(i=1;i<=n;i++){
      U[q][i] = R[i];
    }
  }else{
    for(j=1;j<=n;j++){
      R [j] = (G[j]+U[q][j])/2.;
    }
    fR = ObjectiveFunction2(R);
    if(fR < Ux[q]){
      Ux[q] = fR;
      for(i=1;i<=n;i++){
        U[q][i] = R[i];
      }
    }else{
      for(j=1;j<=n;j++){
        z [j] = (double)rand()/(double)RAND_MAX*(paramax[j]-paramin[j]) + paramin[j];
      }
      Ux[q] = ObjectiveFunction2(z);
      for(i=1;i<=n;i++){
        U[q][i] = z[i];
      }
    }
  }
  delete[] R;
  delete[] z;
}

// 子個体により置き換わったq個の個体をAkに戻す
void ReChild(double ***A, double **U, int *selectq, int k){
  int i,j,ii;
  int count=0;
  for(i=1;i<=m;i++){
    for(ii=1;ii<=m;ii++){
     if(selectq[ii]==i){
       count++;
        for(j=1;j<=n;j++){
          A[k][i][j] = U[count][j];
        }
      }
    }
  }
}

// Akに含まれるm個の個体の並び替え
void SortAk(double ***A, int k){
  double *func;
  double **xx;
  int i,j;
  func = new double[m+1];
  xx = new double*[m+1];
  for(i=0;i<=m+1;i++){
    xx[i] = new double[n+1];
  }

  for(i=1;i<=m;i++){
    for(j=1;j<=n;j++){
      xx[i][j] = A[k][i][j];
  }}
  CalObjectiveFunction(func,xx,m);
  SortOrderIncrease(func,xx,m);
  for(i=1;i<=m;i++){
    for(j=1;j<=n;j++){
      A[k][i][j] = xx[i][j];
    }
  }
  for(i=0;i<=m+1;i++){
    delete[] xx[i];
  }
  delete[] xx;
  delete[] func;
}

// k個すべての集団の混合
void MixAk(double **x, double ***A){
  int i,j,k;
  for(k=1;k<=p;k++){
    for(i=1;i<=m;i++){
      for(j=1;j<=n;j++){
        x[k*m-(m-i)][j] = A[k][i][j];
      }
  }}
}

int main (){

  double *paramax;  // n個のパラメータの最大値
  double *paramin;  // n個のパラメータの最小値
  double **x;   // パラメータn個をもつs=p*m個の個体の配列
  double *fx;   // s=p*m個の個体の目的関数値の配列
  double ***A;  // s個の個体をm個の個体を含むp個の集団の配列
  double *prob; // 選択確率の配列
  double **U;   // パラメータn個をもつq個の親個体の配列
  double *Ux;   // q個の親個体の目的関数値の配列
  double *G;    // q-1個の親個体の重心の配列
  int s;        // 個体生成数(p*m)
  int loop;     // ループ回数
  int i,j,k,ii,jj,kk;
  int *selectq; // 抽出個体のフラグ

  // 配列の確保
  A = new double**[p+1];          // s個の個体をm個の個体を含むp個の集団の配列
  for(k=0;k<=p;k++){
    A[k] = new double*[m+1];
    for(i=0;i<=m;i++){
      A[k][i] = new double[n+1];
  }}
  x = new double*[p*m+1];         // パラメータn個をもつs=p*m個の個体の配列
  fx = new double[p*m+1];         // s=p*m個の個体の目的関数値の配列
  for(i=0;i<=p*m;i++){
    x[i] = new double[n+1];
  }
  U = new double*[q+1];           // パラメータn個をもつq個の親個体の配列
  Ux = new double[q+1];           // q個の親個体の目的関数値の配列
  selectq = new int[q+1];         // q個の親個体のフラグの配列
  for(i=0;i<=q;i++){
    U[i] = new double[n+1];
  }
  prob = new double[m+1];         // 選択確率の配列
  G = new double[n+1];            // q-1個の親個体の重心の配列
  paramax = new double[n+1];      // パラメータの最大値
  paramin = new double[n+1];      // パラメータの最小値

  s=p*m;                          // 総個体数の計算
  // パラメータの最大値と最小値の設定
  for(i=1;i<=n;i++){
    //F1,F4,F5用
    //paramax[i]=5.12;
    //paramin[i]=-5.12;
    //F2用
    //paramax[i]=65.536;
    //paramin[i]=-65.536;
    //F3用
    //paramax[i]=2.048;
    //paramin[i]=-2.048;
    //F6用
    paramax[i]=512.;
    paramin[i]=0.;
    //F7,F8用
    //paramax[i]=512.;
    //paramin[i]=-512.;
  }

  // パラメータサンプルのランダム生成
  srand(time(NULL));
  rand();
  rand();
  rand();
  rand();
  rand();
  rand();
  for(i=1;i<=s;i++){
    for(j=1;j<=n;j++){
      x[i][j] = (double)rand()/(double)RAND_MAX*(paramax[j]-paramin[j]) + paramin[j];
  }}

  // 目的関数値の計算
  CalObjectiveFunction(fx,x,s);

  // 値の昇順並べ替え
  SortOrderIncrease(fx,x,s);

  loop = 0;       // ループ回数の初期化
  while(fx[1]>10E-8 && loop<3000){
    loop++;
    GroupDivision(A,x);
    Probably(prob);
    for(k=1;k<=p;k++){
      for(i=1;i<=beta;i++){
        SelectParent(A,U,prob,selectq,k);
        CalObjectiveFunction(Ux,U,q);
        for(j=1;j<=alpha;j++){
          SortOrderIncrease(Ux,U,q);
          CalcGravPoint(G,U);
          CalcRandMituation(Ux,U,G,paramin,paramax);
        }
        ReChild(A,U,selectq,k);
        SortAk(A,k);
      }
    }
    MixAk(x,A);
    CalObjectiveFunction(fx,x,s);
    SortOrderIncrease(fx,x,s);  
    cout<<2*s*loop<<" "<<loop<<" check1 "<<fx[1]<<" "<<fx[s]<<endl;
  }
  cout<<"Finish Optimization "<<loop<<" "<<fx[1]<<endl;
  for(i=1;i<=n;i++){
    cout<<i<<" "<<x[1][i]<<endl;
  }

  return 0;
}


