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

#define p 10   // �W�c�̌�
#define m 21    // �e�W�c�ɂ�����̂̐�(>=n+1)
#define n 10    // �T������p�����[�^�̐�
#define q 11    // �e�̂̐�
#define alpha 1 // ������
#define beta 21 // ������

  void generate_sample();       // �p�����[�^�T���v���̃����_������
  void ObjectiveFunction(double *fx, double **x, int s);     // �ړI�֐��l�̌v�Z
  void SortOrderIncrease(double *fx, double **x, int s);     // �ړI�֐��l�̏����ɂ��e�̂̕��בւ�
  void GroupDivision(double ***A, double **x);         // �W�c�̕���
  void Probably(double *prob);              // �I���m���̌v�Z
  void SelectParent(double ***A, double **U, double *prob, int *selectq, int k);     // �I���m���ɏ]����q�̐e�̂�񕜌����o����D 
  void CalcGravPoint(double *G, double **U);         //�ړI�֐��l�̏d�S�̌v�Z
  void CalObjectiveFunction(double *fx, double **x, int s);  // �ړI�֐��l�̌v�Z
  double ObjectiveFunction2(double *R);     // �ړI�֐�
  void CalcRandMituation(double *Ux, double **U, double *G, double *paramin, double *paramax);     // �q��R=2G-Uq�̌v�Z�ƓˑR�ψ�
  void ReChild(double ***A, double **U, int *selectq, int k);           // �q�̂ɂ��u���������q�̌̂�Ak�ɖ߂�
  void SortAk(double ***A, int k);                // Ak�Ɋ܂܂��m�̌̂̕��ёւ�
void MixAk(double **x, double ***A);          // k���ׂĂ̏W�c�̍���


// �ړI�֐��l�̌v�Z
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


// �l�̏������בւ�
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

// �W�c�̕���
void GroupDivision(double ***A, double **x){
  int k,i,j;
  for(k=1;k<=p;k++){
    for(i=1;i<=m;i++){
      for(j=1;j<=n;j++){
        A[k][i][j] = x[k+p*(i-1)][j];
  }}} 
}

// �I���m���̌v�Z
void Probably(double *prob){
  int i;
  prob[0] = 0.;
  for(i=1;i<=m;i++){
    prob[i] = prob[i-1] + 2.*(m+1.-i)/(m*(m+1.));
  }
}
  
// �I���m���ɏ]����q�̌̂�񕜌����o����D 
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

// �d�SG�̌v�Z
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

// �q��R=2G-Uq�̌v�Z�ƓˑR�ψ�
void CalcRandMituation(double *Ux, double **U, double *G, double *paramin, double *paramax){
  int check = 0;
  double *R;    // �q�̂̔z��
  double *z;    // �V�����q�̂̔z��
  double fR;
  double fUq;
  int i,j;

  R = new double[n+1];            // �q�̂̃p�����[�^�̔z��
  z = new double[n+1];            // �ˑR�ψَq�̂̃p�����[�^�̔z��

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

// �q�̂ɂ��u���������q�̌̂�Ak�ɖ߂�
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

// Ak�Ɋ܂܂��m�̌̂̕��ёւ�
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

// k���ׂĂ̏W�c�̍���
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

  double *paramax;  // n�̃p�����[�^�̍ő�l
  double *paramin;  // n�̃p�����[�^�̍ŏ��l
  double **x;   // �p�����[�^n������s=p*m�̌̂̔z��
  double *fx;   // s=p*m�̌̖̂ړI�֐��l�̔z��
  double ***A;  // s�̌̂�m�̌̂��܂�p�̏W�c�̔z��
  double *prob; // �I���m���̔z��
  double **U;   // �p�����[�^n������q�̐e�̂̔z��
  double *Ux;   // q�̐e�̖̂ړI�֐��l�̔z��
  double *G;    // q-1�̐e�̂̏d�S�̔z��
  int s;        // �̐�����(p*m)
  int loop;     // ���[�v��
  int i,j,k,ii,jj,kk;
  int *selectq; // ���o�̂̃t���O

  // �z��̊m��
  A = new double**[p+1];          // s�̌̂�m�̌̂��܂�p�̏W�c�̔z��
  for(k=0;k<=p;k++){
    A[k] = new double*[m+1];
    for(i=0;i<=m;i++){
      A[k][i] = new double[n+1];
  }}
  x = new double*[p*m+1];         // �p�����[�^n������s=p*m�̌̂̔z��
  fx = new double[p*m+1];         // s=p*m�̌̖̂ړI�֐��l�̔z��
  for(i=0;i<=p*m;i++){
    x[i] = new double[n+1];
  }
  U = new double*[q+1];           // �p�����[�^n������q�̐e�̂̔z��
  Ux = new double[q+1];           // q�̐e�̖̂ړI�֐��l�̔z��
  selectq = new int[q+1];         // q�̐e�̂̃t���O�̔z��
  for(i=0;i<=q;i++){
    U[i] = new double[n+1];
  }
  prob = new double[m+1];         // �I���m���̔z��
  G = new double[n+1];            // q-1�̐e�̂̏d�S�̔z��
  paramax = new double[n+1];      // �p�����[�^�̍ő�l
  paramin = new double[n+1];      // �p�����[�^�̍ŏ��l

  s=p*m;                          // ���̐��̌v�Z
  // �p�����[�^�̍ő�l�ƍŏ��l�̐ݒ�
  for(i=1;i<=n;i++){
    //F1,F4,F5�p
    //paramax[i]=5.12;
    //paramin[i]=-5.12;
    //F2�p
    //paramax[i]=65.536;
    //paramin[i]=-65.536;
    //F3�p
    //paramax[i]=2.048;
    //paramin[i]=-2.048;
    //F6�p
    paramax[i]=512.;
    paramin[i]=0.;
    //F7,F8�p
    //paramax[i]=512.;
    //paramin[i]=-512.;
  }

  // �p�����[�^�T���v���̃����_������
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

  // �ړI�֐��l�̌v�Z
  CalObjectiveFunction(fx,x,s);

  // �l�̏������בւ�
  SortOrderIncrease(fx,x,s);

  loop = 0;       // ���[�v�񐔂̏�����
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


