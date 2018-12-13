/**
*
* Copyright (C) 2004, 2007, 2011 by the Gascoigne 3D authors
*
* This file is part of Gascoigne 3D
*
* Gascoigne 3D is free software: you can redistribute it
* and/or modify it under the terms of the GNU General Public
* License as published by the Free Software Foundation, either
* version 3 of the License, or (at your option) any later
* version.
*
* Gascoigne 3D is distributed in the hope that it will be
* useful, but WITHOUT ANY WARRANTY; without even the implied
* warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
* PURPOSE.  See the GNU General Public License for more
* details.
*
* Please refer to the file LICENSE.TXT for further information
* on this license.
*
**/


#include  "simpleilu.h"


using namespace std;

/* ----------------------------------------- */

namespace Gascoigne
{
void SimpleIlu::ReInit(int n, int nentries)
{
  SimpleMatrix::ReInit(n,nentries);
  p.reservesize(n,-1);
  q.reservesize(n,-1);
}

/*-----------------------------------------*/

void SimpleIlu::solve(DoubleVector& x) const
{
  hin(x);
  forward ();
  backward();
  her(x);
}

/*-----------------------------------------*/

void SimpleIlu::solve_transpose(DoubleVector& x) const
{
  hin(x);
  forward_transpose ();
  backward_transpose();
  her(x);
}

/*-------------------------------------------------------------*/

void SimpleIlu::hin(const DoubleVector& x) const
{
  int n = x.size();
  yp.reservesize(n);
  assert(n==ST.n());
  for(int i=0;i<n;i++)  yp[i] = x[p[i]];
}

/*-------------------------------------------------------------*/

void SimpleIlu::her(DoubleVector& x) const
{
  for(int i=0;i<ST.n();i++)  x[i] = yp[q[i]];
}

/*-------------------------------------------------------------*/

void SimpleIlu::forward() const
{
  for(int i=1; i<ST.n(); i++)
    {
      int ende = ST.diag(i); 
      for(int pos=ST.start(i); pos<ende; pos++)
	{
	  int j = ST.col(pos);
	  yp[i] -= value[pos]*yp[j];
	}
    }
}

/*-------------------------------------------------------------*/

void SimpleIlu::backward() const
{
  for(int i=ST.n()-1; i>=0; i--)
    {
      int ende = ST.diag(i); 
      for(int pos=ST.stop(i)-1; pos>ende; pos--)
        {
          int j = ST.col(pos);
          yp[i] -= value[pos]*yp[j];
        }
      yp[i]  *= value[ende];
    }
}

/*-------------------------------------------------------------*/

void SimpleIlu::forward_transpose() const
{
  for(int i=0; i<ST.n(); i++)
    {
      int ende = ST.diag(i); 
      for(int pos=ST.start(i); pos<ende; pos++)
        {
          int j = ST.col(pos);
          int pos2 = ST.Find(j,i);
          yp[i] -= value[pos2]*yp[j];
        }
      yp[i]  *= value[ende];
    }
}

/*-------------------------------------------------------------*/

void SimpleIlu::backward_transpose() const
{
  for(int i=ST.n()-1; i>=0; i--)
    {
      int ende = ST.diag(i); 
      for(int pos=ST.stop(i)-1; pos>ende; pos--)
        {
          int j = ST.col(pos);
          int pos2 = ST.Find(j,i);
          yp[i] -= value[pos2]*yp[j];
        }
    }
}

/* ----------------------------------------- */
  
void SimpleIlu::compute_ilu()
{
  // // original
  // for(int i=0; i<ST.n(); i++)
  //   {
  //     for (int pk=ST.start(i); pk<ST.diag(i); pk++)
  //       {
  //         int k = ST.col(pk);
  // 
  //         value[pk] *= value[ST.diag(k)];
  // 
  //         for (int pj=ST.diag(k)+1; pj<ST.stop(k); pj++)
  //           {
  //             int j  = ST.col(pj);
  //             // suche ph
  //             for (int ph=ST.start(i); ph<ST.stop(i); ph++)
  //               {
  //                 if (ST.col(ph)==j)
  //                   {
  //                     value[ph] -= value[pk]*value[pj];
  //                     break;
  //                   }
  //               }
  //           }
  //       }
  //     double d = value[ST.diag(i)];
  //     value[ST.diag(i)] = 1./d;
  //   }
  // ca. 10% Verbesserung der insg. Geschwindigkeit durch das Anlegen von 
  // lokalen Referenzen und Werten 
  // und die Vermeidung von redundanten Aufrufen, tom.
  // (ein "inline"-ing der Methoden von ST (Class ColumnStencil) bringt nichts)
  int ST_n     = ST.n();
  IntVector & ST_start = ST.start();
  IntVector & ST_diag = ST.diag();
  IntVector & ST_col  = ST.col();
  int ST_start_i;
  int ST_stop_i;
  int ST_diag_i;
  int ST_diag_k;
  int ST_stop_k;
  double value_pk;
  for(int i=0; i<ST_n; i++)
    {
      assert(i<ST_start.size());
      ST_start_i = ST_start[i];

      assert(i+1<ST_start.size());
      ST_stop_i = ST_start[i+1]; // stop(i) == start(i+1)

      assert(i<ST_diag.size());
      ST_diag_i = ST_diag[i];

      for (int pk=ST_start_i; pk<ST_diag_i; pk++)
        {
          assert(pk<ST_col.size());
          int k = ST_col[pk];

          assert(k<ST_diag.size());
          ST_diag_k = ST_diag[k];

          value_pk   = value[pk];
          value_pk  *= value[ST_diag_k];

          assert(k+1<ST_diag.size());
          ST_stop_k = ST_start[k+1]; // stop(i) == start(i+1) 

          for (int pj=ST_diag_k+1; pj<ST_stop_k; pj++)
            {
              assert(pj<ST_col.size());
              int j  = ST_col[pj];

              // suche ph
              for (int ph=ST_start_i; ph<ST_stop_i; ph++)
                {
                  if (ST_col[ph]==j)
                    {
                      value[ph] -= value_pk *value[pj];
                      break;
                    }
                }
            }
          value[pk]  = value_pk;
        }
      double d = value[ST_diag_i];
      value[ST_diag_i] = 1./d;
    }
}

/*-------------------------------------------------*/

void SimpleIlu::copy_entries(const MatrixInterface*  A)
{
  const SimpleMatrix* AP = dynamic_cast<const SimpleMatrix*>(A);
  assert(AP);

  const ColumnDiagStencil* AS = dynamic_cast<const ColumnDiagStencil*>(AP->GetStencil());
  assert(AS);

  for(int i=0;i<ST.n();i++)
    {
      int pi = p[i];

      for(int posA=AS->start(pi); posA<AS->stop(pi); posA++)
        {
          int j   = AS->col(posA);
          int pj  = q[j];
          bool found=0;
          for(int pos=ST.start(i); pos<ST.stop(i); pos++)
            {
              int k = ST.col(pos);
              if(k==pj)	
                {
                  value[pos] += AP->GetValue(posA);
                  found=1;
                  break;
                }
            }
          if(!found)
            {
              cout << "not found " << endl;
              abort();
            }
        }
    }
}
}
