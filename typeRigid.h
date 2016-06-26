#ifndef TYPERIGID_H_
#define TYPERIGID_H_

#include <armadillo>

using namespace std;
using namespace arma;

typedef struct _typeA
{
	mat a;
	mat b;
	mat c;
	mat d;

}typeA;

typedef struct _typeRigid
{
	vector <_typeA> A;
	mat normof_v_Pstar;

}typeRigid;

#endif
