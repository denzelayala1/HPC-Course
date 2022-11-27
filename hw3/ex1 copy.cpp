static char help[] = "Takes a set of vectors and outputs a set of orthonormal vectors with the same span using the Gram-Schmidt proccess.\n\n";

#include <petscksp.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>

PetscErrorCode PrintThis(PetscInt i, Vec a[], Vec **u);


int main(int argc, char *argv[]){

   /*  setting up needed data types*/

   PetscInt nVec = 3 , nElem = 5; 
   Vec v[nVec], u[nVec], * u_ptr = u;
   PetscMPIInt size;
   PetscRandom r;


   PetscFunctionBeginUser;
   PetscCall(PetscInitialize(&argc, &argv, (char *)0, help));
   PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));
   //PetscCheck(size == 1, PETSC_COMM_WORLD, PETSC_ERR_WRONG_MPI_SIZE, "This is a uniprocessor example only!");

   PetscCall(PetscOptionsGetInt(NULL, NULL, "-n_Vec", &nVec, NULL));
   PetscCall(PetscOptionsGetInt(NULL, NULL, "-n_Elem", &nElem, NULL));
   

   PetscCall(PetscRandomCreate(PETSC_COMM_WORLD, &r));
   PetscCall(PetscRandomSetType(r,PETSCRAND48));  
   

   /* Creating the vectors that will be passed to the GS-algorithm*/
   PetscCall(VecCreate(PETSC_COMM_WORLD, &v[0]));
   PetscCall(VecSetSizes(v[0], PETSC_DECIDE, nElem));
   PetscCall(VecSetType(v[0],VECSTANDARD));
   PetscCall(VecDuplicate(v[0], &u[0]));
   PetscCall(VecSetRandom(v[0], r));
    //PetscCall(VecView(v[0],PETSC_VIEWER_STDOUT_WORLD));

   for(int k = 1; k < nVec; k++){
      PetscCall(VecDuplicate(v[0], &v[k]));
      PetscCall(VecDuplicate(v[0], &u[k]));
      PetscCall(VecSetRandom(v[k], r));
   }

   /*for(int k = 0; k < nVec; k++){
      PetscCall(VecSetRandom(v[k], r));
   }*/

   for(int k = 0; k < nVec; k++){
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,"---------Vector v%d---------\n",k+1));
      PetscCall(VecView(v[k],PETSC_VIEWER_STDOUT_WORLD));
   }

   for(int k = 0; k < nVec; k++){
      PrintThis(k, v, &u_ptr);
   }

   
   
   for(int k = 0; k < nVec; k++){
      PetscCall(VecDestroy(&v[k]));
   } 
   //PetscCall(VecDestroy(&v));
   PetscCall(PetscRandomDestroy(&r));
   PetscCall(PetscFinalize());


  return 0;
}


PetscErrorCode PrintThis(PetscInt i, Vec a[], Vec **u){
   
   PetscCall(VecCopy(a[i], (*u)[i]));
   PetscCall(PetscPrintf(PETSC_COMM_WORLD,"---------Vector *u%d---------\n",i+1));
   PetscCall(VecView((*u)[i],PETSC_VIEWER_STDOUT_WORLD));

   return 0;
}