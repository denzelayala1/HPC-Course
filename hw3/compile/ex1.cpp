static char help[] = "Takes a set of vectors and outputs a set of orthonormal vectors with the same span using the Gram-Schmidt proccess.\n\n";

#include <petscksp.h>
#include <string.h>
#include <iostream>
#include <stdlib.h>


PetscErrorCode GramSchmidt(PetscInt n, Vec a[], Vec **q);

int main(int argc, char *argv[]){

   /*  setting up needed data types*/

   PetscInt nVec = 4 , nElem = 4; 
   Vec v[nVec], /* Randomized Input vectors */ 
       u[nVec], /* Orthonormalized output vectors */ 
       * u_ptr = u; /* Pointer to the address of the first vector in the output vector set */ 
   PetscMPIInt size; /* number of processes */
   PetscRandom r; /* random number used to generate input vectors*/


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
   PetscCall(VecSetFromOptions(v[0]));
   PetscCall(VecSetRandom(v[0], r)); 
   PetscCall(VecDuplicate(v[0], &u[0]));


   /* Allocating for all the input and output vectors*/
   for(int k = 1; k < nVec; k++){
      PetscCall(VecDuplicate(v[0], &u[k]));
      PetscCall(VecDuplicate(v[0], &v[k]));
      PetscCall(VecSetRandom(v[k], r)); /* instantiating all the input vectors*/
   }
   
   
   GramSchmidt(nVec, v, &u_ptr);
   //GramSchmidt(nVec, v, &&(*u));

   for(int k = 0; k < nVec; k++){
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,"---------Vector v%d---------\n",k+1));
      PetscCall(VecView(v[k],PETSC_VIEWER_STDOUT_WORLD));
   }

   for(int k = 0; k < nVec; k++){
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,"---------Vector u%d---------\n",k+1));
      PetscCall(VecView(u[k],PETSC_VIEWER_STDOUT_WORLD));
   }

   for(int k = 0; k < nVec; k++){
      PetscCall(VecDestroy(&v[k]));
      PetscCall(VecDestroy(&u[k]));
   }
   
   PetscCall(VecDestroy(u_ptr));
   PetscCall(VecDestroy(u));
   PetscCall(VecDestroy(v));
   PetscCall(PetscRandomDestroy(&r));
   PetscCall(PetscFinalize());

  return 0;
}

/*
* Takes in a set of Linearly independent vectors and outputs a 
* orthonormalized set of vectors with the same span
*
*@param n the number of vectors that are being orthonormalized
*@param a[ ] the input vector set
*@param **q a vector pointer that will contain the orthonormalized vectors
*@return Error Codes
*/
PetscErrorCode GramSchmidt(PetscInt n, Vec a[], Vec **q){

   Vec v[n];
   PetscInt i, j;
   PetscScalar Rij;
      
   for(int k = 0; k < n; k++){
      PetscCall(VecDuplicate(a[0], &v[k]));
   }

   PetscCall(VecCopy(a[0], v[0]));

   PetscCall(VecCopy(v[0], (*q)[0]));
   PetscCall(VecNorm(v[0],NORM_2, &Rij));
   PetscCall(VecScale((*q)[0], 1./Rij));

   for(j = 1; j < n; j++){

      PetscCall(VecCopy(a[j], v[j]));

      for(i=0; i < j; i++){

         PetscCall(VecDot((*q)[i],a[j], &Rij));
         PetscCall(VecAXPY(v[j], -Rij, (*q)[i]));
      }

      PetscCall(VecNorm(v[j],NORM_2, &Rij));
      PetscCall(VecCopy(v[j], (*q)[j]));
      PetscCall(VecScale((*q)[j], 1./Rij));

   }
   
   return 0;
}