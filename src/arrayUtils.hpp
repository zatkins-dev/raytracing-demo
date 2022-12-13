#pragma once
// -------------------------------------------
// Double 1D Array
// -------------------------------------------
template <class T> T *array(const int n, const T &init = T{}) {
  T *myArray;
  myArray = new T[n];
  for (int i = 0; i < n; ++i)
    myArray[i] = init;
  return myArray;
}

// uses delete to free an allocated array
template <class T> void freeArray(T *array) { delete[] array; }

// uses delete to free an allocated 2D array
// note: rather than freeing each pointer in outer array, we only free the
// first, as that was the memory allocated with new
template <class T> void freeArray(T **array) {
  delete[] array[0];
  delete[] array;
}

// -------------------------------------------
// Double 2D Array
// -------------------------------------------
template <class T>
T **array2d(const int nRows, const int nCols, const T &init = T{}) {
  T *myArray;
  myArray = new T[nRows * nCols];

  // Create a pointer that points to the beginning of each new row

  T **myArray_ptr;
  myArray_ptr = new T *[nRows];

  for (int row = 0; row < nRows; ++row) {
    myArray_ptr[row] = &myArray[row * nCols];
  }

  for (int i = 0; i < nRows; ++i)
    for (int j = 0; j < nCols; ++j)
      myArray_ptr[i][j] = init;

  // Return that pointer

  return myArray_ptr;
}

// -------------------------------------------
// Standard Ax = b
// -------------------------------------------
template <class T>
void matVec(const T *const *A, const T *x, T *b, const int n) {
  for (int i = 1; i <= n; ++i) {
    b[i] = 0;
    for (int j = 1; j <= n; ++j)
      b[i] += A[i][j] * x[j];
  }
}

// -------------------------------------------
// Sparse Ax = b
// -------------------------------------------
template <class T>
void matVec(const T *const *A, const int *const *J, const T *x, T *b,
            const int n) {
  for (int i = 1; i <= n; ++i) {
    b[i] = 0;
    for (int j = 1; j <= n; ++j)
      b[i] += A[i][j] * x[J[i][j]];
  }
}
