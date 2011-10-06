#include <stdio.h>

int main(void) {
  #define ROWS 2
  #define COLUMNS 3
  int i;
  double data[ROWS * COLUMNS];
  #define M(row, column) data[(column - 1) * ROWS + (row - 1)]
  
  M(1, 1) = 1;
  M(1, 2) = 2;
  M(1, 3) = 3;
  M(2, 1) = 4;
  M(2, 2) = 5;
  M(2, 3) = 6;
  
  for (i = 0; i < ROWS * COLUMNS; i++) {
    printf("%f\n", data[i]);
  }
}