#include <stdio.h>
#include <string.h>
#include <math.h>

int pairedIn(int i, int j, int k, char *structure) {
  if (j - i < 4) {
    return 0;
  } else {
    int counter, stackPointer;
    int pairs[j + 1];
    int stack[j - i + 1];

    stackPointer = 0;
    
    for (counter = i; counter <= j; ++counter) {
      pairs[counter] = 0;
    }

    for (counter = i; counter <= j; ++counter) {
      if (structure[counter] == '(') {
        stack[stackPointer++] = counter;
      } else if (structure[counter] == ')' && stackPointer > 0) {
        stackPointer -= 1;
        pairs[stack[stackPointer]] = counter;
        pairs[counter]             = stackPointer;
      }
    }
        
    return pairs[k];
  }
}

void populatePairings(int i, int j, char *structure, int *pairs) {
  int structurePointer, stackPointer;
  int stack[j - i + 1];
  
  stackPointer = 0;
  
  for (structurePointer = i; structurePointer <= j; ++structurePointer) {
    if (structure[structurePointer] == '(') {
      stack[stackPointer++] = structurePointer;
    } else if (structure[structurePointer] == ')' && stackPointer > 0) {
      pairs[stack[--stackPointer]] = structurePointer;
    }
  }
}

int comparePairings(int i, int j, int *pairs, int *comparisonPairs) {
  int counter, pairings;
  
  pairings = 0;
  
  for (counter = i; counter <= j; ++counter) {
    if (pairs[counter] && pairs[counter] != comparisonPairs[counter]) {
      pairings += 1;
    }
  }
  
  for (counter = i; counter <= j; ++counter) {
    if (comparisonPairs[counter] && comparisonPairs[counter] != pairs[counter]) {
      pairings += 1;
    }
  }
  
  return pairings;
}

int bpDifference(int i, int j, int k, int l, char *structure, int bounds[][2], int numBounds) {
  int counter;
  int pairs[j + 1];
  int comparisonPairs[j + 1];
  
  for (counter = i; counter <= j; ++counter) {
    pairs[counter]           = 0;
    comparisonPairs[counter] = 0;
  }
  
  populatePairings(i, j, structure, pairs);
  
  for (counter = 0; counter < numBounds; ++counter) {
    populatePairings(bounds[counter][0], bounds[counter][1], structure, comparisonPairs);
  }
  
  if (k && l) {
    comparisonPairs[k] = l;
  }
  
  return comparePairings(i, j, pairs, comparisonPairs);
}

int main(int argc, char *argv[]) {
  printf("%d\n", pairedIn(2, 8, 6, "012345(7)"));
  
  int i, array[1][2];
  
  array[0][0] = 4;
  array[0][1] = 7;
  
  printf("%d\n", bpDifference(0, 9, 0, 0, "(()()())())", array, 1));
}