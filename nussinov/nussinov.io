Object do(
  memoize  := method(slotName, value, objectOr(getSlot(slotName), setSlot(slotName, value)))
  objectOr := method(if(call evalArgAt(0), call evalArgAt(0), call evalArgAt(1)))
  tap      := method(block, 
  	block call(self)
  	self)  
)

Stack := List clone do(
  push := method(i, insertAt(i, 0))
  pop  := method(removeFirst)
)

if(System args size != 3, Exception raise("io nussinov.io SEQUENCE STRUCTURE"))
  
rna := Object clone do(
  sequence    := System args at(-2)
  structure   := System args at(-1)
  minBpDist   := 3
  length      := method(sequence size)
  seqAt       := method(i, sequence at(i - 1) asCharacter)
  strAt       := method(i, structure at(i - 1) asCharacter)
  canPair     := method(i, j, j - i > minBpDist and basePairs at(seqAt(i)) at(seqAt(j)))
  basePairs   := Map clone lexicalDo(
    atPut("A", Map clone atPut("U", true))
    atPut("U", Map clone atPut("A", true) atPut("G", true))
    atPut("G", Map clone atPut("C", true) atPut("U", true))
    atPut("C", Map clone atPut("G", true)))
  bpList := memoize("calcBpList", List clone lexicalDo(
    stack := Stack clone
    Range clone setRange(1, length) foreach(i, 
      append(
        if(strAt(i) == ".", 0, if(
          strAt(i) == "(", stack push(i), stack pop(i) tap(block(j, atPut(j - 1, i)))))))))
  numStr := memoize("calcNumStr", Range clone setRange(0, length) map(i, 
    0 to(length) asList reduce(list, j, 
      list append(if(i > 0 and j >= i and j - i <= minBpDist, 1, 0)), List clone)) do(
    at2D     := method(i, j, at(i) at(j))
    atPut2D  := method(i, j, value, at(i) atPut(j, value))
    plusEq2D := method(i, j, add, 
      atPut2D(i, j, at2D(i, j) + add))
    ) lexicalDo (
      (minBpDist + 1) to(length - 1) foreach(d, 
        1 to(length - d) foreach(i,
          j := i + d
          atPut2D(i, j, at2D(i, j - 1))
    
          i to(j - minBpDist - 1) foreach(k, 
            if(canPair(k, j), if(k == i,
              plusEq2D(i, j, at2D(k + 1, j - 1)),
              plusEq2D(i, j, at2D(i, k - 1) * at2D(k + 1, j - 1)))))))
              
    ) at2D(1, length))
)

rna numStr println