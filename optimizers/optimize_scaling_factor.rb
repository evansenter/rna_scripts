# GGGGGCCCCC
# ..........
# k Density   NumStr  Energy    MFE(k)
# 0 0.001613  1 0.000000  ..........
# 1 0.000784  19  1.700000  (......)..
# 2 0.061878  35  -1.400000 ((....))..
# 3 0.935725  7 -3.900000 (((....)))
# 4 0.000000  0 0.000000  
# 5 0.000000  0 0.000000  
# 6 0.000000  0 0.000000  
# 7 0.000000  0 0.000000  
# 8 0.000000  0 0.000000

# CACUUCAACCGAUCGCGGAA
# ....................
# k Density   NumStr  Energy    MFE(k)
# 0 0.006919  1 0.000000  ....................
# 1 0.002220  42  1.100000  ........(........)..
# 2 0.084929  309 -1.500000 ........((......))..
# 3 0.902458  494 -3.000000 ........(((....)))..
# 4 0.002893  167 1.100000  ...(((..(.....)..)))
# 5 0.000579  6 1.800000  ...((...(((....)))))
# 6 0.000000  0 0.000000  
# 7 0.000000  0 0.000000  
# 8 0.000000  0 0.000000  
# 9 0.000000  0 0.000000  
# 10  0.000000  0 0.000000

# GGGGGCCCCCGGGGGCCCCC
# ....................
# k Density   NumStr  Energy    MFE(k)
# 0 0.000000  1 0.000000  ....................
# 1 0.000000  82  0.800000  ..(......)..........
# 2 0.000000  1511  -2.500000 ((......))..........
# 3 0.000000  7173  -5.600000 (((....)))..........
# 4 0.000000  11221 -8.300000 ...((((......))))...
# 5 0.000001  6781  -11.600000  ..(((((......)))))..
# 6 0.000226  1800  -14.900000  .((((((......)))))).
# 7 0.023589  212 -18.000000  .(((((((....))))))).
# 8 0.976183  7 -20.500000  ((((((((....))))))))
# 9 0.000000  0 0.000000  
# 10  0.000000  0 0.000000  
# 11  0.000000  0 0.000000  
# 12  0.000000  0 0.000000  
# 13  0.000000  0 0.000000

# GGGGGCCCCCGGGGGCCCCCGGGGGCCCCC
# ..............................
# k Density   NumStr  Energy    MFE(k)
# 0 0.000000  1 0.000000  ..............................
# 1 0.000000  195 0.700000  ..........(........)..........
# 2 0.000000  9806  -2.800000 ..........((......))..........
# 3 0.000000  179779  -5.900000 ..........(((....)))..........
# 4 0.000000  1290383 -8.300000 ...((((......)))).............
# 5 0.000000  4046555 -11.600000  ..(((((......)))))............
# 6 0.000000  6094496 -15.600000  ((((((........))))))..........
# 7 0.000000  4747469 -19.100000  (((((((......)))))))..........
# 8 0.000000  2056623 -22.200000  ((((((((....))))))))..........
# 9 0.000000  531318  -23.900000  ...(((((((((......)))))))))...
# 10  0.000001  84703 -27.200000  ..((((((((((......))))))))))..
# 11  0.000220  8038  -30.500000  .(((((((((((......))))))))))).
# 12  0.023705  387 -33.600000  .((((((((((((....)))))))))))).
# 13  0.976074  7 -36.100000  (((((((((((((....)))))))))))))
# 14  0.000000  0 0.000000  
# 15  0.000000  0 0.000000  
# 16  0.000000  0 0.000000  
# 17  0.000000  0 0.000000  
# 18  0.000000  0 0.000000

require "./nussinov.rb"

# (->(length) { Array.new(3, length).zip([10.361989255770736, -0.22104738360459786, -0.0009107162603511926]).each_with_index.inject(0) { |sum, (x_coefficient, index)| sum + x_coefficient.last * x_coefficient.first ** index } })[100]

# ruby-1.9.2-p290 :195 > Lagrange.new([10, 8.06044379368964], [20, (5.2333516011638945 + 5.920158557912711) / 2], [30, 2.9109231133167293]).coefficients
#  => [10.361989255770736, -0.22104738360459786, -0.0009107162603511926]

# best_constant = 5.2333516011638945 # Best for "cacuucaaccgaucgcggaa" 6..-1 for 250 consecutive tries (1965.5x improvement)
# best_constant   = 5.920158557912711 # Best for "GGGGGCCCCCGGGGGCCCCC" 9..-1 for 250 consecutive tries
# best_constant   = 8.06044379368964 # Best for "GGGGGCCCCC" 4..-1 for 1000 consecutive tries
# best_constant   = 2.9109231133167293 # Best for "GGGGGCCCCCGGGGGCCCCCGGGGGCCCCC" 14..-1 for 35 consecutive tries

best_constant   = 5.920158557912711
held_for, count = 0, 0
rmsd            = ->(values) { Math.sqrt(values.inject(0) { |sum, value| sum + value ** 2 } / values.count) }
best_rmsd       = rmsd[Rnabor::Nussinov.new(best_constant, "GGGGGCCCCCGGGGGCCCCC").partition_function[9..-1]]

while count < 250
  scaling_factor = (rand(2).zero? ? 1 : -1) * rand * 1.5 ** (-1 * held_for + 1)
  new_constant   = best_constant + scaling_factor
  new_rmsd       = rmsd[Rnabor::Nussinov.new(new_constant, "GGGGGCCCCCGGGGGCCCCC").partition_function[9..-1]]
  
  if new_rmsd < best_rmsd
    best_constant = new_constant
    best_rmsd     = new_rmsd
    held_for      = 0
    count         = 0
  else
    count += 1
    
    if rand < 0.1
      held_for = 0
    else
      held_for += 1
    end
  end
  
  puts "*" * 50
  puts "count:          #{count}"
  puts "held_for:       #{held_for}"
  puts "scaling_factor: #{scaling_factor}"
  puts "new_constant:   #{new_constant}"
  puts "best_constant:  #{best_constant}"
  puts "best_rmsd:      #{best_rmsd}"
  puts "*" * 50
end