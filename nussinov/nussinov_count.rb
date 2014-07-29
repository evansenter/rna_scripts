module Rnabor
  class NussinovCount
    MIN_LOOP_SIZE = 3
    INT_CODES     = {
      "A" => 2,
      "U" => 3,
      "G" => 4,
      "C" => 5
    }
    
    BASE_PAIRINGS = %w|au ua cg gc gu ug|.map(&:upcase).map { |bp| bp.split(//).map { |i| INT_CODES[i] } }.map(&:sort).uniq.map { |i, j| (2 ** i) * (3 ** j) }

    attr_reader :length, :sequence, :table

    def initialize(sequence)
      @length   = sequence.length
      @sequence = [-1].concat(sequence.upcase.split(//).map { |i| INT_CODES[i] })
      @table    = (0..length).map { |i| (0..length).map { |j| i <= j ? 1 : 0 } }
    end
    
    def structure_count
      solve_recurrences[1][-1]
    end

    def solve_recurrences
      ((MIN_LOOP_SIZE + 1)..(length - 1)).each do |d|
        (1..(length - d)).each do |i|
          j = i + d
          
          table[i][j]  = table[i][j - 1]
          table[i][j] += table[i + 1][j - 1] if can_pair?(i, j)

          ((i + 1)..(j - MIN_LOOP_SIZE - 1)).each do |k|
            table[i][j] += table[i][k - 1] * table[k + 1][j - 1] if can_pair?(k, j)
          end
        end
      end
        
      table
    end

    def can_pair?(i, j)
      sequence[i] ^ sequence[j] != 0x0 && BASE_PAIRINGS.include?(sequence[i] < sequence[j] ? (2 ** sequence[i]) * (3 ** sequence[j]) : (2 ** sequence[j]) * (3 ** sequence[i]))
    end
  end
end

if ARGV.empty?
  puts "Call: ruby ./nussinov_count.rb [SEQUENCE]"
else
  puts Rnabor::NussinovCount.new(ARGV.first).structure_count
end