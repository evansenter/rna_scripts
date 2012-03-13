module Rnabor
  class NussinovCount
    MIN_LOOP_SIZE = 3
    BASE_PAIRINGS = {
      "a" => %w[u],
      "u" => %w[a g],
      "g" => %w[c u],
      "c" => %w[g]
    }

    attr_reader :length, :sequence, :table

    def initialize(sequence)
      @length   = sequence.length
      @sequence = " " + sequence.downcase
      @table    = (0..length).map { |i| (0..length).map { |j| i <= j ? 1 : 0 } }
    end
    
    def structure_count
      solve_recurrences[1][-1]
    end

    def solve_recurrences
      ((MIN_LOOP_SIZE + 1)..(length - 1)).each do |d|
        (1..(length - d)).each do |i|
          j = i + d
          
          table[i][j] = table[i][j - 1]

          (i..(j - MIN_LOOP_SIZE - 1)).each do |k|
            if can_pair?(k, j)
              if k == i
                table[i][j] += table[k + 1][j - 1]
              else
                table[i][j] += table[i][k - 1] * table[k + 1][j - 1]
              end 
            end
          end
        end
      end
        
      table
    end

    def can_pair?(i, j)
      j > i + MIN_LOOP_SIZE && BASE_PAIRINGS[sequence[i]].include?(sequence[j])
    end
  end
end

if ARGV.empty?
  puts "Call: ruby ./nussinov_count.rb [SEQUENCE]"
else
  puts Rnabor::NussinovCount.new(ARGV.first).structure_count
end