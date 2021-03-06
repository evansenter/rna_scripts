require "set"
require "awesome_print"

module Rnabor
  class Neighbors
    MIN_LOOP_SIZE = 3
    BASE_PAIRINGS = {
      "a" => %w[u],
      "u" => %w[a g],
      "g" => %w[c u],
      "c" => %w[g]
    }

    attr_reader :length, :sequence, :structure, :table

    def initialize(sequence, structure = nil)
      @length    = sequence.length
      @sequence  = (" " + sequence.downcase).freeze
      @structure = (" " + (structure || "." * length)).freeze
      @table     = generate_table
    end

    def count_neighbors
      flush_table
      
      (1..length).each do |delta|
        ((MIN_LOOP_SIZE + 1)..(length - 1)).each do |distance|
          (1..(length - distance)).each do |i|
            j = i + distance

            table[delta][i][j] = table[delta - (end_base_paired?(i, j) ? 1 : 0)][i][j - 1]
            
            (i..(j - MIN_LOOP_SIZE - 1)).select { |k| can_pair?(k, j) }.each do |k|              
              base_pair_distance = pair_distance(i, k, j)
              if k == i
                table[delta][i][j] += table[delta - base_pair_distance][k + 1][j - 1] || 0
              else
                (0..(delta - base_pair_distance)).each do |upstream_distance|
                  table[delta][i][j] += table[upstream_distance][i][k - 1] * table[delta - base_pair_distance - upstream_distance][k + 1][j - 1]
                end
              end
            end
          end
        end
      end
      
      puts sequence.strip
      puts structure.strip
      
      table.map { |delta_table| delta_table[1].last }.each_with_index do |neighbors, index|
        puts "k = %-10.10sneighbors = %s" % [index, neighbors]
      end
    end
    
    def pair_distance(i, k, j)
      if (@memoized_pair_distance ||= {})[[i, k, j]]
        @memoized_pair_distance[[i, k, j]]
      else
        reference_structure  = match_pairs(structure[i..j]).map { |from, to| [from + i - 1, to + i - 1] }.to_set
        upstream             = match_pairs(k - 1 < i ? "" : structure[i..(k - 1)]).map { |from, to| [from + i - 1, to + i - 1] }
        downstream           = match_pairs(structure[(k + 1)..(j - 1)]).map { |from, to| [from + k, to + k] }
        comparitive_pairings = (upstream + downstream + [[k, j]]).to_set

        @memoized_pair_distance[[i, k, j]] = ((reference_structure - comparitive_pairings) + (comparitive_pairings - reference_structure)).size
      end
    end

    def match_pairs(structure_to_match = structure)
      structure_to_match = " " + structure_to_match unless structure_to_match[0] == " "
      
      if structure_to_match.length >= MIN_LOOP_SIZE + 3 # Padded space + minimum loop size + base pair on either side of loop
        get_pairings(structure_to_match)
      else
        {}
      end
    end
    
    def get_pairings(structure)
      structure.split(//).each_with_index.inject({}) do |hash, (symbol, index)|
        hash.tap do      
          case symbol
          when "(" then hash[index] = nil
          when ")" then hash[hash.select { |from, to| to.nil? }.keys.max] = index
          end
        end
      end
    end
    
    def end_base_paired?(i, j)
      closed_pairs(match_pairs(structure[i..j])).values.compact.include?(j - i + 1)
    end

    def closed_pairs(pair_hash)
      pair_hash.reject { |key, value| key.nil? || value.nil? }
    end

    def can_pair?(i, j)
      BASE_PAIRINGS[sequence[i]].include?(sequence[j])
    end

    def generate_table
      (0..length).map { (0..length).map { Array.new(length + 1) } }
    end
    
    def flush_table
      (0..length).each do |delta|
        (1..length).each do |i|
          (i..length).each do |j|
            if delta == 0
              table[delta][i][j] = 1
            else
              table[delta][i][j] = 0 if j <= i + MIN_LOOP_SIZE
            end
          end
        end
      end
    end
  end
end

# Rnabor::Neighbors.new("gggcc").count_neighbors
# Rnabor::Neighbors.new("gggggccccc").count_neighbors
# Rnabor::Neighbors.new("gggggcccccgggggccccc").count_neighbors
# Rnabor::Neighbors.new("cacuucaaccgaucgcggaa").count_neighbors