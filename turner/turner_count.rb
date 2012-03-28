module Rnabor
  class TurnerCount
    MIN_LOOP_SIZE = 3
    MAX_LOOP_SIZE = 30 # Restrictions not in effect on multiloops
    BASE_PAIRINGS = {
      "a" => %w[u],
      "u" => %w[a g],
      "g" => %w[c u],
      "c" => %w[g]
    }

    attr_reader :length, :sequence, :table_n, :table_b, :table_m, :table_mrc

    def initialize(sequence)
      @length    = sequence.length
      @sequence  = " " + sequence.downcase
      @table_n   = initialize_table { |i, j| i <= j ? (j - i <= MIN_LOOP_SIZE ? 1 : 0) : nil }
      @table_b   = initialize_table { |i, j| 0 if i <= j }
      @table_m   = initialize_table { |i, j| 0 if i <= j }
      @table_mrc = initialize_table { |i, j| 0 if i <= j }
    end
    
    def structure_count
      solve_recurrences[1][-1]
    end

    def solve_recurrences
      ((MIN_LOOP_SIZE + 1)..(length - 1)).each do |d|
        (1..(length - d)).each do |i|
          j = i + d
          
          solve_table_b(i, j)
          solve_table_m(i, j)
          solve_table_mrc(i, j)
          solve_table_n(i, j)
        end
      end
        
      table_n
    end
    
    def solve_table_n(i, j)
      # Main table
      
      table_n[i][j] = table_n[i][j - 1]
      
      i.upto(j - MIN_LOOP_SIZE - 1) do |k|
        if can_pair?(k, j)
          
          if k == i
            set(:n, i, j, table_b[k][j])
          else
            set(:n, i, j, table_n[i][k - 1] * table_b[k][j])
          end 
        end
      end
    end
    
    def solve_table_b(i, j)
      if can_pair?(i, j)
        # Closing hairpin
        set(:b, i, j, 1) if j - i - 1 <= MAX_LOOP_SIZE
        
        # Stacked base pair
        set(:b, i, j, table_b[i + 1][j - 1])
        
        # Left bulge
        (i + 2).upto(j - 1 - MIN_LOOP_SIZE - 1) do |k|
          if can_pair?(k, j - 1) && k - i - 1 <= MAX_LOOP_SIZE
            # k base pairs with j - 1 to form a left bulge where i + 1 < k < j - 1 - MIN_LOOP_SIZE
            set(:b, i, j, table_b[k][j - 1])
          end
        end
        
        # Right bulge
        (j - 2).downto(i + 1 + MIN_LOOP_SIZE + 1) do |k|
          if can_pair?(k, j - 1) && j - k - 1 <= MAX_LOOP_SIZE
            # k base pairs with i + 1 to form a right bulge where i + 1 + MIN_LOOP_SIZE < k < j - 1
            set(:b, i, j, table_b[i + 1][k])
          end
        end
        
        # Interior loop
        (i + 2).upto(j - 1 - MIN_LOOP_SIZE - 1) do |k|
          if k - i - 1 <= MAX_LOOP_SIZE - 1
            (k + MIN_LOOP_SIZE + 1).upto(j - 2) do |l|
              if can_pair?(k, l) && j - l - 1 + k - i - 1 <= MAX_LOOP_SIZE
                # k base pairs with l where i + 1 < k < j - 1 - MIN_LOOP_SIZE and k + MIN_LOOP_SIZE < l < j - 1
                set(:b, i, j, table_b[k][l])
              end
            end
          end
        end
        
        # Multiloop
        (i + 1 + MIN_LOOP_SIZE + 2).upto(j - 1 - MIN_LOOP_SIZE - 1) do |k|
          # i and j close a multiloop with at least 2 components, k marks left closing base of rightmost component
          # with i + 2 + MIN_LOOP_SIZE < k < j - 1 - MIN_LOOP_SIZE
          set(:b, i, j, table_m[i + 1][k - 1] * table_mrc[k][j - 1])
        end
      end
    end
    
    def solve_table_m(i, j)
      # Table assuming i, j are *in* a multiloop (they don't close a multiloop, that's [i - 1, j + 1])
      i.upto(j - MIN_LOOP_SIZE - 1) do |k|
        # There is only one remaining component in [i, j], and it has a left closing base k: i <= k < j - MIN_LOOP_SIZE
        set(:m, i, j, table_mrc[k][j])
      end
      
      (i + MIN_LOOP_SIZE + 2).upto(j - MIN_LOOP_SIZE - 1) do |k|
        # There is more than one remaining component in [i, j], and the rightmost component has a left closing base k: 
        # i <= k < j - MIN_LOOP_SIZE
        set(:m, i, j, table_m[i][k - 1] * table_mrc[k][j])
      end
    end
    
    def solve_table_mrc(i, j)
      # Table assuming i is the left closing base of the rightmost component of a multiloop and the component is closed in [i, j]
      (i + MIN_LOOP_SIZE + 1).upto(j) do |k|
        # i is the leftmost closing base of the only multiloop component in [i, j], and closed within i + MIN_LOOP_SIZE < k <= j
        set(:mrc, i, j, table_b[i][k])
      end
    end
    
    def can_pair?(i, j)
      j - i - 1 >= MIN_LOOP_SIZE && BASE_PAIRINGS[sequence[i]].include?(sequence[j])
    end
    
    def initialize_table
      (0..length).map { |i| (0..length).map { |j| yield i, j } }
    end
    
    def set(table_name, i, j, value)
      raise ArgumentError.new("#{table_name}: (#{i} > #{j})") if i > j
      
      instance_variable_get(:"@table_#{table_name}")[i][j] += value
    end
  end
end

if ARGV.empty?
  puts "Call: ruby ./turner_count.rb [SEQUENCE]"
else
  puts Rnabor::TurnerCount.new(ARGV.first).structure_count
end