require "test/unit"
require "mocha"
require "./nussinov.rb"

class NussinovTest < Test::Unit::TestCase
  def setup
    @rna = Rnabor::Nussinov.new("gggcc", "(...)")
    @rna.table.size.times { |i| @rna.table[i].length.times { |j| @rna.table[i][j] = rand(100) } }
  end
  
  def test_match_pairs__no_pseudoknots
    assert_equal({
      nil => 1, 
      2   => 9, 
      3   => 4, 
      5   => 8, 
      6   => 7, 
      10  => nil
    }, @rna.match_pairs(")(()(()))("))
  end
  
  def test_match_pairs__with_global_structure
    rna = Rnabor::Nussinov.new("ggggccc", "((...))")
    
    assert_equal({ 1 => 7, 2 => 6 }, rna.match_pairs)
  end
  
  def test_match_pairs__with_global_structure__caches_pairings
    rna = Rnabor::Nussinov.new("ggggccc", "((...))")
    
    rna.expects(:get_pairings).once.with(rna.structure).returns({ 1 => 7, 2 => 6 })
    
    assert_equal({ 1 => 7, 2 => 6 }, rna.match_pairs)
    assert_equal({ 1 => 7, 2 => 6 }, rna.match_pairs)
  end
  
  def test_match_pairs__with_substructure__does_not_cache_pairings
    @rna.expects(:get_pairings).twice.with(" ()").returns({ 1 => 2 })
    
    assert_equal({ 1 => 2 }, @rna.match_pairs("()"))
    assert_equal({ 1 => 2 }, @rna.match_pairs("()"))
  end
  
  def test_match_pairs__short_strings_do_not_call_get_pairings
    @rna.expects(:get_pairings).never
    
    assert_equal({}, @rna.match_pairs(""))
    assert_equal({}, @rna.match_pairs("("))
  end
  
  def test_match_pairs__open_base_pairs_marked_as_nil
    assert_equal({ 1 => nil }, @rna.match_pairs("(."))
    assert_equal({ nil => 2 }, @rna.match_pairs(".)"))
  end
  
  def test_base_paired_in
    assert  Rnabor::Nussinov.new("ggggccc", "((...))").base_paired_in?(7, 1, 7)
    assert !Rnabor::Nussinov.new("ggggccc", "((...))").base_paired_in?(7, 1, 6)
    assert  Rnabor::Nussinov.new("ggggccc", "(((..))").base_paired_in?(6, 3, 7)
    assert !Rnabor::Nussinov.new("ggggccc", "(((..))").base_paired_in?(7, 3, 7)
    assert !Rnabor::Nussinov.new("ggggccc", "(((..))").base_paired_in?(7, 3, 6)
  end
  
  def test_closed_pairs
    assert_equal({
      3 => 4,
      5 => 6
    }, @rna.closed_pairs({
      1   => nil,
      nil => 2,
      3   => 4,
      5   => 6
    }))
  end
  
  def test_count_pairs
    assert_equal(2, @rna.count_pairs({
      1   => nil,
      nil => 2,
      3   => 4,
      5   => 6
    }))
  end
  
  def test_can_pair
    rna = Rnabor::Nussinov.new("acgu", "....")
    
    assert !rna.can_pair?(1, 1)
    assert !rna.can_pair?(1, 2)
    assert !rna.can_pair?(1, 3)
    assert  rna.can_pair?(1, 4)
    
    rna = Rnabor::Nussinov.new("gacu", "....")
    
    assert !rna.can_pair?(1, 1)
    assert !rna.can_pair?(1, 2)
    assert  rna.can_pair?(1, 3)
    assert  rna.can_pair?(1, 4)
    
    assert_equal({
      "a" => %w[u],
      "u" => %w[a g],
      "g" => %w[c u],
      "c" => %w[g]
    }, Rnabor::Nussinov::BASE_PAIRINGS)
  end
  
  def test_paired
    assert !@rna.paired?(0, 100)
    assert !@rna.paired?(1, 4)
    assert  @rna.paired?(1, 5)
  end
  
  def test_generate_table
    assert_equal 5.times.map { Array.new(5) }, Rnabor::Nussinov.new("gggcc", "(...)").table
  end
  
  def test_flush_table
    @rna.flush_table
    
    assert_equal [
      [0.0, 0.0, 0.0, 1.0, 1.0],
      [0.0, 0.0, 0.0, 0.0, 1.0],
      [0.0, 0.0, 0.0, 0.0, 0.0],
      [0.0, 0.0, 0.0, 0.0, 0.0],
      [0.0, 0.0, 0.0, 0.0, 0.0]
    ], @rna.table
  end
  
  def test_table_at
    assert_equal @rna.table.first.last, @rna.table_at(1, 5)
    
    @rna.table_at(1, 5, :value)
    
    assert_equal :value, @rna.table.first.last
  end
end