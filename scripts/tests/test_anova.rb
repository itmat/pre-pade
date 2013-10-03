require "test/unit"
require '../anova.rb'

class TestAnova < Test::Unit::TestCase

  def test_something()
    k = 1
    assert_equal(k,k)
  end

  def test_rsruby()
    l = R.runAnova("df",123)
    assert_equal(l,"df")
  end

  def test_calc_length()
    starts = [1,5,9]
    ends = [2,7,14]
    length = calc_length(starts,ends)
    assert_equal(length,8)
  end

  def test_get_genes()
    genes = get_genes("test_data/short.gtf")
    assert_equal(genes[0].name, "CG12449-RB")
    assert_equal(genes[1].starts, [2505272, 2505874, 2517270])
    assert_equal(genes[1].length, 913)
  end

  def test_run()

  end

end