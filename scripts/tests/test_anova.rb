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

  def test_HTSseq()
    htseq_obj = HTSseq.new("test_data/short.htseq")
    htseq_obj.read()
    assert_equal(htseq_obj.number_of_fragments,2767)
  end

  def test_all_fpkm()
    genes = get_genes("test_data/dm3.gtf")
    all_fpkm_values = all_fpkm(["test_data/short_67.htseq",
      "test_data/short.htseq"],genes)
    assert_equal(all_fpkm_values["CG9999-RA"],[0.02536311658242118, 0.03320909286790296])
  end

  def test_format_groups()
    feature, rep = format_groups("3,4,5")
    assert_equal(feature[0],"F0")
    assert_equal(rep[0],3)
  end

  def test_run()

  end

end