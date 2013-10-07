require "test/unit"
require '../anova.rb'

class TestAnova < Test::Unit::TestCase

  def test_something()
    k = 1
    assert_equal(k,k)
  end

  def test_rsruby()
    l = R.runAnovaTest("df",123)
    assert_equal(l,"df")
  end

  #def test_rsruby_anova()
  #   initialize_anova(groups)
  #  test = run_anova([1,2,1,10,4,12],"rep(\"L\",3),rep(\"K\",3)",6)
  #  assert_equal(test,0.03908894533720031)
  #end

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
    feature, rep, val = format_groups("3,4,5")
    assert_equal(val,3)
    assert_equal(rep,"3,4,5")
    assert_equal(feature,"rep(\"F0\",3),rep(\"F1\",4),rep(\"F2\",5)")
  end

  def test_p_values()
    genes = get_genes("test_data/dm3.gtf")
    all_fpkm_values = all_fpkm(["test_data/short_67.htseq","test_data/short_67.htseq","test_data/short_67.htseq",
      "test_data/short.htseq","test_data/short.htseq","test_data/short.htseq"],genes)
    feature, rep, val = format_groups("3,3")
    l = p_values(all_fpkm_values,feature,val)
    assert_equal(l["CG9996-RB"],6.822293201375393e-63)
  end

  #def test_q_values()
  #  genes = get_genes("test_data/dm3.gtf")
  #  all_fpkm_values = all_fpkm(["test_data/short_67.htseq","test_data/short_67.htseq","test_data/short_67.htseq",
  #    "test_data/short.htseq","test_data/short.htseq","test_data/short.htseq"],genes)
  #  feature, rep, val = format_groups("3,3")
  #  all_p_values = p_values(all_fpkm_values,feature,val)
  #  q_val = q_values(all_p_values)
  #  assert_equal(q_val,10)
  #end

  def test_run()

  end

end