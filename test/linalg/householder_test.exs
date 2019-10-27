defmodule HouseHolderTest do
  use ExUnit.Case
  require MatrexUtils.HouseHolder
  doctest MatrexUtils.HouseHolder

  test "create_householder returns housholder vector and norm" do
    assert MatrexUtils.HouseHolder.convert2householder_vec(Matrex.new([[1, 0, 0]]))     == {Matrex.new([[2, 0, 0]]), 4}
    assert MatrexUtils.HouseHolder.convert2householder_vec(Matrex.new([[1, 2, 3]]))     == {Matrex.new([[15, 2, 3]]), 238}
    assert MatrexUtils.HouseHolder.convert2householder_vec(Matrex.new([[-1, -2, -3]]))  == {Matrex.new([[-15, -2, -3]]), 238}
  end

  test "ev22 returns expected value" do
    assert MatrexUtils.HouseHolder.ev22(0, 0, 0, 0) == 0
    assert MatrexUtils.HouseHolder.ev22(2, 1, 1, 2) == 1
    assert MatrexUtils.HouseHolder.ev22(2, 2, 2, 4) == 5.23606797749979
  end
end
