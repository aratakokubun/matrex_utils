defmodule MatrexUtilsTest do
  use ExUnit.Case
  doctest MatrexUtils

  test "greets the world" do
    assert MatrexUtils.hello() == :world
  end
end
