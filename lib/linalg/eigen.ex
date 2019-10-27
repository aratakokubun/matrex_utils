defmodule Eigen do
  @moduledoc """
  Module to calculate eigen vector of matrix.
  """

  require Matrex

  @doc """
  Calculate eigen vector of Hermitian matrix.
  """
  def eigh(%Matrex{} = x, upper_triangle \\ true) do
    # TODO: Support false option cases.
    _eigh(x, upper_triangle)
  end

  defp _eigh(%Matrex{} = x, true) do
    # TODO
  end

  def qr(%Matrex{} = x) do
    # TODO
  end
end
