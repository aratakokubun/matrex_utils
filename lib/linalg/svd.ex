defmodule Svd do
  @moduledoc """
    Module to calculate SVD, singular value decomposition of matrix.
  """

  require Matrex

  def svd(%Matrex{} = x, full_matrices \\ true, compute_uv \\ true, hermitian \\ true) do
    # TODO: Support false option cases.
    _svd(x, full_matrices, compute_uv, hermitian)
  end

  defp _svd(%Matrex{} = x, true, true, true) do
    # TODO
  end
end