defmodule HouseHolder do
  require Matrex
  require MatrexUtils

  @doc """
  Convert vector to mirrored image inverted vector(householder) and squared norm.
  Note that this function must select the non-zero vector.
  @param v: Target vector, must be Matrex of 1 row.

  iex(1)> HouseHolder.convert2householder_vec(Matrex.new([[1, 0, 0]]))
  {Matrex.new([[2, 0, 0]]), 4.0}
  iex(2)> HouseHolder.convert2householder_vec(Matrex.new([[1, 2, 3]]))
  {Matrex.new([[15, 2, 3]]), 238.0}
  iex(3)> HouseHolder.convert2householder_vec(Matrex.new([[-1, -2, -3]]))
  {Matrex.new([[-15, -2, -3]]), 238.0}
  """
  def convert2householder_vec(
        %Matrex{data: <<1::unsigned-integer-little-32, _::unsigned-integer-little-32, _::binary>>} = v) do
    with norm = MatrexUtils.sq_norm(v) do
      householder_vec = _create_householder_vec(v, norm)
      {householder_vec, MatrexUtils.sq_norm(householder_vec)}
    end
  end
  defp _create_householder_vec(%Matrex{} = v, norm) do
    with front_elem = v[1] do
      case front_elem > 0 do
        true  -> Matrex.set(v, 1, 1, front_elem + norm)
        false -> Matrex.set(v, 1, 1, front_elem - norm)
      end
    end
  end

  @doc """
  Execute house holder conversion for rows after already converted.
  @param x: Target matrix
  @param th_row: threshold row, rows below this value are not converted (not changed).
  @param householder_vec: house holder vector
  @param householder_norm: norm of house holder vector

  """
  def convert_with_householder(%Matrex{} = x, th_row, householder_vec, householder_norm) do
    # TODO
  end

  @doc """
  Calculate similar transformed matrix with house holder matrix.
  """
  def st(%Matrex{} = x, row_below) do
    # TODO
  end

  @doc """
  Calculate eigen values of 2x2 matrix for householder method.
  This calculation done by quadratic formula.

  iex(1)> HouseHolder.ev22(0, 0, 0, 0)
  0.0
  iex(2)> HouseHolder.ev22(2, 1, 1, 2)
  1.0
  iex(3)> HouseHolder.ev22(2, 2, 2, 4)
  5.23606797749979
  """
  def ev22(x11, x12, x21, x22) do
    d1 = x11 + x22
    d2 = x11 * x22 - x12 * x21
    a1 = (d1 + :math.sqrt(d1 * d1 - 4 * d2)) / 2
    a2 = (d1 - :math.sqrt(d1 * d1 - 4 * d2)) / 2
    case abs(x22 - a1) < abs(x22 - a2) do
      true  -> a1
      false -> a2
    end
  end
end
