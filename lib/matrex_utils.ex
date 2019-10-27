defmodule MatrexUtils do
  @moduledoc """
  Utility functions to supplement Matrex.
  """

  require Matrex
  @binary_per_data 4

  @doc """
  Adds two matrices which have either matched rows or columns
  """
  def add(%Matrex{} = first, %Matrex{} = second, alpha \\ 1.0, beta \\ 1.0), do: _add(first, second, alpha, beta)

  def _add(
        %Matrex{
          data:
            <<rows::unsigned-integer-little-32, columns::unsigned-integer-little-32, _data1::binary>>
        } = first,
        %Matrex{
          data:
            <<rows::unsigned-integer-little-32, columns::unsigned-integer-little-32, _data2::binary>>
        } = second,
        alpha,
        beta
      ), do: Matrex.add(first, second, alpha, beta)
  def _add(
        %Matrex{
          data:
            <<rows::unsigned-integer-little-32, columns1::unsigned-integer-little-32, _data1::binary>>
        } = first,
        %Matrex{
          data:
            <<rows::unsigned-integer-little-32, columns2::unsigned-integer-little-32, _data2::binary>>
        } = second,
        alpha,
        beta
      ) when columns1 > columns2 and rem(columns1, columns2) == 0 do
    second
    |> broad_cast(:columns, columns1)
    |> Matrex.add(first, alpha, beta)
  end
  def _add(
        %Matrex{
          data:
            <<rows::unsigned-integer-little-32, columns1::unsigned-integer-little-32, _data1::binary>>
        } = first,
        %Matrex{
          data:
            <<rows::unsigned-integer-little-32, columns2::unsigned-integer-little-32, _data2::binary>>
        } = second,
        alpha,
        beta
      ) when columns2 > columns1 and rem(columns2, columns1) == 0 do
    first
    |> broad_cast(:columns, columns2)
    |> Matrex.add(second, alpha, beta)
  end
  def _add(
        %Matrex{
          data:
            <<rows1::unsigned-integer-little-32, columns::unsigned-integer-little-32, _data1::binary>>
        } = first,
        %Matrex{
          data:
            <<rows2::unsigned-integer-little-32, columns::unsigned-integer-little-32, _data2::binary>>
        } = second,
        alpha,
        beta
      ) when rows1 > rows2 and rem(rows1, rows2) == 0 do
    second
    |> broad_cast(:rows, rows1)
    |> Matrex.add(first, alpha, beta)
  end
  def _add(
        %Matrex{
          data:
            <<rows1::unsigned-integer-little-32, columns::unsigned-integer-little-32, _data1::binary>>
        } = first,
        %Matrex{
          data:
            <<rows2::unsigned-integer-little-32, columns::unsigned-integer-little-32, _data2::binary>>
        } = second,
        alpha,
        beta
      ) when rows2 > rows1 and rem(rows2, rows1) == 0 do
    first
    |> broad_cast(:rows, rows2)
    |> Matrex.add(second, alpha, beta)
  end

  @doc """
  Broadcast matrix rows or columns to align to target shape.
  Only accepts broadcasting rows or columns at one time.
  """
  def broad_cast(
        %Matrex{
          data:
            <<rows::unsigned-integer-little-32, columns::unsigned-integer-little-32, body::binary>>},
        :columns,
        target_columns
      ) when rem(target_columns, columns) == 0 do
    new_body = 1..rows
               |> Enum.map(
                    fn index -> _parse_binary(columns, body, index)
                                |> List.duplicate(Kernel.div(target_columns, columns))
                                |> :erlang.list_to_binary()
                    end)
               |> :erlang.list_to_binary()
    %Matrex{data: <<rows::unsigned-integer-little-32,
              target_columns::unsigned-integer-little-32,
              new_body::binary>>}
  end
  def broad_cast(
        %Matrex{
          data:
            <<rows::unsigned-integer-little-32, columns::unsigned-integer-little-32, body::binary>>},
        :rows,
        target_rows
      ) when rem(target_rows, rows) == 0 do
    new_body = body
               |> List.duplicate(Kernel.div(target_rows, rows))
               |> :erlang.list_to_binary()
    %Matrex{data: <<target_rows::unsigned-integer-little-32,
              columns::unsigned-integer-little-32,
              new_body::binary>>}
  end

  @doc """
  Get sum of rows or cols and compose them to Matrex.
  """
  def sum(%Matrex{} = x, :rows) do
    <<rows::unsigned-integer-little-32,
      columns::unsigned-integer-little-32,
      body::binary>> = x.data
    new_body = 1..columns
               |> Enum.map(
                    fn col_index -> 1..rows
                                    |> Stream.map(fn row_index -> _parse_binary(columns, body, row_index, col_index) end)
                                    |> Stream.map(fn <<val::float-little-32>> -> val end)
                                    |> Enum.sum()
                                    |> (&<<&1::float-little-32>>).()
                    end)
               |> :erlang.list_to_binary()
    %Matrex{data: <<1::unsigned-integer-little-32,
              columns::unsigned-integer-little-32,
              new_body::binary>>}
  end
  def sum(%Matrex{} = x, :columns) do
    <<rows::unsigned-integer-little-32,
      columns::unsigned-integer-little-32,
      body::binary>> = x.data
    new_body = 1..rows
               |> Enum.map(
                    fn row_index -> 1..columns
                                    |> Stream.map(fn col_index -> _parse_binary(columns, body, row_index, col_index) end)
                                    |> Stream.map(fn <<val::float-little-32>> -> val end)
                                    |> Enum.sum()
                                    |> (&<<&1::float-little-32>>).()
                    end)
               |> :erlang.list_to_binary()
    %Matrex{data: <<rows::unsigned-integer-little-32,
              1::unsigned-integer-little-32,
              new_body::binary>>}
  end

  @doc """
  Fetch data of specified list of rows and compose them to Matrex.
  """
  def fetch(%Matrex{} = x, [_| _] = row_indices) do
    <<_::unsigned-integer-little-32,
      columns::unsigned-integer-little-32,
      body::binary>> = x.data
    new_body = row_indices
               |> Enum.map(fn index -> _parse_binary(columns, body, index) end)
               |> :erlang.list_to_binary()
    %Matrex{data: <<length(row_indices)::unsigned-integer-little-32,
              columns::unsigned-integer-little-32,
              new_body::binary>>}
  end

  defp _parse_binary(columns, body, row_index) do
    binary_part(body, (row_index - 1) * columns * @binary_per_data, columns * @binary_per_data)
  end
  defp _parse_binary(columns, body, row_index, col_index) do
    binary_part(body, ((row_index - 1) * columns + col_index - 1) * @binary_per_data, @binary_per_data)
  end

  @doc """
  Create a mesh grid 2d array in the rectangle.
  """
  def meshgrid([_| _] = x_range, [_| _] = y_range) do
    {
      List.duplicate(x_range, length(y_range)),
      y_range |> Enum.map(fn elem -> List.duplicate(elem, length(x_range)) end)
    }
  end

  @doc """
  Create a list in the range with the step.
  """
  def arrange(s, e, step) when s < e do
    Stream.iterate(s, &(&1 + step))
    |> Enum.take(Kernel.trunc((e - s) / step))
  end

  @doc """
  Return argmax for row/column indices.
  """
  def argmax(%Matrex{} = x, :rows) do
    1..x[:columns]
    |> Enum.map(fn column_index -> Matrex.column(x, column_index)[:argmax] end)
  end
  def argmax(%Matrex{} = x, :columns) do
    1..x[:rows]
    |> Enum.map(fn row_index -> x[row_index][:argmax] end)
  end

  @doc """
  Create a matrex of flattened list.
  """
  def flattened([_ | _] = list) do
    list
    |> List.flatten()
    |> new()
    |> Matrex.transpose()
  end

  @doc """
  Create a matrex braced with []
  """
  def new([head | _] = list) when is_number(head) do
    Matrex.new([list])
  end

  @doc """
  Normalize matrex with L2 norm, that is sqrt sum of squared all elements.
  """
  def l2_normalize(%Matrex{} = x, eps \\ 1.0e-8) do
    Matrex.divide(x, _l2norm(x) + eps)
  end

  def _l2norm(%Matrex{} = x) do
    x
    |> Matrex.square()
    |> Matrex.sum()
    |> :math.sqrt()
  end

  @doc """
  Create a matrex that has a same shape of given matrex.
  """
  def zeros_like(
        %Matrex{data: <<rows::unsigned-integer-little-32, columns::unsigned-integer-little-32, _::binary>>}) do
    Matrex.zeros(rows, columns)
  end

  @doc """
  Calculate squared norm of specified vector, which must be 1 row of Matrex.

  iex(1)> Matrex.sq_norm(Matrex.new([0, 0, 0]))
  0
  iex(2)> Matrex.sq_norm(Matrex.new([1, 2, 3]))
  14
  iex(3)> Matrex.sq_norm(Matrex.new([-1, -2, -3]))
  14
  """
  def sq_norm(%Matrex{data: <<1::unsigned-integer-little-32, _::unsigned-integer-little-32, _::binary>>} = v) do
    Matrex.dot_nt(v, v) |> Matrex.scalar()
  end
end
