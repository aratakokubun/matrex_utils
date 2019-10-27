defmodule MatrexUtils.MixProject do
  use Mix.Project

  def project do
    [
      app: :matrex_utils,
      version: "0.0.1",
      elixir: "~> 1.9",
      description: "Library to supplement Matrex",
      start_permanent: Mix.env() == :prod,
      package: [
        maintainers: ["kkbnart"],
        licenses: ["MIT"],
        links: %{"GitHub" => "https://github.com/aratakokubun"}
      ],
      deps: deps()
    ]
  end

  # Run "mix help compile.app" to learn about applications.
  def application do
    [
      extra_applications: [:logger]
    ]
  end

  # Run "mix help deps" to learn about dependencies.
  defp deps do
    [
      {:matrex, "~> 0.6.0"},
      {:ex_doc, "~> 0.21.0", only: :dev, runtime: false},
      {:dialyxir, "~> 1.0.0-rc.4", only: [:dev], runtime: false}
    ]
  end
end
