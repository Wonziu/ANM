using Plots

f(x) = x^2 - 3

default(tickfont = (12, :orange), framestyle = :zerolines)
plot(f, 0, 6, linewidth = 2, label = "f(x) = x^2 - 3", legend = :topleft)
annotate!((1, -2, "a"), (5, -2, "b"), (sqrt(3), -2, "alfa"), (3, -2, "c"))
scatter!([(1, 0), (5, 0), (sqrt(3), 0), (3, 0)], c = :red, label = false)

savefig("bisectionPlot.pdf")