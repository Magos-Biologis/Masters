# const output = IOBuffer()
# using REPL
# const out_terminal = REPL.Terminals.TerminalBuffer(output)
# const basic_repl = REPL.BasicREPL(out_terminal)
# const basic_display = REPL.REPLDisplay(basic_repl)

# Base.pushdisplay(basic_display)

using Plots
gr()


test_plot = plot(t -> sin(t))

savefig(test_plot, "figs/test.png")
