using SearchGadgets, Test, Graphs
using SearchGadgets: detect_logic_gate

# c = ̅a
g = SimpleGraph(Edge.([(1,2)]))
detect_logic_gate(g; inputs=(1,), output=2)

# c = ̅a ∧ ̅b
g = SimpleGraph(Edge.([(1,2), (2,3), (3,4), (4,5), (5,1)]))
detect_logic_gate(g; inputs=(1,3), output=2)

# c = a ∨ b
g = SimpleGraph(Edge.([(1,2), (2,3), (3,4), (4,5), (5,1), (2,6)]))
detect_logic_gate(g; inputs=(1,3), output=6, weights=[1, 2, 1, 1, 1, 1])

# c = a ∧ b
g = SimpleGraph(Edge.([(1,2), (2,3), (3,4), (4,5), (5,1), (1,6), (3,7)]))
detect_logic_gate(g; inputs=(6,7), output=2, weights=[2, 1, 2, 1, 1, 1, 1])

