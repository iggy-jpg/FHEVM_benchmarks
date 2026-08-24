# Integrated Compare implementation

`compsigneval.hpp` is the reusable integration of
`testfiles/Compare/compsigneval.cpp`. It retains the three-layer composite
minimax sign approximation and adds integer `GreaterThan` and `Equal` helpers.

The executable target `openfhe_requested_integer_benchmark` uses these helpers
for 4/8-bit EQ and GT. It also contains the requested 4/8/16/32-bit ADD and
4/8/16-bit MULT measurements.
