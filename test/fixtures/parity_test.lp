Minimize
  x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11 + x12
Subject To
  row1: x1 - 2 y1 = 1
  row2: x1 + x2 - 2 y2 = 0
  row3: x2 + x3 - 2 y3 = 1
  row4: x3 + x4 - 2 y4 = 0
  row5: x4 + x5 - 2 y5 = 1
  row6: x5 + x6 - 2 y6 = 0
  row7: x6 + x7 - 2 y7 = 1
  row8: x7 + x8 - 2 y8 = 0
  row9: x8 + x9 - 2 y9 = 1
  row10: x9 + x10 - 2 y10 = 0
  non_parity1: 3 x11 + x12 >= 1
  non_parity2: x10 + x11 - x12 <= 1
Bounds
  0 <= x1 <= 1
  0 <= x2 <= 1
  0 <= x3 <= 1
  0 <= x4 <= 1
  0 <= x5 <= 1
  0 <= x6 <= 1
  0 <= x7 <= 1
  0 <= x8 <= 1
  0 <= x9 <= 1
  0 <= x10 <= 1
  0 <= x11 <= 1
  0 <= x12 <= 1
  0 <= y1 <= 5
  0 <= y2 <= 5
  0 <= y3 <= 5
  0 <= y4 <= 5
  0 <= y5 <= 5
  0 <= y6 <= 5
  0 <= y7 <= 5
  0 <= y8 <= 5
  0 <= y9 <= 5
  0 <= y10 <= 5
Binary
  x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12
General
  y1 y2 y3 y4 y5 y6 y7 y8 y9 y10
End
