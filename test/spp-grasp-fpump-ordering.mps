NAME          FPUMPORDER
ROWS
 N  OBJ
 E  R1
 E  R2
 E  R3
COLUMNS
    x12       R1        1   R2        1   OBJ       1
    x23       R2        1   R3        1   OBJ       1
    x13       R1        1   R3        1   OBJ       1
    x123      R1        1   R2        1   OBJ       2
    x123      R3        1
RHS
    RHS1      R1        1   R2        1   R3        1
BOUNDS
 BV BND1      x12
 BV BND1      x23
 BV BND1      x13
 BV BND1      x123
ENDATA
