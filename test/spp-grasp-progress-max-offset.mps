NAME          SPPPROGRESS
ROWS
 N  OBJ
 E  R1
 E  R2
 E  R3
COLUMNS
    x12       R1        1   R2        1   OBJ       0
    x23       R2        1   R3        1   OBJ       0
    x13       R1        1   R3        1   OBJ       0
    s1        R1        1   OBJ       1
    s2        R2        1   OBJ       1
    s3        R3        1   OBJ       1
RHS
    RHS1      R1        1   R2        1   R3        1
    RHS1      OBJ      -5
BOUNDS
 BV BND1      x12
 BV BND1      x23
 BV BND1      x13
 BV BND1      s1
 BV BND1      s2
 BV BND1      s3
ENDATA
