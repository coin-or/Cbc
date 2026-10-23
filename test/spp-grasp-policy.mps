NAME          SPP_POLICY
ROWS
 N  OBJ
 E  COVER
COLUMNS
    X1        OBJ       1.0     COVER     1.0
    X2        OBJ       2.0     COVER     1.0
    SLACK     OBJ      10.0     COVER     1.0
RHS
    RHS1      COVER     1.0
BOUNDS
 BV BND1      X1
 BV BND1      X2
 BV BND1      SLACK
ENDATA
